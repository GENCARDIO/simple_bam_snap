#!/usr/bin/env python3
"""LocusSnap: an IGV-like genomic snapshot generator built directly on pysam.

Given a BAM and a region, renders every overlapping alignment from its own
parsed CIGAR (and SA tag) rather than re-parsing `samtools tview` text. Reads
can be packed IGV-style or expanded one-per-row and sorted by a chosen metric
- gap_length (indel length from CIGAR, or the implied gap from a split/
supplementary alignment) chief among them - which is what makes it useful for
answering questions like "does aligner A produce more/longer gapped
alignments than aligner B for this indel".

Pass --bam2 to render both BAMs stacked in one comparison image plus a
console summary table.
"""
import argparse
import logging
import os
import re
import sys

import pysam

from locus_snap.annotations import (
    ANNOTATION_DISPLAY_MODES,
    PRIMARY_ISOFORM_MODES,
    AnnotationSource,
    build_annotation_sources,
    build_baf_sources,
    build_custom_annotation_sources,
)
from locus_snap.config import load_config
from locus_snap.downsample import DEFAULT_MAX_ALIGNMENT_DEPTH
from locus_snap.layout import DISPLAY_MODES, HAPLOTYPE_VIEWS, SORT_KEYS
from locus_snap.mate_window import MATE_WINDOW_SOURCES
from locus_snap.read_model import ONLY_TYPES
from locus_snap.refseq import (
    REFSEQ_LABELS, detect_human_assembly, ensure_refseq, normalize_assembly,
)
from locus_snap.render import DEFAULT_COVERAGE_VAF_THRESHOLD, DEFAULT_MAX_REFERENCE_SPAN
from locus_snap.snapshot import OUTPUT_FORMATS, BamSnapshot, compare_snapshots

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("locus_snap")

REGION_RE = re.compile(r"^(?P<chrom>[\w.\-]+):(?P<start>[\d,]+)-(?P<end>[\d,]+)$")
PREFERENCE_ALIASES = {
    "show_alignments": ("no_alignments", True),
    "show_coverage": ("no_coverage", True),
    "show_ideogram": ("no_ideogram", True),
    "pair_colors": ("no_pair_colors", True),
    "mapq_shading": ("no_mapq_shading", True),
    "annotate_gap": ("no_annotate", True),
    "include_supplementary": ("exclude_supplementary", True),
}


def parse_region(region: str, flank: int = 0):
    match = REGION_RE.match(region.strip())
    if not match:
        raise ValueError(
            f"Invalid --region '{region}'. Expected format chrom:start-end (e.g. chr9:101867500-101867650)."
        )
    chrom = match.group("chrom")
    start = int(match.group("start").replace(",", ""))
    end = int(match.group("end").replace(",", ""))
    if end <= start:
        raise ValueError(f"--region end ({end}) must be greater than start ({start}).")
    # user-facing coordinates are 1-based inclusive; internally we use 0-based half-open
    start0 = max(0, start - 1 - flank)
    end0 = end + flank
    return chrom, start0, end0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "LocusSnap: IGV-like genomic snapshots with sortable alignment layouts."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--bam", action="append", required=True, metavar="BAM",
        help="Indexed BAM input; repeat to stack any number of sample panels.",
    )
    parser.add_argument(
        "--bam2",
        help="Legacy second BAM option; equivalent to repeating --bam.",
    )
    parser.add_argument("--region", required=True, help="Region as chrom:start-end (1-based, inclusive).")
    parser.add_argument("--fasta", help="Reference FASTA (indexed or indexable). Enables mismatch/base coloring.")
    parser.add_argument(
        "--max_reference_span", type=int, default=DEFAULT_MAX_REFERENCE_SPAN, metavar="BP",
        help=("Show the coloured FASTA reference-base track only for windows at or below "
              "this span; set to 0 to hide the track while retaining mismatch detection."),
    )
    parser.add_argument("--flank", type=int, default=0, help="Extra bp of context padded on each side of --region.")

    parser.add_argument("--output_dir", default=".", help="Output directory.")
    parser.add_argument(
        "--output_name",
        help=("Output file name. A supported extension selects the format; otherwise "
              ".png is appended by default."),
    )
    parser.add_argument(
        "--output_format", choices=OUTPUT_FORMATS,
        help="Output image format. Overrides a supported extension in --output_name.",
    )
    parser.add_argument("--label1", help="Label for --bam in comparison mode (default: file stem).")
    parser.add_argument("--label2", help="Label for --bam2 in comparison mode (default: file stem).")
    parser.add_argument(
        "--sample_label", action="append", metavar="LABEL",
        help="Label for the corresponding repeated --bam, supplied in the same order.",
    )
    parser.add_argument(
        "--vcf_companion", action="append", metavar="VCF",
        help=("VCF companion for each repeated --bam, in the same order. Repeat once per "
              "BAM and use 'none' for a sample without a VCF. Compressed VCFs require an index."),
    )
    parser.add_argument(
        "--config", help=(
            "YAML file defining default preferences, colours, and rendering styles. "
            "Explicit command-line options override YAML preferences."
        )
    )
    parser.add_argument(
        "--metrics_tsv", help="Write per-read computed metrics (gap length, mismatches, SA, ...) to this TSV path."
    )
    parser.add_argument(
        "--metrics_tsv2", help="TSV path for --bam2's per-read metrics (comparison mode only)."
    )

    track_group = parser.add_argument_group("genomic and quantitative tracks")
    track_group.add_argument(
        "--track", action="append", metavar="PATH",
        help=("BED, GFF/GFF3, GTF, VCF, narrowPeak/broadPeak, signal, SEG, "
              "bedGraph, or log2/CNV track; "
              "repeat for multiple tracks. "
              "Compressed tracks require a .tbi or .csi tabix index."),
    )
    track_group.add_argument(
        "--track_label", action="append", metavar="LABEL",
        help="Optional label for the corresponding --track, supplied in the same order.",
    )
    track_group.add_argument(
        "--custom_track", action="append", nargs="+", metavar="SPEC",
        help=("Define one track as the recommended CSV value "
              "'FILE,TYPE,NAME,COLOR[,DISPLAY[,HEIGHT_IN]]'; repeat as needed. "
              "The legacy FILE TYPE NAME COLOR [DISPLAY [HEIGHT_IN]] form remains accepted. "
              "TYPE is bed, gff, gff3, gtf, vcf, narrowpeak, broadpeak, peak, "
              "signal, seg, bedgraph, log2, cnv, or auto. "
              "COLOR accepts hex, R,G,B, or rgb(R,G,B). "
              "DISPLAY optionally overrides --track_display with collapse, pack, expand, "
              "or density; HEIGHT_IN optionally sets this track's physical height. "
              "BGZF files require .tbi/.csi and are fetched by genomic region with tabix."),
    )
    track_group.add_argument(
        "--track_display", choices=ANNOTATION_DISPLAY_MODES, default="pack",
        help=("Default annotation layout: collapse merges transcript isoforms by gene; "
              "pack shares rows between non-overlapping models; expand uses one row "
              "per transcript; density bins interval overlap into a compact histogram."),
    )
    track_group.add_argument(
        "--primary_isoforms", choices=PRIMARY_ISOFORM_MODES, default="all",
        help=("Gene-transcript selection: all keeps every isoform; prefer uses the "
              "best annotated primary tier per gene but retains all isoforms when no "
              "marker exists; only removes genes without a primary marker."),
    )
    track_group.add_argument(
        "--baf_vcf", action="append", metavar="VCF",
        help=("Add a BAF/LOH track from heterozygous biallelic SNVs in a genotype VCF; "
              "repeat for multiple tracks. Compressed VCF/BCF files require an index."),
    )
    track_group.add_argument(
        "--baf_sample", action="append", metavar="SAMPLE",
        help=("VCF sample for the corresponding --baf_vcf; defaults to the first sample "
              "and is supplied in the same order."),
    )
    track_group.add_argument(
        "--baf_track_label", action="append", metavar="LABEL",
        help="Optional label for the corresponding --baf_vcf, supplied in the same order.",
    )

    parser.add_argument(
        "--layout", choices=["pack", "expand"], default="pack",
        help="'pack': IGV-style greedy row packing. 'expand': one row per read, ordered by --sort_by "
             "(use this to rank alignments by gap length).",
    )
    parser.add_argument(
        "--display_mode", choices=DISPLAY_MODES, default="expand",
        help=("Alignment track appearance: collapse overlays all reads in one row; "
              "expand uses normal-height rows; squish uses compact rows."),
    )
    parser.add_argument(
        "--no_alignments", action="store_true",
        help=("Hide the BAM alignment rows and their legend, for track-only figures "
              "such as ChIP-seq signal profiles."),
    )
    parser.add_argument(
        "--view_as_pairs", action="store_true",
        help=("Place visible primary mates on the same row and link them across "
              "their genomic gap, like IGV's View as pairs option."),
    )
    haplotype_group = parser.add_argument_group("haplotype-aware alignments")
    haplotype_group.add_argument(
        "--haplotype_view", choices=HAPLOTYPE_VIEWS, default="none",
        help=("Haplotype display: none keeps ordinary colours; color colours reads by HP; "
              "split also separates HP values into labelled lanes."),
    )
    haplotype_group.add_argument(
        "--haplotype_filter", nargs="+", metavar="HP",
        help=("Keep selected HP values, such as 1 2; use 'untagged' to retain reads "
              "without the haplotype tag."),
    )
    haplotype_group.add_argument(
        "--haplotype_tag", default="HP", metavar="TAG",
        help="Two-character SAM tag containing the read haplotype.",
    )
    haplotype_group.add_argument(
        "--phase_set_tag", default="PS", metavar="TAG",
        help="Two-character SAM tag containing the phase-set identifier.",
    )
    parser.add_argument(
        "--sort_by", choices=sorted(SORT_KEYS), default="gap_length",
        help=("Sort/priority key. gap_length = max(CIGAR indel length, SA-implied "
              "split gap); base groups reads by their nucleotide at "
              "--sort_base_position and prioritises non-reference alleles."),
    )
    parser.add_argument(
        "--sort_base_position", type=int, metavar="POS",
        help=("1-based genomic position used by --sort_by base. Defaults to the "
              "middle base of the requested window."),
    )
    parser.add_argument(
        "--sort_order", choices=["desc", "asc"], default="desc",
        help="desc puts the largest gap_length (or chosen key) first/top.",
    )
    parser.add_argument("--max_rows", type=int, help="Cap the number of rows drawn (kept rows are highest priority).")
    parser.add_argument(
        "--max_alignment_depth", type=int, default=DEFAULT_MAX_ALIGNMENT_DEPTH,
        help=("Downsample alignment spans above this overlapping depth. Coverage, summaries, and TSVs "
              "still use all reads; set to 0 to disable."),
    )
    parser.add_argument("--long_gap_threshold", type=int, default=10,
                         help="bp threshold for the 'long gap' count in the summary stats.")

    parser.add_argument("--min_mapq", type=int, default=0, help="Skip alignments below this MAPQ.")
    parser.add_argument("--include_secondary", action="store_true", help="Include secondary alignments.")
    parser.add_argument("--exclude_supplementary", action="store_true",
                         help="Exclude supplementary alignments (included by default - needed for SA-gap evidence).")
    parser.add_argument("--include_duplicates", action="store_true", help="Include reads flagged as PCR/optical duplicates.")

    parser.add_argument(
        "--only", nargs="+", choices=sorted(ONLY_TYPES), metavar="TYPE",
        help=f"Isolate only reads matching (OR of) these categories: {', '.join(sorted(ONLY_TYPES))}. "
             "discordant = abnormal pair orientation/insert size/inter-chromosomal.",
    )
    parser.add_argument("--min_softclip", type=int, default=1,
                         help="Minimum soft-clip length (bp) counted as 'soft-clipped' for --only softclip and summary stats.")
    parser.add_argument("--insert_size_sigma", type=float, default=3.0,
                         help="FR-pair insert size beyond median +/- N*robust-stdev (estimated from the window's own "
                              "reads) is flagged small_insert/large_insert (discordant).")

    parser.add_argument("--no_pair_colors", action="store_true",
                         help="Disable IGV-style discordant-pair coloring; use the normal alignment fill for every read.")
    parser.add_argument("--no_mapq_shading", action="store_true",
                         help="Disable lightening read fill color for lower MAPQ.")
    parser.add_argument("--mapq_cap", type=int, default=60,
                         help="MAPQ at which a read's fill reaches full opacity under MAPQ shading.")

    mate_group = parser.add_argument_group("mate view")
    mate_group.add_argument(
        "--mate_view", action="store_true",
        help="Draw the requested region and an automatically selected mate locus as two adjacent panels."
    )
    mate_group.add_argument(
        "--mate_window_source", choices=MATE_WINDOW_SOURCES, default="discordant",
        help=("Evidence used to select the mate locus: mapped mates of discordant reads; "
              "SA entries from split reads; or mapped mates of soft-clipped reads."),
    )
    mate_group.add_argument(
        "--mate_window_size", type=int,
        help="Mate-panel span in bp (default: the final --region span, including --flank).",
    )

    parser.add_argument("--no_ideogram", action="store_true",
                        help="Hide the chromosome overview and red current-window marker.")
    parser.add_argument(
        "--center_guide", action="store_true",
        help="Show a vertical guide through the center of each genomic panel.",
    )
    parser.add_argument(
        "--genome", choices=["auto", "hg19", "grch37", "hg38", "grch38", "none"],
        default="auto",
        help=("Cytoband assembly for the ideogram. auto selects bundled UCSC hg19/hg38 "
              "only from exact BAM contig-length matches; none uses a plain chromosome bar."),
    )
    parser.add_argument(
        "--cytoband_file",
        help=("Custom UCSC five-column cytoBand text file, optionally gzip-compressed. "
              "Overrides --genome."),
    )
    parser.add_argument(
        "--refseq",
        choices=["auto", "hg19", "grch37", "hg38", "grch38", "none"],
        default="auto",
        help=("Automatically load the indexed NCBI RefSeq isoform track for "
              "the BAM assembly; choose an assembly explicitly or none to disable it."),
    )
    parser.add_argument(
        "--refseq_dir", metavar="DIR",
        help=("Cache directory for automatically downloaded RefSeq GFF and tabix index. "
              "Defaults to data/ucsc/refseq in this project."),
    )
    parser.add_argument("--no_coverage", action="store_true", help="Hide the coverage track.")
    rna_group = parser.add_argument_group("RNA-seq / sashimi")
    rna_group.add_argument(
        "--sashimi", action="store_true",
        help="Draw count-labelled splice-junction arcs from CIGAR N operations.",
    )
    rna_group.add_argument(
        "--min_junction_reads", type=int, default=1, metavar="N",
        help="Minimum supporting alignments required to draw a splice junction.",
    )
    rna_group.add_argument(
        "--sashimi_strand", choices=("combined", "split"), default="combined",
        help="Combine strands above the baseline or mirror plus/minus junction arcs.",
    )
    parser.add_argument(
        "--coverage_vaf_threshold", type=float, default=DEFAULT_COVERAGE_VAF_THRESHOLD,
        metavar="FRACTION",
        help=("Colour SNV allele counts in the coverage track when their observed allele "
              "fraction is greater than this value; requires --fasta."),
    )
    parser.add_argument(
        "--min_baseq", type=int, default=0, metavar="Q",
        help="Minimum base quality used for coverage-track SNV depth and evidence.",
    )
    parser.add_argument(
        "--min_variant_mapq", type=int, default=0, metavar="Q",
        help="Minimum read MAPQ used for coverage-track SNV depth and evidence.",
    )
    parser.add_argument(
        "--show_variant_counts", action="store_true",
        help=("Label qualifying coverage SNVs with ALT/depth, VAF, strand counts, "
              "mean base quality, and mean MAPQ when base spacing permits."),
    )
    parser.add_argument("--no_annotate", action="store_true", help="Hide the per-row gap-length annotation (expand layout).")
    parser.add_argument(
        "--show_indel_lengths", action="store_true",
        help="Label CIGAR insertions and deletions with their lengths; hidden by default.",
    )
    parser.add_argument("--fig_width", type=float, default=14.0, help="Figure width in inches.")
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="Raster output resolution in dots per inch (also applies to rasterized vector elements).",
    )

    return parser


def apply_config_preferences(parser: argparse.ArgumentParser, preferences: dict) -> None:
    """Validate YAML preferences against argparse and install them as defaults."""
    actions = {
        action.dest: action for action in parser._actions
        if action.dest not in ("help", "config")
    }
    unknown = sorted(set(preferences) - set(actions) - set(PREFERENCE_ALIASES))
    if unknown:
        raise ValueError(f"Unknown preference key(s): {', '.join(unknown)}.")
    validated = {}
    for key, value in preferences.items():
        destination, invert = PREFERENCE_ALIASES.get(key, (key, False))
        if destination in validated:
            raise ValueError(
                f"Preference {key} duplicates another setting for {destination}."
            )
        action = actions[destination]
        if action.nargs == 0:
            if not isinstance(value, bool):
                raise ValueError(f"Preference {key} must be true or false.")
            validated[destination] = not value if invert else value
            continue
        values = value if isinstance(value, list) else [value]
        converted = []
        for item in values:
            try:
                parsed = action.type(item) if action.type else item
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid value for preference {key}: {item!r}.") from exc
            if action.choices is not None and parsed not in action.choices:
                choices = ", ".join(str(choice) for choice in action.choices)
                raise ValueError(
                    f"Invalid value for preference {key}: {parsed!r}; choose from {choices}."
                )
            converted.append(parsed)
        expects_list = action.nargs in ("+", "*") or isinstance(value, list)
        validated[destination] = converted if expects_list else converted[0]
    parser.set_defaults(**validated)


def main(argv=None) -> int:
    command_args = list(argv) if argv is not None else sys.argv[1:]
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config")
    preliminary, remaining = config_parser.parse_known_args(command_args)
    del remaining
    try:
        visual_config = load_config(preliminary.config)
        parser = build_parser()
        apply_config_preferences(parser, visual_config["preferences"])
    except ValueError as exc:
        log.error(str(exc))
        return 1
    args = parser.parse_args(command_args)

    bam_paths = list(args.bam)
    if args.bam2:
        bam_paths.append(args.bam2)
    sample_labels = list(args.sample_label or [])
    if sample_labels and (args.label1 or args.label2):
        log.error("--sample_label cannot be combined with legacy --label1/--label2.")
        return 1
    if len(sample_labels) > len(bam_paths):
        log.error("More --sample_label values were supplied than --bam files.")
        return 1
    while len(sample_labels) < len(bam_paths):
        sample_labels.append(None)
    if not args.sample_label:
        if args.label1:
            sample_labels[0] = args.label1
        if args.label2 and len(sample_labels) > 1:
            sample_labels[1] = args.label2
    for index, bam_path in enumerate(bam_paths):
        if not sample_labels[index]:
            sample_labels[index] = os.path.splitext(os.path.basename(bam_path))[0]

    companion_vcfs = []
    for value in args.vcf_companion or []:
        companion_vcfs.append(None if value.lower() in ("none", ".", "na") else value)
    if companion_vcfs and len(companion_vcfs) != len(bam_paths):
        log.error(
            "Supply exactly one --vcf_companion per --bam; use 'none' for missing companions."
        )
        return 1
    if not companion_vcfs:
        companion_vcfs = [None] * len(bam_paths)

    try:
        chrom, start, end = parse_region(args.region, flank=args.flank)
    except ValueError as exc:
        log.error(str(exc))
        return 1

    if args.mate_view and len(bam_paths) > 1:
        log.error("--mate_view cannot be combined with multiple BAM panels.")
        return 1
    if args.no_alignments and (args.mate_view or len(bam_paths) > 1):
        log.error("--no_alignments currently supports single-locus, single-BAM figures only.")
        return 1
    if args.mate_window_size is not None and args.mate_window_size <= 0:
        log.error("--mate_window_size must be greater than zero.")
        return 1
    if args.max_alignment_depth < 0:
        log.error("--max_alignment_depth cannot be negative (use 0 to disable downsampling).")
        return 1
    if args.min_junction_reads < 1:
        log.error("--min_junction_reads must be at least one.")
        return 1
    if args.max_reference_span < 0:
        log.error("--max_reference_span cannot be negative (use 0 to hide the reference track).")
        return 1
    if args.fig_width <= 0:
        log.error("--fig_width must be greater than zero.")
        return 1
    if args.dpi <= 0:
        log.error("--dpi must be greater than zero.")
        return 1
    sort_base_position = None
    if args.sort_base_position is not None:
        if args.sort_base_position < 1:
            log.error("--sort_base_position must be a positive 1-based coordinate.")
            return 1
        if args.sort_by != "base":
            log.error("--sort_base_position requires --sort_by base.")
            return 1
        sort_base_position = args.sort_base_position - 1
        if not start <= sort_base_position < end:
            log.error("--sort_base_position must fall inside the requested region.")
            return 1
    if not 0 <= args.coverage_vaf_threshold <= 1:
        log.error("--coverage_vaf_threshold must be between 0 and 1.")
        return 1
    if args.min_baseq < 0 or args.min_variant_mapq < 0:
        log.error("--min_baseq and --min_variant_mapq cannot be negative.")
        return 1
    for option, tag in (("--haplotype_tag", args.haplotype_tag), ("--phase_set_tag", args.phase_set_tag)):
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", tag):
            log.error("%s must be a two-character SAM tag.", option)
            return 1

    haplotype_filter = []
    for value in args.haplotype_filter or []:
        haplotype_filter.append("untagged" if value.lower() == "untagged" else value)

    alignment_colors = visual_config["alignment_colors"]

    try:
        annotation_sources = build_annotation_sources(
            args.track, args.track_label, display_mode=args.track_display,
            primary_isoforms=args.primary_isoforms,
            track_colors=visual_config["track_colors"],
        )
        annotation_sources.extend(build_custom_annotation_sources(
            args.custom_track, default_display_mode=args.track_display,
            primary_isoforms=args.primary_isoforms,
            track_colors=visual_config["track_colors"],
        ))
        annotation_sources.extend(build_baf_sources(
            args.baf_vcf, labels=args.baf_track_label, samples=args.baf_sample,
            track_colors=visual_config["track_colors"],
        ))
    except (OSError, ValueError) as exc:
        log.error(str(exc))
        return 1

    for bam_path in bam_paths:
        if not os.path.isfile(bam_path):
            log.error("Cannot find --bam file: %s", bam_path)
            return 1
    if args.fasta and not os.path.isfile(args.fasta):
        log.error("Cannot find --fasta file: %s", args.fasta)
        return 1

    if args.refseq != "none":
        selected_refseq = args.refseq
        detected_assemblies = set()
        try:
            for bam_path in bam_paths:
                with pysam.AlignmentFile(bam_path, "rb") as bam_file:
                    contig_lengths = dict(zip(bam_file.references, bam_file.lengths))
                detected = detect_human_assembly(contig_lengths)
                if detected:
                    detected_assemblies.add(detected)
            if len(detected_assemblies) > 1:
                assemblies = ", ".join(sorted(detected_assemblies))
                raise ValueError(
                    f"BAM inputs use different detected assemblies ({assemblies}); "
                    "automatic RefSeq cannot be shared."
                )
            detected = next(iter(detected_assemblies), None)
            if selected_refseq == "auto":
                selected_refseq = detected
            else:
                selected_refseq = normalize_assembly(selected_refseq)
                if detected and selected_refseq != detected:
                    raise ValueError(
                        f"Requested RefSeq {selected_refseq} does not match the BAM's "
                        f"detected {detected} assembly."
                    )
            if selected_refseq:
                refseq_path = ensure_refseq(selected_refseq, args.refseq_dir)
                annotation_sources.insert(0, AnnotationSource(
                    str(refseq_path), label=REFSEQ_LABELS[selected_refseq], kind="gff",
                    display_mode=args.track_display,
                    primary_isoforms=args.primary_isoforms,
                    track_colors=visual_config["track_colors"],
                ))
                log.info("Loaded default gene track: %s", REFSEQ_LABELS[selected_refseq])
            else:
                log.warning(
                    "Could not identify hg19/GRCh37 or hg38/GRCh38 from BAM contig "
                    "lengths; the default RefSeq track was not loaded."
                )
        except (OSError, ValueError) as exc:
            log.error(str(exc))
            return 1

    common_kwargs = dict(
        min_mapq=args.min_mapq,
        include_secondary=args.include_secondary,
        include_supplementary=not args.exclude_supplementary,
        include_duplicates=args.include_duplicates,
        max_rows=args.max_rows,
        show_ideogram=not args.no_ideogram,
        show_center_guide=args.center_guide,
        genome=args.genome,
        cytoband_file=args.cytoband_file,
        max_reference_span=args.max_reference_span,
        show_coverage=not args.no_coverage,
        show_sashimi=args.sashimi,
        min_junction_reads=args.min_junction_reads,
        sashimi_strand=args.sashimi_strand,
        coverage_vaf_threshold=args.coverage_vaf_threshold,
        min_baseq=args.min_baseq,
        min_variant_mapq=args.min_variant_mapq,
        show_variant_counts=args.show_variant_counts,
        haplotype_view=args.haplotype_view,
        haplotype_filter=haplotype_filter,
        haplotype_tag=args.haplotype_tag,
        phase_set_tag=args.phase_set_tag,
        annotate_gap=not args.no_annotate,
        show_indel_lengths=args.show_indel_lengths,
        fig_width=args.fig_width,
        dpi=args.dpi,
        long_gap_threshold=args.long_gap_threshold,
        layout=args.layout,
        display_mode=args.display_mode,
        sort_by=args.sort_by,
        sort_base_position=sort_base_position,
        descending=(args.sort_order == "desc"),
        only_types=args.only,
        min_softclip=args.min_softclip,
        insert_size_sigma=args.insert_size_sigma,
        pair_colors=not args.no_pair_colors,
        shade_by_mapq=not args.no_mapq_shading,
        mapq_cap=args.mapq_cap,
        alignment_colors=alignment_colors,
        visual_config=visual_config,
        max_alignment_depth=args.max_alignment_depth,
        annotation_sources=annotation_sources,
        view_as_pairs=args.view_as_pairs,
    )

    if len(bam_paths) > 1:
        try:
            out_path, table = compare_snapshots(
                bam1=bam_paths[0], bam2=bam_paths[1], chrom=chrom, start=start, end=end,
                fasta=args.fasta, output_dir=args.output_dir, output_name=args.output_name,
                output_format=args.output_format,
                label1=sample_labels[0], label2=sample_labels[1],
                additional_bams=bam_paths[2:], additional_labels=sample_labels[2:],
                companion_vcfs=companion_vcfs,
                metrics_tsv_1=args.metrics_tsv, metrics_tsv_2=args.metrics_tsv2,
                **common_kwargs,
            )
        except (OSError, ValueError) as exc:
            log.error(str(exc))
            return 1
        log.info("Wrote comparison snapshot: %s", out_path)
        print(table)
        return 0

    if companion_vcfs[0]:
        try:
            annotation_sources.append(AnnotationSource(
                companion_vcfs[0], label=f"{sample_labels[0]} variants", kind="vcf",
                display_mode="collapse", track_colors=visual_config["track_colors"],
            ))
        except (OSError, ValueError) as exc:
            log.error(str(exc))
            return 1

    snap = BamSnapshot(
        bam=bam_paths[0], chrom=chrom, start=start, end=end, fasta=args.fasta,
        output_dir=args.output_dir, output_name=args.output_name,
        output_format=args.output_format,
        label=sample_labels[0], mate_view=args.mate_view,
        mate_window_source=args.mate_window_source,
        mate_window_size=args.mate_window_size,
        show_alignments=not args.no_alignments,
        **common_kwargs,
    )
    try:
        summary = snap.snap(metrics_tsv=args.metrics_tsv)
    except (OSError, ValueError) as exc:
        log.error(str(exc))
        return 1
    log.info("Wrote snapshot: %s", snap.output_path)
    if snap.downsampled_reads:
        log.info(
            "Downsampled %d alignment(s) above the %dx display-depth cap; full coverage and metrics retained.",
            snap.downsampled_reads, args.max_alignment_depth,
        )
    if snap.mate_window:
        mate = snap.mate_window
        log.info(
            "Selected mate window from %d %s candidate(s): %s:%d-%d",
            mate.candidate_count, mate.source, mate.chrom, mate.start + 1, mate.end,
        )
    print(
        f"{summary.n_reads} reads | gapped: {summary.n_gapped} ({summary.pct_gapped:.1f}%) | "
        f"max gap: {summary.max_gap}bp | split (SA): {summary.n_with_sa} | "
        f"discordant: {summary.n_discordant} ({summary.pct_discordant:.1f}%) | "
        f"soft-clipped: {summary.n_softclipped}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
