#!/usr/bin/env python3
"""simple_bam_snap: an IGV-like snapshot generator built directly on pysam.

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

from src.layout import SORT_KEYS
from src.read_model import ONLY_TYPES
from src.snapshot import BamSnapshot, compare_snapshots

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("simple_bam_snap")

REGION_RE = re.compile(r"^(?P<chrom>[\w.\-]+):(?P<start>[\d,]+)-(?P<end>[\d,]+)$")


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
        description="IGV-like genomic snapshot generator with sortable alignment layout.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", required=True, help="Input BAM file (indexed).")
    parser.add_argument("--bam2", help="Optional second BAM for a side-by-side comparison snapshot.")
    parser.add_argument("--region", required=True, help="Region as chrom:start-end (1-based, inclusive).")
    parser.add_argument("--fasta", help="Reference FASTA (indexed or indexable). Enables mismatch/base coloring.")
    parser.add_argument("--flank", type=int, default=0, help="Extra bp of context padded on each side of --region.")

    parser.add_argument("--output_dir", default=".", help="Output directory.")
    parser.add_argument("--output_name", help="Output file name (.png appended if missing).")
    parser.add_argument("--label1", help="Label for --bam in comparison mode (default: file stem).")
    parser.add_argument("--label2", help="Label for --bam2 in comparison mode (default: file stem).")
    parser.add_argument(
        "--metrics_tsv", help="Write per-read computed metrics (gap length, mismatches, SA, ...) to this TSV path."
    )
    parser.add_argument(
        "--metrics_tsv2", help="TSV path for --bam2's per-read metrics (comparison mode only)."
    )

    parser.add_argument(
        "--layout", choices=["pack", "expand"], default="pack",
        help="'pack': IGV-style greedy row packing. 'expand': one row per read, ordered by --sort_by "
             "(use this to rank alignments by gap length).",
    )
    parser.add_argument(
        "--sort_by", choices=sorted(SORT_KEYS), default="gap_length",
        help="Sort/priority key. gap_length = max(CIGAR indel length, SA-implied split gap).",
    )
    parser.add_argument(
        "--sort_order", choices=["desc", "asc"], default="desc",
        help="desc puts the largest gap_length (or chosen key) first/top.",
    )
    parser.add_argument("--max_rows", type=int, help="Cap the number of rows drawn (kept rows are highest priority).")
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
                         help="Disable IGV-style discordant-pair coloring; fall back to plain forward/reverse strand fill.")
    parser.add_argument("--no_mapq_shading", action="store_true",
                         help="Disable lightening read fill color for lower MAPQ.")
    parser.add_argument("--mapq_cap", type=int, default=60,
                         help="MAPQ at which a read's fill reaches full opacity under MAPQ shading.")

    parser.add_argument("--no_coverage", action="store_true", help="Hide the coverage track.")
    parser.add_argument("--no_annotate", action="store_true", help="Hide the per-row gap-length annotation (expand layout).")
    parser.add_argument("--fig_width", type=float, default=14.0, help="Figure width in inches.")
    parser.add_argument("--dpi", type=int, default=150, help="Output resolution.")

    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    try:
        chrom, start, end = parse_region(args.region, flank=args.flank)
    except ValueError as exc:
        log.error(str(exc))
        return 1

    if not os.path.isfile(args.bam):
        log.error("Cannot find --bam file: %s", args.bam)
        return 1
    if args.fasta and not os.path.isfile(args.fasta):
        log.error("Cannot find --fasta file: %s", args.fasta)
        return 1

    common_kwargs = dict(
        min_mapq=args.min_mapq,
        include_secondary=args.include_secondary,
        include_supplementary=not args.exclude_supplementary,
        include_duplicates=args.include_duplicates,
        max_rows=args.max_rows,
        show_coverage=not args.no_coverage,
        annotate_gap=not args.no_annotate,
        fig_width=args.fig_width,
        dpi=args.dpi,
        long_gap_threshold=args.long_gap_threshold,
        layout=args.layout,
        sort_by=args.sort_by,
        descending=(args.sort_order == "desc"),
        only_types=args.only,
        min_softclip=args.min_softclip,
        insert_size_sigma=args.insert_size_sigma,
        pair_colors=not args.no_pair_colors,
        shade_by_mapq=not args.no_mapq_shading,
        mapq_cap=args.mapq_cap,
    )

    if args.bam2:
        if not os.path.isfile(args.bam2):
            log.error("Cannot find --bam2 file: %s", args.bam2)
            return 1
        out_path, table = compare_snapshots(
            bam1=args.bam, bam2=args.bam2, chrom=chrom, start=start, end=end,
            fasta=args.fasta, output_dir=args.output_dir, output_name=args.output_name,
            label1=args.label1, label2=args.label2,
            metrics_tsv_1=args.metrics_tsv, metrics_tsv_2=args.metrics_tsv2,
            **common_kwargs,
        )
        log.info("Wrote comparison snapshot: %s", out_path)
        print(table)
        return 0

    snap = BamSnapshot(
        bam=args.bam, chrom=chrom, start=start, end=end, fasta=args.fasta,
        output_dir=args.output_dir, output_name=args.output_name,
        label=args.label1, **common_kwargs,
    )
    summary = snap.snap(metrics_tsv=args.metrics_tsv)
    log.info("Wrote snapshot: %s", snap.output_png)
    print(
        f"{summary.n_reads} reads | gapped: {summary.n_gapped} ({summary.pct_gapped:.1f}%) | "
        f"max gap: {summary.max_gap}bp | split (SA): {summary.n_with_sa} | "
        f"discordant: {summary.n_discordant} ({summary.pct_discordant:.1f}%) | "
        f"soft-clipped: {summary.n_softclipped}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
