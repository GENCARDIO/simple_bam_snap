"""Top-level orchestration: region -> reads -> rows -> image (+ optional TSV).

Two entry points:

- ``BamSnapshot``: one BAM, one image.
- ``compare_snapshots``: two BAMs (e.g. bwa vs minibwa) over the same region,
  rendered as one stacked image plus a printable summary table
  answering "which one produced more/longer gapped alignments here".
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import pysam

from src.annotations import AnnotationSource
from src.cytobands import Cytoband, bands_for_chrom, resolve_cytobands
from src.downsample import DEFAULT_MAX_ALIGNMENT_DEPTH, downsample_reads
from src.layout import build_rows, infer_reference_base, truncate_rows
from src.mate_window import MateWindow, choose_mate_window, supporting_query_names
from src.metrics import RegionSummary, format_summary_table, summarize, write_tsv
from src.read_model import AlignedRead, fetch_reads, matches_only
from src.reference import ReferenceWindow
from src.render import (
    AlignmentRenderer,
    DEFAULT_COVERAGE_VAF_THRESHOLD,
    DEFAULT_MAX_REFERENCE_SPAN,
)

OUTPUT_FORMATS = ("png", "svg", "svgz", "pdf", "jpg", "jpeg", "tif", "tiff", "webp")


def resolve_output_path(
    output_dir: str,
    output_name: Optional[str],
    default_stem: str,
    output_format: Optional[str] = None,
) -> str:
    """Resolve an output filename, inferring its format from a known suffix."""
    selected_format = output_format.lower().lstrip(".") if output_format else None
    if selected_format is not None and selected_format not in OUTPUT_FORMATS:
        choices = ", ".join(OUTPUT_FORMATS)
        raise ValueError(f"Unsupported output format {output_format!r}; choose from {choices}.")

    name = output_name or default_stem
    path = Path(name)
    suffix = path.suffix.lower().lstrip(".")
    if selected_format is None:
        if suffix:
            if suffix not in OUTPUT_FORMATS:
                choices = ", ".join(OUTPUT_FORMATS)
                raise ValueError(
                    f"Unsupported output filename extension '.{suffix}'; choose from {choices}."
                )
            selected_format = suffix
        else:
            selected_format = "png"

    if suffix in OUTPUT_FORMATS:
        path = path.with_suffix(f".{selected_format}")
    elif not suffix:
        path = path.with_suffix(f".{selected_format}")
    else:
        path = Path(f"{path}.{selected_format}")
    return str(Path(output_dir) / path)


class BamSnapshot:
    def __init__(
        self,
        bam: str,
        chrom: str,
        start: int,
        end: int,
        fasta: Optional[str] = None,
        output_dir: str = ".",
        output_name: Optional[str] = None,
        layout: str = "pack",
        sort_by: str = "gap_length",
        sort_base_position: Optional[int] = None,
        descending: bool = True,
        min_mapq: int = 0,
        include_secondary: bool = False,
        include_supplementary: bool = True,
        include_duplicates: bool = False,
        max_rows: Optional[int] = None,
        show_alignments: bool = True,
        show_coverage: bool = True,
        annotate_gap: bool = True,
        fig_width: float = 14.0,
        dpi: int = 150,
        long_gap_threshold: int = 10,
        label: Optional[str] = None,
        only_types: Optional[List[str]] = None,
        min_softclip: int = 1,
        insert_size_sigma: float = 3.0,
        pair_colors: bool = True,
        shade_by_mapq: bool = True,
        mapq_cap: int = 60,
        alignment_colors: Optional[Dict[str, Optional[str]]] = None,
        visual_config: Optional[Dict[str, Any]] = None,
        mate_view: bool = False,
        mate_window_source: str = "discordant",
        mate_window_size: Optional[int] = None,
        display_mode: str = "expand",
        max_alignment_depth: int = DEFAULT_MAX_ALIGNMENT_DEPTH,
        annotation_sources: Optional[List[AnnotationSource]] = None,
        show_ideogram: bool = True,
        show_center_guide: bool = False,
        show_sashimi: bool = False,
        min_junction_reads: int = 1,
        sashimi_strand: str = "combined",
        genome: str = "auto",
        cytoband_file: Optional[str] = None,
        max_reference_span: int = DEFAULT_MAX_REFERENCE_SPAN,
        view_as_pairs: bool = False,
        coverage_vaf_threshold: float = DEFAULT_COVERAGE_VAF_THRESHOLD,
        min_baseq: int = 0,
        min_variant_mapq: int = 0,
        show_variant_counts: bool = False,
        show_indel_lengths: bool = False,
        haplotype_view: str = "none",
        haplotype_filter: Optional[List[str]] = None,
        haplotype_tag: str = "HP",
        phase_set_tag: str = "PS",
        output_format: Optional[str] = None,
    ):
        self.bam = bam
        self.chrom = chrom
        self.start = start
        self.end = end
        self.fasta = fasta
        self.output_dir = output_dir
        self.output_name = output_name
        self.output_format = output_format
        self.layout = layout
        self.sort_by = sort_by
        self.sort_base_position = sort_base_position
        self.descending = descending
        self.min_mapq = min_mapq
        self.include_secondary = include_secondary
        self.include_supplementary = include_supplementary
        self.include_duplicates = include_duplicates
        self.max_rows = max_rows
        self.show_alignments = show_alignments
        self.show_coverage = show_coverage
        self.annotate_gap = annotate_gap
        self.fig_width = fig_width
        self.dpi = dpi
        self.long_gap_threshold = long_gap_threshold
        self.label = label or Path(bam).stem
        self.only_types = only_types
        self.min_softclip = min_softclip
        self.insert_size_sigma = insert_size_sigma
        self.pair_colors = pair_colors
        self.shade_by_mapq = shade_by_mapq
        self.mapq_cap = mapq_cap
        self.alignment_colors = alignment_colors
        self.visual_config = visual_config
        self.mate_view = mate_view
        self.mate_window_source = mate_window_source
        self.mate_window_size = mate_window_size
        self.display_mode = display_mode
        self.max_alignment_depth = max_alignment_depth
        self.annotation_sources = list(annotation_sources or [])
        self.show_ideogram = show_ideogram
        self.show_center_guide = show_center_guide
        self.show_sashimi = show_sashimi
        self.min_junction_reads = min_junction_reads
        self.sashimi_strand = sashimi_strand
        self.genome = genome
        self.cytoband_file = cytoband_file
        self.max_reference_span = max_reference_span
        self.view_as_pairs = view_as_pairs
        self.coverage_vaf_threshold = coverage_vaf_threshold
        self.min_baseq = min_baseq
        self.min_variant_mapq = min_variant_mapq
        self.show_variant_counts = show_variant_counts
        self.show_indel_lengths = show_indel_lengths
        self.haplotype_view = haplotype_view
        self.haplotype_filter = list(haplotype_filter or [])
        self.haplotype_tag = haplotype_tag
        self.phase_set_tag = phase_set_tag

        os.makedirs(self.output_dir, exist_ok=True)

        default_prefix = "mate_" if mate_view else ""
        self.output_path = resolve_output_path(
            self.output_dir,
            self.output_name,
            f"{default_prefix}{chrom}_{start}_{end}",
            self.output_format,
        )
        # Kept as a compatibility alias for callers written before multi-format export.
        self.output_png = self.output_path

        self.reads: List[AlignedRead] = []
        self.source_reads: List[AlignedRead] = []
        self.summary: Optional[RegionSummary] = None
        self.mate_window: Optional[MateWindow] = None
        self.downsampled_reads = 0
        self.contig_lengths: Dict[str, int] = {}
        self.cytobands: Dict[str, List[Cytoband]] = {}
        self.cytoband_label: Optional[str] = None
        self.reads_loaded = False

    def load_reads(self) -> List[AlignedRead]:
        with pysam.AlignmentFile(self.bam, "rb") as bam_file:
            self.contig_lengths = dict(zip(bam_file.references, bam_file.lengths))
        if self.show_ideogram:
            self.cytobands, self.cytoband_label = resolve_cytobands(
                self.contig_lengths, genome=self.genome, custom_path=self.cytoband_file
            )
        reference = ReferenceWindow(self.fasta, self.chrom, self.start, self.end)
        self.source_reads = fetch_reads(
            self.bam, self.chrom, self.start, self.end,
            reference=reference,
            min_mapq=self.min_mapq,
            include_secondary=self.include_secondary,
            include_supplementary=self.include_supplementary,
            include_duplicates=self.include_duplicates,
            insert_size_sigma=self.insert_size_sigma,
            only_types=None,
            min_softclip=self.min_softclip,
            haplotype_tag=self.haplotype_tag,
            phase_set_tag=self.phase_set_tag,
            haplotype_filter=self.haplotype_filter,
        )
        self.reads = []
        for read in self.source_reads:
            if matches_only(read, self.only_types, self.min_softclip):
                self.reads.append(read)
        self.reference_window = reference
        self.reads_loaded = True
        return self.reads

    def snap(self, metrics_tsv: Optional[str] = None) -> RegionSummary:
        reads = self.load_reads() if not self.reads_loaded else self.reads
        base_position = self.sort_base_position
        if self.sort_by == "base" and base_position is None:
            base_position = self.start + (self.end - self.start) // 2
        reference_base = (
            self.reference_window.base_at(base_position)
            if base_position is not None else None
        )
        if self.sort_by == "base" and reference_base not in ("A", "C", "G", "T"):
            reference_base = infer_reference_base(reads, base_position)
        base_priority_names = set()
        if self.sort_by == "base":
            for read in reads:
                observed = read.base_at(base_position)
                if observed in ("A", "C", "G", "T") and (
                    reference_base not in ("A", "C", "G", "T")
                    or observed != reference_base
                ):
                    base_priority_names.add(read.query_name)
        genomic_tracks = []
        for source in self.annotation_sources:
            genomic_tracks.append(source.fetch(self.chrom, self.start, self.end))
        display_reads, self.downsampled_reads = downsample_reads(
            reads, max_depth=self.max_alignment_depth,
            priority_names=base_priority_names,
            preserve_pairs=self.view_as_pairs,
        )
        rows = build_rows(
            display_reads, layout=self.layout, sort_by=self.sort_by,
            descending=self.descending, display_mode=self.display_mode,
            view_as_pairs=self.view_as_pairs,
            haplotype_view=self.haplotype_view,
            base_position=base_position, reference_base=reference_base,
        )
        rows, dropped = truncate_rows(rows, self.max_rows)

        renderer = AlignmentRenderer(
            fig_width=self.fig_width, dpi=self.dpi,
            show_alignments=self.show_alignments,
            show_coverage=self.show_coverage, annotate_gap=self.annotate_gap,
            pair_colors=self.pair_colors, shade_by_mapq=self.shade_by_mapq, mapq_cap=self.mapq_cap,
            alignment_colors=self.alignment_colors,
            visual_config=self.visual_config,
            display_mode=self.display_mode,
            show_ideogram=self.show_ideogram,
            max_reference_span=self.max_reference_span,
            view_as_pairs=self.view_as_pairs,
            coverage_vaf_threshold=self.coverage_vaf_threshold,
            min_baseq=self.min_baseq,
            min_variant_mapq=self.min_variant_mapq,
            show_variant_counts=self.show_variant_counts,
            show_indel_lengths=self.show_indel_lengths,
            haplotype_view=self.haplotype_view,
            sort_base_position=base_position if self.sort_by == "base" else None,
            sort_reference_base=reference_base,
            show_center_guide=self.show_center_guide,
            show_sashimi=self.show_sashimi,
            min_junction_reads=self.min_junction_reads,
            sashimi_strand=self.sashimi_strand,
        )
        sort_label = self.sort_by
        if self.sort_by == "base":
            sort_label += f"@{base_position + 1:,}"
            if reference_base:
                sort_label += f" ref={reference_base}"
        title = self.label
        if self.show_alignments:
            title = (
                f"{self.label} -- {len(reads)} reads, display={self.display_mode}, layout={self.layout}, "
                f"view={'pairs' if self.view_as_pairs else 'alignments'}, "
                f"haplotypes={self.haplotype_view}, "
                f"sort_by={sort_label} ({'desc' if self.descending else 'asc'})"
            )
        if self.mate_view:
            self.mate_window = choose_mate_window(
                self.source_reads,
                source=self.mate_window_source,
                window_size=self.mate_window_size or (self.end - self.start),
                contig_lengths=self.contig_lengths,
                min_softclip=self.min_softclip,
            )
            mate = self.mate_window
            mate_reference = ReferenceWindow(self.fasta, mate.chrom, mate.start, mate.end)
            mate_source_reads = fetch_reads(
                self.bam, mate.chrom, mate.start, mate.end,
                reference=mate_reference,
                min_mapq=self.min_mapq,
                include_secondary=self.include_secondary,
                include_supplementary=self.include_supplementary,
                include_duplicates=self.include_duplicates,
                insert_size_sigma=self.insert_size_sigma,
                only_types=None,
                min_softclip=self.min_softclip,
                haplotype_tag=self.haplotype_tag,
                phase_set_tag=self.phase_set_tag,
                haplotype_filter=self.haplotype_filter,
            )
            support_names = supporting_query_names(
                self.source_reads, self.mate_window_source, mate.chrom, self.min_softclip
            )
            mate_reads = []
            for read in mate_source_reads:
                if (
                    matches_only(read, self.only_types, self.min_softclip)
                    or read.query_name in support_names
                ):
                    mate_reads.append(read)
            mate_base_position = base_position
            if (
                self.sort_by == "base"
                and (mate.chrom != self.chrom or not mate.start <= base_position < mate.end)
            ):
                mate_base_position = mate.start + (mate.end - mate.start) // 2
            mate_reference_base = (
                mate_reference.base_at(mate_base_position)
                if mate_base_position is not None else None
            )
            if (
                self.sort_by == "base"
                and mate_reference_base not in ("A", "C", "G", "T")
            ):
                mate_reference_base = infer_reference_base(
                    mate_reads, mate_base_position
                )
            mate_priority_names = set(support_names)
            if self.sort_by == "base":
                for read in mate_reads:
                    observed = read.base_at(mate_base_position)
                    if observed in ("A", "C", "G", "T") and (
                        mate_reference_base not in ("A", "C", "G", "T")
                        or observed != mate_reference_base
                    ):
                        mate_priority_names.add(read.query_name)
            mate_display_reads, mate_downsampled = downsample_reads(
                mate_reads, max_depth=self.max_alignment_depth,
                priority_names=mate_priority_names,
                preserve_pairs=self.view_as_pairs,
            )
            mate_rows = build_rows(
                mate_display_reads, layout=self.layout, sort_by=self.sort_by,
                descending=self.descending, display_mode=self.display_mode,
                view_as_pairs=self.view_as_pairs,
                haplotype_view=self.haplotype_view,
                base_position=mate_base_position,
                reference_base=mate_reference_base,
            )
            mate_rows, mate_dropped = truncate_rows(mate_rows, self.max_rows)
            mate_genomic_tracks = []
            for source in self.annotation_sources:
                mate_genomic_tracks.append(
                    source.fetch(mate.chrom, mate.start, mate.end)
                )
            renderer.render_loci(
                panels=[
                    {
                        "label": f"Primary · {self.label} (n={len(reads)})",
                        "chrom": self.chrom, "start": self.start, "end": self.end,
                        "reference": self.reference_window, "rows": rows,
                        "all_reads_for_coverage": reads, "layout": self.layout,
                        "dropped_reads": dropped,
                        "downsampled_reads": self.downsampled_reads,
                        "sort_base_position": base_position,
                        "sort_reference_base": reference_base,
                        "genomic_tracks": genomic_tracks,
                        "contig_length": self.contig_lengths.get(self.chrom),
                        "cytobands": bands_for_chrom(self.cytobands, self.chrom),
                    },
                    {
                        "label": (
                            f"Mate · {self.mate_window_source}; "
                            f"{mate.candidate_count} candidate(s) (n={len(mate_reads)})"
                        ),
                        "chrom": mate.chrom, "start": mate.start, "end": mate.end,
                        "reference": mate_reference, "rows": mate_rows,
                        "all_reads_for_coverage": mate_reads, "layout": self.layout,
                        "dropped_reads": mate_dropped,
                        "downsampled_reads": mate_downsampled,
                        "sort_base_position": mate_base_position,
                        "sort_reference_base": mate_reference_base,
                        "genomic_tracks": mate_genomic_tracks,
                        "contig_length": self.contig_lengths.get(mate.chrom),
                        "cytobands": bands_for_chrom(self.cytobands, mate.chrom),
                    },
                ],
                out_path=self.output_path,
                suptitle=(
                    f"{self.label} mate view -- display={self.display_mode}, layout={self.layout}, "
                    f"view={'pairs' if self.view_as_pairs else 'alignments'}, "
                    f"haplotypes={self.haplotype_view}, "
                    f"sort_by={sort_label} ({'desc' if self.descending else 'asc'})"
                ),
            )
        else:
            renderer.render(
                rows=rows, chrom=self.chrom, window_start=self.start, window_end=self.end,
                reference=self.reference_window, out_path=self.output_path, title=title,
                layout=self.layout, dropped_reads=dropped, all_reads_for_coverage=reads,
                downsampled_reads=self.downsampled_reads,
                genomic_tracks=genomic_tracks,
                contig_length=self.contig_lengths.get(self.chrom),
                cytobands=bands_for_chrom(self.cytobands, self.chrom),
            )

        self.summary = summarize(
            reads, label=self.label, long_gap_threshold=self.long_gap_threshold, min_softclip=self.min_softclip
        )
        if metrics_tsv:
            write_tsv(reads, metrics_tsv)
        return self.summary


def compare_snapshots(
    bam1: str,
    bam2: str,
    chrom: str,
    start: int,
    end: int,
    fasta: Optional[str] = None,
    output_dir: str = ".",
    output_name: Optional[str] = None,
    label1: Optional[str] = None,
    label2: Optional[str] = None,
    layout: str = "expand",
    sort_by: str = "gap_length",
    sort_base_position: Optional[int] = None,
    descending: bool = True,
    min_mapq: int = 0,
    include_secondary: bool = False,
    include_supplementary: bool = True,
    include_duplicates: bool = False,
    max_rows: Optional[int] = None,
    show_coverage: bool = True,
    annotate_gap: bool = True,
    fig_width: float = 14.0,
    dpi: int = 150,
    long_gap_threshold: int = 10,
    metrics_tsv_1: Optional[str] = None,
    metrics_tsv_2: Optional[str] = None,
    only_types: Optional[List[str]] = None,
    min_softclip: int = 1,
    insert_size_sigma: float = 3.0,
    pair_colors: bool = True,
    shade_by_mapq: bool = True,
    mapq_cap: int = 60,
    alignment_colors: Optional[Dict[str, Optional[str]]] = None,
    visual_config: Optional[Dict[str, Any]] = None,
    display_mode: str = "expand",
    max_alignment_depth: int = DEFAULT_MAX_ALIGNMENT_DEPTH,
    annotation_sources: Optional[List[AnnotationSource]] = None,
    show_ideogram: bool = True,
    show_center_guide: bool = False,
    show_sashimi: bool = False,
    min_junction_reads: int = 1,
    sashimi_strand: str = "combined",
    genome: str = "auto",
    cytoband_file: Optional[str] = None,
    max_reference_span: int = DEFAULT_MAX_REFERENCE_SPAN,
    view_as_pairs: bool = False,
    coverage_vaf_threshold: float = DEFAULT_COVERAGE_VAF_THRESHOLD,
    min_baseq: int = 0,
    min_variant_mapq: int = 0,
    show_variant_counts: bool = False,
    show_indel_lengths: bool = False,
    haplotype_view: str = "none",
    haplotype_filter: Optional[List[str]] = None,
    haplotype_tag: str = "HP",
    phase_set_tag: str = "PS",
    output_format: Optional[str] = None,
    additional_bams: Optional[List[str]] = None,
    additional_labels: Optional[List[str]] = None,
    companion_vcfs: Optional[List[Optional[str]]] = None,
) -> tuple[str, str]:
    """Render two or more BAMs as sample panels sharing one genomic x-axis."""
    os.makedirs(output_dir, exist_ok=True)
    bam_paths = [bam1, bam2]
    bam_paths.extend(additional_bams or [])
    labels = [label1 or Path(bam1).stem, label2 or Path(bam2).stem]
    extra_labels = list(additional_labels or [])
    if len(extra_labels) > len(bam_paths) - 2:
        raise ValueError("More additional labels were supplied than additional BAMs.")
    for index, bam_path in enumerate(bam_paths[2:]):
        label = extra_labels[index] if index < len(extra_labels) else None
        labels.append(label or Path(bam_path).stem)
    tsv_paths = [metrics_tsv_1, metrics_tsv_2]
    while len(tsv_paths) < len(bam_paths):
        tsv_paths.append(None)
    companions = list(companion_vcfs or [])
    if companions and len(companions) != len(bam_paths):
        raise ValueError("Exactly one VCF companion is required per BAM panel.")
    if not companions:
        companions = [None] * len(bam_paths)

    reference = ReferenceWindow(fasta, chrom, start, end)
    base_position = sort_base_position
    if sort_by == "base" and base_position is None:
        base_position = start + (end - start) // 2
    reference_base = reference.base_at(base_position) if base_position is not None else None
    with pysam.AlignmentFile(bam1, "rb") as bam_file:
        contig_lengths = dict(zip(bam_file.references, bam_file.lengths))
    cytobands = (
        resolve_cytobands(
            contig_lengths, genome=genome, custom_path=cytoband_file
        )[0]
        if show_ideogram else {}
    )
    genomic_tracks = []
    for source in annotation_sources or []:
        genomic_tracks.append(source.fetch(chrom, start, end))

    panels = []
    summaries = []
    for panel_index, bam_path in enumerate(bam_paths):
        label = labels[panel_index]
        tsv_path = tsv_paths[panel_index]
        reads = fetch_reads(
            bam_path, chrom, start, end, reference=reference, min_mapq=min_mapq,
            include_secondary=include_secondary, include_supplementary=include_supplementary,
            include_duplicates=include_duplicates, insert_size_sigma=insert_size_sigma,
            only_types=only_types, min_softclip=min_softclip,
            haplotype_tag=haplotype_tag, phase_set_tag=phase_set_tag,
            haplotype_filter=haplotype_filter,
        )
        if sort_by == "base" and reference_base not in ("A", "C", "G", "T"):
            reference_base = infer_reference_base(reads, base_position)
        base_priority_names = set()
        if sort_by == "base":
            for read in reads:
                observed = read.base_at(base_position)
                if observed in ("A", "C", "G", "T") and (
                    reference_base not in ("A", "C", "G", "T")
                    or observed != reference_base
                ):
                    base_priority_names.add(read.query_name)
        display_reads, downsampled = downsample_reads(
            reads, max_depth=max_alignment_depth, preserve_pairs=view_as_pairs,
            priority_names=base_priority_names,
        )
        rows = build_rows(
            display_reads, layout=layout, sort_by=sort_by,
            descending=descending, display_mode=display_mode,
            view_as_pairs=view_as_pairs,
            haplotype_view=haplotype_view,
            base_position=base_position, reference_base=reference_base,
        )
        rows, dropped = truncate_rows(rows, max_rows)
        summary = summarize(reads, label=label, long_gap_threshold=long_gap_threshold, min_softclip=min_softclip)
        summaries.append(summary)
        if tsv_path:
            write_tsv(reads, tsv_path)
        companion_tracks = []
        companion_path = companions[panel_index]
        if companion_path:
            theme_track_colors = None
            if visual_config:
                theme_track_colors = visual_config.get("track_colors")
            companion_source = AnnotationSource(
                companion_path, label=f"{label} variants", kind="vcf",
                display_mode="collapse", track_colors=theme_track_colors,
            )
            companion_tracks.append(companion_source.fetch(chrom, start, end))
        panels.append({
            "label": f"{label}  (n={len(reads)}, gapped={summary.n_gapped}, max_gap={summary.max_gap}bp)",
            "rows": rows,
            "all_reads_for_coverage": reads,
            "layout": layout,
            "dropped_reads": dropped,
            "downsampled_reads": downsampled,
            "companion_tracks": companion_tracks,
        })

    out_path = resolve_output_path(
        output_dir, output_name, f"compare_{chrom}_{start}_{end}", output_format
    )

    renderer = AlignmentRenderer(
        fig_width=fig_width, dpi=dpi, show_coverage=show_coverage, annotate_gap=annotate_gap,
        pair_colors=pair_colors, shade_by_mapq=shade_by_mapq, mapq_cap=mapq_cap,
        alignment_colors=alignment_colors,
        visual_config=visual_config,
        display_mode=display_mode,
        show_ideogram=show_ideogram,
        max_reference_span=max_reference_span,
        view_as_pairs=view_as_pairs,
        coverage_vaf_threshold=coverage_vaf_threshold,
        min_baseq=min_baseq,
        min_variant_mapq=min_variant_mapq,
        show_variant_counts=show_variant_counts,
        show_indel_lengths=show_indel_lengths,
        haplotype_view=haplotype_view,
        sort_base_position=base_position if sort_by == "base" else None,
        sort_reference_base=reference_base,
        show_center_guide=show_center_guide,
        show_sashimi=show_sashimi,
        min_junction_reads=min_junction_reads,
        sashimi_strand=sashimi_strand,
    )
    renderer.render_multi(
        panels=panels, chrom=chrom, window_start=start, window_end=end,
        reference=reference, out_path=out_path,
        suptitle=(
            f"display={display_mode}, layout={layout}, "
            f"view={'pairs' if view_as_pairs else 'alignments'}, "
            f"haplotypes={haplotype_view}, "
            f"sort_by={sort_by}"
            f"{'@' + format(base_position + 1, ',') if sort_by == 'base' else ''} "
            f"({'desc' if descending else 'asc'})"
        ),
        genomic_tracks=genomic_tracks,
        contig_length=contig_lengths.get(chrom),
        cytobands=bands_for_chrom(cytobands, chrom),
    )

    return out_path, format_summary_table(summaries)
