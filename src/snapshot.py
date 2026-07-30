"""Top-level orchestration: region -> reads -> rows -> PNG (+ optional TSV).

Two entry points:

- ``BamSnapshot``: one BAM, one image.
- ``compare_snapshots``: two BAMs (e.g. bwa vs minibwa) over the same region,
  rendered as one stacked side-by-side image plus a printable summary table
  answering "which one produced more/longer gapped alignments here".
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional

from src.layout import build_rows, truncate_rows
from src.metrics import RegionSummary, format_summary_table, summarize, write_tsv
from src.read_model import AlignedRead, fetch_reads
from src.reference import ReferenceWindow
from src.render import AlignmentRenderer


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
        label: Optional[str] = None,
        only_types: Optional[List[str]] = None,
        min_softclip: int = 1,
        insert_size_sigma: float = 3.0,
        pair_colors: bool = True,
        shade_by_mapq: bool = True,
        mapq_cap: int = 60,
    ):
        self.bam = bam
        self.chrom = chrom
        self.start = start
        self.end = end
        self.fasta = fasta
        self.output_dir = output_dir
        self.output_name = output_name
        self.layout = layout
        self.sort_by = sort_by
        self.descending = descending
        self.min_mapq = min_mapq
        self.include_secondary = include_secondary
        self.include_supplementary = include_supplementary
        self.include_duplicates = include_duplicates
        self.max_rows = max_rows
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

        os.makedirs(self.output_dir, exist_ok=True)

        name = self.output_name or f"{chrom}_{start}_{end}.png"
        if not name.endswith(".png"):
            name += ".png"
        self.output_png = str(Path(self.output_dir) / name)

        self.reads: List[AlignedRead] = []
        self.summary: Optional[RegionSummary] = None

    def load_reads(self) -> List[AlignedRead]:
        reference = ReferenceWindow(self.fasta, self.chrom, self.start, self.end)
        self.reads = fetch_reads(
            self.bam, self.chrom, self.start, self.end,
            reference=reference,
            min_mapq=self.min_mapq,
            include_secondary=self.include_secondary,
            include_supplementary=self.include_supplementary,
            include_duplicates=self.include_duplicates,
            insert_size_sigma=self.insert_size_sigma,
            only_types=self.only_types,
            min_softclip=self.min_softclip,
        )
        self._reference = reference
        return self.reads

    def snap(self, metrics_tsv: Optional[str] = None) -> RegionSummary:
        reads = self.load_reads() if not self.reads else self.reads
        rows = build_rows(reads, layout=self.layout, sort_by=self.sort_by, descending=self.descending)
        rows, dropped = truncate_rows(rows, self.max_rows)

        renderer = AlignmentRenderer(
            fig_width=self.fig_width, dpi=self.dpi,
            show_coverage=self.show_coverage, annotate_gap=self.annotate_gap,
            pair_colors=self.pair_colors, shade_by_mapq=self.shade_by_mapq, mapq_cap=self.mapq_cap,
        )
        title = (
            f"{self.label} -- {len(reads)} reads, layout={self.layout}, "
            f"sort_by={self.sort_by} ({'desc' if self.descending else 'asc'})"
        )
        renderer.render(
            rows=rows, chrom=self.chrom, window_start=self.start, window_end=self.end,
            reference=self._reference, out_path=self.output_png, title=title,
            layout=self.layout, dropped_reads=dropped, all_reads_for_coverage=reads,
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
) -> str:
    """Renders both BAMs stacked in one PNG, sharing a genomic x-axis, and
    returns a plain-text comparison table (also handy to print/log)."""
    os.makedirs(output_dir, exist_ok=True)
    label1 = label1 or Path(bam1).stem
    label2 = label2 or Path(bam2).stem

    reference = ReferenceWindow(fasta, chrom, start, end)

    panels = []
    summaries = []
    for bam_path, label, tsv_path in ((bam1, label1, metrics_tsv_1), (bam2, label2, metrics_tsv_2)):
        reads = fetch_reads(
            bam_path, chrom, start, end, reference=reference, min_mapq=min_mapq,
            include_secondary=include_secondary, include_supplementary=include_supplementary,
            include_duplicates=include_duplicates, insert_size_sigma=insert_size_sigma,
            only_types=only_types, min_softclip=min_softclip,
        )
        rows = build_rows(reads, layout=layout, sort_by=sort_by, descending=descending)
        rows, dropped = truncate_rows(rows, max_rows)
        summary = summarize(reads, label=label, long_gap_threshold=long_gap_threshold, min_softclip=min_softclip)
        summaries.append(summary)
        if tsv_path:
            write_tsv(reads, tsv_path)
        panels.append({
            "label": f"{label}  (n={len(reads)}, gapped={summary.n_gapped}, max_gap={summary.max_gap}bp)",
            "rows": rows,
            "all_reads_for_coverage": reads,
            "layout": layout,
            "dropped_reads": dropped,
        })

    name = output_name or f"compare_{chrom}_{start}_{end}.png"
    if not name.endswith(".png"):
        name += ".png"
    out_path = str(Path(output_dir) / name)

    renderer = AlignmentRenderer(
        fig_width=fig_width, dpi=dpi, show_coverage=show_coverage, annotate_gap=annotate_gap,
        pair_colors=pair_colors, shade_by_mapq=shade_by_mapq, mapq_cap=mapq_cap,
    )
    renderer.render_multi(
        panels=panels, chrom=chrom, window_start=start, window_end=end,
        reference=reference, out_path=out_path,
        suptitle=f"layout={layout}, sort_by={sort_by} ({'desc' if descending else 'asc'})",
    )

    return out_path, format_summary_table(summaries)
