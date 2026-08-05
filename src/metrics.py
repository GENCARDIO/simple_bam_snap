"""Turns the per-read features into numbers a human (or a script) can use to
actually answer "does aligner A produce more/longer gapped alignments than
aligner B", instead of just eyeballing a picture.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from statistics import mean
from typing import List, Optional

from src.read_model import AlignedRead

TSV_FIELDS = [
    "read_name",
    "chrom",
    "start",
    "end",
    "strand",
    "mapq",
    "cigar_gap_len",
    "sa_gap_len",
    "gap_length",
    "sa_count",
    "has_cross_chrom_sa",
    "soft_clip_left",
    "soft_clip_right",
    "mismatch_count",
    "insert_size",
    "pair_orientation",
    "pair_category",
    "haplotype",
    "phase_set",
    "mate_chrom",
    "is_secondary",
    "is_supplementary",
    "is_duplicate",
]


def read_to_row(read: AlignedRead) -> dict:
    return {
        "read_name": read.query_name,
        "chrom": read.reference_name,
        "start": read.ref_start,
        "end": read.ref_end,
        "strand": read.strand,
        "mapq": read.mapq,
        "cigar_gap_len": read.cigar_gap_len,
        "sa_gap_len": read.sa_gap_len,
        "gap_length": read.gap_length,
        "sa_count": read.sa_count,
        "has_cross_chrom_sa": read.has_cross_chrom_sa,
        "soft_clip_left": read.soft_clip_left,
        "soft_clip_right": read.soft_clip_right,
        "mismatch_count": read.mismatch_count,
        "insert_size": read.insert_size,
        "pair_orientation": read.pair_orientation,
        "pair_category": read.pair_category,
        "haplotype": read.haplotype,
        "phase_set": read.phase_set,
        "mate_chrom": read.mate_chrom,
        "is_secondary": read.is_secondary,
        "is_supplementary": read.is_supplementary,
        "is_duplicate": read.is_duplicate,
    }


def write_tsv(reads: List[AlignedRead], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        for read in reads:
            writer.writerow(read_to_row(read))


@dataclass
class RegionSummary:
    label: str
    n_reads: int
    n_gapped: int
    n_long_gap: int
    long_gap_threshold: int
    max_gap: int
    mean_gap_of_gapped: float
    total_gap_bp: int
    n_with_sa: int
    n_cross_chrom_sa: int
    n_discordant: int
    n_interchrom: int
    n_softclipped: int
    mean_mapq: float

    @property
    def pct_gapped(self) -> float:
        return 100.0 * self.n_gapped / self.n_reads if self.n_reads else 0.0

    @property
    def pct_long_gap(self) -> float:
        return 100.0 * self.n_long_gap / self.n_reads if self.n_reads else 0.0

    @property
    def pct_discordant(self) -> float:
        return 100.0 * self.n_discordant / self.n_reads if self.n_reads else 0.0


def summarize(
    reads: List[AlignedRead], label: str = "", long_gap_threshold: int = 10, min_softclip: int = 1
) -> RegionSummary:
    n_reads = len(reads)
    gapped = []
    n_long_gapped = 0
    n_with_sa = 0
    n_cross_chrom = 0
    n_discordant = 0
    n_interchrom = 0
    n_softclipped = 0
    for read in reads:
        if read.gap_length > 0:
            gapped.append(read)
        if read.gap_length >= long_gap_threshold:
            n_long_gapped += 1
        if read.sa_count > 0:
            n_with_sa += 1
        if read.has_cross_chrom_sa:
            n_cross_chrom += 1
        if read.is_discordant:
            n_discordant += 1
        if read.pair_category == "interchrom":
            n_interchrom += 1
        if read.soft_clip_total >= min_softclip:
            n_softclipped += 1
    return RegionSummary(
        label=label,
        n_reads=n_reads,
        n_gapped=len(gapped),
        n_long_gap=n_long_gapped,
        long_gap_threshold=long_gap_threshold,
        max_gap=max((r.gap_length for r in reads), default=0),
        mean_gap_of_gapped=mean((r.gap_length for r in gapped)) if gapped else 0.0,
        total_gap_bp=sum(r.gap_length for r in reads),
        n_with_sa=n_with_sa,
        n_cross_chrom_sa=n_cross_chrom,
        n_discordant=n_discordant,
        n_interchrom=n_interchrom,
        n_softclipped=n_softclipped,
        mean_mapq=mean((r.mapq for r in reads)) if reads else 0.0,
    )


def format_summary_table(summaries: List[RegionSummary]) -> str:
    """Render one or more RegionSummary objects as a plain-text side-by-side
    table for terminal output - the quick "who has more gapped alignments"
    answer."""
    rows = [
        ("reads", "{s.n_reads}"),
        ("gapped reads", "{s.n_gapped} ({s.pct_gapped:.1f}%)"),
        (">= {thr}bp gap".format(thr=summaries[0].long_gap_threshold) if summaries else "long gap",
         "{s.n_long_gap} ({s.pct_long_gap:.1f}%)"),
        ("max gap (bp)", "{s.max_gap}"),
        ("mean gap of gapped (bp)", "{s.mean_gap_of_gapped:.1f}"),
        ("total gap bp", "{s.total_gap_bp}"),
        ("reads with SA (split)", "{s.n_with_sa}"),
        ("cross-chrom SA", "{s.n_cross_chrom_sa}"),
        ("discordant pairs", "{s.n_discordant} ({s.pct_discordant:.1f}%)"),
        ("inter-chromosomal pairs", "{s.n_interchrom}"),
        ("soft-clipped reads", "{s.n_softclipped}"),
        ("mean MAPQ", "{s.mean_mapq:.1f}"),
    ]

    labels = []
    label_widths = [len("metric")]
    for index, summary in enumerate(summaries):
        label = summary.label or f"bam{index + 1}"
        labels.append(label)
        label_widths.append(len(label))
    col_width = max(label_widths) + 2
    metric_width = max(len(row[0]) for row in rows) + 2

    lines = []
    header = "metric".ljust(metric_width) + "".join(l.ljust(col_width) for l in labels)
    lines.append(header)
    lines.append("-" * len(header))
    for name, fmt in rows:
        cells = []
        for summary in summaries:
            cells.append(fmt.format(s=summary))
        lines.append(name.ljust(metric_width) + "".join(c.ljust(col_width) for c in cells))
    return "\n".join(lines)
