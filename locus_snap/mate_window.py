"""Select a second genomic window from pair or split-alignment evidence."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from statistics import fmean
from typing import Iterable, List, Mapping, Optional, Tuple

from locus_snap.read_model import AlignedRead


MATE_WINDOW_SOURCES = ("discordant", "split", "softclip")


@dataclass(frozen=True)
class MateWindow:
    chrom: str
    start: int
    end: int
    source: str
    candidate_count: int

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


def mate_candidates(
    reads: Iterable[AlignedRead], source: str, min_softclip: int = 1
) -> List[Tuple[str, int]]:
    """Return ``(chromosome, 0-based coordinate)`` candidates for a source."""
    if source not in MATE_WINDOW_SOURCES:
        raise ValueError(
            f"Unknown mate-window source '{source}'; choose from {', '.join(MATE_WINDOW_SOURCES)}."
        )

    candidates: List[Tuple[str, int]] = []
    for read in reads:
        if source == "split":
            for entry in read.sa_entries:
                candidates.append((entry.rname, (entry.start + entry.end) // 2))
            continue

        eligible = read.is_discordant if source == "discordant" else read.soft_clip_total >= min_softclip
        if eligible and read.mate_chrom is not None and read.mate_start is not None and read.mate_start >= 0:
            candidates.append((read.mate_chrom, read.mate_start))
    return candidates


def supporting_query_names(
    reads: Iterable[AlignedRead], source: str, chrom: str, min_softclip: int = 1
) -> set[str]:
    """Names of primary-window reads that contributed evidence on ``chrom``."""
    if source not in MATE_WINDOW_SOURCES:
        raise ValueError(
            f"Unknown mate-window source '{source}'; choose from {', '.join(MATE_WINDOW_SOURCES)}."
        )

    names = set()
    for read in reads:
        if source == "split":
            matches_chrom = False
            for entry in read.sa_entries:
                if entry.rname == chrom:
                    matches_chrom = True
                    break
            if matches_chrom:
                names.add(read.query_name)
            continue
        eligible = read.is_discordant if source == "discordant" else read.soft_clip_total >= min_softclip
        if eligible and read.mate_chrom == chrom:
            names.add(read.query_name)
    return names


def choose_mate_window(
    reads: Iterable[AlignedRead],
    source: str,
    window_size: int,
    contig_lengths: Optional[Mapping[str, int]] = None,
    min_softclip: int = 1,
) -> MateWindow:
    """Choose the most-supported chromosome and center it on the mean locus.

    Candidates on different chromosomes cannot be meaningfully averaged. We
    therefore select the chromosome with the most evidence (lexical order
    breaks ties), then take the rounded arithmetic mean of its coordinates.
    """
    if window_size <= 0:
        raise ValueError("Mate-window size must be greater than zero.")

    candidates = mate_candidates(reads, source=source, min_softclip=min_softclip)
    if not candidates:
        description = {
            "discordant": "mapped mates of discordant reads",
            "split": "SA split-alignment entries",
            "softclip": "mapped mates of soft-clipped reads",
        }[source]
        raise ValueError(f"Cannot create mate view: the primary window has no {description}.")

    counts = Counter()
    for candidate in candidates:
        counts[candidate[0]] += 1
    selected_chrom = None
    for chrom in counts:
        if selected_chrom is None:
            selected_chrom = chrom
            continue
        if counts[chrom] > counts[selected_chrom]:
            selected_chrom = chrom
        elif counts[chrom] == counts[selected_chrom] and chrom < selected_chrom:
            selected_chrom = chrom
    positions = []
    for chrom, position in candidates:
        if chrom == selected_chrom:
            positions.append(position)
    center = int(round(fmean(positions)))

    start = max(0, center - window_size // 2)
    end = start + window_size
    if contig_lengths and selected_chrom in contig_lengths:
        contig_length = contig_lengths[selected_chrom]
        end = min(end, contig_length)
        start = max(0, end - window_size)

    return MateWindow(
        chrom=selected_chrom,
        start=start,
        end=end,
        source=source,
        candidate_count=len(positions),
    )
