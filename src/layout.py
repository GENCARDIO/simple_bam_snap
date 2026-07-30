"""Turns a flat list of AlignedRead objects into rows to draw.

Two layout strategies are supported:

- ``pack``: the classic IGV "collapsed" behaviour - reads are greedily
  packed left-to-right so non-overlapping reads share a row. Compact, good
  for an overview of a region.
- ``expand``: one read per row, rows ordered by a sort key. This is the
  layout that actually answers "does this BAM have more/longer gapped
  alignments than that one": every read is its own line, ranked by
  ``gap_length`` (or whichever key was requested) so the most heavily
  gapped alignments float straight to the top.

Both accept the same sort-key vocabulary so ``pack`` mode can still express
a priority order (which reads win the earliest rows) even though multiple
reads may end up sharing a row.
"""
from __future__ import annotations

from typing import Callable, List, Optional

from src.read_model import AlignedRead

SortKeyFunc = Callable[[AlignedRead], object]

SORT_KEYS: dict[str, SortKeyFunc] = {
    "gap_length": lambda r: r.gap_length,
    "cigar_gap": lambda r: r.cigar_gap_len,
    "sa_gap": lambda r: r.sa_gap_len,
    "sa_count": lambda r: r.sa_count,
    "soft_clip": lambda r: r.soft_clip_total,
    "mismatch": lambda r: r.mismatch_count,
    "mapq": lambda r: r.mapq,
    "insert_size": lambda r: r.insert_size,
    "start": lambda r: r.ref_start,
    "strand": lambda r: r.strand,
    "read_name": lambda r: r.query_name,
    "none": lambda r: 0,
}


def resolve_sort_key(name: str) -> SortKeyFunc:
    try:
        return SORT_KEYS[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown sort key '{name}'. Choose from: {', '.join(sorted(SORT_KEYS))}"
        ) from exc


def _ordered(reads: List[AlignedRead], sort_by: str, descending: bool) -> List[AlignedRead]:
    key_func = resolve_sort_key(sort_by)
    return sorted(
        reads,
        key=lambda r: (key_func(r), -r.ref_start if descending else r.ref_start),
        reverse=descending,
    )


def pack_rows(
    reads: List[AlignedRead],
    sort_by: str = "start",
    descending: bool = False,
    padding: int = 2,
) -> List[List[AlignedRead]]:
    """Greedy interval packing: process reads in priority order, place each
    in the first row whose rightmost read ends (plus padding) before this
    read starts, else open a new row."""
    ordered = _ordered(reads, sort_by, descending)

    rows: List[List[AlignedRead]] = []
    row_ends: List[int] = []
    for read in ordered:
        placed = False
        for i, row_end in enumerate(row_ends):
            if row_end + padding <= read.ref_start:
                rows[i].append(read)
                row_ends[i] = max(row_end, read.ref_end)
                placed = True
                break
        if not placed:
            rows.append([read])
            row_ends.append(read.ref_end)
    return rows


def expand_rows(
    reads: List[AlignedRead],
    sort_by: str = "gap_length",
    descending: bool = True,
) -> List[List[AlignedRead]]:
    """One read per row, rows ordered by sort_by (descending by default so
    the most heavily gapped alignments land at the top)."""
    ordered = _ordered(reads, sort_by, descending)
    return [[read] for read in ordered]


def build_rows(
    reads: List[AlignedRead],
    layout: str,
    sort_by: str,
    descending: bool,
    padding: int = 2,
) -> List[List[AlignedRead]]:
    if layout == "expand":
        return expand_rows(reads, sort_by=sort_by, descending=descending)
    if layout == "pack":
        return pack_rows(reads, sort_by=sort_by, descending=descending, padding=padding)
    raise ValueError(f"Unknown layout '{layout}'. Choose 'pack' or 'expand'.")


def truncate_rows(
    rows: List[List[AlignedRead]], max_rows: Optional[int]
) -> tuple[List[List[AlignedRead]], int]:
    """Cap the number of rows drawn. Since rows are already sorted by
    priority, truncating keeps the most relevant reads and drops the rest -
    returns (kept_rows, n_reads_dropped)."""
    if max_rows is None or len(rows) <= max_rows:
        return rows, 0
    dropped = sum(len(r) for r in rows[max_rows:])
    return rows[:max_rows], dropped
