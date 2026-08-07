"""Turns a flat list of AlignedRead objects into rows to draw.

Two layout strategies and three display modes are supported:

- ``pack``: classic genome-browser row packing - reads are greedily
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

Display mode is orthogonal to layout: ``expand`` uses normal-height rows,
``squish`` uses shorter rows, and ``collapse`` overlays every read in one row
regardless of the selected layout.
"""
from __future__ import annotations

from bisect import bisect_left
from functools import partial
from operator import attrgetter
from typing import Callable, List, Optional, Tuple

from locus_snap.read_model import AlignedRead

SortKeyFunc = Callable[[AlignedRead], object]


def zero_sort_key(read):
    return 0


SORT_KEYS: dict[str, SortKeyFunc] = {
    "base": zero_sort_key,
    "gap_length": attrgetter("gap_length"),
    "cigar_gap": attrgetter("cigar_gap_len"),
    "sa_gap": attrgetter("sa_gap_len"),
    "sa_count": attrgetter("sa_count"),
    "soft_clip": attrgetter("soft_clip_total"),
    "mismatch": attrgetter("mismatch_count"),
    "mapq": attrgetter("mapq"),
    "insert_size": attrgetter("insert_size"),
    "start": attrgetter("ref_start"),
    "strand": attrgetter("strand"),
    "read_name": attrgetter("query_name"),
    "none": zero_sort_key,
}

BASE_SORT_ORDER = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4}

DISPLAY_MODES = ("collapse", "expand", "squish")
HAPLOTYPE_VIEWS = ("none", "color", "split")


def resolve_sort_key(name: str) -> SortKeyFunc:
    try:
        return SORT_KEYS[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown sort key '{name}'. Choose from: {', '.join(sorted(SORT_KEYS))}"
        ) from exc


def base_sort_key(
    read: AlignedRead, position: int, reference_base: Optional[str] = None
) -> tuple:
    """Rank alternate bases before reference, gaps, and uncovered reads."""
    observed = read.base_at(position)
    observed = observed.upper() if observed else None
    if observed and observed in "ACGT":
        reference = (reference_base or "").upper()
        is_variant = bool(reference in "ACGT" and observed != reference)
        tier = 4 if is_variant else 3
        return tier, BASE_SORT_ORDER[observed], read.mapq
    if observed == "-":
        return 2, 0, read.mapq
    if observed == "~":
        return 1, 0, read.mapq
    return 0, 0, read.mapq


def infer_reference_base(reads: List[AlignedRead], position: int) -> Optional[str]:
    """Use the most frequent observed nucleotide when FASTA is unavailable."""
    counts = {}
    for read in reads:
        observed = read.base_at(position)
        if observed in ("A", "C", "G", "T"):
            counts[observed] = counts.get(observed, 0) + 1
    if not counts:
        return None
    best_base = None
    for base in counts:
        if best_base is None or (counts[base], base) > (counts[best_base], best_base):
            best_base = base
    return best_base


def order_reads(
    reads: List[AlignedRead], sort_by: str, descending: bool,
    base_position: Optional[int] = None, reference_base: Optional[str] = None,
) -> List[AlignedRead]:
    if sort_by == "base":
        if base_position is None:
            raise ValueError("Base sorting requires a genomic base position.")
        key_func = partial(
            base_sort_key, position=base_position, reference_base=reference_base
        )
    else:
        key_func = resolve_sort_key(sort_by)

    def full_sort_key(r):
        return (key_func(r), -r.ref_start if descending else r.ref_start)

    return sorted(reads, key=full_sort_key, reverse=descending)


def group_reads(
    reads: List[AlignedRead], sort_by: str, descending: bool,
    view_as_pairs: bool = False,
    base_position: Optional[int] = None, reference_base: Optional[str] = None,
) -> List[List[AlignedRead]]:
    """Order reads and optionally keep visible primary mates as one unit."""
    ordered = order_reads(
        reads, sort_by, descending,
        base_position=base_position, reference_base=reference_base,
    )
    if not view_as_pairs:
        groups = []
        for read in ordered:
            groups.append([read])
        return groups

    pair_members = {}
    for read in ordered:
        if (
            getattr(read, "is_paired", False)
            and not getattr(read, "is_secondary", False)
            and not getattr(read, "is_supplementary", False)
            and not getattr(read, "mate_is_unmapped", False)
            and getattr(read, "mate_chrom", None) == getattr(read, "reference_name", None)
        ):
            key = (read.query_name, getattr(read, "reference_name", None))
            pair_members.setdefault(key, []).append(read)

    groups = []
    emitted_pairs = set()
    for read in ordered:
        key = (read.query_name, getattr(read, "reference_name", None))
        members = pair_members.get(key, [])
        if read in members and len(members) >= 2:
            if key not in emitted_pairs:
                groups.append(sorted(members, key=attrgetter("ref_start")))
                emitted_pairs.add(key)
        else:
            groups.append([read])
    return groups


def pack_rows(
    reads: List[AlignedRead],
    sort_by: str = "start",
    descending: bool = False,
    padding: int = 2,
    view_as_pairs: bool = False,
    base_position: Optional[int] = None, reference_base: Optional[str] = None,
) -> List[List[AlignedRead]]:
    """Greedily pack priority-ordered groups into non-overlapping rows.

    A row can accept a group before, after, or between intervals already in
    that row. This matters when the priority sort is not genomic: a late
    discordant read must not permanently block the empty left side of a row.
    """
    groups = group_reads(
        reads, sort_by, descending, view_as_pairs,
        base_position=base_position, reference_base=reference_base,
    )

    rows: List[List[AlignedRead]] = []
    row_intervals: List[List[Tuple[int, int]]] = []
    for group in groups:
        group_start = None
        group_end = None
        for read in group:
            if group_start is None or read.ref_start < group_start:
                group_start = read.ref_start
            if group_end is None or read.ref_end > group_end:
                group_end = read.ref_end
        placed = False
        for row_index, intervals in enumerate(row_intervals):
            insert_at = bisect_left(intervals, (group_start, group_end))
            overlaps_left = (
                insert_at > 0
                and intervals[insert_at - 1][1] + padding > group_start
            )
            overlaps_right = (
                insert_at < len(intervals)
                and group_end + padding > intervals[insert_at][0]
            )
            if overlaps_left or overlaps_right:
                continue
            intervals.insert(insert_at, (group_start, group_end))
            rows[row_index].extend(group)
            rows[row_index].sort(key=attrgetter("ref_start"))
            placed = True
            break
        if not placed:
            rows.append(list(group))
            row_intervals.append([(group_start, group_end)])
    return rows


def expand_rows(
    reads: List[AlignedRead],
    sort_by: str = "gap_length",
    descending: bool = True,
    view_as_pairs: bool = False,
    base_position: Optional[int] = None, reference_base: Optional[str] = None,
) -> List[List[AlignedRead]]:
    """One read per row, rows ordered by sort_by (descending by default so
    the most heavily gapped alignments land at the top)."""
    return group_reads(
        reads, sort_by, descending, view_as_pairs,
        base_position=base_position, reference_base=reference_base,
    )


def build_rows(
    reads: List[AlignedRead],
    layout: str,
    sort_by: str,
    descending: bool,
    padding: int = 2,
    display_mode: str = "expand",
    view_as_pairs: bool = False,
    haplotype_view: str = "none",
    base_position: Optional[int] = None, reference_base: Optional[str] = None,
) -> List[List[AlignedRead]]:
    if display_mode not in DISPLAY_MODES:
        raise ValueError(
            f"Unknown display mode '{display_mode}'. Choose from: {', '.join(DISPLAY_MODES)}."
        )
    if haplotype_view not in HAPLOTYPE_VIEWS:
        raise ValueError(
            f"Unknown haplotype view '{haplotype_view}'. Choose from: "
            f"{', '.join(HAPLOTYPE_VIEWS)}."
        )
    if haplotype_view == "split":
        grouped_reads = {}
        for read in reads:
            grouped_reads.setdefault(getattr(read, "haplotype", None), []).append(read)

        def haplotype_label_key(label):
            if label is None:
                return 2, ""
            if str(label).isdigit():
                return 0, int(label)
            return 1, str(label)

        labels = sorted(grouped_reads, key=haplotype_label_key)
        rows = []
        for label in labels:
            rows.extend(build_rows(
                grouped_reads[label], layout=layout, sort_by=sort_by,
                descending=descending, padding=padding,
                display_mode=display_mode, view_as_pairs=view_as_pairs,
                haplotype_view="none",
                base_position=base_position, reference_base=reference_base,
            ))
        return rows
    if display_mode == "collapse":
        groups = group_reads(
            reads, sort_by, descending, view_as_pairs,
            base_position=base_position, reference_base=reference_base,
        )
        ordered = []
        for group in groups:
            ordered.extend(group)
        return [ordered] if ordered else []
    if layout == "expand":
        return expand_rows(
            reads, sort_by=sort_by, descending=descending,
            view_as_pairs=view_as_pairs,
            base_position=base_position, reference_base=reference_base,
        )
    if layout == "pack":
        return pack_rows(
            reads, sort_by=sort_by, descending=descending, padding=padding,
            view_as_pairs=view_as_pairs,
            base_position=base_position, reference_base=reference_base,
        )
    raise ValueError(f"Unknown layout '{layout}'. Choose 'pack' or 'expand'.")


def truncate_rows(
    rows: List[List[AlignedRead]], max_rows: Optional[int]
) -> tuple[List[List[AlignedRead]], int]:
    """Cap the number of rows drawn. Since rows are already sorted by
    priority, truncating keeps the most relevant reads and drops the rest -
    returns (kept_rows, n_reads_dropped)."""
    if max_rows is None or len(rows) <= max_rows:
        return rows, 0
    dropped = 0
    for r in rows[max_rows:]:
        dropped += len(r)
    return rows[:max_rows], dropped
