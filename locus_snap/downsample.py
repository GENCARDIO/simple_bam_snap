"""Deterministic alignment-track downsampling for high-depth regions."""
from __future__ import annotations

import hashlib
from typing import Iterable, List, Optional, Sequence, Set, Tuple

from locus_snap.read_model import AlignedRead


DEFAULT_MAX_ALIGNMENT_DEPTH = 100


def max_alignment_depth(reads: Iterable[AlignedRead]) -> int:
    """Maximum number of alignment spans overlapping any genomic position."""
    events = []
    for read in reads:
        start = read.ref_start
        end = max(read.ref_end, start + 1)
        events.append((start, 1))
        events.append((end, -1))
    depth = maximum = 0
    # End events (-1) sort before start events (+1), matching half-open spans.
    for event in sorted(events):
        depth += event[1]
        maximum = max(maximum, depth)
    return maximum


def retention_priority(read: AlignedRead, priority_names: Set[str]) -> tuple:
    """Higher tuples are retained when overlapping reads compete for space."""
    evidence = sum(
        (
            bool(getattr(read, "is_discordant", False)),
            getattr(read, "gap_length", 0) > 0,
            getattr(read, "sa_count", 0) > 0,
            getattr(read, "soft_clip_total", 0) > 0,
        )
    )
    return (
        read.query_name in priority_names,
        evidence,
        getattr(read, "mapq", 0),
        int.from_bytes(hashlib.blake2b(
            (
                f"{read.query_name}\0{read.ref_start}\0{read.ref_end}\0"
                f"{getattr(read, 'flag', 0)}"
            ).encode("utf-8", errors="replace"),
            digest_size=8,
        ).digest(), "big"),
    )


def downsample_reads(
    reads: Sequence[AlignedRead],
    max_depth: Optional[int] = DEFAULT_MAX_ALIGNMENT_DEPTH,
    priority_names: Optional[Set[str]] = None,
    preserve_pairs: bool = False,
) -> Tuple[List[AlignedRead], int]:
    """Cap concurrent alignment spans while preserving useful evidence.

    The selection is deterministic. When more than ``max_depth`` reads
    overlap, exact mate-view supporters, SV-evidence reads, higher MAPQ reads,
    and finally a stable read-name hash are preferred in that order. A cap of
    zero or ``None`` disables downsampling.
    """
    if max_depth is None or max_depth == 0:
        return list(reads), 0
    if max_depth < 0:
        raise ValueError("Maximum alignment depth cannot be negative.")

    protected = priority_names or set()
    kept = [False] * len(reads)
    active: List[int] = []
    interval_ends = []
    for read in reads:
        interval_ends.append(max(read.ref_end, read.ref_start + 1))
    ordered_indices = sorted(
        range(len(reads)),
        key=lambda i: (reads[i].ref_start, interval_ends[i], reads[i].query_name),
    )

    for index in ordered_indices:
        start = reads[index].ref_start
        still_active = []
        for active_index in active:
            if kept[active_index] and interval_ends[active_index] > start:
                still_active.append(active_index)
        active = still_active
        kept[index] = True
        active.append(index)
        if len(active) > max_depth:
            victim = min(
                active,
                key=lambda active_index: retention_priority(reads[active_index], protected),
            )
            kept[victim] = False
            active.remove(victim)

    selected = []
    for index, read in enumerate(reads):
        if kept[index]:
            selected.append(read)
    if preserve_pairs:
        available_counts = {}
        selected_counts = {}
        for collection, counts in ((reads, available_counts), (selected, selected_counts)):
            for read in collection:
                if (
                    getattr(read, "is_paired", False)
                    and not getattr(read, "is_secondary", False)
                    and not getattr(read, "is_supplementary", False)
                    and not getattr(read, "mate_is_unmapped", False)
                    and getattr(read, "mate_chrom", None) == getattr(read, "reference_name", None)
                ):
                    key = (read.query_name, getattr(read, "reference_name", None))
                    counts[key] = counts.get(key, 0) + 1
        incomplete_names = set()
        for key, count in available_counts.items():
            if count >= 2 and selected_counts.get(key, 0) < 2:
                incomplete_names.add(key)
        complete_pairs = []
        for read in selected:
            key = (read.query_name, getattr(read, "reference_name", None))
            if key not in incomplete_names:
                complete_pairs.append(read)
        selected = complete_pairs
    return selected, len(reads) - len(selected)
