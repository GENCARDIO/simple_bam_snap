import os
import sys
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from src.layout import (
    build_rows, expand_rows, infer_reference_base, pack_rows, truncate_rows,
)


@dataclass
class FakeRead:
    ref_start: int
    ref_end: int
    gap_length: int = 0
    cigar_gap_len: int = 0
    sa_gap_len: int = 0
    sa_count: int = 0
    soft_clip_total: int = 0
    mismatch_count: int = 0
    mapq: int = 60
    insert_size: int = 0
    strand: str = "+"
    query_name: str = "r"
    is_paired: bool = False
    is_secondary: bool = False
    is_supplementary: bool = False
    mate_is_unmapped: bool = False
    mate_chrom: str = "chr1"
    reference_name: str = "chr1"
    haplotype: str = None
    phase_set: str = None
    bases: dict = None

    def base_at(self, position):
        return (self.bases or {}).get(position)


def test_pack_rows_never_overlaps_within_a_row():
    reads = [FakeRead(0, 10), FakeRead(5, 15), FakeRead(20, 30), FakeRead(9, 21)]
    rows = pack_rows(reads, sort_by="start", descending=False, padding=0)
    assert sum(len(r) for r in rows) == len(reads)
    for row in rows:
        row_sorted = sorted(row, key=lambda r: r.ref_start)
        for a, b in zip(row_sorted, row_sorted[1:]):
            assert a.ref_end <= b.ref_start


def test_pack_rows_compacts_non_overlapping_reads_into_one_row():
    reads = [FakeRead(0, 10), FakeRead(20, 30), FakeRead(40, 50)]
    rows = pack_rows(reads, sort_by="start", descending=False, padding=0)
    assert len(rows) == 1
    assert len(rows[0]) == 3


def test_pack_rows_can_insert_before_a_priority_read():
    reads = [
        FakeRead(100, 110, gap_length=10, query_name="priority-right"),
        FakeRead(0, 10, gap_length=0, query_name="left"),
    ]
    rows = pack_rows(reads, sort_by="gap_length", descending=True, padding=0)

    assert len(rows) == 1
    assert [read.query_name for read in rows[0]] == ["left", "priority-right"]


def test_expand_rows_one_read_per_row_sorted_desc_by_gap_length():
    reads = [
        FakeRead(0, 10, gap_length=1, query_name="a"),
        FakeRead(20, 30, gap_length=9, query_name="b"),
        FakeRead(40, 50, gap_length=5, query_name="c"),
    ]
    rows = expand_rows(reads, sort_by="gap_length", descending=True)
    assert [row[0].query_name for row in rows] == ["b", "c", "a"]
    assert all(len(row) == 1 for row in rows)


def test_expand_rows_ties_broken_by_position():
    reads = [
        FakeRead(100, 110, gap_length=0, query_name="late"),
        FakeRead(10, 20, gap_length=0, query_name="early"),
    ]
    rows = expand_rows(reads, sort_by="gap_length", descending=True)
    assert [row[0].query_name for row in rows] == ["early", "late"]


def test_base_sort_prioritises_non_reference_alleles_then_reference_and_gaps():
    reads = [
        FakeRead(0, 20, query_name="reference", bases={9: "G"}),
        FakeRead(0, 20, query_name="alt-a", bases={9: "A"}, mapq=50),
        FakeRead(0, 20, query_name="alt-t", bases={9: "T"}, mapq=60),
        FakeRead(0, 20, query_name="deletion", bases={9: "-"}),
        FakeRead(20, 30, query_name="uncovered"),
    ]

    rows = expand_rows(
        reads, sort_by="base", descending=True,
        base_position=9, reference_base="G",
    )
    names = []
    for row in rows:
        names.append(row[0].query_name)

    assert names == ["alt-t", "alt-a", "reference", "deletion", "uncovered"]


def test_base_sort_requires_a_position():
    with pytest.raises(ValueError, match="requires a genomic base position"):
        expand_rows([FakeRead(0, 10)], sort_by="base", descending=True)


def test_reference_base_can_be_inferred_from_observed_consensus():
    reads = [
        FakeRead(0, 20, bases={9: "C"}),
        FakeRead(0, 20, bases={9: "C"}),
        FakeRead(0, 20, bases={9: "A"}),
        FakeRead(20, 30),
    ]
    assert infer_reference_base(reads, 9) == "C"


def test_build_rows_dispatches_on_layout():
    reads = [FakeRead(0, 10, query_name="a"), FakeRead(20, 30, query_name="b")]
    assert len(build_rows(reads, layout="expand", sort_by="start", descending=False)) == 2
    assert len(build_rows(reads, layout="pack", sort_by="start", descending=False)) == 1


def test_view_as_pairs_keeps_mates_together_as_one_layout_unit():
    reads = [
        FakeRead(0, 10, query_name="pair", is_paired=True),
        FakeRead(90, 100, query_name="pair", is_paired=True),
        FakeRead(40, 50, query_name="between"),
    ]
    expanded = build_rows(
        reads, layout="expand", sort_by="start", descending=False,
        view_as_pairs=True,
    )
    assert [read.query_name for read in expanded[0]] == ["pair", "pair"]
    assert [read.query_name for read in expanded[1]] == ["between"]

    packed = build_rows(
        reads, layout="pack", sort_by="start", descending=False,
        view_as_pairs=True, padding=0,
    )
    assert len(packed) == 2
    assert [read.query_name for read in packed[0]] == ["pair", "pair"]


def test_collapse_display_mode_overlays_every_alignment_in_one_row():
    reads = [FakeRead(0, 10, query_name="a"), FakeRead(5, 15, query_name="b")]
    rows = build_rows(
        reads, layout="expand", sort_by="start", descending=False,
        display_mode="collapse",
    )
    assert len(rows) == 1
    assert [read.query_name for read in rows[0]] == ["a", "b"]


def test_unknown_display_mode_is_rejected():
    with pytest.raises(ValueError, match="Unknown display mode"):
        build_rows([], layout="pack", sort_by="start", descending=False, display_mode="tiny")


def test_split_haplotype_view_groups_and_orders_hp_lanes_even_when_collapsed():
    reads = [
        FakeRead(20, 30, query_name="hp2", haplotype="2", phase_set="200"),
        FakeRead(40, 50, query_name="untagged"),
        FakeRead(0, 10, query_name="hp1", haplotype="1", phase_set="100"),
    ]

    rows = build_rows(
        reads, layout="pack", sort_by="start", descending=False,
        display_mode="collapse", haplotype_view="split",
    )

    assert [[read.query_name for read in row] for row in rows] == [
        ["hp1"], ["hp2"], ["untagged"],
    ]


def test_unknown_haplotype_view_is_rejected():
    with pytest.raises(ValueError, match="Unknown haplotype view"):
        build_rows(
            [], layout="pack", sort_by="start", descending=False,
            haplotype_view="rainbow",
        )


def test_truncate_rows_keeps_highest_priority_rows():
    rows = [[FakeRead(0, 10)], [FakeRead(20, 30)], [FakeRead(40, 50), FakeRead(60, 70)]]
    kept, dropped = truncate_rows(rows, max_rows=2)
    assert len(kept) == 2
    assert dropped == 2  # the two reads in the third row

    kept, dropped = truncate_rows(rows, max_rows=None)
    assert kept == rows and dropped == 0
