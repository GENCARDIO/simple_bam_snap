import os
import sys
from dataclasses import dataclass

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.downsample import downsample_reads, max_alignment_depth


@dataclass
class FakeRead:
    query_name: str
    ref_start: int
    ref_end: int
    mapq: int = 60
    flag: int = 0
    is_discordant: bool = False
    gap_length: int = 0
    sa_count: int = 0
    soft_clip_total: int = 0
    is_paired: bool = False
    is_secondary: bool = False
    is_supplementary: bool = False
    mate_is_unmapped: bool = False
    mate_chrom: str = "chr1"
    reference_name: str = "chr1"


def test_downsampling_caps_overlapping_alignment_depth_at_100():
    reads = [FakeRead(f"read-{i}", 100, 200) for i in range(150)]
    selected, dropped = downsample_reads(reads)
    assert len(selected) == 100
    assert dropped == 50
    assert max_alignment_depth(selected) == 100


def test_downsampling_does_not_remove_non_overlapping_reads():
    reads = [FakeRead(f"read-{i}", i * 20, i * 20 + 10) for i in range(150)]
    selected, dropped = downsample_reads(reads, max_depth=3)
    assert selected == reads
    assert dropped == 0


def test_zero_depth_cap_disables_downsampling():
    reads = [FakeRead(f"read-{i}", 100, 200) for i in range(150)]
    selected, dropped = downsample_reads(reads, max_depth=0)
    assert selected == reads
    assert dropped == 0


def test_negative_depth_cap_is_rejected():
    with pytest.raises(ValueError, match="cannot be negative"):
        downsample_reads([], max_depth=-1)


def test_sv_evidence_and_priority_names_win_competing_slots():
    reads = [
        FakeRead("ordinary", 100, 200, mapq=60),
        FakeRead("gapped", 100, 200, mapq=10, gap_length=30),
        FakeRead("mate-supporter", 100, 200, mapq=0),
    ]
    selected, dropped = downsample_reads(
        reads, max_depth=2, priority_names={"mate-supporter"}
    )
    assert {read.query_name for read in selected} == {"gapped", "mate-supporter"}
    assert dropped == 1


def test_downsampling_is_deterministic():
    reads = [FakeRead(f"read-{i}", 100, 200, mapq=i % 10) for i in range(30)]
    first, first_dropped = downsample_reads(reads, max_depth=8)
    second, second_dropped = downsample_reads(reads, max_depth=8)
    assert [read.query_name for read in first] == [read.query_name for read in second]
    assert first_dropped == second_dropped


def test_pair_preserving_downsampling_does_not_leave_one_visible_mate():
    reads = [
        FakeRead("pair", 100, 200, mapq=60, is_paired=True),
        FakeRead("pair", 100, 200, mapq=0, is_paired=True),
        FakeRead("ordinary", 100, 200, mapq=50),
    ]
    selected, dropped = downsample_reads(
        reads, max_depth=2, preserve_pairs=True
    )

    selected_names = [read.query_name for read in selected]
    assert selected_names.count("pair") in (0, 2)
    assert dropped == len(reads) - len(selected)
