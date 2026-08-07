import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.read_model import (
    AlignedRead,
    SAEntry,
    cigar_reference_length,
    closest_same_chrom_gap,
    fetch_reads,
    iter_cigar_blocks,
)

TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


class FakeReference:
    """Minimal stand-in for ReferenceWindow, just base_at()/available."""

    def __init__(self, start, seq):
        self.start = start
        self.seq = seq
        self.available = True

    def base_at(self, pos):
        idx = pos - self.start
        if 0 <= idx < len(self.seq):
            return self.seq[idx]
        return None


def test_iter_cigar_blocks_walks_ref_and_query_correctly():
    # 3S 5M 2I 3M 4D 6M
    cigartuples = [(4, 3), (0, 5), (1, 2), (0, 3), (2, 4), (0, 6)]
    blocks = list(iter_cigar_blocks(cigartuples, ref_start=100))
    ops = [b.op for b in blocks]
    assert ops == ["S", "M", "I", "M", "D", "M"]

    assert (blocks[0].ref_pos, blocks[0].query_pos) == (100, 0)   # S: consumes query only
    assert (blocks[1].ref_pos, blocks[1].query_pos) == (100, 3)   # M
    assert (blocks[2].ref_pos, blocks[2].query_pos) == (105, 8)   # I: ref unchanged
    assert (blocks[3].ref_pos, blocks[3].query_pos) == (105, 10)  # M
    assert (blocks[4].ref_pos, blocks[4].query_pos) == (108, 13)  # D: query unchanged
    assert (blocks[5].ref_pos, blocks[5].query_pos) == (112, 13)  # M


def test_cigar_reference_length():
    assert cigar_reference_length("45M9D31M") == 85
    assert cigar_reference_length("76M") == 76
    assert cigar_reference_length("5M2I5M") == 10  # insertion doesn't consume reference


def test_closest_same_chrom_gap_adjacent_deletion():
    sa = [SAEntry(rname="chr1", start=200, end=250, strand="+", cigar="50M", mapq=60, nm=0)]
    gap, cross = closest_same_chrom_gap(ref_start=100, ref_end=180, rname="chr1", sa_entries=sa)
    assert gap == 20  # 200 - 180
    assert cross is False


def test_closest_same_chrom_gap_overlap_is_zero():
    sa = [SAEntry(rname="chr1", start=150, end=250, strand="+", cigar="100M", mapq=60, nm=0)]
    gap, cross = closest_same_chrom_gap(100, 180, "chr1", sa)
    assert gap == 0
    assert cross is False


def test_closest_same_chrom_gap_cross_chrom_flagged_not_counted():
    sa = [SAEntry(rname="chr2", start=10, end=60, strand="+", cigar="50M", mapq=60, nm=0)]
    gap, cross = closest_same_chrom_gap(100, 180, "chr1", sa)
    assert gap == 0
    assert cross is True


def test_closest_same_chrom_gap_picks_nearest_of_several():
    sa = [
        SAEntry(rname="chr1", start=500, end=550, strand="+", cigar="50M", mapq=60, nm=0),
        SAEntry(rname="chr1", start=190, end=230, strand="+", cigar="40M", mapq=60, nm=0),
    ]
    gap, cross = closest_same_chrom_gap(100, 180, "chr1", sa)
    assert gap == 10  # nearest partner (190-180), not the far one
    assert cross is False


def test_fetch_reads_finds_the_known_9bp_deletion():
    reads = fetch_reads(TEST_BAM, "chr9", 101867490, 101867600)
    gapped = [r for r in reads if r.cigar_gap_len > 0]
    assert len(gapped) == 9
    assert all(r.cigar_gap_len == 9 for r in gapped)
    assert all(r.gap_length == 9 for r in gapped)
    assert all(r.gap_label() == "9bp" for r in gapped)
    ungapped = [r for r in reads if r.cigar_gap_len == 0]
    assert all(r.gap_length == 0 and r.gap_label() == "" for r in ungapped)


def test_fetch_reads_filters_by_mapq_and_duplicates():
    all_reads = fetch_reads(TEST_BAM, "chr9", 101867480, 101867620)
    strict = fetch_reads(TEST_BAM, "chr9", 101867480, 101867620, min_mapq=61)
    assert len(strict) == 0
    assert len(all_reads) > 0


def test_mismatch_detection_against_fake_reference():
    with pysam.AlignmentFile(TEST_BAM, "rb") as f:
        segment = next(f.fetch("chr9", 101867425, 101867426))
    seq = segment.query_sequence
    mutated_base = "A" if seq[0] != "A" else "C"
    mutated_ref = mutated_base + seq[1:]
    ref = FakeReference(segment.reference_start, mutated_ref)

    read = AlignedRead(segment, ref)
    assert read.mismatch_count == 1
    assert read.mismatches[0] == (segment.reference_start, seq[0].upper())
    assert read.mismatch_details[0][:2] == (segment.reference_start, seq[0].upper())
    assert read.mismatch_details[0][2] == segment.query_qualities[0]
    assert read.base_at(segment.reference_start) == seq[0].upper()
    assert read.base_at(segment.reference_start - 1) is None


def test_base_at_reports_deletions():
    reads = fetch_reads(TEST_BAM, "chr9", 101867490, 101867600)
    deletion_read = next(read for read in reads if read.deletions)
    deletion_start = deletion_read.deletions[0][0]
    assert deletion_read.base_at(deletion_start) == "-"


def test_no_reference_means_no_mismatch_computation():
    with pysam.AlignmentFile(TEST_BAM, "rb") as f:
        segment = next(f.fetch("chr9", 101867425, 101867426))
    read = AlignedRead(segment, reference=None)
    assert read.mismatch_count == 0
    assert read.mismatches == []
