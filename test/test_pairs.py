import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.read_model import (
    AlignedRead,
    classify_insert_sizes,
    compute_pair_orientation,
    matches_only,
)

HEADER = pysam.AlignmentHeader.from_dict(
    {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}],
    }
)


def make_segment(
    qname, chrom, pos, is_reverse, mate_chrom=None, mate_pos=None, mate_reverse=False,
    tlen=0, cigar="50M", mate_unmapped=False, paired=True,
):
    seg = pysam.AlignedSegment(HEADER)
    seg.query_name = qname
    seg.query_sequence = "A" * 50
    seg.query_qualities = pysam.qualitystring_to_array("I" * 50)
    seg.reference_id = HEADER.get_tid(chrom)
    seg.reference_start = pos
    seg.mapping_quality = 60
    seg.cigarstring = cigar
    seg.is_paired = paired
    seg.is_reverse = is_reverse
    if paired:
        seg.mate_is_unmapped = mate_unmapped
        if not mate_unmapped:
            seg.next_reference_id = HEADER.get_tid(mate_chrom or chrom)
            seg.next_reference_start = mate_pos if mate_pos is not None else 0
            seg.mate_is_reverse = mate_reverse
        seg.template_length = tlen
    return seg


def test_orientation_fr_normal_innie():
    seg = make_segment("r1", "chr1", 100, is_reverse=False, mate_pos=200, mate_reverse=True)
    assert compute_pair_orientation(seg) == "FR"


def test_orientation_fr_symmetric_when_this_read_is_rightmost():
    # this read is the rightmost (reverse) member of an otherwise-normal FR pair
    seg = make_segment("r1", "chr1", 200, is_reverse=True, mate_pos=100, mate_reverse=False)
    assert compute_pair_orientation(seg) == "FR"


def test_orientation_rf_everted():
    seg = make_segment("r1", "chr1", 100, is_reverse=True, mate_pos=200, mate_reverse=False)
    assert compute_pair_orientation(seg) == "RF"


def test_orientation_ff_same_strand():
    seg = make_segment("r1", "chr1", 100, is_reverse=False, mate_pos=200, mate_reverse=False)
    assert compute_pair_orientation(seg) == "FF"


def test_orientation_rr_same_strand():
    seg = make_segment("r1", "chr1", 100, is_reverse=True, mate_pos=200, mate_reverse=True)
    assert compute_pair_orientation(seg) == "RR"


def test_orientation_none_for_interchrom():
    seg = make_segment("r1", "chr1", 100, is_reverse=False, mate_chrom="chr2", mate_pos=200, mate_reverse=True)
    assert compute_pair_orientation(seg) is None


def test_orientation_none_for_mate_unmapped():
    seg = make_segment("r1", "chr1", 100, is_reverse=False, mate_unmapped=True)
    assert compute_pair_orientation(seg) is None


def test_orientation_none_for_unpaired():
    seg = make_segment("r1", "chr1", 100, is_reverse=False, paired=False)
    assert compute_pair_orientation(seg) is None


def test_pair_category_on_aligned_read():
    interchrom = AlignedRead(make_segment("a", "chr1", 100, False, mate_chrom="chr2", mate_pos=50, mate_reverse=True))
    assert interchrom.pair_category == "interchrom"
    assert interchrom.is_discordant

    everted = AlignedRead(make_segment("b", "chr1", 100, True, mate_pos=200, mate_reverse=False))
    assert everted.pair_category == "everted"

    ff = AlignedRead(make_segment("c", "chr1", 100, False, mate_pos=200, mate_reverse=False))
    assert ff.pair_category == "ff"
    assert ff.is_discordant

    rr = AlignedRead(make_segment("f", "chr1", 100, True, mate_pos=200, mate_reverse=True))
    assert rr.pair_category == "rr"
    assert rr.is_discordant

    normal = AlignedRead(make_segment("d", "chr1", 100, False, mate_pos=200, mate_reverse=True, tlen=176))
    assert normal.pair_category == "normal"
    assert not normal.is_discordant

    unpaired = AlignedRead(make_segment("e", "chr1", 100, False, paired=False))
    assert unpaired.pair_category == "normal"
    assert not unpaired.is_discordant


def test_classify_insert_sizes_flags_outliers_not_the_cluster():
    reads = []
    # a tight cluster of "normal" FR insert sizes around 200bp
    for i, size in enumerate([190, 195, 198, 200, 200, 202, 205, 208, 210, 195]):
        reads.append(AlignedRead(make_segment(f"norm{i}", "chr1", 1000 + i, False, mate_pos=1000 + i + size,
                                               mate_reverse=True, tlen=size)))
    large = AlignedRead(make_segment("big", "chr1", 5000, False, mate_pos=5000 + 5000, mate_reverse=True, tlen=5000))
    small = AlignedRead(make_segment("tiny", "chr1", 6000, False, mate_pos=6000 + 20, mate_reverse=True, tlen=20))
    reads += [large, small]

    classify_insert_sizes(reads, sigma=3.0)

    assert large.pair_category == "large_insert"
    assert large.is_large_insert
    assert small.pair_category == "small_insert"
    assert small.is_small_insert
    for r in reads[:10]:
        assert r.pair_category == "normal"


def test_classify_insert_sizes_noop_with_too_few_pairs():
    reads = [AlignedRead(make_segment("only", "chr1", 100, False, mate_pos=100000, mate_reverse=True, tlen=99900))]
    classify_insert_sizes(reads, sigma=3.0, min_pairs=10)
    assert reads[0].pair_category == "normal"  # not enough data to call it an outlier


def test_matches_only_discordant():
    interchrom = AlignedRead(make_segment("a", "chr1", 100, False, mate_chrom="chr2", mate_pos=50, mate_reverse=True))
    normal = AlignedRead(make_segment("b", "chr1", 100, False, mate_pos=200, mate_reverse=True, tlen=176))
    assert matches_only(interchrom, ["discordant"]) is True
    assert matches_only(normal, ["discordant"]) is False
    assert matches_only(normal, None) is True  # no filter => everything matches


def test_matches_only_softclip_respects_threshold():
    seg = make_segment("a", "chr1", 100, False, mate_pos=200, mate_reverse=True, tlen=176, cigar="3S47M")
    read = AlignedRead(seg)
    assert read.soft_clip_total == 3
    assert matches_only(read, ["softclip"], min_softclip=3) is True
    assert matches_only(read, ["softclip"], min_softclip=4) is False


def test_matches_only_gapped_and_split():
    gapped = AlignedRead(make_segment("g", "chr1", 100, False, mate_pos=300, mate_reverse=True, tlen=250, cigar="20M5D25M"))
    plain = AlignedRead(make_segment("p", "chr1", 100, False, mate_pos=300, mate_reverse=True, tlen=250, cigar="50M"))
    assert matches_only(gapped, ["gapped"]) is True
    assert matches_only(plain, ["gapped"]) is False
    assert matches_only(plain, ["split"]) is False  # no SA tag set on either fixture
