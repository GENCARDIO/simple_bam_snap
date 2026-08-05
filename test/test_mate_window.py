import os
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.mate_window import choose_mate_window, mate_candidates, supporting_query_names
from src.annotations import AnnotationSource
from src.read_model import SAEntry
from src.snapshot import BamSnapshot

TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def fake_read(
    *, name="read", discordant=False, mate_chrom=None, mate_start=None, softclip=0, sa_entries=None
):
    return SimpleNamespace(
        query_name=name,
        is_discordant=discordant,
        mate_chrom=mate_chrom,
        mate_start=mate_start,
        soft_clip_total=softclip,
        sa_entries=sa_entries or [],
    )


def test_discordant_window_uses_most_supported_chromosome_and_mean_position():
    reads = [
        fake_read(discordant=True, mate_chrom="chr2", mate_start=100),
        fake_read(discordant=True, mate_chrom="chr2", mate_start=300),
        fake_read(discordant=True, mate_chrom="chr1", mate_start=900),
        fake_read(discordant=False, mate_chrom="chr2", mate_start=500),
    ]

    window = choose_mate_window(reads, "discordant", window_size=100)

    assert (window.chrom, window.start, window.end) == ("chr2", 150, 250)
    assert window.candidate_count == 2


def test_mate_window_is_clamped_to_contig_edges_without_changing_size():
    reads = [fake_read(discordant=True, mate_chrom="chr1", mate_start=995)]
    window = choose_mate_window(
        reads, "discordant", window_size=100, contig_lengths={"chr1": 1000}
    )
    assert (window.start, window.end) == (900, 1000)


def test_split_source_uses_sa_entry_centres():
    entries = [
        SAEntry("chr3", 100, 150, "+", "50M", 60, 0),
        SAEntry("chr3", 300, 400, "-", "100M", 60, 0),
    ]
    window = choose_mate_window(
        [fake_read(sa_entries=entries)], "split", window_size=50
    )
    # Mean of SA centres 125 and 350 is 237.5; Python rounds ties to even.
    assert (window.chrom, window.start, window.end) == ("chr3", 213, 263)


def test_softclip_source_uses_only_reads_meeting_threshold():
    reads = [
        fake_read(mate_chrom="chr4", mate_start=100, softclip=4),
        fake_read(mate_chrom="chr4", mate_start=300, softclip=8),
    ]
    assert mate_candidates(reads, "softclip", min_softclip=5) == [("chr4", 300)]


def test_supporting_names_identify_reads_behind_selected_mate_chromosome():
    reads = [
        fake_read(name="keep", discordant=True, mate_chrom="chr2", mate_start=100),
        fake_read(name="other-chrom", discordant=True, mate_chrom="chr3", mate_start=100),
        fake_read(name="normal", discordant=False, mate_chrom="chr2", mate_start=100),
    ]
    assert supporting_query_names(reads, "discordant", "chr2") == {"keep"}


def test_no_mate_candidates_has_actionable_error():
    with pytest.raises(ValueError, match="no mapped mates of discordant reads"):
        choose_mate_window([], "discordant", window_size=100)


def test_bam_snapshot_renders_two_panel_softclip_mate_view(tmp_path):
    annotation = tmp_path / "regions.bed"
    cnv = tmp_path / "copy-number.seg"
    annotation.write_text(
        "chr9\t101867400\t101867460\tFEATURE_A\n"
        "chr9\t101867500\t101867590\tFEATURE_B\n",
        encoding="utf-8",
    )
    cnv.write_text(
        "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
        "Tumour\tchr9\t1\t141213431\t100\t0.35\n",
        encoding="utf-8",
    )
    snap = BamSnapshot(
        bam=TEST_BAM,
        chrom="chr9",
        start=101867480,
        end=101867620,
        output_dir=str(tmp_path),
        output_name="mate-test.png",
        mate_view=True,
        mate_window_source="softclip",
        min_softclip=3,
        max_rows=3,
        display_mode="squish",
        haplotype_view="split",
        show_coverage=False,
        fig_width=6,
        dpi=40,
        annotation_sources=[
            AnnotationSource(str(annotation), "Regions"),
            AnnotationSource(str(cnv), "Copy number"),
        ],
    )

    summary = snap.snap()

    assert summary.n_reads > 0
    assert snap.mate_window is not None
    assert snap.mate_window.chrom == "chr9"
    assert os.path.isfile(snap.output_png)
    assert os.path.getsize(snap.output_png) > 0
