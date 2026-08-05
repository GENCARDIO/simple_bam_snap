import os
import sys
from types import SimpleNamespace

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.annotations import AnnotationItem, LoadedAnnotationTrack
from src.cytobands import Cytoband
from src.read_model import CigarBlock
from src.render import (
    BASE_COLORS,
    CNV_GAIN_COLOR,
    CNV_LOSS_COLOR,
    IDEOGRAM_WINDOW_COLOR,
    AlignmentRenderer,
    compute_coverage,
    compute_snv_evidence,
    compute_snv_counts,
    ellipsize,
    haplotype_color,
)


def test_compute_coverage_counts_match_blocks_but_not_deletions():
    reads = [
        SimpleNamespace(blocks=[
            CigarBlock("M", 100, 0, 3),
            CigarBlock("D", 103, 3, 2),
            CigarBlock("M", 105, 3, 2),
        ]),
        SimpleNamespace(blocks=[CigarBlock("M", 101, 0, 5)]),
    ]
    assert compute_coverage(reads, 100, 107) == [1, 2, 2, 1, 1, 2, 1]


def test_coverage_colors_only_snvs_above_vaf_threshold():
    reads = []
    for index in range(5):
        mismatches = []
        if index < 2:
            mismatches.append((100, "A"))
        if index == 2:
            mismatches.append((100, "C"))
        reads.append(SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], mismatches=mismatches,
            mismatch_details=[(position, base, 35) for position, base in mismatches],
            query_sequence=mismatches[0][1] if mismatches else "G",
            query_qualities=[35], mapq=60, is_reverse=index % 2 == 1,
        ))

    assert compute_snv_counts(reads, 100, 101) == {100: {"A": 2, "C": 1}}

    renderer = AlignmentRenderer(
        coverage_vaf_threshold=0.20, show_variant_counts=True
    )
    fig, ax = plt.subplots()
    renderer.draw_coverage_track(ax, reads, 100, 101)

    assert len(ax.patches) == 2  # grey depth plus the 40% A allele; 20% C is hidden
    assert ax.patches[1].get_height() == 2
    assert to_hex(ax.patches[1].get_facecolor()) == BASE_COLORS["A"]
    assert any(
        text.get_text() == "A 2/5 40% F1/R1 BQ35 MQ60"
        for text in ax.texts
    )
    plt.close(fig)


def test_snv_evidence_applies_quality_filters_and_tracks_strand_means():
    reads = [
        SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], query_sequence="G",
            query_qualities=[35], mapq=60, is_reverse=False,
            mismatch_details=[(100, "G", 35)],
        ),
        SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], query_sequence="G",
            query_qualities=[25], mapq=40, is_reverse=True,
            mismatch_details=[(100, "G", 25)],
        ),
        SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], query_sequence="G",
            query_qualities=[10], mapq=60, is_reverse=False,
            mismatch_details=[(100, "G", 10)],
        ),
    ]

    depth, evidence = compute_snv_evidence(
        reads, 100, 101, min_baseq=20, min_mapq=30
    )
    allele = evidence[100]["G"]

    assert depth == [2]
    assert (allele.count, allele.forward, allele.reverse) == (2, 1, 1)
    assert allele.base_quality_sum / allele.count == 30
    assert allele.mapq_sum / allele.count == 50


def test_squish_rows_are_shorter_than_expanded_rows():
    expanded = AlignmentRenderer(display_mode="expand")
    squished = AlignmentRenderer(display_mode="squish")
    assert squished.row_height_in < expanded.row_height_in


def test_long_labels_are_ellipsized_to_the_available_lane():
    assert ellipsize("Candidate regions", 10) == "Candidate…"
    assert ellipsize("short", 10) == "short"


def test_legend_clusters_related_terms_by_topic():
    renderer = AlignmentRenderer(fig_width=14)
    fig = plt.figure(figsize=(14, 2))
    legends = renderer.draw_legends(fig, fig_height=2)

    assert [legend.get_title().get_text() for legend in legends] == [
        "Alignment events", "Pair evidence", "Base identity",
    ]
    assert [text.get_text() for text in legends[0].get_texts()] == [
        "Normal / concordant", "Insertion", "Deletion",
    ]
    assert [text.get_text() for text in legends[2].get_texts()] == list("ACGT")
    legend_ax = fig.axes[-1]
    assert len(legend_ax.patches) == 2  # alternating compartment fill + outer rail
    assert len(legend_ax.lines) == 2  # internal compartment dividers
    plt.close(fig)


def test_haplotype_view_colours_reads_and_replaces_pair_legend_compartment():
    renderer = AlignmentRenderer(haplotype_view="color")
    read = SimpleNamespace(
        haplotype="2", pair_category="large_insert", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
    )
    color, alpha = renderer.read_style(read)
    assert color == haplotype_color("2")
    assert alpha == pytest.approx(0.9)

    fig = plt.figure(figsize=(14, 2))
    legends = renderer.draw_legends(fig, fig_height=2)
    assert [legend.get_title().get_text() for legend in legends] == [
        "Alignment events", "Haplotype", "Base identity",
    ]
    assert [text.get_text() for text in legends[1].get_texts()] == [
        "HP 1", "HP 2", "Other HP", "Untagged",
    ]
    plt.close(fig)


def test_split_haplotype_lanes_show_hp_and_phase_set_labels():
    renderer = AlignmentRenderer(haplotype_view="split")
    rows = [
        [SimpleNamespace(haplotype="1", phase_set="100")],
        [SimpleNamespace(haplotype="1", phase_set="100")],
        [SimpleNamespace(haplotype="2", phase_set="200")],
        [SimpleNamespace(haplotype=None, phase_set=None)],
    ]
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_ylim(4, 0)

    renderer.draw_haplotype_lanes(ax, rows)

    assert [text.get_text() for text in ax.texts] == [
        "HP 1 · PS 100", "HP 2 · PS 200", "Untagged",
    ]
    assert len(ax.lines) == 2
    assert len(ax.patches) == 2
    plt.close(fig)


@pytest.mark.parametrize("fig_width", [5, 8, 14])
def test_legend_keeps_rendered_text_clear_of_plot_and_coordinate_labels(fig_width):
    renderer = AlignmentRenderer(fig_width=fig_width, view_as_pairs=True)
    fig_height = 6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(
        left=0.08, right=0.92, top=0.95,
        bottom=renderer.legend_margin_in / fig_height,
    )
    legends = renderer.draw_legends(fig, fig_height, 0.08, 0.92)
    renderer.separate_legend_from_plots(fig, [ax], legends)

    canvas_renderer = fig.canvas.get_renderer()
    legend_top = max(
        legends[0].axes.get_window_extent(canvas_renderer).y1,
        max(legend.get_window_extent(canvas_renderer).y1 for legend in legends),
    )
    plot_bottom = ax.get_tightbbox(canvas_renderer).y0
    assert plot_bottom >= legend_top + renderer.legend_plot_gap_in * fig.dpi - 1
    plt.close(fig)


def test_small_window_reference_track_draws_coloured_base_cells_and_letters():
    renderer = AlignmentRenderer(max_reference_span=4)
    sequence = "ACGT"
    reference = SimpleNamespace(
        available=True,
        base_at=lambda position: sequence[position],
    )
    fig, ax = plt.subplots(figsize=(4, 1))
    ax.set_xlim(0, 4)

    renderer.draw_reference_track(ax, reference, 0, 4, available_width_in=4)

    assert len(ax.patches) == 4
    assert [text.get_text() for text in ax.texts] == ["A", "C", "G", "T", "reference"]
    plt.close(fig)


def test_view_as_pairs_draws_a_link_between_visible_primary_mates():
    renderer = AlignmentRenderer(view_as_pairs=True, shade_by_mapq=False)
    common = dict(
        query_name="pair", reference_name="chr1", mate_chrom="chr1",
        is_paired=True, is_secondary=False, is_supplementary=False,
        mate_is_unmapped=False, is_duplicate=False, pair_category="normal",
        mapq=60, blocks=[], gap_label=lambda: "",
    )
    left = SimpleNamespace(ref_start=10, ref_end=20, **common)
    right = SimpleNamespace(ref_start=40, ref_end=50, **common)
    fig, ax = plt.subplots()

    renderer.draw_alignment_row(
        ax, [left, right], y0=0.1, h=0.8,
        render_base_detail=False, layout="expand",
    )

    assert len(ax.lines) == 1
    assert list(ax.lines[0].get_xdata()) == [20, 40]
    assert len(ax.patches) == 2
    plt.close(fig)


def test_base_sort_highlights_the_alternative_allele_cell():
    renderer = AlignmentRenderer(
        sort_base_position=10, sort_reference_base="C", shade_by_mapq=False
    )
    read = SimpleNamespace(
        ref_start=5, ref_end=15, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[CigarBlock("M", 5, 0, 10)], mismatches=[],
        query_sequence="CCCCCACCCC",
        base_at=lambda position: "A" if position == 10 else "C",
    )
    fig, ax = plt.subplots()

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=True)

    assert len(ax.patches) == 2
    assert ax.patches[-1].get_x() == 10
    assert to_hex(ax.patches[-1].get_facecolor()) == BASE_COLORS["A"]
    plt.close(fig)


def test_center_guide_is_hidden_by_default_and_centered_when_enabled():
    fig, (hidden_ax, visible_ax) = plt.subplots(nrows=2)
    AlignmentRenderer().draw_center_guide(hidden_ax, 100, 200)
    AlignmentRenderer(show_center_guide=True).draw_center_guide(
        visible_ax, 100, 200
    )

    assert len(hidden_ax.lines) == 0
    assert len(visible_ax.lines) == 1
    assert list(visible_ax.lines[0].get_xdata()) == [150, 150]
    assert visible_ax.lines[0].get_linestyle() == "--"
    plt.close(fig)


def test_ideogram_marks_the_window_in_red():
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()
    renderer.draw_ideogram(ax, "chr1", 100, 200, 1_000)
    fig.canvas.draw()

    assert len(ax.patches) == 2
    assert ax.patches[0].get_x() == 100
    assert ax.patches[0].get_width() == 100
    chromosome_bounds = ax.patches[0].get_window_extent(fig.canvas.get_renderer())
    axes_bounds = ax.get_window_extent()
    assert chromosome_bounds.x0 == pytest.approx(axes_bounds.x0)
    assert chromosome_bounds.x1 == pytest.approx(axes_bounds.x1)
    assert to_hex(ax.patches[1].get_facecolor()) == IDEOGRAM_WINDOW_COLOR
    assert {text.get_text() for text in ax.texts} == {"chr1", "0.0 Mb", "window"}
    plt.close(fig)


def test_ideogram_draws_cytobands_and_centromere():
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()
    bands = [
        Cytoband("chr1", 0, 400, "p11", "gneg"),
        Cytoband("chr1", 400, 500, "p10", "acen"),
        Cytoband("chr1", 500, 600, "q10", "acen"),
        Cytoband("chr1", 600, 1_000, "q11", "gpos100"),
    ]
    renderer.draw_ideogram(ax, "chr1", 100, 200, 1_000, bands)

    assert len(ax.patches) == 8  # base, four bands, bridge, outline, marker
    assert sum(isinstance(patch, Polygon) for patch in ax.patches) == 2
    chromosome_vertices = ax.patches[0].get_xy()
    neck_y = [
        y for x, y in chromosome_vertices if x == pytest.approx(150)
    ]
    assert len(neck_y) == 2
    assert min(neck_y) > 0.28
    assert max(neck_y) < 0.72
    assert min(neck_y) < 0.5 < max(neck_y)
    assert to_hex(ax.patches[5].get_facecolor()) == "#b84b4b"
    assert to_hex(ax.patches[-1].get_facecolor()) == IDEOGRAM_WINDOW_COLOR
    plt.close(fig)


def test_gene_track_repeats_orientation_arrows_across_introns():
    renderer = AlignmentRenderer()
    item = AnnotationItem(
        100, 200, "TX1", "+", blocks=[(100, 120), (180, 200)]
    )
    track = LoadedAnnotationTrack("Genes", "gtf", "#17217a", [item], [[item]])
    fig, ax = plt.subplots(figsize=(8, 1))
    ax.set_xlim(90, 210)

    renderer.draw_annotation_track(ax, track, 90, 210)

    assert sum(line.get_marker() == ">" for line in ax.lines) > 1
    assert not any(line.get_marker() == "<" for line in ax.lines)
    plt.close(fig)


def test_primary_isoform_label_has_visible_marker():
    renderer = AlignmentRenderer()
    item = AnnotationItem(
        100, 160, "GENE1", "+", blocks=[(100, 160)],
        transcript_label="TX1", primary_rank=1, primary_label="MANE Select",
    )
    track = LoadedAnnotationTrack(
        "Genes", "gtf", "#17217a", [item], [[item]], display_mode="expand"
    )
    fig, ax = plt.subplots(figsize=(8, 1))
    ax.set_xlim(90, 170)

    renderer.draw_annotation_track(ax, track, 90, 170)

    assert any(text.get_text() == "TX1 ★" for text in ax.texts)
    plt.close(fig)


def test_cnv_track_draws_log2_segments_around_zero_with_gain_loss_colours():
    renderer = AlignmentRenderer()
    loss = AnnotationItem(100, 130, value=-0.6, sample="Tumour")
    gain = AnnotationItem(130, 170, value=0.8, sample="Tumour")
    track = LoadedAnnotationTrack(
        "Tumour CNV", "seg", "#555555", [loss, gain], [[loss, gain]],
        color_by_sign=True,
    )
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_xlim(90, 180)

    renderer.draw_annotation_track(ax, track, 90, 180)

    assert len(ax.patches) == 2
    assert [to_hex(patch.get_facecolor()) for patch in ax.patches] == [
        CNV_LOSS_COLOR, CNV_GAIN_COLOR,
    ]
    assert [line.get_ydata()[0] for line in ax.lines] == [0, -0.6, 0.8]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["-1", "0", "1"]
    assert any(text.get_text() == "Tumour" for text in ax.texts)
    plt.close(fig)


def test_baf_track_draws_heterozygous_snvs_around_half_baseline():
    renderer = AlignmentRenderer()
    low = AnnotationItem(100, 101, "rs1", value=0.18, sample="Tumour")
    balanced = AnnotationItem(130, 131, "rs2", value=0.52, sample="Tumour")
    high = AnnotationItem(160, 161, "rs3", value=0.84, sample="Tumour")
    track = LoadedAnnotationTrack(
        "Tumour BAF", "baf", "#7a1f5c",
        [low, balanced, high], [[low, balanced, high]],
    )
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_xlim(90, 180)

    renderer.draw_annotation_track(ax, track, 90, 180)

    assert len(ax.collections) == 1
    assert ax.collections[0].get_offsets().tolist() == [
        [100.5, 0.18], [130.5, 0.52], [160.5, 0.84],
    ]
    assert list(ax.lines[0].get_ydata()) == [0.5, 0.5]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["0.0", "0.5", "1.0"]
    assert ax.get_ylabel() == "BAF"
    plt.close(fig)
