import os
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from simple_bam_snap import apply_config_preferences, build_parser
from src.config import DEFAULT_ALIGNMENT_COLORS, load_alignment_colors, load_config
from src.render import AlignmentRenderer, INSERTION_COLOR


def test_partial_alignment_color_config_inherits_defaults(tmp_path):
    config = tmp_path / "colors.yaml"
    config.write_text(
        "alignment_colors:\n"
        "  normal: '#123456'\n"
        "  ff: gold\n"
        "  interchrom: '#abcdef'\n",
        encoding="utf-8",
    )

    colors = load_alignment_colors(str(config))

    assert colors["normal"] == "#123456"
    assert colors["ff"] == "gold"
    assert colors["interchrom"] == "#abcdef"
    assert colors["large_insert"] == DEFAULT_ALIGNMENT_COLORS["large_insert"]


@pytest.mark.parametrize(
    "yaml_text, message",
    [
        ("alignment_colors:\n  typo: red\n", "Unknown alignment color"),
        ("alignment_colors:\n  normal: definitely-not-a-color\n", "Invalid color"),
        ("alignment_colors: blue\n", "must be a YAML mapping"),
    ],
)
def test_invalid_alignment_color_config_is_rejected(tmp_path, yaml_text, message):
    config = tmp_path / "bad.yaml"
    config.write_text(yaml_text, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_alignment_colors(str(config))


def test_renderer_uses_configured_pair_category_color():
    renderer = AlignmentRenderer(
        alignment_colors={"ff": "#010203"}, shade_by_mapq=False
    )
    read = SimpleNamespace(
        pair_category="ff", mate_chrom="chr1", is_secondary=False,
        is_duplicate=False, mapq=60,
    )
    color, alpha = renderer.read_style(read)
    assert color == "#010203"
    assert alpha == 0.9


def test_default_small_insert_colour_matches_cigar_insertions():
    assert DEFAULT_ALIGNMENT_COLORS["small_insert"] == INSERTION_COLOR


def test_interchrom_null_keeps_chromosome_palette():
    renderer = AlignmentRenderer(
        alignment_colors={"interchrom": None}, shade_by_mapq=False
    )
    read = SimpleNamespace(
        pair_category="interchrom", mate_chrom="chr7", is_secondary=False,
        is_duplicate=False, mapq=60,
    )
    color = renderer.read_style(read)[0]
    assert color.startswith("#")


def test_complete_theme_sections_are_merged_and_used(tmp_path):
    config = tmp_path / "theme.yaml"
    config.write_text(
        "base_colors:\n"
        "  A: '#112233'\n"
        "track_colors:\n"
        "  bed: '#223344'\n"
        "visual_colors:\n"
        "  coverage: '#334455'\n"
        "haplotype_colors:\n"
        "  '1': '#445566'\n"
        "cytoband_colors:\n"
        "  acen: '#556677'\n"
        "chromosome_palette: ['#667788', '#778899']\n"
        "styles:\n"
        "  row_height_in: 0.31\n"
        "  alignment_alpha: 0.73\n",
        encoding="utf-8",
    )

    theme = load_config(str(config))
    renderer = AlignmentRenderer(visual_config=theme, shade_by_mapq=False)
    read = SimpleNamespace(
        pair_category="normal", mate_chrom="chr1", is_secondary=False,
        is_duplicate=False, mapq=60,
    )

    assert theme["base_colors"]["A"] == "#112233"
    assert theme["base_colors"]["C"] != "#112233"  # inherited default
    assert theme["track_colors"]["bed"] == "#223344"
    assert renderer.visual_colors["coverage"] == "#334455"
    assert renderer.haplotype_colors["1"] == "#445566"
    assert renderer.cytoband_colors["acen"] == "#556677"
    assert renderer.chromosome_palette == ["#667788", "#778899"]
    assert renderer.row_height_in == pytest.approx(0.31)
    assert renderer.read_style(read)[1] == pytest.approx(0.73)


def test_yaml_preferences_become_defaults_but_cli_still_wins():
    parser = build_parser()
    apply_config_preferences(parser, {
        "display_mode": "squish",
        "max_alignment_depth": 175,
        "view_as_pairs": True,
        "show_coverage": False,
        "include_supplementary": False,
    })

    configured = parser.parse_args([
        "--bam", "reads.bam", "--region", "chr1:1-10",
    ])
    overridden = parser.parse_args([
        "--bam", "reads.bam", "--region", "chr1:1-10",
        "--display_mode", "expand", "--max_alignment_depth", "80",
    ])

    assert configured.display_mode == "squish"
    assert configured.max_alignment_depth == 175
    assert configured.view_as_pairs
    assert configured.no_coverage
    assert configured.exclude_supplementary
    assert overridden.display_mode == "expand"
    assert overridden.max_alignment_depth == 80


@pytest.mark.parametrize(
    "yaml_text, message",
    [
        ("visual_colors:\n  typo: red\n", "Unknown visual_colors key"),
        ("styles:\n  alignment_alpha: 1.5\n", "between 0 and 1"),
        ("styles:\n  row_height_in: 0\n", "greater than zero"),
        ("styles:\n  center_guide_line_style: zigzag\n", "must be one of"),
        ("chromosome_palette: []\n", "non-empty YAML list"),
        ("mystery_section: {}\n", "Unknown config section"),
    ],
)
def test_invalid_theme_configuration_is_rejected(tmp_path, yaml_text, message):
    config = tmp_path / "bad-theme.yaml"
    config.write_text(yaml_text, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_config(str(config))


def test_invalid_yaml_preference_is_rejected():
    parser = build_parser()
    with pytest.raises(ValueError, match="Invalid value for preference display_mode"):
        apply_config_preferences(parser, {"display_mode": "microscopic"})
    with pytest.raises(ValueError, match="Unknown preference"):
        apply_config_preferences(parser, {"not_an_option": 3})
