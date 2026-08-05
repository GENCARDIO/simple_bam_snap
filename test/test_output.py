import os
import struct
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.snapshot import BamSnapshot, resolve_output_path


TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def test_output_path_defaults_infers_and_overrides_formats(tmp_path):
    output_dir = str(tmp_path)

    assert resolve_output_path(output_dir, None, "locus").endswith("locus.png")
    assert resolve_output_path(output_dir, "locus.svg", "unused").endswith("locus.svg")
    assert resolve_output_path(
        output_dir, "locus.png", "unused", "pdf"
    ).endswith("locus.pdf")
    assert resolve_output_path(
        output_dir, "locus.v1", "unused", "svg"
    ).endswith("locus.v1.svg")


def test_unknown_output_extension_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="Unsupported output filename extension"):
        resolve_output_path(str(tmp_path), "locus.pmg", "unused")


def test_snapshot_writes_svg_selected_by_filename(tmp_path):
    snap = BamSnapshot(
        bam=TEST_BAM,
        chrom="chr9",
        start=101867480,
        end=101867620,
        output_dir=str(tmp_path),
        output_name="vector.svg",
        show_ideogram=False,
        show_coverage=False,
        max_rows=2,
        fig_width=4,
        dpi=40,
    )

    snap.snap()

    assert snap.output_path.endswith("vector.svg")
    assert snap.output_png == snap.output_path
    assert "<svg" in (tmp_path / "vector.svg").read_text(encoding="utf-8")[:500]


def test_png_pixel_width_follows_figure_width_and_dpi(tmp_path):
    snap = BamSnapshot(
        bam=TEST_BAM,
        chrom="chr9",
        start=101867480,
        end=101867620,
        output_dir=str(tmp_path),
        output_name="raster",
        output_format="png",
        show_ideogram=False,
        show_coverage=False,
        max_rows=1,
        fig_width=4,
        dpi=50,
    )

    snap.snap()

    with open(snap.output_path, "rb") as image_file:
        header = image_file.read(24)
    width = struct.unpack(">I", header[16:20])[0]
    assert width == 200
