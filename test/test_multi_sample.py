import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.snapshot import compare_snapshots


TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def write_companion(path, position, identifier):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr9\t{position}\t{identifier}\tA\tG\t.\tPASS\tEND={position + 20}\n",
        encoding="utf-8",
    )


def test_three_bams_render_with_matched_vcf_companion_tracks(tmp_path):
    first_vcf = tmp_path / "first.vcf"
    second_vcf = tmp_path / "second.vcf"
    third_vcf = tmp_path / "third.vcf"
    write_companion(first_vcf, 101867500, "first-variant")
    write_companion(second_vcf, 101867540, "second-variant")
    write_companion(third_vcf, 101867590, "third-variant")

    output, summary = compare_snapshots(
        TEST_BAM,
        TEST_BAM,
        "chr9",
        101867480,
        101867620,
        additional_bams=[TEST_BAM],
        label1="Tumour",
        label2="Normal",
        additional_labels=["Relapse"],
        companion_vcfs=[str(first_vcf), str(second_vcf), str(third_vcf)],
        output_dir=str(tmp_path),
        output_name="multi-sample.svg",
        show_ideogram=False,
        show_coverage=False,
        max_rows=1,
        fig_width=7,
        dpi=40,
    )

    svg = (tmp_path / "multi-sample.svg").read_text(encoding="utf-8")
    assert output.endswith("multi-sample.svg")
    assert "Tumour variants" in svg
    assert "Normal variants" in svg
    assert "Relapse variants" in svg
    assert "first-variant" in svg
    assert "second-variant" in svg
    assert "third-variant" in svg
    assert "Tumour" in summary
    assert "Normal" in summary
    assert "Relapse" in summary
