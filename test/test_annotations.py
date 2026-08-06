import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.annotations import (
    AnnotationSource,
    BafSource,
    build_annotation_sources,
    build_baf_sources,
    build_custom_annotation_sources,
    infer_track_format,
    normalize_track_color,
)
from src.snapshot import compare_snapshots

TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def test_format_inference_handles_plain_and_compressed_names():
    assert infer_track_format("genes.gtf") == "gtf"
    assert infer_track_format("genes.gff3.gz") == "gff3"
    assert infer_track_format("regions.bed.bgz") == "bed"
    assert infer_track_format("variants.vcf.gz") == "vcf"
    assert infer_track_format("tumour.seg") == "seg"
    assert infer_track_format("bins.bedgraph.gz") == "bedgraph"
    assert infer_track_format("ratios.bdg") == "bedgraph"
    assert infer_track_format("sample.log2") == "log2"
    assert infer_track_format("H3K27ac.narrowPeak.gz") == "narrowpeak"
    assert infer_track_format("H3K27me3.broadPeak") == "broadpeak"
    assert infer_track_format("DNase.signal.bgz") == "signal"
    with pytest.raises(ValueError, match="Cannot infer"):
        infer_track_format("track.txt")


def test_bed12_builds_thick_coding_blocks_and_thin_utrs(tmp_path):
    bed = tmp_path / "genes.bed"
    bed.write_text(
        # Two exons: [100,150), [200,260); coding region [120,230).
        "chr1\t100\t260\tTX1\t0\t+\t120\t230\t0\t2\t50,60\t0,100\n",
        encoding="utf-8",
    )
    track = AnnotationSource(str(bed), "Genes").fetch("chr1", 90, 270)
    item = track.items[0]
    assert track.label == "Genes"
    assert track.color == "#000000"
    assert item.blocks == [(120, 150), (200, 230)]
    assert item.utrs == [(100, 120), (230, 260)]
    assert item.strand == "+"


def test_gtf_groups_features_by_transcript_and_infers_utrs(tmp_path):
    gtf = tmp_path / "genes.gtf"
    attributes = 'gene_id "g1"; transcript_id "tx1"; gene_name "GENE1";'
    gtf.write_text(
        "chr1\ttest\texon\t101\t150\t.\t-\t.\t" + attributes + "\n"
        "chr1\ttest\texon\t201\t260\t.\t-\t.\t" + attributes + "\n"
        "chr1\ttest\tCDS\t121\t150\t.\t-\t0\t" + attributes + "\n"
        "chr1\ttest\tCDS\t201\t230\t.\t-\t0\t" + attributes + "\n",
        encoding="utf-8",
    )
    item = AnnotationSource(str(gtf)).fetch("chr1", 90, 270).items[0]
    assert item.name == "GENE1"
    assert item.strand == "-"
    assert item.blocks == [(120, 150), (200, 230)]
    assert item.utrs == [(100, 120), (230, 260)]


def test_annotation_default_colours_depend_on_track_type(tmp_path):
    bed = tmp_path / "regions.bed"
    gtf = tmp_path / "genes.gtf"
    vcf = tmp_path / "variants.vcf"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")
    gtf.write_text(
        'chr1\ttest\texon\t11\t20\t.\t+\t.\tgene_id "g"; transcript_id "t";\n',
        encoding="utf-8",
    )
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t15\t.\tA\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )

    assert AnnotationSource(str(bed)).fetch("chr1", 0, 30).color == "#000000"
    assert AnnotationSource(str(gtf)).fetch("chr1", 0, 30).color == "#17217a"
    assert AnnotationSource(str(vcf)).fetch("chr1", 0, 30).color == "#7a1f5c"
    assert AnnotationSource(str(bed), color="#abcdef").fetch("chr1", 0, 30).color == "#abcdef"


def test_yaml_track_palette_changes_type_defaults(tmp_path):
    bed = tmp_path / "regions.bed"
    gtf = tmp_path / "genes.gtf"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")
    gtf.write_text(
        'chr1\ttest\texon\t11\t20\t.\t+\t.\tgene_id "g"; transcript_id "t";\n',
        encoding="utf-8",
    )
    colors = {
        "bed": "#112233", "gene": "#223344", "vcf": "#334455",
        "cnv": "#445566", "baf": "#556677",
    }

    assert AnnotationSource(
        str(bed), track_colors=colors
    ).fetch("chr1", 0, 30).color == "#112233"
    assert AnnotationSource(
        str(gtf), track_colors=colors
    ).fetch("chr1", 0, 30).color == "#223344"


def test_custom_track_accepts_explicit_type_name_and_rgb_colour(tmp_path):
    track_file = tmp_path / "regions.data"
    track_file.write_text("chr1\t10\t20\tregion\n", encoding="utf-8")

    sources = build_custom_annotation_sources([
        [str(track_file), "bed", "Important regions", "rgb(12, 34, 56)"],
    ])
    track = sources[0].fetch("chr1", 0, 30)

    assert track.kind == "bed"
    assert track.label == "Important regions"
    assert track.color == "#0c2238"
    assert len(track.items) == 1


def test_custom_track_accepts_per_track_display_override(tmp_path):
    track_file = tmp_path / "regions.data"
    track_file.write_text(
        "chr1\t10\t20\ta\nchr1\t30\t40\tb\n", encoding="utf-8"
    )
    source = build_custom_annotation_sources([
        [str(track_file), "bed", "Regions", "#000000", "expand"],
    ])[0]

    track = source.fetch("chr1", 0, 50)
    assert source.display_mode == "expand"
    assert len(track.rows) == 2


def test_custom_track_accepts_per_track_height_override(tmp_path):
    track_file = tmp_path / "regions.data"
    track_file.write_text("chr1\t10\t20\tregion\n", encoding="utf-8")
    specification = f"{track_file},bed,Regions,rgb(12,34,56),pack,0.65"

    source = build_custom_annotation_sources([[specification]])[0]
    track = source.fetch("chr1", 0, 30)

    assert source.height_in == pytest.approx(0.65)
    assert track.height_in == pytest.approx(0.65)


@pytest.mark.parametrize("height", ["0", "-1", "large"])
def test_custom_track_rejects_invalid_height(tmp_path, height):
    track_file = tmp_path / "regions.data"
    track_file.write_text("chr1\t10\t20\tregion\n", encoding="utf-8")

    with pytest.raises(ValueError, match="HEIGHT_IN"):
        build_custom_annotation_sources([[
            str(track_file), "bed", "Regions", "#000000", "pack", height,
        ]])


def test_narrowpeak_preserves_signal_and_summit(tmp_path):
    peak_file = tmp_path / "H3K27ac.narrowPeak"
    peak_file.write_text(
        "chr1\t100\t180\tpeak-1\t500\t.\t12.5\t8\t7\t30\n"
        "chr1\t190\t240\tpeak-2\t250\t.\t-1\t4\t3\t-1\n",
        encoding="utf-8",
    )

    track = AnnotationSource(str(peak_file), display_mode="collapse").fetch(
        "chr1", 90, 250
    )

    assert track.kind == "narrowpeak"
    assert track.items[0].value == pytest.approx(12.5)
    assert track.items[0].summit == 130
    assert track.items[1].value == pytest.approx(250)
    assert track.items[1].summit is None


def test_broadpeak_and_signal_tracks_use_quantitative_rows(tmp_path):
    broad_file = tmp_path / "H3K27me3.broadPeak"
    signal_file = tmp_path / "DNase.signal"
    broad_file.write_text(
        "chr1\t100\t220\tbroad-1\t400\t.\t18.0\t6\t5\n",
        encoding="utf-8",
    )
    signal_file.write_text("chr1\t100\t120\t3.5\n", encoding="utf-8")

    broad = AnnotationSource(str(broad_file)).fetch("chr1", 90, 230)
    signal = AnnotationSource(str(signal_file)).fetch("chr1", 90, 130)

    assert broad.kind == "broadpeak"
    assert broad.items[0].summit is None
    assert broad.items[0].value == pytest.approx(18.0)
    assert signal.kind == "signal"
    assert signal.items[0].value == pytest.approx(3.5)


def test_epigenomic_signal_rejects_negative_cnv_values(tmp_path):
    signal_file = tmp_path / "signed.signal"
    signal_file.write_text("chr1\t100\t120\t-0.5\n", encoding="utf-8")

    with pytest.raises(ValueError, match="require non-negative"):
        AnnotationSource(str(signal_file)).fetch("chr1", 90, 130)


def test_peak_density_display_is_accepted(tmp_path):
    peak_file = tmp_path / "DNase.narrowPeak"
    peak_file.write_text(
        "chr1\t100\t150\ta\t100\t.\t5\t3\t2\t20\n",
        encoding="utf-8",
    )

    source = AnnotationSource(str(peak_file), display_mode="density")
    track = source.fetch("chr1", 90, 160)

    assert source.display_mode == "density"
    assert track.display_mode == "density"
    assert len(track.rows) == 1


def test_bgzip_narrowpeak_is_fetched_with_tabix(tmp_path):
    plain = tmp_path / "DNase.narrowPeak"
    compressed = tmp_path / "DNase.narrowPeak.gz"
    plain.write_text(
        "chr1\t100\t150\tinside\t100\t.\t5\t3\t2\t20\n"
        "chr1\t300\t350\toutside\t100\t.\t6\t3\t2\t25\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="bed", force=True)

    source = AnnotationSource(str(compressed))
    track = source.fetch("chr1", 90, 160)

    assert source.compressed
    assert track.kind == "narrowpeak"
    assert track.items[0].name == "inside"
    assert track.items[0].summit == 120


def test_custom_track_accepts_one_comma_separated_argument(tmp_path):
    track_file = tmp_path / "regions.data"
    track_file.write_text("chr1\t10\t20\tregion\n", encoding="utf-8")
    specification = f'{track_file},bed,"Candidate, somatic regions",#123456,collapse'

    source = build_custom_annotation_sources([[specification]])[0]
    track = source.fetch("chr1", 0, 30)

    assert source.display_mode == "collapse"
    assert track.label == "Candidate, somatic regions"
    assert track.color == "#123456"


@pytest.mark.parametrize("color", ["12,34,56", "rgb(12,34,56)"])
def test_comma_separated_custom_track_preserves_rgb_colour(tmp_path, color):
    track_file = tmp_path / "regions.data"
    track_file.write_text("chr1\t10\t20\tregion\n", encoding="utf-8")
    specification = f"{track_file},bed,Regions,{color},expand"

    source = build_custom_annotation_sources([[specification]])[0]

    assert source.color == "#0c2238"
    assert source.display_mode == "expand"


def test_comma_separated_custom_track_rejects_missing_fields(tmp_path):
    with pytest.raises(ValueError, match="requires one CSV value"):
        build_custom_annotation_sources([[f"{tmp_path}/regions.bed,bed,Regions"]])


def test_gene_track_collapse_pack_and_expand_transcript_isoforms(tmp_path):
    gtf = tmp_path / "isoforms.gtf"
    gtf.write_text(
        'chr1\ttest\ttranscript\t101\t260\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; gene_name "GENE1";\n'
        'chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; gene_name "GENE1";\n'
        'chr1\ttest\texon\t201\t260\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; gene_name "GENE1";\n'
        'chr1\ttest\ttranscript\t101\t260\t.\t+\t.\tgene_id "g1"; transcript_id "tx2"; gene_name "GENE1";\n'
        'chr1\ttest\texon\t101\t130\t.\t+\t.\tgene_id "g1"; transcript_id "tx2"; gene_name "GENE1";\n'
        'chr1\ttest\texon\t181\t260\t.\t+\t.\tgene_id "g1"; transcript_id "tx2"; gene_name "GENE1";\n',
        encoding="utf-8",
    )

    collapsed = AnnotationSource(str(gtf), display_mode="collapse").fetch("chr1", 90, 270)
    packed = AnnotationSource(str(gtf), display_mode="pack").fetch("chr1", 90, 270)
    expanded = AnnotationSource(str(gtf), display_mode="expand").fetch("chr1", 90, 270)

    assert len(collapsed.items) == 2  # source transcript models remain available
    assert len(collapsed.rows) == 1
    assert len(collapsed.rows[0]) == 1
    assert collapsed.rows[0][0].name == "GENE1"
    assert collapsed.rows[0][0].blocks == [(100, 150), (180, 260)]
    assert len(packed.rows) == 2
    assert len(expanded.rows) == 2
    assert all(len(row) == 1 for row in expanded.rows)
    assert [row[0].transcript_label for row in expanded.rows] == ["tx1", "tx2"]


def test_ncbi_gff_cds_parent_merges_with_refseq_transcript(tmp_path):
    gff = tmp_path / "refseq.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tBestRefSeq\tgene\t101\t250\t.\t+\t.\tID=gene-G1;Name=GENE1;gene=GENE1\n"
        "chr1\tBestRefSeq\tmRNA\t101\t250\t.\t+\t.\tID=rna-NM_1.2;Parent=gene-G1;"
        "Name=NM_1.2;gene=GENE1;transcript_id=NM_1.2;tag=RefSeq Select\n"
        "chr1\tBestRefSeq\texon\t101\t250\t.\t+\t.\tParent=rna-NM_1.2;"
        "gene=GENE1;transcript_id=NM_1.2\n"
        "chr1\tBestRefSeq\tCDS\t151\t220\t.\t+\t0\tID=cds-NP_1.1;"
        "Parent=rna-NM_1.2;Name=NP_1.1;gene=GENE1;protein_id=NP_1.1\n",
        encoding="utf-8",
    )

    track = AnnotationSource(str(gff), display_mode="expand").fetch("chr1", 90, 260)

    assert len(track.items) == 1
    assert track.items[0].transcript_label == "NM_1.2"
    assert track.items[0].group_label == "GENE1"
    assert track.items[0].blocks == [(150, 220)]
    assert track.items[0].primary_label == "RefSeq Select"


def test_primary_isoform_selection_uses_marker_priority_and_gene_fallback(tmp_path):
    gtf = tmp_path / "primary.gtf"
    records = [
        (101, 150, 'gene_id "g1"; transcript_id "tx_mane"; gene_name "G1"; tag "MANE_Select";'),
        (101, 180, 'gene_id "g1"; transcript_id "tx_canonical"; gene_name "G1"; tag "Ensembl_canonical";'),
        (101, 220, 'gene_id "g1"; transcript_id "tx_other"; gene_name "G1";'),
        (301, 350, 'gene_id "g2"; transcript_id "tx_unmarked_a"; gene_name "G2";'),
        (301, 390, 'gene_id "g2"; transcript_id "tx_unmarked_b"; gene_name "G2";'),
        (501, 560, 'gene_id "g3"; transcript_id "tx_appris"; gene_name "G3"; tag "appris_principal_1";'),
    ]
    gtf.write_text("".join(
        f"chr1\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attributes}\n"
        f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attributes}\n"
        for start, end, attributes in records
    ), encoding="utf-8")

    all_track = AnnotationSource(
        str(gtf), display_mode="expand", primary_isoforms="all"
    ).fetch("chr1", 90, 600)
    preferred = AnnotationSource(
        str(gtf), display_mode="expand", primary_isoforms="prefer"
    ).fetch("chr1", 90, 600)
    strict = AnnotationSource(
        str(gtf), display_mode="expand", primary_isoforms="only"
    ).fetch("chr1", 90, 600)

    assert {item.transcript_label: item.primary_rank for item in all_track.items} == {
        "tx_mane": 1,
        "tx_canonical": 2,
        "tx_other": None,
        "tx_unmarked_a": None,
        "tx_unmarked_b": None,
        "tx_appris": 3,
    }
    assert [item.transcript_label for item in preferred.items] == [
        "tx_mane", "tx_unmarked_a", "tx_unmarked_b", "tx_appris",
    ]
    assert [item.transcript_label for item in strict.items] == [
        "tx_mane", "tx_appris",
    ]


def test_primary_isoform_preference_is_applied_before_gene_collapse(tmp_path):
    gtf = tmp_path / "collapsed-primary.gtf"
    gtf.write_text(
        'chr1\ttest\ttranscript\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "primary"; tag "MANE_Select";\n'
        'chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "primary"; tag "MANE_Select";\n'
        'chr1\ttest\ttranscript\t101\t250\t.\t+\t.\tgene_id "g1"; transcript_id "long";\n'
        'chr1\ttest\texon\t101\t250\t.\t+\t.\tgene_id "g1"; transcript_id "long";\n',
        encoding="utf-8",
    )

    track = AnnotationSource(
        str(gtf), display_mode="collapse", primary_isoforms="prefer"
    ).fetch("chr1", 90, 270)

    assert len(track.items) == 1
    assert track.rows[0][0].blocks == [(100, 150)]
    assert track.rows[0][0].primary_label == "MANE Select"


def test_invalid_gene_track_display_mode_is_rejected(tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="display mode"):
        AnnotationSource(str(bed), display_mode="stack-everything")
    with pytest.raises(ValueError, match="primary-isoform mode"):
        AnnotationSource(str(bed), primary_isoforms="guess")


def test_compressed_custom_track_fetches_with_tabix(tmp_path):
    plain = tmp_path / "regions.data"
    compressed = tmp_path / "regions.data.gz"
    plain.write_text(
        "chr1\t10\t20\tinside\nchr1\t100\t120\toutside\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="bed", force=True)

    source = build_custom_annotation_sources([
        [str(compressed), "bed", "Indexed regions", "#000000"],
    ])[0]
    track = source.fetch("chr1", 0, 30)

    assert source.compressed
    assert [item.name for item in track.items] == ["inside"]


def test_track_colour_accepts_hex_and_rgb_but_rejects_invalid_values():
    assert normalize_track_color("#17217a") == "#17217a"
    assert normalize_track_color("23,33,122") == "#17217a"
    assert normalize_track_color("rgb(23, 33, 122)") == "#17217a"
    with pytest.raises(ValueError, match="between 0 and 255"):
        normalize_track_color("256,0,0")
    with pytest.raises(ValueError, match="Invalid track colour"):
        normalize_track_color("not-a-colour")


def test_gff3_transcript_id_is_used_instead_of_gene_parent(tmp_path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t260\t.\t+\t.\tID=g1;Name=GENE1\n"
        "chr1\ttest\tmRNA\t101\t260\t.\t+\t.\tID=tx1;Parent=g1;Name=TX1\n"
        "chr1\ttest\texon\t101\t150\t.\t+\t.\tParent=tx1\n"
        "chr1\ttest\texon\t201\t260\t.\t+\t.\tParent=tx1\n",
        encoding="utf-8",
    )
    track = AnnotationSource(str(gff)).fetch("chr1", 90, 270)
    assert len(track.items) == 1
    assert track.items[0].name == "TX1"
    assert track.items[0].blocks == [(100, 150), (200, 260)]


def test_compressed_track_requires_tabix_index(tmp_path):
    compressed = tmp_path / "regions.bed.gz"
    compressed.write_bytes(b"not-bgzip")
    with pytest.raises(ValueError, match="requires a tabix index"):
        AnnotationSource(str(compressed))


def test_bgzip_bed_is_fetched_through_tabix_index(tmp_path):
    plain = tmp_path / "regions.bed"
    compressed = tmp_path / "regions.bed.gz"
    plain.write_text("chr1\t100\t120\ta\nchr2\t50\t80\tb\n", encoding="utf-8")
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="bed", force=True)

    track = AnnotationSource(str(compressed)).fetch("chr1", 90, 130)

    assert [item.name for item in track.items] == ["a"]
    assert os.path.isfile(str(compressed) + ".tbi")


def test_bgzip_gff3_is_grouped_after_tabix_fetch(tmp_path):
    plain = tmp_path / "genes.gff3"
    compressed = tmp_path / "genes.gff3.gz"
    plain.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tmRNA\t101\t260\t.\t+\t.\tID=tx1;Name=TX1\n"
        "chr1\ttest\texon\t101\t150\t.\t+\t.\tParent=tx1\n"
        "chr1\ttest\texon\t201\t260\t.\t+\t.\tParent=tx1\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="gff", force=True)

    track = AnnotationSource(str(compressed)).fetch("chr1", 120, 220)

    assert len(track.items) == 1
    assert track.items[0].name == "TX1"
    assert track.items[0].blocks == [(100, 150), (200, 260)]


def test_plain_vcf_uses_zero_based_spans_ids_and_info_end(tmp_path):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\trs1\tA\tG\t.\tPASS\t.\n"
        "chr1\t121\t.\tATC\tA\t.\tPASS\t.\n"
        "chr1\t151\tsv1\tN\t<DEL>\t.\tPASS\tEND=180\n"
        "chr1\t301\toutside\tC\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )

    track = AnnotationSource(str(vcf), "Variants").fetch("chr1", 90, 200)

    assert [item.name for item in track.items] == ["rs1", "ATC>A", "sv1"]
    assert [(item.start, item.end) for item in track.items] == [
        (100, 101), (120, 123), (150, 180),
    ]


def test_bgzip_vcf_is_fetched_through_tabix_index(tmp_path):
    plain = tmp_path / "variants.vcf"
    compressed = tmp_path / "variants.vcf.gz"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\tinside\tA\tG\t.\tPASS\t.\n"
        "chr1\t301\toutside\tC\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)

    track = AnnotationSource(str(compressed)).fetch("chr1", 90, 130)

    assert track.kind == "vcf"
    assert [item.name for item in track.items] == ["inside"]
    assert os.path.isfile(str(compressed) + ".tbi")


def test_seg_parses_header_samples_values_and_one_based_coordinates(tmp_path):
    seg = tmp_path / "tumour.seg"
    seg.write_text(
        "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
        "Tumour\tchr1\t101\t150\t12\t-0.65\n"
        "Tumour\tchr1\t151\t220\t18\t0.42\n"
        "Normal\tchr1\t221\t260\t9\t0.03\n"
        "Tumour\tchr1\t261\t280\t4\tNA\n"
        "Tumour\tchr2\t101\t150\t8\t1.2\n",
        encoding="utf-8",
    )

    track = AnnotationSource(str(seg), "Copy number").fetch("chr1", 90, 240)

    assert track.kind == "seg"
    assert track.color == "#555555"
    assert track.color_by_sign
    assert len(track.rows) == 1
    assert [(item.start, item.end, item.value, item.sample) for item in track.items] == [
        (100, 150, -0.65, "Tumour"),
        (150, 220, 0.42, "Tumour"),
        (220, 260, 0.03, "Normal"),
    ]
    custom = AnnotationSource(str(seg), color="#663399").fetch("chr1", 90, 240)
    assert custom.color == "#663399"
    assert not custom.color_by_sign


def test_bedgraph_and_log2_use_zero_based_intervals(tmp_path):
    ratios = tmp_path / "bins.log2"
    ratios.write_text(
        "chr1\t100\t120\t-0.3\n"
        "chr1\t120\t140\t0.8\n"
        "chr1\t140\t160\tnan\n"
        "chr2\t100\t120\t1.0\n",
        encoding="utf-8",
    )

    track = AnnotationSource(str(ratios)).fetch("chr1", 110, 130)

    assert track.kind == "log2"
    assert [(item.start, item.end, item.value) for item in track.items] == [
        (100, 120, -0.3), (120, 140, 0.8),
    ]


def test_bgzip_seg_is_fetched_through_custom_tabix_columns(tmp_path):
    plain = tmp_path / "tumour.seg"
    compressed = tmp_path / "tumour.seg.gz"
    plain.write_text(
        "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
        "Tumour\tchr1\t101\t150\t12\t-0.65\n"
        "Tumour\tchr1\t301\t350\t10\t0.7\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(
        str(compressed), seq_col=1, start_col=2, end_col=3,
        line_skip=1, zerobased=False, force=True,
    )

    track = AnnotationSource(str(compressed)).fetch("chr1", 90, 180)

    assert track.kind == "seg"
    assert [(item.start, item.end, item.value) for item in track.items] == [
        (100, 150, -0.65),
    ]


def test_baf_vcf_selects_sample_heterozygous_snvs_and_uses_ad(tmp_path):
    vcf = tmp_path / "genotypes.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNormal\tTumour\n"
        "chr1\t101\trs1\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:16,4\n"
        "chr1\t121\thom\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:9,11\t1/1:0,20\n"
        "chr1\t141\tindel\tA\tAT\t.\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:10,10\n"
        "chr1\t401\tout\tG\tC\t.\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:10,10\n",
        encoding="utf-8",
    )

    track = BafSource(str(vcf), sample="Tumour").fetch("chr1", 90, 200)

    assert track.kind == "baf"
    assert track.label == "genotypes BAF"
    assert [(item.start, item.name, item.value, item.sample) for item in track.items] == [
        (100, "rs1", 0.2, "Tumour"),
    ]


def test_bgzip_baf_vcf_uses_tabix_and_builder_pairs_options(tmp_path):
    plain = tmp_path / "genotypes.vcf"
    compressed = tmp_path / "genotypes.vcf.gz"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumour\n"
        "chr1\t101\trs1\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:7,13\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)

    sources = build_baf_sources(
        [str(compressed)], labels=["Tumour BAF"], samples=["Tumour"]
    )
    track = sources[0].fetch("chr1", 90, 120)

    assert track.label == "Tumour BAF"
    assert track.items[0].value == 0.65
    with pytest.raises(ValueError, match="not present"):
        BafSource(str(compressed), sample="Missing")


def test_track_labels_follow_track_order(tmp_path):
    first = tmp_path / "first.bed"
    second = tmp_path / "second.gtf"
    first.write_text("chr1\t0\t1\n", encoding="utf-8")
    second.write_text("", encoding="utf-8")
    sources = build_annotation_sources([str(first), str(second)], ["First"])
    assert [source.label for source in sources] == ["First", "second"]
    with pytest.raises(ValueError, match="More --track_label"):
        build_annotation_sources([str(first)], ["one", "two"])


def test_annotation_track_renders_in_bam_comparison(tmp_path):
    bed = tmp_path / "regions.bed"
    seg = tmp_path / "tumour.seg"
    bed.write_text(
        "chr9\t101867490\t101867540\tLEFT\n"
        "chr9\t101867550\t101867600\tRIGHT\n",
        encoding="utf-8",
    )
    seg.write_text(
        "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
        "Tumour\tchr9\t101867481\t101867550\t20\t-0.5\n"
        "Tumour\tchr9\t101867551\t101867620\t22\t0.7\n",
        encoding="utf-8",
    )
    output, summary = compare_snapshots(
        TEST_BAM, TEST_BAM, "chr9", 101867480, 101867620,
        output_dir=str(tmp_path), output_name="compare-tracks.png",
        display_mode="collapse", haplotype_view="split",
        show_coverage=False, dpi=40,
        annotation_sources=[
            AnnotationSource(str(bed), "Regions"),
            AnnotationSource(str(seg), "Copy number"),
        ],
    )
    assert os.path.isfile(output)
    assert os.path.getsize(output) > 0
    assert "metric" in summary
