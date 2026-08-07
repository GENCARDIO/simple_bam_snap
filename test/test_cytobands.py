import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.cytobands import bands_for_chrom, load_cytobands, resolve_cytobands


def test_loads_ucsc_cytoband_format_and_resolves_chr_alias(tmp_path):
    path = tmp_path / "bands.txt"
    path.write_text(
        "chr1\t0\t40\tp11\tgneg\n"
        "chr1\t40\t50\tp10\tacen\n"
        "chr1\t50\t100\tq11\tgpos100\n",
        encoding="utf-8",
    )
    bands = load_cytobands(path)

    assert [band.stain for band in bands_for_chrom(bands, "1")] == [
        "gneg", "acen", "gpos100"
    ]


def test_auto_detects_hg19_from_exact_contig_length():
    bands, assembly = resolve_cytobands({"chr9": 141_213_431})

    assert assembly == "hg19"


def test_grch_aliases_select_corresponding_ucsc_cytobands():
    grch37, assembly37 = resolve_cytobands({}, genome="grch37")
    grch38, assembly38 = resolve_cytobands({}, genome="grch38")

    assert assembly37 == "hg19"
    assert assembly38 == "hg38"
    assert bands_for_chrom(grch37, "chr9")
    assert bands_for_chrom(grch38, "chr9")


def test_auto_detection_does_not_guess_unknown_assembly():
    bands, assembly = resolve_cytobands({"chrDemo": 1_000})

    assert bands == {}
    assert assembly is None
