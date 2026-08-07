import gzip
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.refseq import detect_human_assembly, ensure_refseq, normalize_assembly


def test_detects_supported_human_assemblies_from_exact_contig_lengths():
    assert detect_human_assembly({"chr9": 141_213_431}) == "hg19"
    assert detect_human_assembly({"9": 138_394_717}) == "hg38"
    assert detect_human_assembly({"chrDemo": 1_000}) is None


def test_grch_names_normalize_to_ucsc_assemblies():
    assert normalize_assembly("GRCh37") == "hg19"
    assert normalize_assembly("GRCh38") == "hg38"


def test_refseq_download_is_bgzipped_indexed_and_reused(tmp_path):
    source = tmp_path / "remote.gtf.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            'chr1\tRefSeq\texon\t201\t250\t.\t+\t.\tgene_id "G1"; '
            'transcript_id "NM_1"; gene_name "GENE1";\n'
            'chr1\tRefSeq\texon\t101\t150\t.\t+\t.\tgene_id "G2"; '
            'transcript_id "NM_2"; gene_name "GENE2";\n'
        )

    cache = tmp_path / "cache"
    first = ensure_refseq("grch37", cache, source.as_uri())
    source.unlink()
    second = ensure_refseq("hg19", cache, "file:///does/not/exist")

    assert first == second
    assert first.name == "hg19.ncbiRefSeq.gff.gz"
    assert os.path.isfile(f"{first}.tbi")
    with pysam.TabixFile(str(first)) as tabix:
        records = list(tabix.fetch("chr1", 100, 250))
    assert len(records) == 2
    assert records[0].split("\t")[3] == "101"


def test_unknown_refseq_assembly_is_rejected():
    with pytest.raises(ValueError, match="Unknown RefSeq assembly"):
        normalize_assembly("hg18")
