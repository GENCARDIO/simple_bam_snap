"""Provision indexed UCSC NCBI RefSeq GTF tracks for human assemblies."""
from __future__ import annotations

import gzip
import heapq
import logging
import os
from pathlib import Path
import shutil
import tempfile
from typing import Mapping, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pysam

from src.cytobands import resolve_cytobands


log = logging.getLogger("simple_bam_snap.refseq")

DEFAULT_REFSEQ_DIR = Path(__file__).resolve().parent.parent / "data" / "ucsc" / "refseq"
REFSEQ_SOURCES = {
    "hg19": {
        "annotation": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/"
            "105.20220307/GCF_000001405.25_GRCh37.p13/"
            "GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
        ),
        "assembly_report": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/"
            "105.20220307/GCF_000001405.25_GRCh37.p13/"
            "GCF_000001405.25_GRCh37.p13_assembly_report.txt"
        ),
    },
    "hg38": {
        "annotation": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/"
            "GCF_000001405.40_GRCh38.p14/"
            "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
        ),
        "assembly_report": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/"
            "GCF_000001405.40_GRCh38.p14/"
            "GCF_000001405.40_GRCh38.p14_assembly_report.txt"
        ),
    },
}
REFSEQ_LABELS = {
    "hg19": "NCBI RefSeq (hg19/GRCh37)",
    "hg38": "NCBI RefSeq (hg38/GRCh38)",
}
ASSEMBLY_ALIASES = {
    "hg19": "hg19", "grch37": "hg19",
    "hg38": "hg38", "grch38": "hg38",
}


def normalize_assembly(value: str) -> str:
    normalized = value.lower()
    if normalized not in ASSEMBLY_ALIASES:
        choices = ", ".join(sorted(ASSEMBLY_ALIASES))
        raise ValueError(f"Unknown RefSeq assembly {value!r}; choose {choices}.")
    return ASSEMBLY_ALIASES[normalized]


def detect_human_assembly(contig_lengths: Mapping[str, int]) -> Optional[str]:
    """Identify hg19 or hg38 using exact chromosome lengths."""
    return resolve_cytobands(contig_lengths, genome="auto")[1]


def validate_refseq(path: str | Path) -> None:
    """Confirm that a cached GTF is BGZF-compressed and tabix-readable."""
    source = Path(path)
    index_path = Path(f"{source}.tbi")
    if not source.is_file() or not index_path.is_file():
        raise ValueError(f"RefSeq cache is incomplete: {source}")
    try:
        with pysam.TabixFile(str(source)) as tabix:
            if not tuple(tabix.contigs):
                raise ValueError(f"Indexed RefSeq track contains no contigs: {source}")
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read indexed RefSeq track {source}: {exc}") from exc


def write_gff_chunk(records: list, path: Path) -> None:
    """Write one coordinate-sorted temporary chunk with mergeable sort keys."""
    records.sort(key=lambda item: (item[0], item[1], item[2]))
    with path.open("wt", encoding="utf-8") as handle:
        for chrom, start, order, line in records:
            handle.write(f"{chrom}\t{start:012d}\t{order:012d}\t{line}")


def ensure_refseq(
    assembly: str, cache_dir: Optional[str | Path] = None,
    source_url: Optional[str] = None,
) -> Path:
    """Return a local BGZF/tabix RefSeq GTF, downloading it once if absent."""
    selected = normalize_assembly(assembly)
    directory = Path(cache_dir) if cache_dir else DEFAULT_REFSEQ_DIR
    output_path = directory / f"{selected}.ncbiRefSeq.gff.gz"
    index_path = Path(f"{output_path}.tbi")
    if output_path.is_file() and index_path.is_file():
        validate_refseq(output_path)
        return output_path

    try:
        directory.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"Cannot create RefSeq cache directory {directory}: {exc}") from exc

    annotation_url = source_url or REFSEQ_SOURCES[selected]["annotation"]
    report_url = None if source_url else REFSEQ_SOURCES[selected]["assembly_report"]
    log.info("Downloading %s from NCBI", REFSEQ_LABELS[selected])
    try:
        with tempfile.TemporaryDirectory(prefix=f"{selected}-", dir=directory) as temp_name:
            temp_dir = Path(temp_name)
            archive_path = temp_dir / "source.gff.gz"
            report_path = temp_dir / "assembly_report.txt"
            header_path = temp_dir / "headers.gff"
            plain_path = temp_dir / "source.gff"
            bgzf_path = temp_dir / output_path.name
            try:
                request = Request(
                    annotation_url, headers={"User-Agent": "simple-bam-snap/1"}
                )
                with urlopen(request, timeout=180) as response, archive_path.open("wb") as target:
                    shutil.copyfileobj(response, target, length=1024 * 1024)
                if report_url:
                    request = Request(
                        report_url, headers={"User-Agent": "simple-bam-snap/1"}
                    )
                    with urlopen(request, timeout=60) as response, report_path.open("wb") as target:
                        shutil.copyfileobj(response, target, length=1024 * 1024)
            except (HTTPError, URLError, TimeoutError, OSError) as exc:
                raise ValueError(
                    f"Could not download {REFSEQ_LABELS[selected]} from NCBI: {exc}. "
                    "Retry with network access or use --refseq none."
                ) from exc

            aliases = {}
            if report_path.is_file():
                with report_path.open("rt", encoding="utf-8") as report:
                    for line in report:
                        if line.startswith("#"):
                            continue
                        fields = line.rstrip().split("\t")
                        if len(fields) < 10 or fields[9] == "na":
                            continue
                        ucsc_name = fields[9]
                        for alias_index in (0, 4, 6):
                            alias = fields[alias_index]
                            if alias != "na":
                                aliases[alias] = ucsc_name

            try:
                chunks = []
                records = []
                record_order = 0
                with gzip.open(archive_path, "rt", encoding="utf-8") as compressed, header_path.open(
                    "wt", encoding="utf-8"
                ) as headers:
                    for line in compressed:
                        if line.startswith("##sequence-region "):
                            parts = line.rstrip().split(" ")
                            if len(parts) >= 2:
                                parts[1] = aliases.get(parts[1], parts[1])
                            headers.write(" ".join(parts) + "\n")
                        elif line.startswith("#"):
                            if line.strip() != "###":
                                headers.write(line)
                        else:
                            chrom, separator, remainder = line.partition("\t")
                            normalized_chrom = aliases.get(chrom, chrom)
                            normalized_line = normalized_chrom + separator + remainder
                            fields = normalized_line.split("\t", 4)
                            if len(fields) < 4:
                                continue
                            records.append(
                                (normalized_chrom, int(fields[3]), record_order, normalized_line)
                            )
                            record_order += 1
                            if len(records) >= 50_000:
                                chunk_path = temp_dir / f"chunk-{len(chunks):05d}.gff"
                                write_gff_chunk(records, chunk_path)
                                chunks.append(chunk_path)
                                records = []
                if records:
                    chunk_path = temp_dir / f"chunk-{len(chunks):05d}.gff"
                    write_gff_chunk(records, chunk_path)
                    chunks.append(chunk_path)

                handles = []
                try:
                    with plain_path.open("wt", encoding="utf-8") as plain, header_path.open(
                        "rt", encoding="utf-8"
                    ) as headers:
                        shutil.copyfileobj(headers, plain, length=1024 * 1024)
                        for chunk_path in chunks:
                            handles.append(chunk_path.open("rt", encoding="utf-8"))
                        for keyed_line in heapq.merge(*handles):
                            plain.write(keyed_line.split("\t", 3)[3])
                finally:
                    for handle in handles:
                        handle.close()
                pysam.tabix_compress(str(plain_path), str(bgzf_path), force=True)
                pysam.tabix_index(str(bgzf_path), preset="gff", force=True)
            except (EOFError, OSError, ValueError) as exc:
                raise ValueError(
                    f"Downloaded RefSeq GFF for {selected} could not be indexed: {exc}"
                ) from exc

            os.replace(bgzf_path, output_path)
            os.replace(Path(f"{bgzf_path}.tbi"), index_path)
    except OSError as exc:
        raise ValueError(f"Cannot prepare RefSeq cache in {directory}: {exc}") from exc

    validate_refseq(output_path)
    log.info("Cached indexed RefSeq track: %s", output_path)
    return output_path
