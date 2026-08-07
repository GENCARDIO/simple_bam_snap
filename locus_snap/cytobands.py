"""Load UCSC cytoband tables and select bundled human assemblies safely."""
from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple


@dataclass(frozen=True)
class Cytoband:
    chrom: str
    start: int
    end: int
    name: str
    stain: str


DATA_DIR = Path(__file__).resolve().parent / "data" / "ucsc"
BUNDLED_CYTOBANDS = {
    "hg19": DATA_DIR / "hg19.cytoBand.txt.gz",
    "hg38": DATA_DIR / "hg38.cytoBand.txt.gz",
}
GENOME_ALIASES = {"grch37": "hg19", "grch38": "hg38"}


def load_cytobands(path: str | Path) -> Dict[str, List[Cytoband]]:
    """Read UCSC's five-column cytoBand text format."""
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"Cytoband file does not exist: {source}")

    by_chrom: Dict[str, List[Cytoband]] = {}
    handle_context = (
        gzip.open(source, "rt", encoding="utf-8")
        if source.suffix.lower() in {".gz", ".bgz", ".bgzf"}
        else source.open("rt", encoding="utf-8")
    )
    with handle_context as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"{source}:{line_number}: expected UCSC cytoBand columns "
                    "chrom, start, end, name, gieStain"
                )
            chrom, start_text, end_text, name, stain = fields[:5]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(
                    f"{source}:{line_number}: cytoband coordinates must be integers"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{source}:{line_number}: invalid cytoband interval {start}-{end}"
                )
            by_chrom.setdefault(chrom, []).append(
                Cytoband(chrom, start, end, name, stain)
            )

    if not by_chrom:
        raise ValueError(f"Cytoband file contains no records: {source}")

    def band_sort_key(band: Cytoband):
        return (band.start, band.end, band.name)

    for bands in by_chrom.values():
        bands.sort(key=band_sort_key)
    return by_chrom


def bands_for_chrom(
    cytobands: Mapping[str, List[Cytoband]], chrom: str
) -> List[Cytoband]:
    """Resolve the common ``chr9``/``9`` naming difference."""
    if chrom in cytobands:
        return cytobands[chrom]
    alternate = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    return cytobands.get(alternate, [])


def resolve_cytobands(
    contig_lengths: Mapping[str, int],
    genome: str = "auto",
    custom_path: Optional[str] = None,
) -> Tuple[Dict[str, List[Cytoband]], Optional[str]]:
    """Resolve custom, explicit, or safely auto-detected cytoband data.

    Auto-detection requires at least one exact chromosome-length match and a
    unique best bundled assembly. It therefore never guesses from chromosome
    names alone.
    """
    genome = GENOME_ALIASES.get(genome.lower(), genome.lower())
    if custom_path:
        return load_cytobands(custom_path), Path(custom_path).name
    if genome == "none":
        return {}, None
    if genome in BUNDLED_CYTOBANDS:
        return load_cytobands(BUNDLED_CYTOBANDS[genome]), genome
    if genome != "auto":
        raise ValueError(
            f"Unknown genome '{genome}'. Choose auto, hg19/GRCh37, hg38/GRCh38, or none."
        )

    candidates = []
    for assembly, path in BUNDLED_CYTOBANDS.items():
        bands = load_cytobands(path)
        match_count = 0
        for chrom, length in contig_lengths.items():
            chrom_bands = bands_for_chrom(bands, chrom)
            if not chrom_bands:
                continue
            max_end = 0
            for band in chrom_bands:
                if band.end > max_end:
                    max_end = band.end
            if max_end == length:
                match_count += 1
        candidates.append((match_count, assembly, bands))

    def candidate_score(item) -> int:
        return item[0]

    candidates.sort(key=candidate_score, reverse=True)
    best_score = candidates[0][0]
    if best_score == 0 or (len(candidates) > 1 and candidates[1][0] == best_score):
        return {}, None
    assembly, bands = candidates[0][1:]
    return bands, assembly
