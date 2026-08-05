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


DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "ucsc"
BUNDLED_CYTOBANDS = {
    "hg19": DATA_DIR / "hg19.cytoBand.txt.gz",
    "hg38": DATA_DIR / "hg38.cytoBand.txt.gz",
}


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
    for bands in by_chrom.values():
        bands.sort(key=lambda band: (band.start, band.end, band.name))
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
    if custom_path:
        return load_cytobands(custom_path), Path(custom_path).name
    if genome == "none":
        return {}, None
    if genome in BUNDLED_CYTOBANDS:
        return load_cytobands(BUNDLED_CYTOBANDS[genome]), genome
    if genome != "auto":
        raise ValueError(f"Unknown genome '{genome}'. Choose auto, hg19, hg38, or none.")

    candidates = []
    for assembly, path in BUNDLED_CYTOBANDS.items():
        bands = load_cytobands(path)
        match_count = sum(
            bool(chrom_bands) and max(band.end for band in chrom_bands) == length
            for chrom, length in contig_lengths.items()
            for chrom_bands in [bands_for_chrom(bands, chrom)]
        )
        candidates.append((match_count, assembly, bands))
    candidates.sort(key=lambda item: item[0], reverse=True)
    best_score = candidates[0][0]
    if best_score == 0 or (len(candidates) > 1 and candidates[1][0] == best_score):
        return {}, None
    assembly, bands = candidates[0][1:]
    return bands, assembly
