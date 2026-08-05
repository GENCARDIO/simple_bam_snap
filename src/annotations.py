"""Read interval, gene, variant, and copy-number annotation tracks."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pysam
from matplotlib.colors import is_color_like, to_hex

from src.config import DEFAULT_TRACK_COLORS


# Track-type defaults: simple regions are black; transcript structures use a
# deep UCSC-like navy.
BED_ANNOTATION_COLOR = DEFAULT_TRACK_COLORS["bed"]
GENE_ANNOTATION_COLOR = DEFAULT_TRACK_COLORS["gene"]
VCF_ANNOTATION_COLOR = DEFAULT_TRACK_COLORS["vcf"]
CNV_ANNOTATION_COLOR = DEFAULT_TRACK_COLORS["cnv"]
BAF_ANNOTATION_COLOR = DEFAULT_TRACK_COLORS["baf"]
ANNOTATION_COLOR = GENE_ANNOTATION_COLOR  # backwards-compatible public name
SUPPORTED_TRACK_FORMATS = (
    "bed", "gff", "gff3", "gtf", "vcf", "seg", "bedgraph", "bdg", "log2", "cnv",
)
CNV_TRACK_FORMATS = ("seg", "bedgraph", "bdg", "log2", "cnv")
BAF_TRACK_FORMATS = ("baf",)
ANNOTATION_DISPLAY_MODES = ("collapse", "pack", "expand")
PRIMARY_ISOFORM_MODES = ("all", "prefer", "only")
COMPRESSED_SUFFIXES = (".gz", ".bgz", ".bgzf")
TRANSCRIPT_TYPES = {"transcript", "mrna", "ncrna", "trna", "rrna"}
UTR_TYPES = {"utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr"}

Interval = Tuple[int, int]


@dataclass
class AnnotationItem:
    """One interval or transcript model drawn on a single annotation row."""

    start: int
    end: int
    name: str = ""
    strand: str = "."
    blocks: List[Interval] = field(default_factory=list)  # full-height exon/CDS blocks
    utrs: List[Interval] = field(default_factory=list)    # thin UTR blocks
    group: str = ""
    group_label: str = ""
    transcript_label: str = ""
    value: Optional[float] = None
    sample: str = ""
    primary_rank: Optional[int] = None
    primary_label: str = ""


@dataclass
class LoadedAnnotationTrack:
    label: str
    kind: str
    color: str
    items: List[AnnotationItem]
    rows: List[List[AnnotationItem]]
    display_mode: str = "pack"
    color_by_sign: bool = False


def default_label(path: str) -> str:
    name = Path(path).name
    lower = name.lower()
    for suffix in COMPRESSED_SUFFIXES:
        if lower.endswith(suffix):
            name = name[:-len(suffix)]
            lower = lower[:-len(suffix)]
            break
    for suffix in (
        ".bedgraph", ".gff3", ".gff", ".gtf", ".bed", ".vcf",
        ".seg", ".bdg", ".log2", ".cnv",
    ):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return name


def infer_track_format(path: str) -> str:
    name = Path(path).name.lower()
    for suffix in COMPRESSED_SUFFIXES:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    suffix = Path(name).suffix.lstrip(".")
    if suffix == "bdg":
        return "bedgraph"
    if suffix not in SUPPORTED_TRACK_FORMATS:
        raise ValueError(
            f"Cannot infer annotation format for '{path}'. Expected .bed, .gff, .gff3, .gtf, "
            ".vcf, .seg, .bedgraph/.bdg, .log2, or .cnv"
            " (optionally followed by .gz/.bgz/.bgzf)."
        )
    return suffix


def normalize_track_color(value: str) -> str:
    """Validate a track colour and normalize hex or 0-255 RGB to hex."""
    text = value.strip()
    rgb_match = re.fullmatch(
        r"(?:rgb\s*\()?\s*(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})\s*\)?",
        text,
        flags=re.IGNORECASE,
    )
    if rgb_match:
        channels = tuple(int(channel) for channel in rgb_match.groups())
        if any(channel > 255 for channel in channels):
            raise ValueError(f"RGB track colour channels must be between 0 and 255: {value!r}")
        return "#{:02x}{:02x}{:02x}".format(*channels)
    if not is_color_like(text):
        raise ValueError(
            f"Invalid track colour {value!r}; use a hex colour or R,G,B / rgb(R,G,B)."
        )
    return to_hex(text)


def merge_intervals(intervals: Iterable[Interval]) -> List[Interval]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    result = []
    for start, end in merged:
        result.append((start, end))
    return result


def subtract_intervals(intervals: Iterable[Interval], masks: Iterable[Interval]) -> List[Interval]:
    masks = merge_intervals(masks)
    result: List[Interval] = []
    for start, end in merge_intervals(intervals):
        cursor = start
        for mask_start, mask_end in masks:
            if mask_end <= cursor:
                continue
            if mask_start >= end:
                break
            if mask_start > cursor:
                result.append((cursor, min(mask_start, end)))
            cursor = max(cursor, mask_end)
            if cursor >= end:
                break
        if cursor < end:
            result.append((cursor, end))
    return result


def pack_annotation_items(
    items: Sequence[AnnotationItem], padding: int = 1
) -> List[List[AnnotationItem]]:
    """Greedily place non-overlapping annotation models on shared rows."""
    rows: List[List[AnnotationItem]] = []
    row_ends: List[int] = []
    for item in sorted(items, key=lambda value: (value.start, value.end, value.name)):
        for row_index, row_end in enumerate(row_ends):
            if row_end + padding <= item.start:
                rows[row_index].append(item)
                row_ends[row_index] = item.end
                break
        else:
            rows.append([item])
            row_ends.append(item.end)
    return rows


def collapse_annotation_items(items: Sequence[AnnotationItem]) -> List[AnnotationItem]:
    """Merge transcript isoforms into gene-level annotation models."""
    grouped = {}
    for index, item in enumerate(items):
        group = item.group or f"item:{index}"
        grouped.setdefault(group, []).append(item)

    collapsed = []
    for group, members in grouped.items():
        blocks = merge_intervals(
            block for member in members for block in member.blocks
        )
        utrs = subtract_intervals(
            (utr for member in members for utr in member.utrs), blocks
        )
        strands = {member.strand for member in members if member.strand in ("+", "-")}
        collapsed.append(AnnotationItem(
            start=min(member.start for member in members),
            end=max(member.end for member in members),
            name=next(
                (member.group_label for member in members if member.group_label),
                next((member.name for member in members if member.name), group),
            ),
            strand=next(iter(strands)) if len(strands) == 1 else ".",
            blocks=blocks,
            utrs=utrs,
            group=group,
            group_label=next(
                (member.group_label for member in members if member.group_label), ""
            ),
            primary_rank=min(
                (member.primary_rank for member in members if member.primary_rank is not None),
                default=None,
            ),
            primary_label=next(
                (member.primary_label for member in members if member.primary_label), ""
            ),
        ))
    return sorted(collapsed, key=lambda item: (item.start, item.end, item.name))


def parse_gff_attributes(raw: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for field in raw.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        else:
            match = re.match(r"(\S+)\s+[\"']?(.*?)[\"']?$", field)
            if not match:
                continue
            key, value = match.groups()
        attributes[key.strip()] = value.strip().strip('"')
    return attributes


def primary_isoform_annotation(raw: str) -> Tuple[Optional[int], str]:
    """Return marker priority and label for common primary-transcript annotations."""
    normalized = re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")
    if "mane_select" in normalized:
        return 1, "MANE Select"
    if (
        "ensembl_canonical" in normalized
        or "transcript_is_canonical_1" in normalized
        or "transcript_is_canonical_true" in normalized
    ):
        return 2, "Ensembl canonical"
    if "appris_principal" in normalized:
        return 3, "APPRIS principal"
    generic_markers = (
        "is_canonical_true", "canonical_true", "canonical_1",
        "is_canonical_yes", "canonical_yes",
        "primary_transcript_true", "primary_transcript_1",
        "primary_transcript_yes", "is_primary_true", "is_primary_1", "is_primary_yes",
    )
    if any(marker in normalized for marker in generic_markers):
        return 4, "Primary transcript"
    return None, ""


def select_primary_isoforms(
    items: Sequence[AnnotationItem], mode: str
) -> List[AnnotationItem]:
    """Select the best available annotated-primary tier independently per gene."""
    if mode == "all":
        return list(items)
    grouped = {}
    for index, item in enumerate(items):
        grouped.setdefault(item.group or f"item:{index}", []).append(item)
    selected = []
    for members in grouped.values():
        ranks = []
        for item in members:
            if item.primary_rank is not None:
                ranks.append(item.primary_rank)
        if ranks:
            best_rank = min(ranks)
            selected.extend(item for item in members if item.primary_rank == best_rank)
        elif mode == "prefer":
            selected.extend(members)
    return sorted(selected, key=lambda item: (item.start, item.end, item.name))


def parse_bed(lines: Iterable[str], chrom: str, start: int, end: int) -> List[AnnotationItem]:
    items: List[AnnotationItem] = []
    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) < 3:
            fields = line.split()
        if len(fields) < 3 or fields[0] != chrom:
            continue
        try:
            item_start, item_end = int(fields[1]), int(fields[2])
        except ValueError as exc:
            raise ValueError(f"Invalid BED coordinates: {line.rstrip()}") from exc
        if item_end <= start or item_start >= end:
            continue

        name = fields[3] if len(fields) > 3 and fields[3] != "." else ""
        strand = fields[5] if len(fields) > 5 and fields[5] in ("+", "-") else "."
        exons = [(item_start, item_end)]
        thick_blocks = list(exons)
        utrs: List[Interval] = []

        if len(fields) >= 12:
            try:
                block_count = int(fields[9])
                sizes = []
                for value in fields[10].rstrip(",").split(","):
                    if value:
                        sizes.append(int(value))
                offsets = []
                for value in fields[11].rstrip(",").split(","):
                    if value:
                        offsets.append(int(value))
                if len(sizes) != block_count or len(offsets) != block_count:
                    raise ValueError("BED blockCount does not match blockSizes/blockStarts")
                exons = []
                for size, offset in zip(sizes, offsets):
                    exons.append(
                        (item_start + offset, item_start + offset + size)
                    )
                thick_start, thick_end = int(fields[6]), int(fields[7])
            except ValueError as exc:
                raise ValueError(f"Invalid BED12 record: {line.rstrip()} ({exc})") from exc

            if thick_end > thick_start:
                coding_range = [(thick_start, thick_end)]
                thick_blocks = []
                for exon_start, exon_end in exons:
                    coding_start = max(exon_start, thick_start)
                    coding_end = min(exon_end, thick_end)
                    if coding_start < coding_end:
                        thick_blocks.append((coding_start, coding_end))
                utrs = subtract_intervals(exons, coding_range)
            else:
                thick_blocks = exons

        items.append(AnnotationItem(
            start=item_start, end=item_end, name=name, strand=strand,
            blocks=merge_intervals(thick_blocks), utrs=merge_intervals(utrs),
        ))
    return items


def parse_vcf(lines: Iterable[str], chrom: str, start: int, end: int) -> List[AnnotationItem]:
    """Convert VCF records into zero-based annotation intervals."""
    items = []
    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) < 5 or fields[0] != chrom:
            continue
        try:
            variant_start = int(fields[1]) - 1
        except ValueError as exc:
            raise ValueError(f"Invalid VCF position: {line.rstrip()}") from exc
        reference = fields[3]
        alternate = fields[4]
        variant_end = variant_start + max(len(reference), 1)
        if len(fields) > 7:
            for info_field in fields[7].split(";"):
                if info_field.startswith("END="):
                    try:
                        variant_end = max(variant_end, int(info_field[4:]))
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid VCF INFO/END value: {line.rstrip()}"
                        ) from exc
                    break
        if variant_end <= start or variant_start >= end:
            continue
        identifier = fields[2] if fields[2] != "." else f"{reference}>{alternate}"
        items.append(AnnotationItem(
            start=variant_start, end=variant_end, name=identifier,
            blocks=[(variant_start, variant_end)],
        ))
    return items


def parse_bedgraph(
    lines: Iterable[str], chrom: str, start: int, end: int
) -> List[AnnotationItem]:
    """Parse zero-based bedGraph-style chrom/start/end/log2 intervals."""
    items = []
    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) < 4:
            fields = line.split()
        if len(fields) < 4 or fields[0] != chrom:
            continue
        try:
            item_start, item_end = int(fields[1]), int(fields[2])
            value = float(fields[3])
        except ValueError as exc:
            if fields[3].strip().lower() in {".", "na", "nan", "null"}:
                continue
            normalized = {field.lower() for field in fields[:4]}
            if normalized & {"chrom", "chromosome", "start", "end", "value", "log2"}:
                continue
            raise ValueError(f"Invalid bedGraph/log2 record: {line.rstrip()}") from exc
        if not isfinite(value):
            continue
        if item_end <= item_start:
            raise ValueError(f"CNV interval end must exceed start: {line.rstrip()}")
        if item_end <= start or item_start >= end:
            continue
        sample = fields[4] if len(fields) > 4 else ""
        items.append(AnnotationItem(
            item_start, item_end, sample, blocks=[(item_start, item_end)],
            value=value, sample=sample,
        ))
    return items


def parse_seg(
    lines: Iterable[str], chrom: str, start: int, end: int
) -> List[AnnotationItem]:
    """Parse standard SEG tables with 1-based inclusive coordinates."""
    items = []
    columns = None
    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) < 4:
            fields = line.split()
        if len(fields) < 4:
            continue

        normalized = []
        for field in fields:
            normalized.append(re.sub(r"[^a-z0-9]", "", field.lower()))
        if columns is None and any(value in ("chrom", "chr", "chromosome") for value in normalized):
            candidates = {
                "sample": ("sample", "sampleid", "id"),
                "chrom": ("chrom", "chr", "chromosome"),
                "start": ("start", "locstart", "startpos"),
                "end": ("end", "locend", "endpos"),
                "value": (
                    "segmentmean", "segmean", "mean", "value", "log2",
                    "log2ratio", "copynumberratio",
                ),
            }
            columns = {}
            for role, names in candidates.items():
                columns[role] = next(
                    (index for index, value in enumerate(normalized) if value in names), None
                )
            required = (columns["chrom"], columns["start"], columns["end"], columns["value"])
            if any(index is None for index in required):
                raise ValueError(
                    "SEG header must contain chromosome, start, end, and segment-mean/log2 columns."
                )
            continue

        if columns:
            chrom_index = columns["chrom"]
            start_index = columns["start"]
            end_index = columns["end"]
            value_index = columns["value"]
            sample_index = columns["sample"]
        elif len(fields) >= 6:
            sample_index, chrom_index, start_index, end_index, value_index = 0, 1, 2, 3, 5
        elif len(fields) == 5:
            sample_index, chrom_index, start_index, end_index, value_index = 0, 1, 2, 3, 4
        else:
            sample_index, chrom_index, start_index, end_index, value_index = None, 0, 1, 2, 3

        required_indexes = (chrom_index, start_index, end_index, value_index)
        if max(required_indexes) >= len(fields):
            raise ValueError(f"Invalid SEG record: {line.rstrip()}")
        if fields[chrom_index] != chrom:
            continue
        if fields[value_index].strip().lower() in {".", "na", "nan", "null"}:
            continue
        try:
            item_start = int(fields[start_index]) - 1
            item_end = int(fields[end_index])
            value = float(fields[value_index])
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Invalid SEG record: {line.rstrip()}") from exc
        if not isfinite(value):
            continue
        if item_end <= item_start:
            raise ValueError(f"CNV interval end must exceed start: {line.rstrip()}")
        if item_end <= start or item_start >= end:
            continue
        sample = fields[sample_index] if sample_index is not None else ""
        items.append(AnnotationItem(
            item_start, item_end, sample, blocks=[(item_start, item_end)],
            value=value, sample=sample,
        ))
    return items


def gff_group_ids(feature_type: str, attributes: Dict[str, str]) -> List[str]:
    transcript_id = attributes.get("transcript_id")
    if transcript_id:
        return [transcript_id]
    if feature_type in TRANSCRIPT_TYPES:
        identifier = attributes.get("ID")
        if identifier:
            return [identifier]
    parent = attributes.get("Parent")
    if parent:
        parents = []
        for value in parent.split(","):
            if value:
                parents.append(value)
        return parents
    fallback = attributes.get("gene_id") or attributes.get("ID")
    return [fallback] if fallback else []


def parse_gff(lines: Iterable[str], chrom: str, start: int, end: int) -> List[AnnotationItem]:
    models: Dict[str, dict] = {}
    genes: List[AnnotationItem] = []
    gene_labels: Dict[str, str] = {}

    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) != 9 or fields[0] != chrom:
            continue
        feature_type = fields[2].lower()
        try:
            feature_start, feature_end = int(fields[3]) - 1, int(fields[4])
        except ValueError as exc:
            raise ValueError(f"Invalid GFF/GTF coordinates: {line.rstrip()}") from exc
        if feature_end <= start or feature_start >= end:
            continue
        attributes = parse_gff_attributes(fields[8])
        primary_rank, primary_label = primary_isoform_annotation(fields[8])
        strand = fields[6] if fields[6] in ("+", "-") else "."
        label = (
            attributes.get("Name") or attributes.get("gene_name") or
            attributes.get("transcript_name") or attributes.get("gene_id") or ""
        )

        if feature_type == "gene":
            gene_id = attributes.get("gene_id") or attributes.get("ID") or label
            if gene_id:
                gene_labels[gene_id] = label or gene_id
            genes.append(AnnotationItem(
                feature_start, feature_end, label, strand,
                blocks=[(feature_start, feature_end)], group=gene_id,
                group_label=label or gene_id,
            ))
            continue
        if feature_type not in TRANSCRIPT_TYPES | {"exon", "cds"} | UTR_TYPES:
            continue

        group_ids = gff_group_ids(feature_type, attributes)
        if not group_ids:
            group_ids = [f"anonymous:{fields[0]}:{feature_start}:{feature_end}:{feature_type}"]
        for group_id in group_ids:
            gene_group = attributes.get("gene_id")
            if not gene_group and feature_type in TRANSCRIPT_TYPES:
                gene_group = attributes.get("Parent", "").split(",")[0]
            gene_group = gene_group or group_id
            gene_label = attributes.get("gene_name") or gene_labels.get(gene_group, "")
            transcript_label = (
                attributes.get("transcript_name")
                or attributes.get("Name")
                or attributes.get("transcript_id")
                or group_id
            )
            model = models.setdefault(group_id, {
                "start": feature_start, "end": feature_end, "strand": strand,
                "name": label or group_id, "group": gene_group,
                "group_label": gene_label,
                "transcript_label": transcript_label,
                "exons": [], "cds": [], "utrs": [],
                "primary_rank": primary_rank,
                "primary_label": primary_label,
            })
            model["start"] = min(model["start"], feature_start)
            model["end"] = max(model["end"], feature_end)
            if label and (not model["name"] or model["name"] == group_id):
                model["name"] = label
            if strand != ".":
                model["strand"] = strand
            if gene_group != group_id or not model["group"]:
                model["group"] = gene_group
            if gene_label:
                model["group_label"] = gene_label
            if transcript_label and transcript_label != group_id:
                model["transcript_label"] = transcript_label
            if primary_rank is not None and (
                model["primary_rank"] is None or primary_rank < model["primary_rank"]
            ):
                model["primary_rank"] = primary_rank
                model["primary_label"] = primary_label
            if feature_type == "exon":
                model["exons"].append((feature_start, feature_end))
            elif feature_type == "cds":
                model["cds"].append((feature_start, feature_end))
            elif feature_type in UTR_TYPES:
                model["utrs"].append((feature_start, feature_end))

    if not models:
        return genes

    items: List[AnnotationItem] = []
    for model in models.values():
        exons = merge_intervals(model["exons"])
        cds = merge_intervals(model["cds"])
        explicit_utrs = merge_intervals(model["utrs"])
        if not exons:
            exons = merge_intervals(cds + explicit_utrs)
        if not exons:
            exons = [(model["start"], model["end"])]

        if cds:
            blocks = cds
            utrs = merge_intervals(explicit_utrs + subtract_intervals(exons, cds))
        elif explicit_utrs:
            blocks = subtract_intervals(exons, explicit_utrs)
            utrs = explicit_utrs
        else:
            blocks = exons
            utrs = []
        items.append(AnnotationItem(
            start=min(model["start"], min(value[0] for value in exons)),
            end=max(model["end"], max(value[1] for value in exons)),
            name=model["name"], strand=model["strand"], blocks=blocks, utrs=utrs,
            group=model["group"],
            group_label=gene_labels.get(model["group"], model["group_label"]),
            transcript_label=model["transcript_label"],
            primary_rank=model["primary_rank"],
            primary_label=model["primary_label"],
        ))
    return items


class AnnotationSource:
    """A validated annotation file that can fetch one genomic window."""

    def __init__(
        self,
        path: str,
        label: Optional[str] = None,
        color: Optional[str] = None,
        kind: Optional[str] = None,
        display_mode: str = "pack",
        primary_isoforms: str = "all",
        track_colors: Optional[Dict[str, str]] = None,
    ):
        self.path = str(path)
        self.label = label or default_label(self.path)
        explicit_kind = kind.lower() if kind else None
        if explicit_kind == "auto":
            explicit_kind = None
        if explicit_kind == "bdg":
            explicit_kind = "bedgraph"
        if explicit_kind is not None and explicit_kind not in SUPPORTED_TRACK_FORMATS:
            raise ValueError(
                f"Unsupported annotation type '{kind}'. Choose bed, gff, gff3, gtf, vcf, "
                "seg, bedgraph, log2, cnv, or auto."
            )
        self.kind = explicit_kind or infer_track_format(self.path)
        if display_mode not in ANNOTATION_DISPLAY_MODES:
            raise ValueError(
                f"Unsupported annotation display mode '{display_mode}'. Choose "
                f"{', '.join(ANNOTATION_DISPLAY_MODES)}."
            )
        self.display_mode = display_mode
        if primary_isoforms not in PRIMARY_ISOFORM_MODES:
            raise ValueError(
                f"Unsupported primary-isoform mode '{primary_isoforms}'. Choose "
                f"{', '.join(PRIMARY_ISOFORM_MODES)}."
            )
        self.primary_isoforms = primary_isoforms
        colors = dict(DEFAULT_TRACK_COLORS)
        colors.update(track_colors or {})
        self.color_by_sign = color is None and self.kind in CNV_TRACK_FORMATS
        if color is not None:
            self.color = normalize_track_color(color)
        elif self.kind == "bed":
            self.color = colors["bed"]
        elif self.kind == "vcf":
            self.color = colors["vcf"]
        elif self.kind in CNV_TRACK_FORMATS:
            self.color = colors["cnv"]
        else:
            self.color = colors["gene"]
        self.compressed = Path(self.path).name.lower().endswith(COMPRESSED_SUFFIXES)
        if not Path(self.path).is_file():
            raise ValueError(f"Annotation track not found: {self.path}")
        if self.compressed:
            if not (
                Path(self.path + ".tbi").is_file() or Path(self.path + ".csi").is_file()
            ):
                raise ValueError(
                    f"Compressed annotation track requires a tabix index (.tbi or .csi): {self.path}"
                )
            # Validate both BGZF data and its index during CLI setup instead
            # of failing later during rendering. Region fetches below always
            # use TabixFile.fetch(), never a full compressed-file scan.
            try:
                with pysam.TabixFile(self.path) as tabix:
                    tuple(tabix.contigs)
            except (OSError, ValueError) as exc:
                raise ValueError(
                    f"Cannot open BGZF/tabix annotation track '{self.path}': {exc}"
                ) from exc

    def iter_lines(self, chrom: str, start: int, end: int) -> Iterator[str]:
        if not self.compressed:
            with open(self.path, "r", encoding="utf-8") as handle:
                for line in handle:
                    if line.strip() and not line.startswith(("#", "track", "browser")):
                        yield line
            return

        with pysam.TabixFile(self.path) as tabix:
            contigs = set(tabix.contigs)
            query_chrom = chrom
            if query_chrom not in contigs:
                alternate = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                if alternate not in contigs:
                    return
                query_chrom = alternate
            for line in tabix.fetch(query_chrom, max(0, start), end):
                if line.strip() and not line.startswith("#"):
                    # Parsers compare against the requested/BAM contig.
                    if query_chrom != chrom:
                        fields = line.split("\t")
                        for index, value in enumerate(fields):
                            if value == query_chrom:
                                fields[index] = chrom
                                break
                        line = "\t".join(fields)
                    yield line

    def fetch(self, chrom: str, start: int, end: int) -> LoadedAnnotationTrack:
        lines = self.iter_lines(chrom, start, end)
        if self.kind == "bed":
            items = parse_bed(lines, chrom, start, end)
        elif self.kind == "vcf":
            items = parse_vcf(lines, chrom, start, end)
        elif self.kind in ("bedgraph", "log2"):
            items = parse_bedgraph(lines, chrom, start, end)
        elif self.kind in ("seg", "cnv"):
            items = parse_seg(lines, chrom, start, end)
        else:
            items = parse_gff(lines, chrom, start, end)
        if self.kind in ("gff", "gff3", "gtf"):
            items = select_primary_isoforms(items, self.primary_isoforms)
        if self.kind in CNV_TRACK_FORMATS:
            rows = [items] if items else []
        elif self.display_mode == "expand":
            rows = []
            ordered_items = sorted(
                items, key=lambda value: (value.start, value.end, value.name)
            )
            for item in ordered_items:
                rows.append([item])
        elif self.display_mode == "collapse" and self.kind in ("gff", "gff3", "gtf"):
            rows = pack_annotation_items(collapse_annotation_items(items))
        elif self.display_mode == "collapse":
            rows = [list(items)] if items else []
        else:
            rows = pack_annotation_items(items)
        return LoadedAnnotationTrack(
            label=self.label, kind=self.kind, color=self.color,
            items=items, rows=rows, display_mode=self.display_mode,
            color_by_sign=self.color_by_sign,
        )


class BafSource:
    """A VCF sample rendered as B-allele fractions at heterozygous SNVs."""

    def __init__(
        self, path: str, label: Optional[str] = None,
        sample: Optional[str] = None, color: Optional[str] = None,
        track_colors: Optional[Dict[str, str]] = None,
    ):
        self.path = str(path)
        self.label = label or f"{default_label(self.path)} BAF"
        colors = dict(DEFAULT_TRACK_COLORS)
        colors.update(track_colors or {})
        self.color = normalize_track_color(color or colors["baf"])
        path_object = Path(self.path)
        if not path_object.is_file():
            raise ValueError(f"BAF VCF not found: {self.path}")
        lower_name = path_object.name.lower()
        self.indexed = lower_name.endswith(COMPRESSED_SUFFIXES) or lower_name.endswith(".bcf")
        if self.indexed and not (
            Path(self.path + ".tbi").is_file() or Path(self.path + ".csi").is_file()
        ):
            raise ValueError(
                f"Compressed BAF VCF requires a tabix/CSI index (.tbi or .csi): {self.path}"
            )
        try:
            with pysam.VariantFile(self.path) as variant_file:
                samples = list(variant_file.header.samples)
        except (OSError, ValueError) as exc:
            raise ValueError(f"Cannot open BAF VCF '{self.path}': {exc}") from exc
        if not samples:
            raise ValueError(f"BAF VCF has no genotype samples: {self.path}")
        if sample and sample not in samples:
            raise ValueError(
                f"BAF sample '{sample}' is not present in {self.path}; choose from {', '.join(samples)}."
            )
        self.sample = sample or samples[0]

    def fetch(self, chrom: str, start: int, end: int) -> LoadedAnnotationTrack:
        items = []
        with pysam.VariantFile(self.path) as variant_file:
            contigs = set(variant_file.header.contigs)
            query_chrom = chrom
            if query_chrom not in contigs:
                alternate = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                if alternate not in contigs:
                    return LoadedAnnotationTrack(
                        self.label, "baf", self.color, [], [], "collapse"
                    )
                query_chrom = alternate
            if self.indexed:
                records = variant_file.fetch(query_chrom, max(0, start), end)
            else:
                records = (
                    record for record in variant_file
                    if record.contig == query_chrom and record.stop > start and record.start < end
                )
            for record in records:
                if (
                    len(record.ref) != 1 or not record.alts or len(record.alts) != 1
                    or len(record.alts[0]) != 1
                ):
                    continue
                sample_data = record.samples[self.sample]
                genotype = sample_data.get("GT")
                if not genotype or set(allele for allele in genotype if allele is not None) != {0, 1}:
                    continue
                allele_depths = sample_data.get("AD")
                baf = None
                if allele_depths and len(allele_depths) >= 2:
                    ref_depth, alt_depth = allele_depths[0], allele_depths[1]
                    if ref_depth is not None and alt_depth is not None and ref_depth + alt_depth > 0:
                        baf = alt_depth / (ref_depth + alt_depth)
                if baf is None:
                    allele_fraction = sample_data.get("AF")
                    if isinstance(allele_fraction, (tuple, list)):
                        allele_fraction = allele_fraction[0] if allele_fraction else None
                    if allele_fraction is not None:
                        baf = float(allele_fraction)
                if baf is None or not isfinite(baf) or not 0 <= baf <= 1:
                    continue
                position = record.start
                name = record.id or f"{record.ref}>{record.alts[0]}"
                items.append(AnnotationItem(
                    position, position + 1, name, blocks=[(position, position + 1)],
                    value=baf, sample=self.sample,
                ))
        return LoadedAnnotationTrack(
            self.label, "baf", self.color, items,
            [items] if items else [], "collapse",
        )


def build_annotation_sources(
    paths: Optional[Sequence[str]], labels: Optional[Sequence[str]] = None,
    display_mode: str = "pack",
    primary_isoforms: str = "all",
    track_colors: Optional[Dict[str, str]] = None,
) -> List[AnnotationSource]:
    paths = list(paths or [])
    labels = list(labels or [])
    if len(labels) > len(paths):
        raise ValueError("More --track_label values were supplied than --track files.")
    labels.extend([None] * (len(paths) - len(labels)))
    sources = []
    for path, label in zip(paths, labels):
        sources.append(AnnotationSource(
            path, label, display_mode=display_mode,
            primary_isoforms=primary_isoforms,
            track_colors=track_colors,
        ))
    return sources


def build_custom_annotation_sources(
    specifications: Optional[Sequence[Sequence[str]]],
    default_display_mode: str = "pack",
    primary_isoforms: str = "all",
    track_colors: Optional[Dict[str, str]] = None,
) -> List[AnnotationSource]:
    """Build repeatable ``FILE TYPE NAME COLOR [DISPLAY]`` track definitions."""
    sources = []
    for specification in specifications or []:
        if len(specification) not in (4, 5):
            raise ValueError(
                "Each --custom_track requires FILE TYPE NAME COLOR and an optional DISPLAY."
            )
        path, kind, name, color = specification[:4]
        display_mode = specification[4] if len(specification) == 5 else default_display_mode
        sources.append(AnnotationSource(
            path=path, kind=kind, label=name, color=color,
            display_mode=display_mode, primary_isoforms=primary_isoforms,
            track_colors=track_colors,
        ))
    return sources


def build_baf_sources(
    paths: Optional[Sequence[str]], labels: Optional[Sequence[str]] = None,
    samples: Optional[Sequence[str]] = None,
    track_colors: Optional[Dict[str, str]] = None,
) -> List[BafSource]:
    paths = list(paths or [])
    labels = list(labels or [])
    samples = list(samples or [])
    if len(labels) > len(paths):
        raise ValueError("More --baf_track_label values were supplied than --baf_vcf files.")
    if len(samples) > len(paths):
        raise ValueError("More --baf_sample values were supplied than --baf_vcf files.")
    labels.extend([None] * (len(paths) - len(labels)))
    samples.extend([None] * (len(paths) - len(samples)))
    sources = []
    for path, label, sample in zip(paths, labels, samples):
        sources.append(BafSource(
            path, label=label, sample=sample, track_colors=track_colors
        ))
    return sources
