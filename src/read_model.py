"""Parses raw pysam records into the per-read features the rest of the tool
sorts and draws by.

The whole point of this redesign is to stop treating an alignment as a blob
of `samtools tview` text and instead treat it as structured data: CIGAR
blocks, soft/hard clips, mismatches against the reference, and the SA tag of
split (supplementary) alignments. Two "gap" signals fall out of that:

- ``cigar_gap_len``: total length of I/D operations in the read's own CIGAR
  - the classic way an aligner represents an indel inline.
- ``sa_gap_len``: the reference distance between this alignment and its
  closest same-chromosome supplementary partner - the way an aligner
  represents an indel/SV it could *not* fit into one CIGAR string, by
  splitting the read instead.

``gap_length = max(cigar_gap_len, sa_gap_len)`` is the unified "how gapped is
this alignment" metric used for sorting: a read carrying a long indel shows
up under either representation.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import pysam

from src.reference import ReferenceWindow

CIGAR_OPS = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}
REF_CONSUMING = ("M", "D", "N", "=", "X")
QUERY_CONSUMING = ("M", "I", "S", "=", "X")
SA_FIELD_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class CigarBlock:
    op: str
    ref_pos: int    # 0-based genomic coordinate where this block starts
    query_pos: int  # index into query_sequence where this block starts
    length: int


@dataclass
class SAEntry:
    rname: str
    start: int  # 0-based
    end: int    # 0-based, exclusive
    strand: str
    cigar: str
    mapq: int
    nm: int


def iter_cigar_blocks(cigartuples, ref_start: int):
    """Walk a pysam cigartuples list, yielding one CigarBlock per operation."""
    ref_pos = ref_start
    query_pos = 0
    for op_code, length in cigartuples:
        op = CIGAR_OPS[op_code]
        yield CigarBlock(op=op, ref_pos=ref_pos, query_pos=query_pos, length=length)
        if op in REF_CONSUMING:
            ref_pos += length
        if op in QUERY_CONSUMING:
            query_pos += length


def cigar_reference_length(cigar_str: str) -> int:
    """Reference span implied by a CIGAR string (e.g. an SA tag's CIGAR field)."""
    total = 0
    for num, op in SA_FIELD_PATTERN.findall(cigar_str):
        if op in REF_CONSUMING:
            total += int(num)
    return total


def parse_sa_tag(raw: str) -> List[SAEntry]:
    entries = []
    for record in raw.strip().split(";"):
        if not record:
            continue
        parts = record.split(",")
        if len(parts) < 6:
            continue
        rname, pos, strand, cigar, mapq, nm = parts[:6]
        start = int(pos) - 1
        entries.append(
            SAEntry(
                rname=rname,
                start=start,
                end=start + cigar_reference_length(cigar),
                strand=strand,
                cigar=cigar,
                mapq=int(mapq) if mapq.lstrip("-").isdigit() else 0,
                nm=int(nm) if nm.strip().lstrip("-").isdigit() else 0,
            )
        )
    return entries


def compute_pair_orientation(seg: pysam.AlignedSegment) -> Optional[str]:
    """Classify a same-chromosome mapped pair's relative orientation the way
    IGV does, using only this read's own FLAG/PNEXT fields (no mate fetch):

    - "FR": normal Illumina "innie" pair.
    - "RF": everted / mate-pair-style "outie".
    - "FF" / "RR": both mates on the same strand.

    Returns None for unpaired reads, reads with an unmapped mate, or
    inter-chromosomal pairs (where "orientation" isn't a same-axis concept).
    """
    if not seg.is_paired or seg.is_unmapped or seg.mate_is_unmapped:
        return None
    if seg.next_reference_name != seg.reference_name:
        return None

    this_pos, mate_pos = seg.reference_start, seg.next_reference_start
    this_rev, mate_rev = seg.is_reverse, seg.mate_is_reverse
    if mate_pos > this_pos or (mate_pos == this_pos and not this_rev):
        left_rev, right_rev = this_rev, mate_rev
    else:
        left_rev, right_rev = mate_rev, this_rev

    if not left_rev and right_rev:
        return "FR"
    if left_rev and not right_rev:
        return "RF"
    if not left_rev and not right_rev:
        return "FF"
    return "RR"


def classify_insert_sizes(reads: List["AlignedRead"], sigma: float = 3.0, min_pairs: int = 10) -> None:
    """Flags FR pairs with an outlier insert size as small_insert/large_insert,
    mutating each read's `pair_category`/`is_small_insert`/`is_large_insert`
    in place.

    The "expected" insert size range is estimated from this same read set
    (median +/- sigma * robust MAD-based stdev) rather than a fixed constant,
    since expected fragment size varies by library/protocol. With too few FR
    pairs to estimate a range, nothing is flagged.
    """
    samples = []
    for r in reads:
        if r.pair_orientation == "FR" and r.insert_size > 0:
            samples.append(r.insert_size)
    samples.sort()
    if len(samples) < min_pairs:
        return

    n = len(samples)
    median = samples[n // 2] if n % 2 else (samples[n // 2 - 1] + samples[n // 2]) / 2
    deviations = []
    for s in samples:
        deviations.append(abs(s - median))
    deviations.sort()
    mad = deviations[n // 2] if n % 2 else (deviations[n // 2 - 1] + deviations[n // 2]) / 2
    robust_std = mad * 1.4826 or 1.0
    lo, hi = median - sigma * robust_std, median + sigma * robust_std

    for r in reads:
        if r.pair_orientation != "FR" or r.insert_size <= 0:
            continue
        if r.insert_size < lo:
            r.is_small_insert = True
            r.pair_category = "small_insert"
        elif r.insert_size > hi:
            r.is_large_insert = True
            r.pair_category = "large_insert"


def closest_same_chrom_gap(
    ref_start: int, ref_end: int, rname: str, sa_entries: List[SAEntry]
) -> Tuple[int, bool]:
    """Reference-space distance to the nearest same-chromosome SA partner.

    Returns (gap_length, has_cross_chrom_sa). Overlapping partners (e.g. a
    duplication rather than a deletion) report a gap of 0 rather than a
    negative number.
    """
    same_chrom = []
    for entry in sa_entries:
        if entry.rname == rname:
            same_chrom.append(entry)
    has_cross_chrom = len(same_chrom) < len(sa_entries)
    if not same_chrom:
        return 0, has_cross_chrom

    best = None
    for e in same_chrom:
        if e.end <= ref_start:
            gap = ref_start - e.end
        elif e.start >= ref_end:
            gap = e.start - ref_end
        else:
            gap = 0
        if best is None or gap < best:
            best = gap
    return best, has_cross_chrom


class AlignedRead:
    """A single alignment record plus every feature the layout/render/metrics
    layers need, computed once up front."""

    def __init__(
        self, segment: pysam.AlignedSegment,
        reference: Optional[ReferenceWindow] = None,
        haplotype_tag: str = "HP", phase_set_tag: str = "PS",
    ):
        self.segment = segment
        seg = segment

        self.query_name = seg.query_name
        self.query_sequence = seg.query_sequence
        self.query_qualities = list(seg.query_qualities or [])
        self.ref_start = seg.reference_start
        self.ref_end = seg.reference_end if seg.reference_end is not None else seg.reference_start
        self.reference_name = seg.reference_name
        self.mapq = seg.mapping_quality
        self.strand = "-" if seg.is_reverse else "+"
        self.is_reverse = seg.is_reverse
        self.is_duplicate = seg.is_duplicate
        self.is_secondary = seg.is_secondary
        self.is_supplementary = seg.is_supplementary
        self.is_paired = seg.is_paired
        self.mate_is_unmapped = bool(seg.is_paired and seg.mate_is_unmapped)
        self.is_proper_pair = bool(seg.is_paired and seg.is_proper_pair)
        self.insert_size = abs(seg.template_length) if seg.is_paired else 0
        self.flag = seg.flag
        self.haplotype = str(seg.get_tag(haplotype_tag)) if seg.has_tag(haplotype_tag) else None
        self.phase_set = str(seg.get_tag(phase_set_tag)) if seg.has_tag(phase_set_tag) else None

        # --- pair orientation / discordance (see compute_pair_orientation) ---
        has_mapped_mate = seg.is_paired and not seg.mate_is_unmapped
        self.mate_chrom = seg.next_reference_name if has_mapped_mate else None
        self.mate_start = seg.next_reference_start if has_mapped_mate else None
        self.is_interchrom = bool(has_mapped_mate and self.mate_chrom != self.reference_name)
        self.pair_orientation = compute_pair_orientation(seg)
        self.is_small_insert = False
        self.is_large_insert = False  # refined by classify_insert_sizes() over the whole read set

        if self.is_interchrom:
            self.pair_category = "interchrom"
        elif self.pair_orientation == "RF":
            self.pair_category = "everted"
        elif self.pair_orientation == "FF":
            self.pair_category = "ff"
        elif self.pair_orientation == "RR":
            self.pair_category = "rr"
        else:
            self.pair_category = "normal"

        self.blocks: List[CigarBlock] = []
        self.insertions: List[Tuple[int, int]] = []
        self.deletions: List[Tuple[int, int, bool]] = []  # (ref_pos, length, is_skip)
        self.soft_clip_left = 0
        self.soft_clip_right = 0
        self.hard_clip_left = 0
        self.hard_clip_right = 0
        self.mismatches: List[Tuple[int, str]] = []
        self.mismatch_details: List[Tuple[int, str, int]] = []

        cigartuples = seg.cigartuples
        if cigartuples:
            query_seq = seg.query_sequence
            n_ops = len(cigartuples)
            has_ref = reference is not None and reference.available
            for i, blk in enumerate(iter_cigar_blocks(cigartuples, self.ref_start)):
                self.blocks.append(blk)
                if blk.op == "I":
                    self.insertions.append((blk.ref_pos, blk.length))
                elif blk.op == "D":
                    self.deletions.append((blk.ref_pos, blk.length, False))
                elif blk.op == "N":
                    self.deletions.append((blk.ref_pos, blk.length, True))
                elif blk.op == "S":
                    if i == 0:
                        self.soft_clip_left = blk.length
                    if i == n_ops - 1:
                        self.soft_clip_right = blk.length
                elif blk.op == "H":
                    if i == 0:
                        self.hard_clip_left = blk.length
                    if i == n_ops - 1:
                        self.hard_clip_right = blk.length
                elif blk.op in ("M", "=", "X") and query_seq and has_ref:
                    for offset in range(blk.length):
                        rpos = blk.ref_pos + offset
                        qbase = query_seq[blk.query_pos + offset].upper()
                        rbase = reference.base_at(rpos)
                        if rbase and rbase != "N" and qbase != rbase:
                            self.mismatches.append((rpos, qbase))
                            base_quality = (
                                self.query_qualities[blk.query_pos + offset]
                                if blk.query_pos + offset < len(self.query_qualities) else 0
                            )
                            self.mismatch_details.append((rpos, qbase, base_quality))

        deletion_len = 0
        for deletion in self.deletions:
            if not deletion[2]:
                deletion_len += deletion[1]
        insertion_len = 0
        for insertion in self.insertions:
            insertion_len += insertion[1]
        self.cigar_gap_len = deletion_len + insertion_len
        self.soft_clip_total = self.soft_clip_left + self.soft_clip_right
        self.mismatch_count = len(self.mismatches)

        self.sa_entries = parse_sa_tag(seg.get_tag("SA")) if seg.has_tag("SA") else []
        self.sa_count = len(self.sa_entries)
        self.sa_gap_len, self.has_cross_chrom_sa = closest_same_chrom_gap(
            self.ref_start, self.ref_end, self.reference_name, self.sa_entries
        )

        self.gap_length = max(self.cigar_gap_len, self.sa_gap_len)

    @property
    def is_discordant(self) -> bool:
        return self.pair_category != "normal"

    def base_at(self, ref_position: int) -> Optional[str]:
        """Return the aligned base, deletion, skip, or no-call at a locus."""
        sequence = self.query_sequence or ""
        for block in self.blocks:
            block_end = block.ref_pos + block.length
            if not block.ref_pos <= ref_position < block_end:
                continue
            if block.op in ("M", "=", "X"):
                query_index = block.query_pos + ref_position - block.ref_pos
                if query_index < len(sequence):
                    return sequence[query_index].upper()
                return "N"
            if block.op == "D":
                return "-"
            if block.op == "N":
                return "~"
        return None

    def gap_label(self) -> str:
        """Short human-readable summary of the dominant gap signal, for annotating rows."""
        if self.cigar_gap_len == 0 and self.sa_gap_len == 0:
            return ""
        if self.cigar_gap_len >= self.sa_gap_len:
            return f"{self.cigar_gap_len}bp"
        return f"~{self.sa_gap_len}bp SA"

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"AlignedRead({self.query_name!r}, {self.reference_name}:{self.ref_start}-{self.ref_end}, "
            f"gap={self.gap_length}, sa={self.sa_count})"
        )


ONLY_TYPES = ("discordant", "gapped", "split", "softclip")


def matches_only(read: AlignedRead, only_types: Optional[List[str]], min_softclip: int = 1) -> bool:
    """OR filter: True if `read` matches any of the requested --only categories
    (or if no filter was requested at all)."""
    if not only_types:
        return True
    if "discordant" in only_types and read.is_discordant:
        return True
    if "gapped" in only_types and read.gap_length > 0:
        return True
    if "split" in only_types and read.sa_count > 0:
        return True
    if "softclip" in only_types and read.soft_clip_total >= min_softclip:
        return True
    return False


def fetch_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    reference: Optional[ReferenceWindow] = None,
    min_mapq: int = 0,
    include_secondary: bool = False,
    include_supplementary: bool = True,
    include_duplicates: bool = False,
    include_unmapped: bool = False,
    insert_size_sigma: float = 3.0,
    only_types: Optional[List[str]] = None,
    min_softclip: int = 1,
    haplotype_tag: str = "HP",
    phase_set_tag: str = "PS",
    haplotype_filter: Optional[List[str]] = None,
) -> List[AlignedRead]:
    """Fetch + filter + featurize every alignment overlapping [start, end).

    Insert-size discordance is classified across the whole fetched cohort
    (it needs a distribution to compare against), so that step - and any
    --only filtering that depends on it - happens after every read in the
    window has been collected, not per-read during the fetch loop.
    """
    reads = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for segment in bam.fetch(chrom, max(0, start), end):
            if segment.is_unmapped and not include_unmapped:
                continue
            if segment.is_duplicate and not include_duplicates:
                continue
            if segment.is_secondary and not include_secondary:
                continue
            if segment.is_supplementary and not include_supplementary:
                continue
            if segment.mapping_quality < min_mapq:
                continue
            read = AlignedRead(
                segment, reference,
                haplotype_tag=haplotype_tag, phase_set_tag=phase_set_tag,
            )
            if haplotype_filter:
                haplotype_value = read.haplotype if read.haplotype is not None else "untagged"
                if haplotype_value not in haplotype_filter:
                    continue
            reads.append(read)

    classify_insert_sizes(reads, sigma=insert_size_sigma)

    if only_types:
        matching_reads = []
        for read in reads:
            if matches_only(read, only_types, min_softclip):
                matching_reads.append(read)
        reads = matching_reads
    return reads
