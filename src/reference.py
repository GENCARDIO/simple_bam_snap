"""Thin wrapper around a reference FASTA for a single genomic window.

Reference bases are only needed to call mismatches (everything else -
indels, clips, split-read gaps - comes straight out of the CIGAR/SA tags),
so this stays deliberately small and is entirely optional: tools that call
into here should degrade gracefully when no FASTA is supplied.
"""
from __future__ import annotations

import os
from typing import Optional

import pysam


class ReferenceWindow:
    """Loads (and caches) the reference sequence for one chrom:start-end window."""

    def __init__(self, fasta_path: Optional[str], chrom: str, start: int, end: int):
        self.fasta_path = fasta_path
        self.chrom = chrom
        self.start = start  # 0-based, inclusive
        self.end = end      # 0-based, exclusive
        self.sequence: Optional[str] = None

        if fasta_path:
            self.sequence = self._load(fasta_path)

    def _load(self, fasta_path: str) -> str:
        if not os.path.isfile(fasta_path):
            raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")
        try:
            fasta = pysam.FastaFile(fasta_path)
        except (OSError, ValueError):
            # Missing/stale .fai index - (re)build it once and retry.
            pysam.faidx(fasta_path)
            fasta = pysam.FastaFile(fasta_path)

        try:
            contig = self._resolve_contig(fasta, self.chrom)
            return fasta.fetch(contig, max(0, self.start), self.end).upper()
        finally:
            fasta.close()

    @staticmethod
    def _resolve_contig(fasta: "pysam.FastaFile", chrom: str) -> str:
        """Tolerate chr1 vs 1 naming mismatches between BAM and FASTA."""
        names = set(fasta.references)
        if chrom in names:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        if alt in names:
            return alt
        raise ValueError(
            f"Contig '{chrom}' not found in reference FASTA (tried '{alt}' too)"
        )

    def base_at(self, ref_pos: int) -> Optional[str]:
        """Return the reference base at a 0-based genomic coordinate, or None."""
        if self.sequence is None:
            return None
        idx = ref_pos - self.start
        if 0 <= idx < len(self.sequence):
            return self.sequence[idx]
        return None

    @property
    def available(self) -> bool:
        return self.sequence is not None
