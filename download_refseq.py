#!/usr/bin/env python3
"""Pre-populate both supported human NCBI RefSeq annotation caches."""
from src.refseq import ensure_refseq


for assembly in ("hg19", "hg38"):
    print(ensure_refseq(assembly))
