#!/usr/bin/env python3
"""Backward-compatible launcher for the renamed LocusSnap CLI."""

from locus_snap.cli import *  # noqa: F401,F403
from locus_snap.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
