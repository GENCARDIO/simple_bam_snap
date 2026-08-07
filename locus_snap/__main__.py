"""Enables ``python -m locus_snap``."""
import sys

from locus_snap.cli import main

if __name__ == "__main__":
    sys.exit(main())
