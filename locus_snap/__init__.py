"""LocusSnap: IGV-like genomic snapshots from an indexed BAM, built on pysam."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("locus-snap")
except PackageNotFoundError:
    # Running from a source checkout that was never `pip install`-ed.
    __version__ = "0.0.0+unknown"

from locus_snap.cli import apply_config_preferences, build_parser, main

__all__ = ["__version__", "apply_config_preferences", "build_parser", "main"]
