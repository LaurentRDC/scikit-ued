# -*- coding: utf-8 -*-
"""
scikit-ued command-line utilities
"""
import argparse
import sys
from pathlib import Path
from crystals import Crystal

from . import __version__
from .io import WITH_PYQTGRAPH, diffshow

parser = argparse.ArgumentParser(
    prog="skued", description=f"scikit-ued {__version__} command-line utilities."
)

subparsers = parser.add_subparsers(title="command", dest="command")


DIFFSHOW_HELP = "" if WITH_PYQTGRAPH else "[UNAVAILABLE] "
DIFFSHOW_HELP += "Read a file and show interactive window. Requires PyQtGraph."

DIFFSHOW_FILENAME_HELP = """Path to file. All formats supported by 
``skued.diffread`` are supported here, including TIFFs (*.tif, *.tiff), 
Digital Micrograph 3/4 (*.dm3, *.dm4), Merlin Image Binary (*.mib), and all 
formats supported by scikit-image."""

# Possibility in the future to add more utilities with more subparsers
diffshow_parser = subparsers.add_parser("diffshow", description=DIFFSHOW_HELP)
diffshow_parser.add_argument("filename", type=Path, help=DIFFSHOW_FILENAME_HELP)


def main(args=None):
    if args is None:
        args = parser.parse_args()

    if args.command == "diffshow":
        main_diffshow(args.filename)


def main_diffshow(fname):
    """Display an interactive window"""
    if not WITH_PYQTGRAPH:
        print(
            "PyQtGraph is required for this functionality. You can install PyQtGraph either with \
            conda or pip, depending on what you used to install scikit-ued."
        )
        sys.exit(1)
    diffshow(fname)
    sys.exit(0)


if __name__ == "__main__":
    main()
