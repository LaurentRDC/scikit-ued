# -*- coding: utf-8 -*-
"""
scikit-ued command-line utilities
"""
import argparse
import sys
from pathlib import Path

from . import __version__
from .io import WITH_PYQTGRAPH, diffshow

parser = argparse.ArgumentParser(
    prog="skued",
    description=f"scikit-ued {__version__} command-line utilities.",
)

subparsers = parser.add_subparsers(title="command", dest="command")

# Possibility in the future to add more utilities with more subparsers
diffshow_parser = subparsers.add_parser(
    "diffshow", help="Read a file and show interactive window. Requires PyQtGraph."
)
diffshow_parser.add_argument(
    "filename",
    type=Path,
    help="Path to file. All formats supported by ``skued.diffread`` are supported here, including TIFFs and DM3/DM4",
)

def main(args=None):
    if args is None:
        args = parser.parse_args()
    
    if args.command == 'diffshow':
        main_diffshow(args.filename)

def main_diffshow(fname):
    """ Display an interactive window """
    if not WITH_PYQTGRAPH:
        print("PyQtGraph is required for this functionality.")
        sys.exit(1)
    diffshow(fname)
    sys.exit(0)


if __name__ == "__main__":
    main()
