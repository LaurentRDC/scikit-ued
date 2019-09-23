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
Digital Micrograph 3/4 (*.dm3, *.dm4),Merlin Image Binary (*.mib), and all 
formats supported by scikit-image."""

# Possibility in the future to add more utilities with more subparsers
diffshow_parser = subparsers.add_parser("diffshow", description=DIFFSHOW_HELP)
diffshow_parser.add_argument("filename", type=Path, help=DIFFSHOW_FILENAME_HELP)

CRYSTINFO_HELP = "Display the crystallographic information related to a crystal file"

CRYSTINFO_FILENAME_HELP = """Path to the file. Supported formats are Crystallography 
Information File (*.cif) and Quantum ESPRESSO PWSCF (*.pwscf). """

crystinfo = subparsers.add_parser('crystinfo', description=CRYSTINFO_HELP)
crystinfo.add_argument("filename", type=Path, help=CRYSTINFO_FILENAME_HELP)


def main(args=None):
    if args is None:
        args = parser.parse_args()

    if args.command == "diffshow":
        main_diffshow(args.filename)
    
    if args.command == 'crystinfo':
        main_crystinfo(args.filename)


def main_diffshow(fname):
    """ Display an interactive window """
    if not WITH_PYQTGRAPH:
        print(
            "PyQtGraph is required for this functionality. You can install PyQtGraph either with \
            conda or pip, depending on what you used to install scikit-ued."
        )
        sys.exit(1)
    diffshow(fname)
    sys.exit(0)

def main_crystinfo(fname):
    """ Display information associated with crystal file """
    cryst = None
    for constructor in (Crystal.from_cif, Crystal.from_pwscf,):
        try:
            cryst = constructor(fname)
        except:
            continue
        else:
            break
    
    if cryst is None:
        print(f"{fname} could not be parsed.")
        sys.exit(1)

    unitcell = repr(cryst).replace("< ", "").replace(" >", "") # This is from python's object representation
    sym = f"""Symmetry information:
    International symbol 
                (short) ..... {cryst.international_symbol}
                (full) ...... {cryst.international_full}
    International number .... {cryst.international_number}
    Hermann-Mauguin symbol .. {cryst.hm_symbol}
    Pointgroup .............. {cryst.pointgroup}
    Hall Number ............. {cryst.hall_number}
    Centering ............... {cryst.centering}"""

    print(unitcell)
    print(sym)

if __name__ == "__main__":
    main()
