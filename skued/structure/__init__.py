# -*- coding: utf-8 -*-
"""
Structure package
-----------------
This package allows for manipulation and modelling of atomic structures, especially
in crystalline form.
"""
class ParseError(IOError):
    pass
    
from .atom import Atom, real_coords, frac_coords
from .atom_data import ELEM_TO_MAGMOM, ELEM_TO_MASS, ELEM_TO_NAME, ELEM_TO_NUM, NUM_TO_ELEM
from .base import AtomicStructure, Base
from .lattice import Lattice, lattice_vectors_from_parameters
from .pdb_parser import PDBParser
from .spg_data import Hall2Number
from .cif_parser import CIFParser
from .crystal import Crystal, symmetry_expansion