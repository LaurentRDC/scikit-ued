# -*- coding: utf-8 -*-
"""
Structure package
-----------------
This package allows for manipulation and modelling of atomic structures, especially
in crystalline form.
"""


from crystals import (
    Atom,
    AtomicStructure,
    Base,
    Crystal,
    Lattice,
    LatticeSystem,
    lattice_system,
    CIFParser,
    PDBParser,
    CODParser,
    symmetry_expansion,

    NUM_TO_ELEM,
    ELEM_TO_NUM,
    ELEM_TO_NAME

)
