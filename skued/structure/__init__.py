# -*- coding: utf-8 -*-
"""
Structure package
-----------------
This package allows for manipulation and modelling of atomic structures, especially
in crystalline form.
"""
from .transformable import Transformable
from .atom import Atom, real_coords, frac_coords, atomic_number
from .lattice import Lattice, lattice_vectors_from_parameters
from .structure import AtomicStructure, CompositeStructure, Structure
from .pdb_parser import PDBParser
from .cif_parser import CIFParser
from .crystal import Crystal

###########################
# Graphite built-in Crystal
import numpy as np
a = 1.42; c = 6.714
a1 = a * np.array([3/2, np.sqrt(3)/2, 0])
a2 = a * np.array([-3/2,np.sqrt(3)/2, 0])
a3 = c*np.array([0,0,1])
lattice_vectors = [a1,a2,a3]

#Basis atom positions
r1 = np.array([0,0,0])
r2 = a * np.array([1/2, np.sqrt(3)/2, 0])
r3 = c * np.array([0,0,1/2])
r4 = np.array([-a/2, a*np.sqrt(3)/2, c/2])
            
#Generate unit cell as a list of Atoms
unitcell = list()
for coordinates in (r1,r2,r3,r4):
    unitcell.append(Atom(element = 'C', coords = coordinates))

graphite = Crystal(atoms = unitcell, lattice_vectors = lattice_vectors)