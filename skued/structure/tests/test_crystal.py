
from copy import deepcopy, copy
from itertools import permutations
import numpy as np
from .. import Crystal, Atom, Lattice, graphite
from ... import rotation_matrix, transform
import unittest

# large test crystal
br = Crystal.from_pdb('1fbb')

def vect_angle(a1, a2):
    """ Angle (in rads) between vectors """
    return np.arccos(np.vdot(a1,a2)/(np.linalg.noam(a1)*np.linalg.norm(a2)))

def equal_crystal_transformation(cryst1, cryst2):
    """
    Returns True if lattice vectors, symmetry operators, and atomic positions are 
    the same. 
    
    Does not check atomic composition.
    """
    lv_equal = np.allclose(np.array(cryst1.lattice_vectors), np.array(cryst2.lattice_vectors))
    sym_ops_equal = all([np.allclose(sym1, sym2) for sym1, sym2 in zip(cryst1.symmetry_operators, cryst2.symmetry_operators)])
    atoms_places_equal = all([np.allclose(atm1.coords, atm2.coords) for atm1, atm2 in zip(cryst1.atoms, cryst2.atoms)]) 
    unitcell_places_equal = all([np.allclose(atm1.coords, atm2.coords) for atm1, atm2 in zip(cryst1.unitcell, cryst2.unitcell)])

    return lv_equal and sym_ops_equal and atoms_places_equal and unitcell_places_equal

class TestBoundedReflections(unittest.TestCase):

	def setUp(self):
		self.crystal = deepcopy(graphite)

	def test_bounded_reflections_negative(self):
		""" Test that negative reflection bounds raise an Exception.
		Otherwise, an infinite number of reflections will be generated """
		with self.assertRaises(ValueError):
			hkl = list(self.crystal.bounded_reflections(-1))
	
	def test_bounded_reflections_zero(self):
		""" Check that bounded_reflections returns (000) for a zero bound """
		h, k, l = self.crystal.bounded_reflections(0)
		[self.assertEqual(len(i), 1) for i in (h, k, l)]
		[self.assertEqual(i[0], 0) for i in (h, k, l)]
	
	def test_bounded_reflections_all_within_bounds(self):
		""" Check that every reflection is within the bound """
		bound = 10
		Gx, Gy, Gz = self.crystal.scattering_vector(*self.crystal.bounded_reflections(nG = bound))
		norm_G = np.sqrt(Gx**2 + Gy**2 + Gz**2)
		self.assertTrue(np.all(norm_G <= bound))

class TestCrystalRotations(unittest.TestCase):

    def setUp(self):
        self.crystal = deepcopy(graphite)
    
    def test_crystal_equality(self):
        """ Tests that the function 'equal_crystal_transformation' is working properly """
        self.assertTrue(equal_crystal_transformation(self.crystal, self.crystal))

        cryst2 = deepcopy(self.crystal)
        cryst2.transform(2*np.eye(3)) # This stretches lattice vectors, symmetry operators
        self.assertFalse(self.crystal is cryst2)
        self.assertFalse(equal_crystal_transformation(self.crystal, cryst2))

        cryst2.transform(0.5*np.eye(3))
        self.assertTrue(equal_crystal_transformation(self.crystal, cryst2))
    
    def test_trivial_rotation(self):
        """ Test rotation by 360 deg around all axes. """
        unrotated = deepcopy(self.crystal)
        self.crystal.rotate(360, [1,0,0])
        self.assertTrue(equal_crystal_transformation(self.crystal, unrotated))
    
    def test_identity_transform(self):
        """ Tests the trivial identity transform """
        transf = deepcopy(self.crystal)
        transf.transform(np.eye(3))
        self.assertTrue(equal_crystal_transformation(self.crystal, transf))
    
    def test_one_axis_rotation(self):
        """ Tests the crystal orientation after rotations. """
        orient = self.crystal.unitcell[1].coords - self.crystal.unitcell[0].coords
        self.crystal.rotate(45, [0,1,0])
        r_orient = self.crystal.unitcell[1].coords - self.crystal.unitcell[0].coords
        self.assertFalse(np.allclose(orient, r_orient))
        self.crystal.rotate(-45, [0,1,0])
        r_orient_2 = self.crystal.unitcell[1].coords - self.crystal.unitcell[0].coords
        self.assertTrue(np.allclose(orient, r_orient_2))

    def test_wraparound_rotation(self):
        cryst1 = deepcopy(self.crystal)
        cryst2 = deepcopy(self.crystal)

        cryst1.rotate(22.3, [0,0,1])
        cryst2.rotate(-(360 - 22.3), [0,0,1])
        self.assertTrue(equal_crystal_transformation(cryst1, cryst2))

if __name__ == '__main__':
    unittest.main()