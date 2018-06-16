# -*- coding: utf-8 -*-
import numpy as np
import unittest
from .. import lorentz, electron_wavelength, interaction_parameter

class TestLorentz(unittest.TestCase):
    
    def test_trivial(self):
        """ Test that the lorentz factor for 0 kV is unity """
        self.assertEqual(lorentz(0), 1)
    
    def test_vectorized(self):
        """ Test lorentz() on an array of energies """
        kV = np.zeros((128,), dtype = np.float)
        self.assertTrue(np.allclose(lorentz(kV), np.ones_like(kV)))
    
    def test_range(self):
        """ Test that lorentz factor is always in the range [1, infty) """
        kv = np.linspace(0, 1e6, num = 256)
        factors = lorentz(kv)
        self.assertTrue(np.all(np.greater_equal(factors, 1)))

class TestElectronWavelength(unittest.TestCase):
    
    def test_trivial(self):
        """ Test that the electron wavelength at zero energy is zero """
        self.assertAlmostEqual(electron_wavelength(10), 0.122, places = 3)
        self.assertAlmostEqual(electron_wavelength(200), 0.0250, places = 3)

class TestInteractionParameter(unittest.TestCase):
    
    def test_100kV(self):
        """ Test that the interaction_parameter(100) is what is expected
        from Kirkland 2010 """
        self.assertAlmostEqual(interaction_parameter(100), 0.924*1e-3, places = 6)
    
    def test_vectorized(self):
        """ Test that interaction_parameters() is vectorized """
        kV = np.full((64,), fill_value = 100.0)
        i = np.full_like(kV, fill_value = 0.924*1e-3)
        self.assertTrue(np.allclose(interaction_parameter(kV), i, atol = 1e-6))
    
if __name__ == '__main__':
    unittest.main()