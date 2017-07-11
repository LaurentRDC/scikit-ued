# -*- coding: utf-8 -*-
import unittest
from .. import spectrum_colors, rgb_sweep

class TestSpectrumColors(unittest.TestCase):

	def test_on_ints(self):
		""" Test spectrum_colors on an int """
		colors = spectrum_colors(10)
		self.assertEqual(len(list(colors)), 10)

	def test_on_sized_iterable(self):
		""" Test on iterable that has a __len__ attribute: list, tuple, etc. """
		colors = spectrum_colors( [1,2,3,4,5] )
		self.assertEqual(len(list(colors)), 5)

	def test_on_unsized_iterable(self):
		""" Test spectrum_colors on unsized_iterable (e.g. generator) """
		colors = spectrum_colors( range(0, 10) )
		self.assertEqual(len(list(colors)), 10)
	
	def test_on_length_1_iterable(self):
		""" Test that spectrum_colors is working on single-length iterables """
		self.assertSequenceEqual(list(spectrum_colors(1)),
								 list(spectrum_colors([0])))
								
class TestRGBSweep(unittest.TestCase):

	def test_on_ints(self):
		""" Test the number of elements yielded from rgb_sweep """
		colors = rgb_sweep(10, source = (1,0,0), dest = (0,1,0))
		self.assertEqual(len(list(colors)), 10)
	
	def test_source_equal_to_dest(self):
		""" Test that rgb_sweep still works if source is equal to destination."""
		colors = rgb_sweep(10, source = (1,0,0), dest = (1,0,0))
		self.assertEqual(len(list(colors)), 10)
	
	def test_hex_colors(self):
		""" Test that rgb_sweep works on hexadecimal strings """
		colors = list(rgb_sweep(10, source = '#ffffff', dest = '#000000'))

		self.assertSequenceEqual(colors[0], (1.0, 1.0, 1.0))
		self.assertSequenceEqual(colors[-1], (0.0, 0.0, 0.0))
		
if __name__ == '__main__':
	unittest.main()