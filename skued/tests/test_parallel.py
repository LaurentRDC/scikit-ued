
from ..parallel import pmap, preduce
from functools import reduce
import numpy as np
from operator import add
import unittest

def identity(obj, *args, **kwargs):
    """ ignores args and kwargs """
    return obj

class TestParallelReduce(unittest.TestCase):
	
	def test_preduce_one_process(self):
		""" Test that preduce reduces to functools.reduce for a single process """
		integers = list(range(0, 10))
		preduce_results = preduce(add, integers, processes = 1)
		reduce_results = reduce(add, integers)

		self.assertEqual(preduce_results, reduce_results)

	def test_preduce_multiple_processes(self):
		""" Test that preduce reduces to functools.reduce for a single process """
		integers = list(range(0, 10))
		preduce_results = preduce(add, integers, processes = 2)
		reduce_results = reduce(add, integers)

		self.assertEqual(preduce_results, reduce_results)

	def test_on_numpy_arrays(self):
		""" Test sum of numpy arrays as parallel reduce"""
		arrays = [np.zeros((32,32)) for _ in range(10)]
		s = preduce(add, arrays, processes = 2)

		self.assertTrue(np.allclose(s, arrays[0]))

	def test_with_kwargs(self):
		""" Test preduce with keyword-arguments """
		pass

class TestParallelMap(unittest.TestCase):

	def test_map_one_process(self):
		""" Tests that pmap reduces to map for a single process."""
		# Since map returns an iterator of type MapObject,
		# we can check that pmap reduces to map by checking
		# the return type.

		integers = list(range(0, 10))
		pmap_result = pmap(identity, integers, processes = 1)
		map_result = map(identity, integers)

		self.assertIsInstance(pmap_result, type(map_result))

	def test_trivial_map_no_args(self):
		""" Test that pmap is working with no positional arguments """
		integers = list(range(0,10))
		result = pmap(identity, integers)
		self.assertEqual(integers, result)
	
	def test_trivial_map_with_args_and_kwargs(self):
		""" Test that pmap is working with args and kwargs """
		integers = list(range(0,10))
		result = pmap(identity, integers, args = (1,), kwargs = {'test' : True})
		self.assertEqual(result, integers)

if __name__ == '__main__':
    unittest.main()