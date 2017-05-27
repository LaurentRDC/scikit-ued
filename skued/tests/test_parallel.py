# -*- coding: utf-8 -*-
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

	def test_trivial_map_no_args(self):
		""" Test that pmap is working with no positional arguments """
		integers = list(range(0,10))
		result = list(pmap(identity, integers))
		self.assertEqual(integers, result)
	
	def test_trivial_map_kwargs(self):
		""" Test that pmap is working with args and kwargs """
		integers = list(range(0,10))
		result = list(pmap(identity, integers, kwargs = {'test' : True}))
		self.assertEqual(result, integers)

if __name__ == '__main__':
    unittest.main()