
from ..parallel import pmap
import unittest

def identity(obj, *args, **kwargs):
    """ ignores args and kwargs """
    return obj

class ParallelMapTest(unittest.TestCase):

	def test_map_one_process(self):
		""" Tests that pmap reduces to map for a single process."""
		# Since map returns an iterator of type MapObject,
		# we can check that pmap reduces to map by checking
		# the return type.

		integers = list(range(0, 10))
		pmap_result = pmap(identity, integers, processes = 1)
		map_result = map(identity, integers)

		self.assertTrue(isinstance(pmap_result, type(map_result)))

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