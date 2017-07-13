# -*- coding: utf-8 -*-
import unittest
from itertools import repeat
from .. import last, chunked, linspace, multilinspace

class TestLast(unittest.TestCase):

    def test_trivial(self):
        """ Test last() on iterable of identical values """
        i = repeat(1, 10)
        self.assertEqual(last(i), 1)
    
    def test_on_empty_iterable(self):
        """ Test that last() return None for empty iterable """
        i = list()
        self.assertIs(last(i), None)

class TestLinspace(unittest.TestCase):

    def test_endpoint(self):
        """ Test that the endpoint is included by linspace() when appropriate"""
        with self.subTest('endpoint = True'):
            space = linspace(0, 1, num = 10, endpoint = True)
            self.assertEqual(last(space), 1)

        with self.subTest('endpoint = False'):
            space = linspace(0, 1, num = 10, endpoint = False)
            self.assertAlmostEqual(last(space), 0.9)
    
    def test_length(self):
        """ Test that linspace() returns an iterable of the correct length """
        with self.subTest('endpoint = True'):
            space = list(linspace(0, 1, num = 13, endpoint = True))
            self.assertEqual(len(space), 13)

        with self.subTest('endpoint = False'):
            space = list(linspace(0, 1, num = 13, endpoint = False))
            self.assertEqual(len(space), 13)

class TestMultilinspace(unittest.TestCase):

    def test_endpoint(self):
        """ Test that the endpoint is included by linspace() when appropriate"""
        with self.subTest('endpoint = True'):
            space = multilinspace((0,0), (1,1), num = 10, endpoint = True)
            self.assertSequenceEqual(last(space), (1,1))

        with self.subTest('endpoint = False'):
            space = multilinspace((0,0), (1,1), num = 10, endpoint = False)
            # Unfortunately there is no assertSequenceAlmostEqual
            self.assertSequenceEqual(last(space), (0.8999999999999999, 0.8999999999999999))
    
    def test_length(self):
        """ Test that linspace() returns an iterable of the correct length """
        with self.subTest('endpoint = True'):
            space = list(multilinspace((0,0), (1,1), num = 13, endpoint = True))
            self.assertEqual(len(space), 13)

        with self.subTest('endpoint = False'):
            space = list(multilinspace((0,0), (1,1), num = 13, endpoint = False))
            self.assertEqual(len(space), 13)

class TestChunked(unittest.TestCase):
    
    def test_larger_chunksize(self):
        """ Test chunked() with a chunksize larger that the iterable itself """
        i = repeat(1, 10)
        chunks = chunked(i, chunksize = 15)
        self.assertEqual(len(list(chunks)), 1)  # One single chunk is returned

    def test_on_infinite_generator(self):
        """ Test chunked() on an infinite iterable """
        i = repeat(1)
        chunks = chunked(i, chunksize = 15)
        for _ in range(10):
            self.assertEqual(len(next(chunks)), 15)
		
if __name__ == '__main__':
	unittest.main()