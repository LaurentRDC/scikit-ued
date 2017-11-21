# -*- coding: utf-8 -*-
from .. import McGillRawDataset, AbstractRawDataset, ExperimentalParameter
import unittest

class TestRawDataset(AbstractRawDataset):

    test = ExperimentalParameter('test', int, default = 0)
    readonly = ExperimentalParameter('readonly', str, default = '', readonly = True)

    def raw_data(self): return None

class TestAbstractRawDataset(unittest.TestCase):

    def test_abstract_methods(self):
        """ Test that instantiation of AbstractRawDataset 
        raises an error """
        with self.assertRaises(TypeError):
            AbstractRawDataset('')
    
    def test_minimal_methods(self):
        """ 
        Test implementing the minimal methods:

        * raw_data
        """
        TestRawDataset()
    
    def test_experimental_parameters(self):
        """ Test the behavior of the ExperimentalParameter descriptor """
        
        test_dataset = TestRawDataset()

        with self.subTest('Default value'):
            self.assertEqual(test_dataset.test, 0)
        
        with self.subTest('Changing value'):
            test_dataset.test = 1
            self.assertEqual(test_dataset.test, 1)
        
        with self.subTest('Parameter type'):
            with self.assertRaises(TypeError):
                test_dataset.test = 'test'
        
        with self.subTest('Read-only'):
            with self.assertRaises(AttributeError):
                test_dataset.readonly = 'something else'
    
    def test_valid_metadata(self):
        """ Test that the class attribute 'valid_metadata' is working as intended """
        
        self.assertIn('test', TestRawDataset.valid_metadata)
        self.assertIn('readonly', TestRawDataset.valid_metadata)
        self.assertLessEqual(AbstractRawDataset.valid_metadata, TestRawDataset.valid_metadata)