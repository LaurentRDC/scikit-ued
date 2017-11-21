# -*- coding: utf-8 -*-
""" Dataset interactions """

from .raw import AbstractRawDataset, ExperimentalParameter
from .mcgill import McGillRawDataset

RAW_DATASET_CLASSES = { McGillRawDataset }