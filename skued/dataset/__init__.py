# -*- coding: utf-8 -*-
""" Dataset interactions """

from .meta import ExperimentalParameter, MetaRawDataset
from .raw import AbstractRawDataset
from .mcgill import McGillRawDataset

RAW_DATASET_CLASSES = { McGillRawDataset }