# -*- coding: utf-8 -*-
"""
I/O package
-----------
This package provdes utility functions for interfacing with less common file formats.
"""
from .io import diffread, diffshow
from .merlin import mibheader, mibread, imibread
from .dm import dmread