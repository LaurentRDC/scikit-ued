# -*- coding: utf-8 -*-
"""
I/O package
-----------
This package provdes utility functions for interfacing with less common file formats.
"""
from .dm import dmread
from .diffshow import diffshow, WITH_PYQTGRAPH
from .io import diffread
from .merlin import imibread, mibheader, mibread
