# -*- coding: utf-8 -*-
"""
Time-series analysis
--------------------
This package allows for exploration of time-series data, especially
in the context of ultrafast diffraction.
"""

from .fitting import biexponential, exponential
from .nfft_routines import nfft, nfftfreq
from .robust import mad
from .selections import (
    ArbitrarySelection,
    DiskSelection,
    RectSelection,
    RingSelection,
    Selection,
)
from .time_zero import register_time_shift, register_time_shifts
