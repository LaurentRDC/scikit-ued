# -*- coding: utf-8 -*-
"""
Baseline-determination algorithms based on the wavelet tranform.

Functions
---------
baseline_dt
    Baseline determination of signals using the dual-tree complex wavelet transform. Available in 1D or along an axis. [1, 2]

baseline_dwt
    Baseline determination of signals using the discrete wavelet transform. Provided for comparison
    with the dual-tree equivalent 'baseline'. Modified algorithm from [3]. For better performance in 1D,
	use baseline_dt.

References
----------
[1] L. P. René de Cotret and B. J. Siwick, 'A general method for baseline-removal in ultrafast electron powder diffraction data 
	using the dual-tree complex wavelet transform,' Struct. Dyn. 4 (2016)

[2] Selesnick, I. W. et al. 'The Dual-tree Complex Wavelet Transform', IEEE Signal Processing Magazine pp. 123 - 151, November 2005.

[3] Galloway et al. 'An Iterative Algorithm for Background Removal in Spectroscopy by Wavelet Transforms', 
	Applied Spectroscopy pp. 1370 - 1376, September 2009.
"""
from .algorithms import baseline_dt, baseline_dwt
from .dtcwt import (ALL_COMPLEX_WAV, ALL_FIRST_STAGE, available_dt_filters,
                    available_first_stage_filters, dt_max_level, dtcwt, idtcwt)
