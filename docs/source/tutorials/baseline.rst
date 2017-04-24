.. include:: ../references.txt

.. _baseline_tutorial:

************************
Baseline Tutorials
************************

Contents
========

* :ref:`dual_tree_baseline`

Due to the high electron cross-section, background signals (or baseline) are
much more of a problem for electron diffraction than equivalent X-ray experiments.

In the case of polycrystalline samples (i.e. 1D diffraction signals),
the definite way of removing background signals is to use an iterative
approach based on the dual-tree complex wavelet transform.

First, we load an example dataset::

	import matplotlib.pyplot as plt
	import numpy as np

	s, intensity = np.load('powder.py')

.. plot :: tutorials/plots/powder_plot.py

We can add a background typical of silicon nitride substrates, as well
as inelastic scattering effects::

	from skued import gaussian
	
	background = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

.. plot:: tutorials/plots/powder_w_background.py

Scikit-ued offers two ways of removing the background.

.. _dual_tree_baseline:

Iterative Baseline Determination using the Dual-Tree Complex Wavelet Transform
==============================================================================

.. plot:: tutorials/plots/powder_dt_baseline_progression.py

:ref:`Return to Top <baseline_tutorial>`