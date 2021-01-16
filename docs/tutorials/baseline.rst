.. include:: ../references.txt

.. _baseline_tutorial:

.. currentmodule:: skued

**********************
Baseline-determination
**********************

Contents
========

* :ref:`dwt_baseline`
* :ref:`dual_tree_baseline`

Due to the high electron cross-section, background signals (or baseline) are
much more of a problem for electron diffraction than equivalent X-ray experiments.

In the case of polycrystalline samples (i.e. 1D diffraction signals),
the definite way of removing background signals is to use an iterative
approach based on the dual-tree complex wavelet transform. Consider the following
example polycrystalline vanadium dioxide pattern:

.. plot::

	import matplotlib.pyplot as plt
	import numpy as np

	s, intensity = np.load('data/powder.npy')
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(s, intensity, 'k')
	ax.set_xlabel('s ($4 \pi / \AA$)')
	ax.set_ylabel('Diffracted intensity (counts)')
	ax.set_title('Background-subtracted diffraction pattern of rutile VO$_2$')
	plt.show()

It would be very difficult to interpolate a background from elastic scattering-free regions
of the diffraction pattern, and this is a moderately We can add a background typical of silicon nitride 
substrates, as well as inelastic scattering effects:

	>>> from skued import gaussian
	>>> import numpy as np
	>>> 
	>>> s, intensity = np.load("docs/tutorials/data/powder.npy")
	>>>
	>>> background = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	>>> substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	>>> substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

.. plot:: 

	import matplotlib.pyplot as plt
	import numpy as np
	from skued import gaussian

	s, intensity = np.load('data/powder.npy')

	# Double exponential inelastic background and substrate effects
	background = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

	plt.plot(s, intensity + background + substrate1 + substrate2, 'k-',
			label = 'Diffraction')
	plt.plot(s, background + substrate1 + substrate2, 'r-', 
			label = 'Baseline')

	plt.legend()

	plt.title('Diffraction pattern of rutile VO$_2$')
	plt.xlabel('s ($4 \pi / \AA$)')
	plt.ylabel('Diffracted intensity (counts)')
	plt.xlim([s.min(), s.max()])
	plt.ylim([0, 100])
	plt.show()

Scikit-ued offers two ways of removing the background.

.. _dwt_baseline:

Iterative Baseline Determination using the Discrete Wavelet Transform
=====================================================================

The procedure and rational for the :func:`baseline_dwt` routine is described in detail in:

Galloway et al. 'An Iterative Algorithm for Background Removal in Spectroscopy by Wavelet 
Transforms', Applied Spectroscopy pp. 1370 - 1376, September 2009.

Here is a usage example for the data presented above:

	>>> import numpy as np
	>>> from skued import gaussian
	>>> from skued import baseline_dwt
	>>> 
	>>> s, intensity = np.load("docs/tutorials/data/powder.npy")
	>>> 
	>>> # Double exponential inelastic background and substrate effects
	>>> diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	>>> substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	>>> substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)
	>>>
	>>> baseline = baseline_dwt(s, level = 6, max_iter = 150, wavelet = 'sym6')

.. plot::

	import matplotlib.pyplot as plt
	import numpy as np
	from skued import gaussian, spectrum_colors
	from skued import baseline_dwt

	s, intensity = np.load('data/powder.npy')

	# Double exponential inelastic background and substrate effects
	diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

	signal = intensity + diffuse + substrate1 + substrate2
	levels = list(range(1,7))
	colors = spectrum_colors(levels)

	fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = (9,3))

	ax1.plot(s, signal, 'k-', label = 'Diffraction')
	ax1.plot(s, diffuse + substrate1 + substrate2, 'r', label = 'Background')
	ax1.set_title('Diffraction pattern of rutile VO$_2$')

	for l, c in zip(levels, colors):
		baseline = baseline_dwt(signal, level = l, max_iter = 150, wavelet = 'sym6')
		ax2.plot(s, baseline, color = c, label = f'Level {l}')
	ax2.set_title('Baseline examples (DWT)')

	for ax in (ax1, ax2):
		ax.set_xlabel('s ($4 \pi / \AA$)')
		ax.set_ylabel('Diffracted intensity (counts)')
		ax.set_xlim([0.2, 0.5])
		ax.set_ylim([30, 80])
		ax.legend()
	plt.show()

.. _dual_tree_baseline:

Iterative Baseline Determination using the Dual-Tree Complex Wavelet Transform
==============================================================================

In the case of 1D data (or along a 1D axis), there is a more performant alternative to :func:`baseline_dwt`. The 
**dual-tree complex wavelet transform** improves on the discrete wavelet transform in many ways.
Therefore, the method presented in this section should be preferred.

For more information on the why and how, please refer to:

L. P. René de Cotret and B. J. Siwick, A general method for baseline-removal in ultrafast 
electron powder diffraction data using the dual-tree complex wavelet transform, Struct. Dyn. 4 (2017)

Here is a usage example for the data presented above:

	>>> import numpy as np
	>>> from skued import gaussian
	>>> from skued import baseline_dt
	>>>
	>>> s, intensity = np.load("docs/tutorials/data/powder.npy")
	>>> 
	>>> # Double exponential inelastic background and substrate effects
	>>> diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	>>> substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	>>> substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)
	>>> 
	>>> baseline = baseline_dt(s, wavelet = 'qshift3', level = 6, max_iter = 150)

.. plot::

	import matplotlib.pyplot as plt
	import numpy as np
	from skued import gaussian, spectrum_colors
	from skued import baseline_dt

	s, intensity = np.load('data/powder.npy')

	# Double exponential inelastic background and substrate effects
	diffuse = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
	substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
	substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

	signal = intensity + diffuse + substrate1 + substrate2
	levels = list(range(1,7))
	colors = spectrum_colors(levels)

	fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = (9,3))

	ax1.plot(s, signal, 'k-', label = 'Diffraction')
	ax1.plot(s, diffuse + substrate1 + substrate2, 'r', label = 'Background')
	ax1.set_title('Diffraction pattern of rutile VO$_2$')

	for l, c in zip(levels, colors):
		baseline = baseline_dt(signal, level = l, max_iter = 150, wavelet = 'qshift3')
		ax2.plot(s, baseline, color = c, label = f'Level {l}')
	ax2.set_title('Baseline examples (DTCWT)')

	for ax in (ax1, ax2):
		ax.set_xlabel('s ($4 \pi / \AA$)')
		ax.set_ylabel('Diffracted intensity (counts)')
		ax.set_xlim([0.2, 0.5])
		ax.set_ylim([30, 80])
		ax.legend()
	plt.show()

The :func:`baseline_dt` routine will usually be more accurate than its :func:`baseline_dwt` counterpart.
However, :func:`baseline_dwt` can be applied to 1D and 2D data.

:ref:`Return to Top <baseline_tutorial>`