.. include:: ../references.txt

.. _image_analysis_tutorial:

.. currentmodule:: skued

*************************
Image Analysis/Processing
*************************

Diffraction patterns analysis is essentially specialized image processing. This tutorial
will show some of the image processing and analysis techniques that are part of the scikit-ued.

.. note::
    Use scikit-ued in combination with `npstreams`_ to process electron diffraction data in parallel.

Contents
========

* :ref:`io`
* :ref:`alignment`
* :ref:`symmetry`
* :ref:`pixel_masks`
* :ref:`powder`
* :ref:`denoising`

.. _io:

Reading Images 
==============

Diffraction patterns can come in a variety of exotic file formats. Scikit-ued has built-in support for the following file formats:

* Gatan's closed source DM3 and DM4 (`*.dm3`, `*.dm4`);
* Merlin Image Binary (`*.mib`);
* TIFF images (`*.tif`, `*.tiff`);
* All other file formats supported by `scikit-image`_.

The :func:`diffread` function will transparently distinguish between those formats and dispatch to the right functions. 

.. _alignment:

Automatic center-finding
========================

Many analyses involving diffraction patterns require the knowledge of the center. For this purpose, ``scikit-ued``
provides :func:`autocenter`. It can be used trivially as follows:

	>>> from skued import diffread, autocenter
	>>> import numpy as np
	>>> 
	>>> im = diffread('docs/tutorials/Cr_1.tif')
	>>>
	>>> # Invalid pixels are masked with a False
	>>> mask = np.ones_like(im, dtype = np.bool)
	>>> mask[0:1250, 975:1225] = False
	>>>
	>>> center = autocenter(im, mask=mask)

Let's take a look at the result. The center is shown with a red dot:

.. plot::

	from skued import diffread, autocenter
	import matplotlib.pyplot as plt

	im = diffread('Cr_1.tif')

	mask = np.ones_like(im, dtype = np.bool)
	mask[0:1250, 975:1225] = False

	# Reduce size of images because of memory usage of ReadTheDocs
	im = im[::3, ::3]
	mask = mask[::3, ::3]

	rc, cc = autocenter(im, mask=mask)

	fig, ax1 = plt.subplots(figsize = (4.5, 4.5))
	ax1.imshow(im, vmin = 0, vmax = 200, cmap='inferno')
	ax1.scatter(cc, rc, color='r')

	ax1.get_xaxis().set_visible(False)
	ax1.get_yaxis().set_visible(False)

	plt.tight_layout()
	plt.show()

We can do the same with single-crystal diffraction patterns:

	>>> from skued import diffread, autocenter
	>>> import numpy as np
	>>> 
	>>> im = diffread('docs/tutorials/graphite.tif')
	>>> mask = diffread('docs/tutorials/graphite_mask.tif').astype(np.bool)
	>>>
	>>> center = autocenter(im, mask=mask)

Let's take a look at the result. The center is shown with a red dot:

.. plot::

	from skued import diffread, autocenter
	import matplotlib.pyplot as plt

	im = diffread('graphite.tif')
	mask = diffread('graphite_mask.tif').astype(np.bool)

	# Reduce size of images because of memory usage of ReadTheDocs
	im = im[::3, ::3]
	mask = mask[::3, ::3]

	rc, cc = autocenter(im, mask=mask)

	fig, ax1 = plt.subplots(figsize = (4.5, 4.5))
	ax1.imshow(im, vmin = 0, vmax = 200, cmap='inferno')
	ax1.scatter(cc, rc, color='r')

	ax1.get_xaxis().set_visible(False)
	ax1.get_yaxis().set_visible(False)

	plt.tight_layout()
	plt.show()

It's worth emphasizing that apart from the mask, there are no parameters to play with; :func:`autocenter` is as automatic as it can be!

Diffraction pattern alignment
=============================

Diffraction patterns can drift over a period of a few minutes, and for reliable data synthesis
it is important to align patterns to a reference.

Diffraction patterns all have a fixed feature: the position of the beam-block. Therefore, some pixels 
in each diffraction pattern must be ignored in the computation of the cross-correlation. 

Setting the 'invalid pixels' to 0 will not work, at those will correlate with the invalid pixels from the reference. One must use
the **masked normalized cross-correlation**.

All of this is taken care of with the :func:`align` function. Let's look at some polycrystalline Chromium:

.. plot::

	from skued import diffread
	import matplotlib.pyplot as plt

	ref = diffread('Cr_1.tif')
	im = diffread('Cr_2.tif')

	fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize = (9,3))
	ax1.imshow(ref, vmin = 0, vmax = 200, cmap='inferno')
	ax2.imshow(im, vmin = 0, vmax = 200, cmap='inferno')
	ax3.imshow(ref - im, cmap = 'RdBu_r')

	for ax in (ax1, ax2, ax3):
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)

	ax1.set_title('Reference')
	ax2.set_title('Data')
	ax3.set_title('Difference')

	plt.tight_layout()
	plt.show()

From the difference pattern, we can see that the 'Data' pattern is shifted from 'Reference' quite a bit, 
but the beamblock **has not moved**. To determine the exact shift, we need to use a mask that obscures the 
beam-block and main beam:

	>>> from skued import diffread, align
	>>> import numpy as np
	>>> 
	>>> ref = diffread('docs/tutorials/Cr_1.tif')
	>>> im = diffread('docs/tutorials/Cr_2.tif')
	>>>
	>>> # Invalid pixels are masked with a False
	>>> mask = np.ones_like(ref, dtype = np.bool)
	>>> mask[0:1250, 975:1225] = False
	>>>
	>>> aligned = align(image=im, reference=ref, mask=mask)

Let's look at the difference:

.. plot::

	from skued import diffread, align
	from pathlib import Path

	ref = diffread("Cr_1.tif")
	im = diffread("Cr_2.tif")

	mask = np.ones_like(ref, dtype=np.bool)
	mask[0:1250, 975:1225] = False

	# Reduce size of images because of memory usage of ReadTheDocs
	mask = mask[::3, ::3]
	ref = ref[::3, ::3]
	im = im[::3, ::3]

	aligned = align(image=im, reference=ref, mask=mask)

	fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(9, 6))
	ax1.imshow(ref, vmin=0, vmax=200, cmap='inferno')
	ax2.imshow(im, vmin=0, vmax=200, cmap='inferno')
	ax3.imshow(ref - im, cmap="RdBu_r")
	ax4.imshow(np.logical_not(mask) * ref, vmin=0, vmax=200, cmap="inferno")
	ax5.imshow(aligned, vmin=0, vmax=200, cmap='inferno')
	ax6.imshow(ref - aligned, cmap="RdBu_r")

	for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)

	ax1.set_title("Reference")
	ax2.set_title("Data")
	ax3.set_title("Difference")
	ax4.set_title("Mask")
	ax5.set_title("Aligned data")
	ax6.set_title("Difference after alignment")

	plt.tight_layout()

.. _symmetry:

Image processing involving symmetry
===================================

Rotational symmetry
-------------------
Diffraction patterns exhibit rotational symmetry based on the crystal structure. We can
take advantage of such symmetry to correct images in case of artifacts or defects. A useful
routine is :func:`nfold`, which averages portions of a diffraction pattern with itself based on
rotational symmetry.

.. plot::

	import matplotlib.pyplot as plt
	from skued import nfold, diffread, autocenter
	import numpy as np

	image = diffread('graphite.tif')

	mask = np.ones_like(image, dtype = np.bool)
	mask[1100::, 442:480] = False # Artifact line
	mask[0:1260, 900:1140] = False # beamblock

	# Reduce size of images because of memory usage of ReadTheDocs
	image = image[::3, ::3]
	mask = mask[::3, ::3]

	yc, xc = autocenter(image, mask=mask)
	av = nfold(image, mod = 6, center = (xc, yc), mask = mask)

	fig , (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (9,3))
	ax1.imshow(image, vmin = 0, vmax = 150, cmap='inferno')
	ax2.imshow(np.logical_not(mask) * image, vmin = 0, vmax = 150, cmap='inferno')
	ax3.imshow(av, vmin = 0, vmax = 150, cmap='inferno')

	for ax in (ax1, ax2, ax3):
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)

	ax1.set_title('Graphite')
	ax2.set_title('Mask')
	ax3.set_title('Averaged')

	plt.tight_layout()
	plt.show()

To use :func:`nfold`, all you need to know is the center of the diffraction pattern:

    >>> from skued import nfold, diffread
	>>> 
    >>> im = diffread('docs/tutorials/graphite.tif')
    >>> av = nfold(im, mod = 6, center = (1024, 1024))    # mask is optional


.. _pixel_masks:

Pixel Masks
===========

Image data can be rejected on a per-pixel basis by using pixel masks. These masks are represented
by boolean arrays that evaluate to ``Falose`` on invalid pixels, and ``True`` otherwise

:mod:`scikit-ued` offers some functions related to creation and manipulation of pixel masks.

Creation of a pixel mask
------------------------

A pixel mask can be created from a set of images sharing the same properties. For example, diffraction patterns
before photoexcitation (i.e. dark runs) form a set of images that should be identical.

Let's imaging a set of such images with filenames `dark_run_*.tif`. We can create a pixel mask with the :func:`mask_from_collection`:

    >>> from glob import iglob
    >>> from skued import mask_from_collection, diffread
	>>> 
    >>> dark_runs = map(diffread, iglob('dark_run_*.tif')) # doctest: +SKIP
    >>> mask = mask_from_collection(dark_runs) # doctest: +SKIP

In the above example, pixel values outside opf the [0, 30000] range will be marked as invalid (default behaviour). Moreover,
the per-pixel standard deviation over the image set is computed; pixels that fluctuate too much are also rejected.

Note that since :func:`mask_from_collection` uses :mod:`npstreams` under the hood, the collection of files used to compute the 
mask can be huge, much larger than available memory.

.. _powder:

Image analysis on polycrystalline diffraction patterns
======================================================

Angular average
---------------

First, let's load an image:

	>>> import numpy as np
	>>> import matplotlib.pyplot as plt
	>>> from skued import diffread, autocenter
	>>> 
	>>> image = diffread("docs/tutorials/Cr_1.tif")
	>>> 
	>>> mask = np.ones_like(image, dtype=np.bool)
	>>> mask[0:1150, 1000:1170] = False
	>>> 
	>>> plt.imshow(image)	# doctest: +SKIP
	>>> plt.show() # doctest: +SKIP

.. plot::

	import numpy as np
	import matplotlib.pyplot as plt
	from skued import diffread

	image = diffread("Cr_1.tif")

	mask = np.ones_like(image, dtype=np.bool)
	mask[0:1150, 1000:1170] = False

	fig, (ax1, ax2) = plt.subplots(1,2)
	ax1.imshow(image, vmin=0, vmax=200, cmap='inferno')
	ax2.imshow(np.logical_not(mask) * image, vmin=0, vmax=200, cmap='inferno')

	ax1.set_title('Chromium')
	ax2.set_title('Masked area')

	for ax in (ax1, ax2):
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)

	plt.tight_layout()
	plt.show()


... and we can easily compute an angular average:
	
	>>> from skued import azimuthal_average, autocenter
	>>> 
	>>> yc, xc = autocenter(image, mask=mask) # doctest: +SKIP
	>>> radius, intensity = azimuthal_average(image, (xc, yc), mask=mask) # doctest: +SKIP
	>>> 
	>>> plt.plot(radius, intensity) # doctest: +SKIP

.. plot::
	
	from skued import azimuthal_average, diffread, autocenter
	import numpy as np
	import matplotlib.pyplot as plt

	image = diffread("Cr_1.tif")

	mask = np.ones_like(image, dtype=np.bool)
	mask[0:1150, 1000:1170] = False

	# Reduce size of images because of memory usage of ReadTheDocs
	image = image[::3, ::3]
	mask = mask[::3, ::3]

	yc, xc = autocenter(image, mask=mask)
	radius, intensity = azimuthal_average(image, (xc, yc), mask=mask)

	plt.plot(radius, intensity)
	plt.xlabel('Radius [px]')
	plt.ylabel('Diffracted intensity [a.u.]')
	plt.tight_layout()
	plt.show()

.. _denoising:

Bonus : Removing Hot Spots 
==========================

An interesting use-case of baseline-removal (described in :ref:`baseline_tutorial`) is the removal of hot spots from images.

Consider the following diffraction pattern:

.. plot::

	import matplotlib.pyplot as plt
	from skued import diffread

	im = diffread('hotspots.tif')
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.imshow(im, vmin = 0, vmax = 2e3, cmap='inferno')
	ax.xaxis.set_visible(False)
	ax.yaxis.set_visible(False)
	plt.tight_layout()
	plt.show()

We can consider the image *without hotspots* as the baseline of the image *with hotspots* :

	>>> from skued import diffread, baseline_dwt
	>>> 
	>>> im = diffread('docs/tutorials/hotspots.tif')
	>>> denoised = baseline_dwt(im, max_iter = 250, level = 1, wavelet = 'sym2', axis = (0, 1))

The result is plotted below:

.. plot::

	import matplotlib.pyplot as plt
	from skued import diffread, baseline_dwt

	im = diffread('hotspots.tif')
	denoised = baseline_dwt(im, max_iter = 250, level = 1, wavelet = 'sym2', axis = (0, 1))

	fig, (ax1, ax2) = plt.subplots(1, 2)
	ax1.imshow(im, vmin = 0, vmax = 2e3, cmap='inferno')
	ax2.imshow(denoised, vmin = 0, vmax = 2e3, cmap='inferno')

	for ax in (ax1, ax2):
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
	plt.tight_layout()
	plt.show()

Try different combinations of wavelets, levels, and number of iterations (``max_iter``).

Notice that the baseline-removal function used in the case of an image is :func:`baseline_dwt`, which works on 2D arrays.
The same is not possible with :func:`baseline_dt`, which only works on 1D arrays at this time.

:ref:`Return to Top <image_analysis_tutorial>`