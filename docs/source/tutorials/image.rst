.. include:: ../references.txt

.. _image_analysis_tutorial:

**********************************
Image Analysis/Processing Tutorial
**********************************

Diffraction patterns analysis is essentially specialized image processing. This tutorial
will show some of the image processing and analysis techniques that are part of the :mod:`skued.image` module.

.. note::
    A lot of functionality has been moved to the package `npstreams`_. Use scikit-ued in combination
    with `npstreams`_ to process electron diffraction data in parallel.

Contents
========

* :ref:`alignment`
* :ref:`symmetry`
* :ref:`pixel_masks`
* :ref:`powder`
* :ref:`denoising`

.. _alignment:

Diffraction pattern alignment
=============================

Diffraction patterns can drift over a period of a few minutes, and for reliable data synthesis
it is important to align patterns to a reference.

The procedure of detecting, or registering, the translation between two similar images is usually
done by measuring the cross-correlation between images. When images are very similar, this procedure
is fine; take a look at scikit-image's :func:`skimage.feature.register_translation` for example. 

However, diffraction patterns all have a fixed feature: the position of the beam-block. Therefore, some pixels 
in each diffraction pattern must be ignored in the computation of the cross-correlation. 

Setting the 'invalid pixels' to 0 will not work, at those will correlate with the invalid pixels from the reference. One must use
the **masked normalized cross-correlation** through scikit-ued's :func:`mnxc2`.

All of this is taken care of in scikit-ued's :func:`diff_register` function. Let's look at some polycrystalline Chromium:

.. plot::

	from skued import diffread
	import matplotlib.pyplot as plt

	ref = diffread('Cr_1.tif')
	im = diffread('Cr_2.tif')

	fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize = (9,3))
	ax1.imshow(ref, vmin = 0, vmax = 200)
	ax2.imshow(im, vmin = 0, vmax = 200)
	ax3.imshow(ref - im, cmap = 'RdBu_r', vmin = -100, vmax = 100)

	for ax in (ax1, ax2, ax3):
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)

	ax1.set_title('Reference')
	ax2.set_title('Data')
	ax3.set_title('Difference')

	plt.tight_layout()
	plt.show()

From the difference pattern, we can see that the 'Data' pattern is shifted from 'Reference' quite a bit.
To determine the exact shift, we need to use a mask that obscures the beam-block and main beam::

	from skued import diff_register, shift_image, diffread
	import numpy as np

	ref = diffread('Cr_1.tif')
	im = diffread('Cr_2.tif')

	mask = np.zeros_like(ref, dtype = np.bool)
	mask[0:1250, 950:1250] = True

	shift = diff_register(im, reference = ref, mask = mask)
	im = shift_image(im, shift)

Let's look at the difference:

.. plot::

	import matplotlib.pyplot as plt
	import numpy as np
	from skued import diff_register, shift_image, diffread

	ref = diffread('Cr_1.tif')
	im = diffread('Cr_2.tif')

	mask = np.zeros_like(ref, dtype = np.bool)
	mask[0:1250, 950:1250] = True

	shift = diff_register(im, ref, mask)
	shifted = shift_image(im, -shift)

	fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize = (9,6))
	ax1.imshow(ref, vmin = 0, vmax = 200)
	ax2.imshow(im, vmin = 0, vmax = 200)
	ax3.imshow(ref - im, cmap = 'RdBu_r', vmin = -100, vmax = 100)
	ax4.imshow(mask, vmin = 0, vmax = 1, cmap = 'binary')
	ax5.imshow(shifted, vmin = 0, vmax = 200)
	ax6.imshow(ref - shifted, cmap = 'RdBu_r', vmin = -100, vmax = 100)

	for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)

	ax1.set_title('Reference')
	ax2.set_title('Data')
	ax3.set_title('Difference')
	ax4.set_title('Mask')
	ax5.set_title('Aligned data')
	ax6.set_title('Diff. after shift')

	plt.tight_layout()
	plt.show()

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
    from skued import nfold, diffread
    import numpy as np

    center = (1010, 1111)

    mask = np.zeros((2048, 2048), dtype = np.bool)
    mask[1100::, 442:480] = True # Artifact line
    mask[0:1260, 900:1140] = True # beamblock

    image = diffread('graphite.tif')
    av = nfold(image, mod = 6, center = center, mask = mask)

    fig , (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (9,3))
    ax1.imshow(image, vmin = 0, vmax = 150)
    ax2.imshow(mask, vmin = 0, vmax = 1)
    ax3.imshow(av, vmin = 0, vmax = 150)

    for ax in (ax1, ax2, ax3):
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

    ax1.set_title('Graphite')
    ax2.set_title('Mask')
    ax3.set_title('Averaged')

    plt.tight_layout()
    plt.show()

To use :func:`nfold`, all you need to know is the center of the diffraction pattern::

    from skued import nfold, diffread

    im = diffread('graphite.tif')
    av = nfold(im, mod = 6, center = center)    # mask is optional


.. _pixel_masks:

Pixel Masks
===========

Image data can be rejected on a per-pixel basis by using pixel masks. These masks are represented
by boolean arrays that evaluate to ``True`` on invalid pixels.

:mod:`scikit-ued` offers some functions related to creation and manipulation of pixel masks.

Creation of a pixel mask
------------------------

A pixel mask can be created from a set of images sharing the same properties. For example, diffraction patterns
before photoexcitation (i.e. dark runs) form a set of images that should be identical.

Let's imaging a set of such images with filenames `dark_run_*.tif`. We can create a pixel mask with the :func:`mask_from_collection`::

    from glob import iglob
    from skued import mask_from_collection, diffread

    dark_runs = map(diffread, iglob('dark_run_*.tif'))    # Can be a huge stack of images
    mask = mask_from_collection(dark_runs)

In the above example, pixel values outside opf the [0, 30000] range will be marked as invalid (default behaviour). Moreover,
the per-pixel standard deviation over the image set is computed; pixels that fluctuate too much are also rejected.

Note that since :func:`mask_from_collection` uses :mod:`npstreams` under the hood, the collection used to compute the 
mask can be huge.

.. _powder:

Image analysis on polycrystalline diffraction patterns
======================================================

Center-finding
--------------
Polycrystalline diffraction patterns display concentric rings, and finding
the center of those concentric rings is important. Let's load a test image:

.. plot::

	from skued import diffread
    import matplotlib.pyplot as plt
    path = 'Cr_1.tif'

    im = diffread(path)
    mask = np.zeros_like(im, dtype = np.bool)
    mask[0:1250, 950:1250] = True

    im[mask] = 0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(im, vmin = 0, vmax = 200)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.show()

This is a noisy diffraction pattern of polycrystalline vanadium dioxide. 
Finding the center of such a symmetry pattern can be done with the 
:func:`powder_center` routine::
	
	from skued import powder_center
	ic, jc = powder_center(im, mask = mask)
	
	# Plotting the center as a black disk
	import numpy as np
	ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),
	                     indexing = 'ij')
	rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
	im[rr < 100] = 0

	plt.imshow(im, vmax = 200)
	plt.show()

.. plot::

    from skued import powder_center, diffread
    import numpy as np
    import matplotlib.pyplot as plt
    path = 'Cr_1.tif'
    im = diffread(path)
    mask = np.zeros_like(im, dtype = np.bool)
    mask[0:1250, 950:1250] = True
    ic, jc = powder_center(im, mask = mask)
    ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),indexing = 'ij')
    rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
    im[rr < 100] = 1e6

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(im, vmin = 0, vmax = 200)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    plt.show()

Angular average
---------------

First, we create a test image::

	import numpy as np
	import matplotlib.pyplot as plt
	from skued import gaussian

	image = np.zeros( (256, 256) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center

	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	rr = np.sqrt((xx - xc)**2 + (yy-yc)**2)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
	image[np.logical_and(rr < 40, rr > 38)] = 1
	image[np.logical_and(rr < 100, rr > 98)] = 0.5
	image /= image.max()	# Normalize max to 1
	image += np.random.random(size = image.shape)

	plt.imshow(image)
	plt.show()

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from skued import gaussian
    image = np.zeros( (256, 256) )
    xc, yc = image.shape[0]/2, image.shape[1]/2	# center
    extent = np.arange(0, image.shape[0])
    xx, yy = np.meshgrid(extent, extent)
    rr = np.sqrt((xx - xc)**2 + (yy-yc)**2)
    image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
    image[np.logical_and(rr < 40, rr > 38)] = 1
    image[np.logical_and(rr < 100, rr > 98)] = 0.5
    image /= image.max()	# Normalize max to 1
    image += np.random.random(size = image.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(image)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.show()


... and we can easily compute an angular average::
	
	from skued import azimuthal_average

	radius, intensity = azimuthal_average(image, (xc, yc))

	plt.plot(radius, intensity)

.. plot::
	
	from skued import azimuthal_average, gaussian
	import numpy as np
	import matplotlib.pyplot as plt
	image = np.zeros( (256, 256) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	rr = np.sqrt((xx - xc)**2 + (yy-yc)**2)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
	image[np.logical_and(rr < 40, rr > 38)] = 1
	image[np.logical_and(rr < 100, rr > 98)] = 0.5
	image /= image.max()	# Normalize max to 1
	image += np.random.random(size = image.shape)
	radius, intensity = azimuthal_average(image, (xc, yc))
	plt.plot(radius, intensity)
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
	ax.imshow(im, vmin = 0, vmax = 2e3)
	ax.xaxis.set_visible(False)
	ax.yaxis.set_visible(False)
	plt.show()

We can consider the image *without hotspots* as the baseline of the image *with hotspots* ::

	from skued import diffread, baseline_dwt

	im = diffread('hotspots.tif')
	denoised = baseline_dwt(im, max_iter = 250, level = 1, wavelet = 'sym2', axis = (0, 1))

The result is plotted below:

.. plot::

	import matplotlib.pyplot as plt
	from skued import diffread, baseline_dwt

	im = diffread('hotspots.tif')
	denoised = baseline_dwt(im, max_iter = 250, level = 1, wavelet = 'sym2', axis = (0, 1))

	fig, (ax1, ax2) = plt.subplots(1, 2)
	ax1.imshow(im, vmin = 0, vmax = 2e3)
	ax2.imshow(denoised, vmin = 0, vmax = 2e3)

	for ax in (ax1, ax2):
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
	plt.show()

Try different combinations of wavelets, levels, and number of iterations (``max_iter``).

Notice that the baseline-removal function used in the case of an image is :func:`baseline_dwt`, which works on 2D arrays.
The same is not possible with :func:`baseline_dt`, which only works on 1D arrays at this time.

:ref:`Return to Top <image_analysis_tutorial>`