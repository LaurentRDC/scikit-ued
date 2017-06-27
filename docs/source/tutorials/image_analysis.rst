.. include:: ../references.txt

.. _image_analysis_tutorial:

***********************
Image Analysis Tutorial
***********************

Due to the high electron cross-section, background signals (or baseline) are
much more of a problem for electron diffraction than equivalent X-ray experiments.

Contents
========

* :ref:`register_translation`
* :ref:`powder`

.. _register_translation:

Image translation registration
==============================



Translation registration is an important problem in image-processing, and solutions
to this are commonplace. What distinguishes diffraction experiments from the typical
translation registration problem is the fact that diffraction images often have invalid
zones.

Let's load a test image::
	
	from skimage import img_as_uint
	from skimage.io import imread
	import matplotlib.pyplot as plt

	path = '\\data\\vo2.tif'
	im = img_as_uint(imread(path, plugin = 'tifffile'))
	plt.imshow(im, vmin = im.min(), vmax = 1200)
	plt.show()

.. plot::

	from skimage import img_as_uint
	from skimage.io import imread
	import matplotlib.pyplot as plt
	path = 'data\\vo2.tif'
	im = img_as_uint(imread(path, plugin = 'tifffile'))
	plt.imshow(im, vmin = im.min(), vmax = 1200)
	plt.show()

During an experiment, the diffraction pattern will shift around; however, the beam-block
and the center beam (of which we see the edges) will not shift with the pattern. Moreover,
detectors can have dead pixels or hot pixels, which will throw off translation registration. 
Therefore, in comparing two diffraction patterns for alignment, we need to declare zones of 
the diffraction pattern as invalid.


.. _powder:

Image analysis on polycrystalline diffraction patterns
======================================================

Center-finding
--------------
Polycrystalline diffraction patterns display concentric rings, and finding
the center of those concentric rings is important.

Consider the diffraction pattern presented in :ref:`register_translation`.
Finding the center of such a symmetric pattern can be done with the 
:code:`powder_center` routine::
	
	from skued.image_analysis import powder_center
	ic, jc = powder_center(im)
	
	# Plotting the center as a black disk
	import numpy as np
	ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),
	                     indexing = 'ij')
	rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
	im[rr < 10] = 30000

	plt.imshow(im, vmin = im.min(), vmax = 1200)
	plt.show()

.. plot::

	from skimage import img_as_uint
	from skimage.io import imread
	import numpy as np
	import matplotlib.pyplot as plt
	path = 'data\\vo2.tif'
	im = img_as_uint(imread(path, plugin = 'tifffile'))
	from skued.image_analysis import powder_center
	ic, jc = powder_center(im)
	ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),indexing = 'ij')
	rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
	im[rr < 10] = 30000
	plt.imshow(im, vmin = im.min(), vmax = 1200)
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
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
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
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 100)
	image /= image.max()
	image += np.random.random(size = image.shape)
	plt.imshow(image)
	plt.show()


... and we can easily compute an angular average::
	
	from skued.image_analysis import angular_average

	radius, intensity = angular_average(image, (xc, yc))

	plt.plot(radius, intensity)

.. plot::
	
	from skued.image_analysis import angular_average
	from skued import gaussian
	import numpy as np
	import matplotlib.pyplot as plt
	from skued import gaussian
	image = np.zeros( (256, 256) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
	image /= image.max()
	image += np.random.random(size = image.shape)
	radius, intensity = angular_average(image, (xc, yc))
	plt.plot(radius, intensity)
	plt.show()

We can also get the standard error in the mean of this average::

	extras = dict()
	radius, intensity = angular_average(image, (xc, yc), extras = extras)

	plt.errorbar(radius, intensity, yerr - extras['error'])

.. plot::
	
	from skued.image_analysis import angular_average
	from skued import gaussian
	import numpy as np
	import matplotlib.pyplot as plt
	from skued import gaussian
	image = np.zeros( (256, 256) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
	image /= image.max()
	image += np.random.random(size = image.shape)
	extras = dict()
	radius, intensity = angular_average(image, (xc, yc), extras = extras)
	plt.errorbar(radius, intensity, yerr = extras['error'])
	plt.show()

As expected, the points close to the center and very far from it display as higher
amount of noise due to the number of points used in the averaging.

:ref:`Return to Top <image_analysis_tutorial>`