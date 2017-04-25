.. include:: ../references.txt

.. _image_analysis_tutorial:

***********************
Image Analysis Tutorial
***********************

Due to the high electron cross-section, background signals (or baseline) are
much more of a problem for electron diffraction than equivalent X-ray experiments.

Contents
========

* :ref:`symmetry`

.. _symmetry:

Angular average of polycrystalline diffraction patterns
=======================================================

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