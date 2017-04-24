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

	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center

	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
	image[np.logical_and(128 < rr,rr < 130)] = 1 # ring

	plt.imshow(image)
	plt.show()

.. plot::

	import numpy as np
	import matplotlib.pyplot as plt
	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
	image[np.logical_and(128 < rr,rr < 130)] = 1 # ring
	plt.imshow(image)
	plt.show()


... and we can easily compute an angular average::
	
	from skued.image_analysis import angular_average

	intensity = angular_average(image, (xc, yc))

	plt.plot(intensity)

.. plot::
	
	from skued.image_analysis import angular_average
	import numpy as np
	import matplotlib.pyplot as plt
	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
	image[np.logical_and(120 < rr,rr < 130)] = 1 # ring
	intensity = angular_average(image, (xc, yc))
	plt.plot(intensity)
	plt.show()

We can extract more information from this averaging by passing an extra
dictionary ::
	
	extras = dict()
	intensity = angular_average(image, (xc, yc), extras = extras)
	radius, error = extras['radius'], extras['error']

	plt.plot(radius, intensity)

.. plot ::
	
	from skued.image_analysis import angular_average
	import numpy as np
	import matplotlib.pyplot as plt
	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
	image[np.logical_and(120 < rr,rr < 130)] = 1 # ring
	extras = dict()
	intensity = angular_average(image, (xc, yc), extras = extras)
	radius, error = extras['radius'], extras['error']
	plt.plot(radius, intensity)
	plt.show()

:ref:`Return to Top <image_analysis_tutorial>`