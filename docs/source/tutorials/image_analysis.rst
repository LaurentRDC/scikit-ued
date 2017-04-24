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

	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center

	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)

	plt.imshow(image)
	plt.show()

.. plot::

	import numpy as np
	import matplotlib.pyplot as plt
	from skued import gaussian
	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
	plt.imshow(image)
	plt.show()


... and we can easily compute an angular average::
	
	from skued.image_analysis import angular_average

	intensity = angular_average(image, (xc, yc))

	plt.plot(intensity)

.. plot::
	
	from skued.image_analysis import angular_average
	from skued import gaussian
	import numpy as np
	import matplotlib.pyplot as plt
	from skued import gaussian
	image = np.zeros( (512, 512) )
	xc, yc = image.shape[0]/2, image.shape[1]/2	# center
	extent = np.arange(0, image.shape[0])
	xx, yy = np.meshgrid(extent, extent)
	image += gaussian([xx, yy], center = [xc, yc], fwhm = 200)
	intensity = angular_average(image, (xc, yc))
	plt.plot(intensity)
	plt.show()

:ref:`Return to Top <image_analysis_tutorial>`