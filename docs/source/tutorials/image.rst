.. include:: ../references.txt

.. _image_analysis_tutorial:

***********************
Image Analysis Tutorial
***********************

Due to the high electron cross-section, background signals (or baseline) are
much more of a problem for electron diffraction than equivalent X-ray experiments.

Contents
========

* :ref:`streaming`
* :ref:`powder`

.. _streaming:

Streaming Image Processing
==========================

Diffraction datasets can be large, much larger than a machine can handle
in-memory. In this case, it makes sense to assemble processing pipelines instead
of working on the data all at once.

Consider the following snippet to combine 50 images 
from an iterable :code:`source`::

	import numpy as np

	images = np.empty( shape = (2048, 2048, 50) )
	from index, im in enumerate(source):
	    images[:,:,index] = im
	
	avg = np.average(images, axis = 2)

If the :code:`source` iterable provided 1000 images, the above routine would
not work on most machines. Moreover, what if we want to transform the images 
one by one before averaging them? What about looking at the average while it 
is being computed?

Scikit-ued provides some functions that can make streaming processing possible. These
function will have an 'i' prefix (for :code:`iterator`). Let's look at an example::

	import numpy as np
	from skued.image import ialign, iaverage
	from skimage.io import imread

	stream = map(imread, list_of_filenames)
	aligned = ialign(stream)
	averaged = iaverage(aligned)

At this point, the generators :code:`map`, :code:`ialign`, and :code:`iaverage` are 'wired'
but will not compute anything until it is requested. We can use the function
:code:`last` to get at the final average, but we could also look at the average
step-by-step by calling :code:`next`::

	avg = next(averaged) # only one images is loaded, aligned and added to the average
	total = last(averaged) # average of the entire stream

An important advantage of processing images in a streaming fashion is the lower
memory usage; this allows the use of multiple processes in parallel::

	from skued import pmap

	def align_and_avg(filenames):
	    stream = map(imread, filenames)
	    aligned = ialign(stream)
	    return last(iaverage(aligned))
	
	batches = [list_of_filenames1, list_of_filenames2, list_of_filenames3]
	
	for avg in pmap(align_and_avg, batches, processes = 3):
	    # write to disk of display
	    pass

Example: averaging with error
------------------------------

It is possible to combine :code:`iaverage` and :code:`isem` into a single stream using :code:`itertools.tee`. 
Here is a recipe for it::

	def iaverage_with_error(images, weights):    
	    """ 
	    Combined streaming average and standard error of diffraction images. 
		
	    Parameters
	    ----------
	    images : iterable of ndarrays
	    	Images to be averaged. This iterable can also a generator.
	    weights : iterable of ndarray, iterable of floats, or None, optional
	    	Array of weights. See `numpy.average` for further information. If None (default), 
	    	total picture intensity of valid pixels is used to weight each picture.
	    
	    Yields
	    ------
	    avg: `~numpy.ndarray`
	    	Weighted average. 
	    sem: `~numpy.ndarray`
	    	Weighted average. 
	    """
	    stream1, stream2 = itertools.tee(images, n = 2)
	    averages = iaverage(stream1, weights = weights)
	    errors = isem(stream2)
	    yield from zip(averages, errors)

.. _powder:

Image analysis on polycrystalline diffraction patterns
======================================================

Center-finding
--------------
Polycrystalline diffraction patterns display concentric rings, and finding
the center of those concentric rings is important.

Let's load a test image::
	
	from skimage import img_as_uint
	from skimage.io import imread
	import matplotlib.pyplot as plt

	path = '\\data\\vo2.tif'
	im = img_as_uint(imread(path, plugin = 'tifffile'))

	mask = np.zeros_like(im, dtype = np.bool)
	mask[0:1250, 700:1100] = True
	im[mask] = 0

	plt.imshow(im, vmin = 1000, vmax = 1200)
	plt.show()

.. plot::

	from skimage import img_as_uint
	from skimage.io import imread
	import matplotlib.pyplot as plt
	path = 'data\\vo2.tif'
	im = img_as_uint(imread(path, plugin = 'tifffile'))
	mask = np.zeros_like(im, dtype = np.bool)
	mask[0:1250, 700:1100] = True
	im[mask] = 0
	plt.imshow(im, vmin = 1000, vmax = 1200)
	plt.show()

This is a noisy diffraction pattern of polycrystalline vanadium dioxide. 
Finding the center of such a symmetry pattern can be done with the 
:code:`powder_center` routine::
	
	from skued.image import powder_center
	ic, jc = powder_center(im, mask = mask)
	
	# Plotting the center as a black disk
	import numpy as np
	ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),
	                     indexing = 'ij')
	rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
	im[rr < 100] = 0

	plt.imshow(im, vmax = 1200)
	plt.show()

.. plot::

	from skimage import img_as_uint
	from skimage.io import imread
	import numpy as np
	import matplotlib.pyplot as plt
	path = 'data\\vo2.tif'
	im = img_as_uint(imread(path, plugin = 'tifffile'))
	from skued.image import powder_center
	mask = np.zeros_like(im, dtype = np.bool)
	mask[0:1250, 700:1100] = True
	ic, jc = powder_center(im, mask = mask)
	ii, jj = np.meshgrid(np.arange(im.shape[0]), np.arange(im.shape[1]),indexing = 'ij')
	rr = np.sqrt((ii - ic)**2 + (jj - jc)**2)
	im[rr < 100] = 0
	plt.imshow(im, vmin = 1000, vmax = 1200)
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
	plt.imshow(image)
	plt.show()


... and we can easily compute an angular average::
	
	from skued.image import angular_average

	radius, intensity = angular_average(image, (xc, yc))

	plt.plot(radius, intensity)

.. plot::
	
	from skued.image import angular_average
	from skued import gaussian
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
	radius, intensity = angular_average(image, (xc, yc))
	plt.plot(radius, intensity)
	plt.show()

:ref:`Return to Top <image_analysis_tutorial>`