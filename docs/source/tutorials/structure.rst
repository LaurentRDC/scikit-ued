.. include:: ../references.txt

.. _structure_tutorial:

******************
Structure Tutorial
******************

Modeling atomic structures

Contents
========

* :ref:`Atom`
* :ref:`Lattice`
* :ref:`Crystal`

.. _atom:

The Atom Class
==============
The basis of structure manipulations is to manipulate atoms. `Atom` objects are in the
category of `Transformable` objects, meaning that their coordinates can be transformed
according to any affine transform.

To create an atom, simply provide its element and coordinates::
	
	from skued.structure import Atom

	copper = Atom(element = 'Cu', coords = [0,0,0])

One important feature of the `Atom` class is the possibility to compute the electrostatic
potential across meshes::

	import numpy as np
	import matplotlib.pyplot as plt

	xx, yy = np.meshgrid(np.linspace(-0.3, 0.3, num = 100), 
						 np.linspace(-0.3, 0.3, num = 100))
	dist = np.sqrt(xx**2 + yy**2)	# distance from the atom in Angstroms

	es_potential = copper.potential(dist)
	plt.imshow(es_potential)

After plot formatting:

.. plot::
	
	import numpy as np
	import matplotlib.pyplot as plt
	from skued.structure import Atom
	copper = Atom(element = 'Cu', coords = [0,0,0])
	xx, yy = np.meshgrid(np.linspace(-0.3, 0.3, num = 100), 
						 np.linspace(-0.3, 0.3, num = 100))
	dist = np.sqrt(xx**2 + yy**2)	# distance from the atom in Angstroms
	es_potential = copper.potential(dist)
	plt.title('Atomic potential of Cu (log-scale)')
	plt.imshow(np.log(1 + es_potential), extent = [xx.min(), xx.max(), yy.min(), yy.max()])
	plt.ylabel('x-direction ($\AA$)')
	plt.xlabel('y-direction ($\AA$)')
	plt.show()

.. _lattice:

The Lattice Class
=================


.. _crystal:

The Crystal Class
=================

:ref:`Return to Top <image_analysis_tutorial>`