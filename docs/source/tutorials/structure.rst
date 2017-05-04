.. include:: ../references.txt

.. _structure_tutorial:

******************
Structure Tutorial
******************

Modeling atomic structures

Contents
========

* :ref:`Atom`
* :ref:`Crystal`

.. _atom:

The :code:`Atom` Class
======================
The basis of structure manipulations is to manipulate atoms. :code:`Atom` objects are in the
category of `Transformable` objects, meaning that their coordinates can be transformed
according to any affine transform.

To create an atom, simply provide its element and coordinates::
	
	from skued.structure import Atom

	copper = Atom(element = 'Cu', coords = [0,0,0])

One important feature of the :code:`Atom` class is the possibility to compute the electrostatic
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

.. _crystal:

The :code:`Crystal` Class
=========================
Diffraction experiments relying on the redundancy of crystals to get good experimental signals;
hence, handling crystal models is the main feature of the :code:`skued.structure` package.

Constructing a :code:`Crystal` object from a file
-------------------------------------------------
Creating a :code:`Crystal` object can be done most easily from a Crystal Information File (CIF, .cif) or 
a Protein DataBank file (PDB, .pdb). For a CIF file::
	
	from skued.structure import Crystal

	TiSe2 = Crystal.from_cif('tise2.cif')

Protein DataBank files are even easier to handle; simply provide the 4-letter identification code
and the structure file will be downloaded, cached, and parsed::
	
	from skued.structure import Crystal

	bacteriorhodopsin = Crystal.from_pdb('1fbb')

Constructing a :code:`Crystal` object by hand
---------------------------------------------
If you don't have a file on hand, or want to create an idealized crystal, consider building a :code:`Crystal`
object by hand.

To do this, you need:
1. iterable of :code:`Atom` objects, with coordinates. These atoms can either be the full unit cell
or the asymmetric unit cell;
2. three lattice vectors;
3. Symmetry operators (optional). These symmetry operators will be applied to the atoms to generate
the full unit cell. Hence, if your iterable of atoms contains the entire unit cell, symmetry operators do
not need to be provided.

As an example, let's create the simplest crystal structure known: 
`alpha-Polonium (simple cubic)<https://en.wikipedia.org/wiki/Polonium#Solid_state_form>`::
	
	from skued.structure import Crystal, Atom
	import numpy as np

	lattice_vectors = 3.35 * np.eye(3)
	atoms = [Atom('Po', coords = [0,0,0])]

	polonium = Crystal(atoms = atoms, lattice_vectors = lattice_vectors)

That's it!

Crystal attributes
------------------
The :code:`Crystal` object provides some interfaces for easy structure manipulation. First, a :code:`Crystal` is an iterable::

	from skued.structure import graphite

	for atm in graphite:	#Loops over atoms in the unit cell
	    print(atm)

	for atm in graphite.unitcell:	#equivalent
	    print(atm)

Note that :code:`iter(graphite)` is a generator, whereas :code:`graphite.unitcell` is a list; this
distinction is important when handling large crystals. 

To access atoms only present in the (user-defined) asymmetric cell::

	for atm in graphite.atoms:
	    print(atm)

:code:`Crystal` objects also provide interoperability with :code:`spglib`::

	import spglib

	spglib.get_symmetry_dataset(cell = graphite.spglib_cell, symprec = 1e-4)

Lattice vectors and reciprocal space
-------------------------------------
Once a :code:`Crystal` object is ready, you can manipulate the lattice parameters via the underlying :code:`Lattice`
super-class. Let's use the built-in example of graphite::

	from skued.structure import graphite
	
	a1, a2, a3 = graphite.lattice_vectors
	b1, b2, b3 = graphite.reciprocal_vectors

The standard `three lengths and angles` description of a lattice is also accessible::

	a, b, c, alpha, beta, gamma = graphite.parameters

Atomic potential
----------------
The crystal electrostatic potential (the scattering potential leading to electron diffraction) can be
computed from a :code:`Crystal`::

	import numpy as np
	import matplotlib.pyplot as plt
	from skued.structure import graphite

	xx, yy = np.meshgrid(np.linspace(-1, 1, num = 100), 
	                     np.linspace(-1, 1, num = 100))
	zz = np.zeros_like(xx)

	potential = graphite.potential(xx, yy, zz)
	plt.imshow(es_potential)

After plot formatting:

.. plot::
	
	import numpy as np
	import matplotlib.pyplot as plt
	from skued.structure import graphite
	xx, yy = np.meshgrid(np.linspace(-5, 5, num = 512), 
						 np.linspace(-5, 5, num = 512))
	zz = np.zeros_like(xx)
	potential = graphite.potential(xx, yy, zz)
	plt.title('Electrostatic potential of graphite (log-scale)')
	plt.imshow(np.log(1 + 0.01*potential), extent = [xx.min(), xx.max(), yy.min(), yy.max()])
	plt.ylabel('x-direction ($\AA$)')
	plt.xlabel('y-direction ($\AA$)')
	plt.show()

Note that while the :code:`graphite` crystal only has four atoms in its unitcell, the :code:`graphite.potential` method
will process the meshes :code:`xx`, :code:`yy`, and :code:`zz` through the :code:`skued.minimum_image_distance` function
to implement the **minimum-image distance convention** for periodic boundary conditions All this to say that the potential 
is computed for a seemingly-infinite crystal.

:ref:`Return to Top <structure_tutorial>`