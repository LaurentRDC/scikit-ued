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

The :class:`Atom` Class
=======================
The basis of structure manipulations is to manipulate atoms. :class:`Atom` objects are in the
category of `Transformable` objects, meaning that their coordinates can be transformed
according to any affine transform.

To create an atom, simply provide its element and coordinates::
	
	from skued import Atom

	copper = Atom(element = 'Cu', coords = [0,0,0])

Since we are most concerned with atoms in crystals, the coordinates here are assumed to be fractional.
The real-space position with respect to a :class:`Crystal` or :class:`Lattice` can be accessed using the 
:meth:`xyz` method::

    from skued.structure import graphite
    
    carbon = list(graphite)[-1]
    fractional = carbon.coords
    real = carbon.xyz(lattice = graphite)

One important feature of the :class:`Atom` class is the possibility to compute the electrostatic
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

The :class:`Crystal` Class
==========================
Diffraction experiments rely on the redundancy of crystals to get good experimental signals;
hence, handling crystal models is the main feature of the :mod:`skued.structure` subpackage.

Constructing a :class:`Crystal` object
--------------------------------------
Creating a :class:`Crystal` object can be done most easily from a Crystal Information File (CIF, .cif)::
	
	from skued import Crystal

	TiSe2 = Crystal.from_cif('tise2.cif')

Scikit-ued also has an internal database of CIF files. Valid names are stored in :attr:`Crystal.builtins` and can be
constructed like so::

	assert 'Au' in Crystal.builtins
	gold = Crystal.from_database('Au')

Protein DataBank files are even easier to handle; simply provide the 4-letter identification code
and the structure file will be taken care of by scikit-ued::
	
	hemoglobin = Crystal.from_pdb('1gxz')

Another convenient way to construct a :class:`Crystal` is through the `Crystallography Open Database <http://www.crystallography.net/cod/>`_::

	# Default is the latest revision
	vo2 = Crystal.from_cod(1521124)

	# Revisions are accessible as well
	old_vo2 = Crystal.from_cod(1521124, revision = 140771)

Constructing a :class:`Crystal` object by hand
----------------------------------------------
If you don't have a file on hand, or want to create an idealized crystal, consider building a :class:`Crystal`
object by hand.

To do this, you need:
1. iterable of :class:`Atom` objects, with coordinates. These atoms can either be the full unit cell
or the asymmetric unit cell;
2. three lattice vectors;
3. Symmetry operators (optional). These symmetry operators will be applied to the atoms to generate
the full unit cell. Hence, if your iterable of atoms contains the entire unit cell, symmetry operators do
not need to be provided. The symmetry operators must be expressed in the reduced (or fractional) basis.

As an example, let's create the simplest crystal structure known: 
`alpha-Polonium (simple cubic) <https://en.wikipedia.org/wiki/Polonium#Solid_state_form>`_::
	
	from skued import Crystal, Atom
	import numpy as np

	lattice_vectors = 3.35 * np.eye(3)
	atoms = [Atom('Po', coords = [0,0,0])]

	polonium = Crystal(atoms = atoms, lattice_vectors = lattice_vectors)

That's it!

Crystal attributes
------------------
The :class:`Crystal` object provides some interfaces for easy structure manipulation. First, a :class:`Crystal` is an iterable::

	from skued.structure import graphite

	for atm in graphite:	#Loops over atoms in the unit cell
	    print(atm.element, atm.coords)

Note that iterating over the :attr:`crystal.atoms` attribute may or may not be equivalent to 
:data:`iter(crystal)`, due to the way symmetry operators are defined.

Lattice vectors and reciprocal space
-------------------------------------
Once a :class:`Crystal` object is ready, you can manipulate the lattice parameters via the underlying :class:`Lattice`
super-class. Let's use the built-in example of graphite::

	from skued.structure import graphite
	
	a1, a2, a3 = graphite.lattice_vectors
	b1, b2, b3 = graphite.reciprocal_vectors

The standard `three lengths and angles` description of a lattice is also accessible::

	a, b, c, alpha, beta, gamma = graphite.lattice_parameters

The unit cell volume (and by extensions, density) is also accessible:

	vol = graphite.volume
	density = vol/len(graphite)

Space-group Information
-----------------------
Thanks to `spglib <http://atztogo.github.io/spglib/>`_, we can get symmetry and space-group information 
from a :class:`Crystal` instance::

	from skued import Crystal
	
	gold = Crystal.from_database('Au')
	spg_info = gold.spacegroup_info()

In the above example, :data:`spg_info` is a dictionary with the following four keys:

* ``'international_symbol'``: International Tables of Crystallography space-group symbol (short)

* ``'hall_symbol'`` : Hall symbol

* ``'international_number'`` : International Tables of Crystallography space-group number (between 1 and 230)

* ``'hall_number'`` : Hall number (between 1 and 531)

You can get even more information by using :mod:`spglib` functions directly::

	from spglib import get_symmetry_dataset

	all_the_info = get_symmetry_dataset(gold.spglib_cell)

The content of the :data:`all_the_info` dictionary is documented `here <http://atztogo.github.io/spglib/python-spglib.html#get-symmetry-dataset>`_.
Many of :mod:`spglib`'s routines can be used with :attr:`Crystal.spglib_cell`.

Scattering utilities
--------------------
:class:`Crystal` objects have a few methods that make life easier when dealing with scattering data and modeling.

The conversion between Miller indices and scattering vectors is available:: 

	from skued.structure import graphite

	G = graphite.scattering_vector(1,0,0)
	h, k, l = graphite.miller_indices(G) #1, 0, 0

Arrays of Miller indices can be generated for all Miller indices that fall below a bound::

	h, k, l = graphite.bounded_reflections(12) 	# All reflections below 12 Angs^-1

In this example, :data:`h`, :data:`k`, and :data:`l` are arrays of integers; each combined row is a reflection.

Static structure factor calculation is also possible, both for a single reflection and arrays of reflections::

	import numpy as np

	# For a single reflection
	SF = graphite.structure_factor_miller(1, 0, 0)

	# For an array of reflections: vectorized calculation
	h, k, l = graphite.bounded_reflections(12)
	SF = graphite.structure_factor_miller(h, k, l)
	SF.shape == h.shape 	# True

:ref:`Return to Top <structure_tutorial>`