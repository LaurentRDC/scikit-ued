.. include:: ../references.txt

.. _structure_tutorial:

**************************
Modeling atomic structures
**************************

Contents
========

* :ref:`Crystal`
* :ref:`Atom`


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
	
	hemoglobin = Crystal.from_pdb('1gzx')

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
3. (optional) symmetry operators that generate the full unit cell from the asymmetric unit cell.

As an example, let's create the simplest crystal structure known: 
`alpha-Polonium (simple cubic) <https://en.wikipedia.org/wiki/Polonium#Solid_state_form>`_::
	
	from skued import Crystal, Atom
	import numpy as np

	lattice_vectors = 3.35 * np.eye(3)
	unitcell = [Atom('Po', coords = [0,0,0])]

	polonium = Crystal(unitcell, lattice_vectors)

In the case where atoms are given as an asymmetric unit cell and a set of symmetry operators, you can use the
:func:`symmetry_expansion` function to generate a set of *unique* atoms (even if some symmetry operators might be redundant).

Crystal attributes
------------------
The :class:`Crystal` object provides some interfaces for easy structure manipulation. First, a :class:`Crystal` is an iterable::

	from skued.structure import graphite

	for atm in graphite:	#Loops over atoms in the unit cell
	    print(atm.element, atm.coords)
    
The :func:`len` of a :class:`Crystal` is the unit cell size (in number of atoms)::

    from skued import Crystal

    c = Crystal.from_pdb('1gzx') # hemoglobin
    len(c) # 17536

:class:`Crystal` instances can be equated to each other::

    gold = Crystal.from_database('Au')
    silver = Crystal.from_database('Ag')

    assert gold == silver # false

If a :class:`Crystal` was generated from a file, the path to its file can be retrieved
from the :attr:`source` attribute::

    c = Crystal.from_pdb('1gzx')
    print(c.source)

:class:`Crystal` instances can be converted to NumPy arrays as well::

    import numpy

    arr = numpy.array(Crystal.from_database('Si'))

:data:`arr` will contain one row per unit cell atom:

.. table::
    :widths: grid

    +--------------+-----------------------+
    |Atomic Number |Fractional coordinates |
    +==============+=======+========+======+
    |      Z1      |  x1   |   y1   |  z1  |
    +--------------+-------+--------+------+
    |      Z2      |  x2   |   y2   |  z2  |
    +--------------+-------+--------+------+
    |      Z3      |  x3   |   y3   |  z3  |
    +--------------+-------+--------+------+
    |            ...                       |
    +--------------------------------------+

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

* ``'international_symbol'``: International Tables of Crystallography space-group symbol (short);

* ``'international_full'``: International Tables of Crystallography space-group full symbol;

* ``'hall_symbol'`` : Hall symbol;

* ``'pointgroup'`` : International Tables of Crystallography point-group;

* ``'international_number'`` : International Tables of Crystallography space-group number (between 1 and 230);

* ``'hall_number'`` : Hall number (between 1 and 531).

Scattering utilities
--------------------
:class:`Lattice` objects have a few methods that make life easier when dealing with scattering data and modeling.

The conversion between Miller indices and scattering vectors is available:: 

	from skued.structure import graphite

    # Behavior inherited from Lattice superclass
	G = graphite.scattering_vector(1,0,0)
	h, k, l = graphite.miller_indices(G) #1, 0, 0

Compatibility with ASE
----------------------
The `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/index.html>`_ is a powerful tool. You can harness its power and convert
between :class:`ase.Atoms` and :class:`skued.Crystal` at will.

To create an :class:`ase.Atoms` object from a :class:`Crystal`, use the :meth:`Crystal.ase_atoms` method::

	from ase.calculators.abinit import Abinit
	from skued import Crystal
	
	gold = Crystal.from_database('Au')
	ase_gold = gold.ase_atoms(calculator = Abinit(...))

All keywords of the :class:`ase.Atoms` constructor are supported. To get back to a :class:`Crystal` instance::

	gold2 = Crystal.from_ase(ase_gold)

.. _atom:

Compatibility with ASE
----------------------
The `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/index.html>`_ is a powerful tool. You can harness its power and convert
between :class:`ase.Atoms` and :class:`skued.Crystal` at will.

To create an :class:`ase.Atoms` object from a :class:`Crystal`, use the :meth:`Crystal.ase_atoms` method::

	from ase.calculators.abinit import Abinit
	from skued import Crystal
	
	gold = Crystal.from_database('Au')
	ase_gold = gold.ase_atoms(calculator = Abinit(...))

All keywords of the :class:`ase.Atoms` constructor are supported. To get back to a :class:`Crystal` instance::

	gold2 = Crystal.from_ase(ase_gold)

.. _atom:

The :class:`Atom` Class
=======================
The basis of structure manipulations is to manipulate atoms. :class:`Atom` objects are in the
category of `Transformable` objects, meaning that their coordinates can be transformed
according to any affine transform.

To create an atom, simply provide its element and coordinates::
	
	from skued import Atom

	copper = Atom(element = 'Cu', coords = [0,0,0])

Optional information can be give, such as magnetic moment and mean-squared displacement. For users of :mod:`ase`, 
another possibility is to instantiate an :class:`Atom` from an :class:`ase.Atom` using the :meth:`Atom.from_ase` 
constructor.

Since we are most concerned with atoms in crystals, the coordinates here are assumed to be fractional.
The real-space position with respect to a :class:`Crystal` or :class:`Lattice` can be accessed using the 
:meth:`xyz` method::

    from skued.structure import graphite
    
    carbon = list(graphite)[-1]
    fractional = carbon.coords
    real = carbon.xyz(lattice = graphite)

:ref:`Return to Top <structure_tutorial>`