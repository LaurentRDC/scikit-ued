.. include:: references.txt

.. _manipulating_structures:

***********************
Manipulating Structures
***********************

.. currentmodule:: skued

Handling crystal structure information is crucial for many data analysis and modelling tasks.
See the :ref:`Structure tutorial <structure_tutorial>` for some examples on how to use the following
classes.

=================
The Crystal class
=================

.. autoclass:: Crystal
    :show-inheritance:
    :members:
    :inherited-members:

==============
The Atom class
==============

.. autoclass:: Atom

============
Base classes
============

The :class:`Lattice` class allows for manipulating lattice information separately from
atomic information.

.. autoclass:: Lattice

The possible lattice systems are stored in the :class:`LatticeSystem` enumeration:

.. autoclass:: LatticeSystem
    :show-inheritance:
    :members:
    :undoc-members:

In order to abstract away the possibility of crystal substructures in the future,
the :class:`Crystal` class derives from the following base class:

.. autoclass:: AtomicStructure
    :members:
    :special-members:

=========
Utilities
=========

To help with fleshing out unit cell atoms from symmetry operators:

.. autofunction:: symmetry_expansion
.. autofunction:: lattice_system