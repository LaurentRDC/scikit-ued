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

==============
The Atom class
==============

.. autoclass:: Atom

=================
The Lattice class
=================

.. autoclass:: Lattice

=========
Utilities
=========

To help with fleshing out unit cell atoms from symmetry operators:

.. autofunction:: symmetry_expansion