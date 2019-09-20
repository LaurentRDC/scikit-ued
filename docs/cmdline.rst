.. include:: references.txt

.. _cmdline:

**********************
Command-line utilities
**********************

.. currentmodule:: skued

Scikit-ued includes command-line utilities for repetitive tasks. To see the available 
commands on your system::

    > skued --help

Interactive image viewer
------------------------

One command-line utility is the interactive image viewer. To read and display
a diffraction pattern/image located at `path/to/file`::

    > skued diffshow path/to/file

Make sure that `PyQtGraph`, an optional dependency, is installed.

Crystal information
-------------------

Another command-line utility is the crystallographic information script, ``crystinfo``. 
To show information about a crystal file (``vo2.cif``, for example)::

    > skued crystinfo vo2.cif
    Crystal object with following unit cell:
        Atom O  @ (0.20, 0.80, 0.50)
        Atom O  @ (0.70, 0.70, 0.00)
        Atom O  @ (0.80, 0.20, 0.50)
        Atom O  @ (0.30, 0.30, 0.00)
        Atom V  @ (0.50, 0.50, 0.50)
        Atom V  @ (0.00, 0.00, 0.00)
    Lattice parameters:
        a=4.517Å, b=4.517Å, c=2.872Å
        α=90.000°, β=90.000°, γ=90.000°
    Chemical composition:
        O: 66.667%
        V: 33.333%
    Source: 
        vo2.cif
    Symmetry information:
        International symbol 
                    (short) ..... P4_2/mnm
                    (full) ...... P 4_2/m 2_1/n 2/m
        International number .... 136
        Hermann-Mauguin symbol .. P42/mnm
        Pointgroup .............. D4h
        Hall Number ............. 419
        Centering ............... CenteringType.primitive
