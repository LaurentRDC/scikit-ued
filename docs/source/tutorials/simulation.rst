.. include:: ../references.txt

.. _simulation_tutorial:

*******************
Simulation Tutorial
*******************

Simulating diffraction patterns

Contents
========

* :ref:`powdersim`
* :ref:`electrostatic`

.. _powdersim:

Polycrystalline Diffraction Simulation
======================================
The simplest diffraction simulation available in scikit-ued is polycrystalline
diffraction simulation.

All you need is a :class:`Crystal` object and a range of scattering length, defined from 
scattering angles `s` as :math:`s = \sin{\theta}/\lambda`::
	
	import matplotlib.pyplot as plt
	import numpy as np
	from skued import powdersim
	from skued import Crystal
	graphite = Crystal.from_database('C')

	s = np.linspace(0.1, 0.8, 1000)
	diff = powdersim(crystal = graphite, scattering_length = s)

	plt.figure()
	plt.plot(s, diff/diff.max())

After plot formatting:

.. plot::
	
	import matplotlib.pyplot as plt
	import numpy as np
	from skued import Crystal
	graphite = Crystal.from_database('C')
	from skued import powdersim
	s = np.linspace(0.1, 0.8, 1000)
	diff = powdersim(crystal = graphite, scattering_length = s)
	plt.figure()
	plt.plot(s, diff/diff.max())
	plt.xlim([s.min(), s.max()])
	plt.xlabel('$s = \sin{\\theta}/\lambda$')
	plt.ylabel('Diffracted intensity (A.u.)')
	plt.title('Polycrystalline graphite diffraction')

.. _electrostatic:

Electrostatic Potential Simulation
==================================
The scattering potential of electrons is the crystal electrostatic potential; hence
computing such potential is a useful tool.

To compute the electrostatic potential for an infinite crystal on an arbitrary 3D mesh,
take a look at :func:`electrostatic`::

    from skued import Crystal
    from skued import electrostatic

    graphite = Crystal.from_database('C')

    extent = np.linspace(-10, 10, 256)
    xx, yy, zz = np.meshgrid(extent, extent, extent)
    potential = electrostatic(graphite, xx, yy, zz)

Another possibility is to calculate the electrostatic potential for an infinite slab in the 
x-y plane, but finite in z-direction, using :func:`pelectrostatic` (p for projected)::

    from skued import pelectrostatic

    extent = np.linspace(-5, 5, 256)
    xx, yy= np.meshgrid(extent, extent)
    potential = pelectrostatic(graphite, xx, yy)

After plot formatting:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from skued import Crystal
    from skued import pelectrostatic

    extent = np.linspace(-5, 5, 256)
    xx, yy = np.meshgrid(extent, extent)

    potential = pelectrostatic(Crystal.from_database('C'), xx, yy)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_title('Electrostatic potential of graphite')
    im = ax.imshow(potential)
    cbar = plt.colorbar(im)
    cbar.set_label('Electrostatic potential ($V \cdot \AA$)')

:ref:`Return to Top <simulation_tutorial>`