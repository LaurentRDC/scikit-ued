.. include:: ../references.txt

.. _simulation_tutorial:

.. currentmodule:: skued

***********
Simulations
***********

Simulating diffraction patterns

.. _powdersim:

Polycrystalline Diffraction Simulation
======================================
The simplest diffraction simulation available in scikit-ued is polycrystalline
diffraction simulation.

All you need is a :class:`crystals.Crystal` object and a range of scattering length, defined from 
scattering angles `s` as :math:`s = \sin{\theta}/\lambda`:
	
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from skued import powdersim
    >>> from crystals import Crystal
    >>> 
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> q = np.linspace(1, 10, 1024)
    >>> diff = powdersim(graphite, q)
    >>> 
    >>> fig, ax = plt.subplots(1,1)
    >>> ax.plot(q, diff/diff.max()) # doctest: +SKIP

After plot formatting:

.. plot::
	
	import matplotlib.pyplot as plt
	import numpy as np
	from crystals import Crystal
	graphite = Crystal.from_database('C')
	from skued import powdersim
	q = np.linspace(1, 10, 1024)
	diff = powdersim(graphite, q)
	plt.figure()
	plt.plot(q, diff/diff.max())
	plt.xlim([q.min(), q.max()])
	plt.xlabel('$q (1/\AA)$')
	plt.ylabel('Diffracted intensity (A.u.)')
	plt.title('Polycrystalline graphite diffraction')


.. _kinematicsim:

Kinematical single-crystal simulation
=====================================
Single-crystal kinematical diffraction simulation is available via the :func:`kinematicsim` function:

    >>> from crystals import Crystal
    >>> from skued import kinematicsim
    >>> import numpy as np
    >>>
    >>> extent = np.linspace(-10, 10, num=1024)
    >>> kx, ky = np.meshgrid(extent, extent)
    >>> kk = np.sqrt(kx**2 + ky**2)
    >>> 
    >>> crystal = Crystal.from_database('C') # graphite
    >>> I = kinematicsim(crystal, kx=kx, ky=ky, energy=50)
    
.. plot::

    from crystals import Crystal
    from skued import kinematicsim
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter

    extent = np.linspace(-10, 10, num=1024)
    kx, ky = np.meshgrid(extent, extent)
    kk = np.sqrt(kx**2 + ky**2)

    I = kinematicsim(Crystal.from_database("C"), kx=kx, ky=ky, energy=50)
    I[kk < 1] = 0.0

    # Clean up the image a bit
    I[:] = gaussian_filter(I, sigma=3)

    fig, ax = plt.subplots(1,1, figsize=(4.5, 4.5))
    ax.imshow(I, vmax=0.9*I.max(), cmap='inferno')
    ax.axis('off')

    plt.show()

.. _electrostatic:

Electrostatic Potential Simulation
==================================
The scattering potential of electrons is the crystal electrostatic potential; hence
computing such potential is a useful tool.

To compute the electrostatic potential for an infinite crystal on an arbitrary 3D mesh,
take a look at :func:`electrostatic`:

    >>> import numpy as np
    >>> from crystals import Crystal
    >>> from skued import electrostatic
    >>>
    >>> graphite = Crystal.from_database('C')
    >>> 
    >>> extent = np.linspace(-10, 10, 128)
    >>> xx, yy, zz = np.meshgrid(extent, extent, extent)
    >>> potential = electrostatic(graphite, xx, yy, zz)

In order to look at a slice in 2D, take a look at the function :func:`plane_mesh`:

    >>> import numpy as np
    >>> from crystals import Crystal
    >>> from skued import electrostatic, plane_mesh
    >>>
    >>> vo2 = Crystal.from_database('vo2-m1')
    >>> a, b, _ = vo2.lattice_vectors
    >>> 
    >>> extent = np.linspace(-10, 10, 128)
    >>> xx, yy, zz = plane_mesh(a, b, extent, extent)
    >>> potential = electrostatic(vo2, xx, yy, zz)

Another possibility is to calculate the electrostatic potential for an infinite slab in the 
x-y plane, but finite in z-direction, using :func:`pelectrostatic` (p for projected):

    >>> import numpy as np
    >>> from skued import pelectrostatic
    >>> 
    >>> extent = np.linspace(-5, 5, 256)
    >>> xx, yy= np.meshgrid(extent, extent)
    >>> potential = pelectrostatic(graphite, xx, yy)

After plot formatting:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from crystals import Crystal
    from skued import pelectrostatic

    extent = np.linspace(-5, 5, 256)
    xx, yy = np.meshgrid(extent, extent)

    potential = pelectrostatic(Crystal.from_database('C'), xx, yy)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_title('Electrostatic potential of graphite')
    im = ax.imshow(potential, cmap='inferno')
    cbar = plt.colorbar(im)
    cbar.set_label('Electrostatic potential ($V \cdot \AA$)')

:ref:`Return to Top <simulation_tutorial>`