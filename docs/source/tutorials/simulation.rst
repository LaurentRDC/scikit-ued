.. include:: ../references.txt

.. _simulation_tutorial:

*******************
Simulation Tutorial
*******************

Simulating diffraction patterns

Contents
========

* :ref:`powdersim`

.. _powdersim:

Polycrystalline Diffraction Simulation
======================================
The simplest diffraction simulation available in scikit-ued is polycrystalline
diffraction simulation.

All you need is a :code:`Crystal` object and a range of scattering length, defined from 
scattering angles `s` as :math:`s = \sin{\theta}/\lambda`::
	
	import matplotlib.pyplot as plt
	import numpy as np
	from skued.structure import graphite
	from skued.simulation import powdersim

	s = np.linspace(0.1, 0.8, 1000)
	diff = powdersim(crystal = graphite, scattering_length = s)

	plt.figure()
	plt.plot(s, diff/diff.max())

After plot formatting:

.. plot::
	
	import matplotlib.pyplot as plt
	import numpy as np
	from skued.structure import graphite
	from skued.simulation import powdersim
	s = np.linspace(0.1, 0.8, 1000)
	diff = powdersim(crystal = graphite, scattering_length = s)
	plt.figure()
	plt.plot(s, diff/diff.max())
	plt.xlim([s.min(), s.max()])
	plt.xlabel('$s = \sin{\theta}/\lambda$')
	plt.ylabel('Diffracted intensity (A.u.)')
	plt.title('Polycrystalline graphite diffraction')

:ref:`Return to Top <simulation_tutorial>`