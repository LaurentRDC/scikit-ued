"""
Polycrystalline diffraction pattern simulation
"""
from ..voigt import pseudo_voigt
import numpy as np

def powdersim(crystal, scattering_length, broadening = True, **kwargs):
	"""
	Simulates polycrystalline diffraction pattern.

	Parameters
	----------
	crystal : `skued.structure.Crystal`
		Crystal from which to diffract.
	scattering_length : `~numpy.ndarray`, shape (N,)
		Range of scattering length over which to compute the diffraction pattern [2pi/Angs].
	broadening: bool, optional
		If True (default), returned pattern will display Voigt-related broadening.
	
	Returns
	-------
	pattern : `~numpy.ndarray`, shape (N,)
		Diffraction pattern
	"""
	Gx, Gy, Gz = crystal.scattering_vector(*crystal.bounded_reflections(4*np.pi*scattering_length.max()))
	scatt_length = np.sqrt(Gx**2 + Gy**2 + Gz**2)/(4*np.pi)
	intensities = np.absolute(crystal.structure_factor((Gx, Gy, Gz)))**2

	pattern = np.zeros_like(scattering_length)
	if broadening:
		for s, I in zip(scatt_length, intensities):
			pattern += I * pseudo_voigt(scattering_length, s, 0.01, 0.02)
	else:
		for s, I in zip(scatt_length, intensities):
			pattern += I * pseudo_voigt(scattering_length, s, 1e-5, 1e-5)

	return pattern