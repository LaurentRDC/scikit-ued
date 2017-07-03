"""
Image correlation and related functions
=======================================
"""
import numpy as np
from scipy.fftpack import fft2, ifft2, ifftshift

# TODO: axes parameter
def masked_xcorr(arr1, arr2, m1, m2 = None):
	"""
	Normalized cross-correlation between two images with invalid pixels. 

	Parameters
	----------
	arr1, arr2 : `~numpy.ndarray`, shape (N,M)
		Images to be cross-correlated
	m1 : `~numpy.ndarray`, shape (N,M)
		Mask of `arr1`. The mask should evaluate to `True`
		(or 1) on invalid pixels.
	m2 : `~numpy.ndarray`, shape (N,M) or None, optional
		Mask of `arr2`. The mask should evaluate to `True`
		(or 1) on invalid pixels. If None (default), `m2` is 
		taken to be the same as `m1`.
		
	Returns
	-------
	out : `~numpy.ndarray`, dtype complex
		Masked, normalized cross-correlation.
		
	References
	----------
	.. [#] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
		IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
	"""
	if m2 is None:
		m2 = m1
		
	m1, m2 = np.asfarray(~m1), np.asfarray(~m2)
	arr1, arr2 = np.asfarray(arr1), np.asfarray(arr2)

	# Contrary to the reference above, we do not compute the fft up to the 'full' size
	# In my experience this has not been a problem.	
	F1 = fft2(arr1 * m1)
	F2s = fft2(arr2 * m2).conj()

	M1 = fft2(m1)
	M2s = fft2(m2).conj()

	# I have noticed no clear performance boost by storing
	# repeated calculation (e.g. ifft2(M1 * M2s)); however, the following
	# is already hard enough to read...
	numerator = ifft2(F1 * F2s)
	numerator -= ifft2(F1 * M2s) * ifft2(M1 * F2s) / ifft2(M1 * M2s)

	denominator = ifft2(fft2(arr1*arr1) * M2s) - (ifft2(F1 * M2s))**2/ifft2(M1 * M2s)
	denominator *= ifft2(M1*fft2(arr2*arr2).conjugate()) - ifft2(M1 * F2s)**2/(ifft2(M1 * M2s))
	np.sqrt(denominator, out = denominator)

	return ifftshift(numerator / denominator)