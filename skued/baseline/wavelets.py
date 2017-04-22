"""
Dual-tree complex wavelet transform module.
"""
from itertools import cycle
import numpy as np
from pywt import dwt, idwt, dwt2, idwt2, dwt_max_level, wavelist, Wavelet
from functools import lru_cache
from os.path import join, dirname

DEFAULT_MODE = 'constant'
DEFAULT_FIRST_STAGE = 'sym4'
DEFAULT_CMP_WAV = 'qshift4'

DATADIR = join(dirname(__file__), 'data')
ALL_QSHIFT = ('qshift1', 'qshift2', 'qshift3', 'qshift4', 'qshift5', 'qshift6')
ALL_COMPLEX_WAV = ('kingsbury99',) + ALL_QSHIFT
ALL_FIRST_STAGE = ('kingsbury99_fs',) + tuple([wav for wav in wavelist(kind = 'discrete') if wav != 'dmey'])

def dualtree(data, first_stage = DEFAULT_FIRST_STAGE, wavelet = DEFAULT_CMP_WAV, level = 'max', mode = DEFAULT_MODE, axis = -1):
    """
    Dual-tree complex wavelet transform, implemented from [1], in 1D, along an axis. 

    Parameters
    ----------
    data: array_like
        Input data. Can be of any shape, but the transform can only be applied in 1D (i.e. along a single axis).
        The length along the axis must be even.
    first_stage : str, optional
        Wavelet to use for the first stage. See dualtree.ALL_FIRST_STAGE for a list of suitable arguments
    wavelet : str, optional
        Wavelet to use in stages > 1. Must be appropriate for the dual-tree complex wavelet transform.
        See dualtree.ALL_COMPLEX_WAV for possible
    level : int or 'max', optional
        Decomposition level (must be >= 0). If level is 'max' (default) then it
        will be calculated using the ``dt_max_level`` function.
    mode : str, optional
        Signal extension mode, see pywt.Modes.
    axis : int, optional
        Axis over which to compute the transform. Default is -1

    Returns
    -------
    [cA_n, cD_n, cD_n-1, ..., cD2, cD1] : list
        Ordered list of coefficients arrays
        where `n` denotes the level of decomposition. The first element
        (`cA_n`) of the result is approximation coefficients array and the
        following elements (`cD_n` - `cD_1`) are details coefficients arrays.
        Arrays have the same number of dimensions as the input.
    
    Raises
    ------
    ValueError
        Raised if axis argument is invalid (e.g. too large).
    
    Notes
    -----
    The implementation uses two tricks presented in [1]:
        `` Different first-stage wavelet ``
            The first level of the dual-tree complex wavelet transform involves a combo of shifted wavelets.
        
        `` Swapping of filters at each stage ``
            At each level > 1, the filters (separated into real and imaginary trees) are swapped.
    
    References
    ----------
    [1] Selesnick, I. W. et al. 'The Dual-tree Complex Wavelet Transform', IEEE Signal Processing Magazine pp. 123 - 151, November 2005.
    """
    data = np.asarray(data, dtype = np.float)/np.sqrt(2)

    if level == 'max':
        level = dt_max_level(data = data, first_stage = first_stage, wavelet = wavelet, axis = axis)
    elif level == 0:
        return [data]
    
    # Check axis bounds
    if axis > data.ndim - 1:
        raise ValueError('Input array has {} dimensions, but input axis is {}'.format(data.ndim, axis))
        
    real_wavelet, imag_wavelet = dualtree_wavelet(wavelet)
    real_first, imag_first = dualtree_first_stage(first_stage)
    
    real_coeffs = _single_tree_analysis_1d(data = data, first_stage = real_first, wavelet = (real_wavelet, imag_wavelet), level = level, mode = mode, axis = axis)
    imag_coeffs = _single_tree_analysis_1d(data = data, first_stage = imag_first, wavelet = (imag_wavelet, real_wavelet), level = level, mode = mode, axis = axis)

    # Combine coefficients into complex form
    # TODO: This is slow
    return [real + 1j*imag for real, imag in zip(real_coeffs, imag_coeffs)]

def idualtree(coeffs, first_stage = DEFAULT_FIRST_STAGE, wavelet = DEFAULT_CMP_WAV, mode = DEFAULT_MODE, axis = -1):
    """
    Inverse dual-tree complex wavelet transform, implemented from [1], along an axis.

    Parameters
    ----------
    coeffs : array_like
        Coefficients list [cAn, cDn, cDn-1, ..., cD2, cD1]
    first_stage : str, optional
        Wavelet to use for the first stage. See dualtree.ALL_FIRST_STAGE for a list of possible arguments.
    wavelet : str, optional
        Wavelet to use in stages > 1. Must be appropriate for the dual-tree complex wavelet transform.
        See dualtree.ALL_COMPLEX_WAV for possible arguments.
    mode : str, optional
        Signal extension mode, see pywt.Modes.
    axis : int, optional
        Axis over which to compute the inverse transform.
        
    Returns
    -------
    reconstructed : ndarray

    Raises
    ------
    ValueError
        If the input coefficients are too few
    
    Notes
    -----
    The implementation uses two tricks presented in [1]:
        `` Different first-stage wavelet ``
            The first level of the dual-tree complex wavelet transform involves a combo of shifted wavelets.
        
        `` Swapping of filters at each stage ``
            At each level > 1, the filters (separated into real and imaginary trees) are swapped.
    
    References
    ----------
    [1] Selesnick, I. W. et al. 'The Dual-tree Complex Wavelet Transform', IEEE Signal Processing Magazine pp. 123 - 151, November 2005.
    """
    if len(coeffs) < 1:
        raise ValueError("Coefficient list too short with {} elements (minimum 1 array required).".format(len(coeffs)))
    elif len(coeffs) == 1: # level 0 inverse transform
        real, imag = coeffs[0], coeffs[0]
    else:
        real_wavelet, imag_wavelet = dualtree_wavelet(wavelet)
        real_first, imag_first = dualtree_first_stage(first_stage)

        real = _single_tree_synthesis_1d(coeffs = [coeff.real for coeff in coeffs], first_stage = real_first, 
                                      wavelet = (real_wavelet, imag_wavelet), mode = mode, axis = axis)
        imag = _single_tree_synthesis_1d(coeffs = [coeff.imag for coeff in coeffs], first_stage = imag_first, 
                                      wavelet = (imag_wavelet, real_wavelet), mode = mode, axis = axis)
    
    return np.sqrt(2)*(real + imag)/2

def dt_max_level(data, first_stage, wavelet, axis = -1):
    """
    Returns the maximum decomposition level from the dual-tree complex wavelet transform.

    Parameters
    ----------
    data : ndarray
        Input data. Can be of any dimension.
    first_stage : str
        Wavelet used in the first stage of the dual-tree cwt. See pywt.wavelist() for suitable arguments.
    wavelet : str
        Dual-tree complex wavelet to use. Argument must be supported by dualtree_wavelet
    axis : int, optional
        Axis over which to compute the transform. Default is -1
        
    Returns
    -------
    max_level : int
    """
    real_wavelet, imag_wavelet = dualtree_wavelet(wavelet)
    return dwt_max_level(data_len = data.shape[axis], filter_len = max([real_wavelet.dec_len, imag_wavelet.dec_len]))

def _normalize_size_axis(approx, detail, axis):
    """ 
    Adjust the approximate coefficients' array size to that of the detail coefficients' array.
    
    Parameters
    ---------- 
    approx : ndarray
    detail: ndarray
    axis : int

    Returns
    -------
    ndarray
        Same shape as detail input.
    """
    if approx.shape[axis] == detail.shape[axis]:
        return approx
    # Swap axes to bring the specific axis to front, truncate array, and re-swap.
    # This is an extension of the 1D case:
    # >>> approx = approx[:-1] 
    return np.swapaxes( np.swapaxes(approx, axis, 0)[:-1] ,0, axis)

def _single_tree_analysis_1d(data, first_stage, wavelet, level, mode, axis):
    """
    Single tree of the forward dual-tree complex wavelet transform.
    
    Parameters
    ----------
    data : ndarray, ndim 1
    first_stage : Wavelet object
    wavelet : 2-tuple of Wavelet object
    level : int
    mode : str
    axis : int

    Returns
    -------
    [cA_n, cD_n, cD_n-1, ..., cD2, cD1] : list
        Ordered list of coefficients arrays
        where `n` denotes the level of decomposition. The first element
        (`cA_n`) of the result is approximation coefficients array and the
        following elements (`cD_n` - `cD_1`) are details coefficients arrays.
    """
    approx, first_detail = dwt(data = data, wavelet = first_stage, mode = mode, axis = axis)    
    coeffs_list = [first_detail]
    for i, wav in zip(range(level - 1), cycle(wavelet)):
        approx, detail = dwt(data = approx, wavelet = wav, mode = mode, axis = axis)
        coeffs_list.append(detail)
    
    # Format list ot be compatible to PyWavelet's format. See pywt.wavedec documentation.
    coeffs_list.append(approx)
    return list(reversed(coeffs_list))

def _single_tree_synthesis_1d(coeffs, first_stage, wavelet, mode, axis):
    """
    Single tree of the inverse dual-tree complex wavelet transform.

    Parameters
    ----------
    coeffs : list
    first_stage : Wavelet object
    wavelet : 2-tuple of Wavelet objects
    mode : str
    axis : int

    Returns
    -------
    reconstructed : ndarray, ndim 1
    """
    # Determine the level except first stage:
    # The order of wavelets depends on whether
    # the level is even or odd.
    level = len(coeffs) - 1
    approx, detail_coeffs, first_stage_detail = coeffs[0], coeffs[1:-1], coeffs[-1]

    # Set the right alternating order between real and imag wavelet
    if level % 2 == 1:
        wavelet = reversed(wavelet)
        
    # In the case of level = 1, coeffs[1:-1] is an empty list. Then, the following
    # loop is not run since zip() iterates as long as the shortest iterable.
    for detail, wav in zip(detail_coeffs, cycle(wavelet)):
        approx = _normalize_size_axis(approx, detail, axis = axis)
        approx = idwt(cA = approx, cD = detail, wavelet = wav, mode = mode, axis = axis)
    
    approx = _normalize_size_axis(approx, first_stage_detail, axis = axis)
    return idwt(cA = approx, cD = first_stage_detail, wavelet = first_stage, mode = mode, axis = axis)

@lru_cache(maxsize = len(ALL_COMPLEX_WAV))
def dualtree_wavelet(name):
    """
    Returns a complex wavelet suitable for dual-tree cwt from a name.

    Parameters
    ----------
    name : str, {'qshift1', 'qshift2', 'qshift3', 'qshift4', 'kingsbury99'}
        Valid arguments can be found in dualtree.ALL_COMPLEX_WAV
    
    Returns
    -------
    real, imag : pywt.Wavelet objects.
    
    Raises
    ------
    ValueError
        If illegal wavelet name.
    """
    if name not in ALL_COMPLEX_WAV:
        raise ValueError('{} is not associated with any implemented complex wavelet. Possible choices are {}'.format(name, ALL_COMPLEX_WAV))
    
    # This function is simply for now
    if name == 'kingsbury99':
        return kingsbury99()
    
    return _qshift(name)

@lru_cache(maxsize = len(ALL_FIRST_STAGE))
def dualtree_first_stage(wavelet = 'kingsbury99_fs'):
    """
    Returns two wavelets to be used in the dual-tree complex wavelet transform, at the first stage.

    Parameters
    ----------
    wavelet : str or Wavelet
        Wavelet to be shifted for first-stage use. Valid arguments can be found in dualtree.ALL_FIRST_STAGE

    Return
    ------
    wav1, wav2 : Wavelet objects

    Raises
    ------
    ValueError
        If invalid first stage wavelet.
    """
    # Special case, preshifted
    if wavelet == 'kingsbury99_fs':
        return kingsbury99_fs()

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)
    
    if wavelet.name not in ALL_FIRST_STAGE:
        raise ValueError('{} is an invalid first stage wavelet.'.format(wavelet.name))
    
    # extend filter bank with zeros
    filter_bank = [np.array(f, copy = True) for f in wavelet.filter_bank]
    for filt in filter_bank:
        extended = np.zeros( shape = (filt.shape[0] + 2,), dtype = np.float)
        extended[1:-1] = filt
        filt = extended

    # Shift deconstruction filters to one side, and reconstruction
    # to the other side
    shifted_fb = [np.array(f, copy = True) for f in wavelet.filter_bank]
    for filt in shifted_fb[::2]:    # Deconstruction filters
        filt = np.roll(filt, 1)
    for filt in shifted_fb[2::]:    # Reconstruction filters
        filt = np.roll(filt, -1)
    
    return Wavelet(name = wavelet.name, filter_bank = filter_bank), Wavelet(name = wavelet.name, filter_bank = shifted_fb)

@lru_cache(maxsize = len(ALL_QSHIFT))
def _qshift(name):
    """
    Returns a complex qshift wavelet by name.

    Parameters
    ----------
    name : str
        Wavelet to use. Valid arguments can be found in dualtree.ALL_QSHIFT
    
    Returns
    -------
    wav1, wav2 : pywt.Wavelet objects
        real and imaginary wavelet
    
    Raises
    ------ 
    ValueError 
        If illegal wavelet family name.
    """
    filters = ('h0a', 'h0b', 'g0a', 'g0b', 'h1a', 'h1b', 'g1a', 'g1b')
    
    filename = join(DATADIR, name + '.npz')
    with np.load(filename) as mat:
        try:
            (dec_real_low, dec_imag_low, rec_real_low, rec_imag_low, 
             dec_real_high, dec_imag_high, rec_real_high, rec_imag_high) = tuple([mat[k].flatten() for k in filters])
        except KeyError:
            raise ValueError('Wavelet does not define ({0}) coefficients'.format(', '.join(filters)))
    
    real_filter_bank = [dec_real_low, dec_real_high, rec_real_low, rec_real_high]
    imag_filter_bank = [dec_imag_low, dec_imag_high, rec_imag_low, rec_imag_high]

    return Wavelet(name = 'real:' + name, filter_bank = real_filter_bank), Wavelet(name = 'imag:' + name, filter_bank = imag_filter_bank)

#############################################################################################
#                           EXAMPLE COMPLEX WAVELETS FROM
#               http://eeweb.poly.edu/iselesni/WaveletSoftware/dt1D.html 
#############################################################################################
@lru_cache(maxsize = 1)
def kingsbury99_fs():
    """
    Returns a first-stage complex wavelet as published in [1]. 

    References
    ----------
    [1] Kingsbury, N. 'Image Processing with Complex Wavelets'. Philocophical Transactions of the Royal Society A pp. 2543-2560, September 1999
    """
    real_dec_lo = np.array([0, -0.08838834764832, 0.08838834764832, 0.69587998903400,0.69587998903400, 0.08838834764832, -0.08838834764832, 0.01122679215254, 0.01122679215254, 0])
    real_dec_hi = np.array([0, -0.01122679215254, 0.01122679215254, 0.08838834764832, 0.08838834764832, -0.69587998903400, 0.69587998903400, -0.08838834764832, -0.08838834764832, 0])
    real_rec_lo, real_rec_hi = real_dec_lo[::-1], real_dec_hi[::-1]

    imag_dec_lo = np.array([0.01122679215254, 0.01122679215254, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 0.08838834764832, -0.08838834764832, 0, 0])
    imag_dec_hi = np.array([0, 0, -0.08838834764832, -0.08838834764832, 0.69587998903400, -0.69587998903400, 0.08838834764832, 0.08838834764832, 0.01122679215254, -0.01122679215254])
    imag_rec_lo, imag_rec_hi = imag_dec_lo[::-1], imag_dec_hi[::-1]

    real_fb = [real_dec_lo, real_dec_hi, real_rec_lo, real_rec_hi]
    imag_fb = [imag_dec_lo, imag_dec_hi, imag_rec_lo, imag_rec_hi]

    return Wavelet(name = 'real:', filter_bank = real_fb), Wavelet(name = 'imag:', filter_bank = imag_fb)

@lru_cache(maxsize = 1)
def kingsbury99():
    """
    Returns a late-stage complex wavelet as published in [1].

    References
    ----------
    [1] Kingsbury, N. 'Image Processing with Complex Wavelets'. Philocophical Transactions of the Royal Society A pp. 2543-2560, September 1999
    """
    real_dec_lo = np.array([ 0.03516384000000, 0, -0.08832942000000, 0.23389032000000, 0.76027237000000, 0.58751830000000, 0, -0.11430184000000, 0, 0])
    real_dec_hi = np.array([0, 0, -0.11430184000000, 0, 0.58751830000000, -0.76027237000000, 0.23389032000000, 0.08832942000000, 0, -0.03516384000000])
    real_rec_lo, real_rec_hi = real_dec_lo[::-1], real_dec_hi[::-1]

    imag_dec_lo = np.array([ 0, 0, -0.11430184000000, 0, 0.58751830000000, 0.76027237000000, 0.23389032000000, -0.08832942000000, 0, 0.03516384000000])
    imag_dec_hi = np.array([-0.03516384000000, 0, 0.08832942000000, 0.23389032000000, -0.76027237000000, 0.58751830000000, 0, -0.11430184000000, 0, 0])
    imag_rec_lo, imag_rec_hi = imag_dec_lo[::-1], imag_dec_hi[::-1]

    real_fb = [real_dec_lo, real_dec_hi, real_rec_lo, real_rec_hi]
    imag_fb = [imag_dec_lo, imag_dec_hi, imag_rec_lo, imag_rec_hi]

    return Wavelet(name = 'real:', filter_bank = real_fb), Wavelet(name = 'imag:', filter_bank = imag_fb)