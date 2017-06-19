# -*- coding: utf-8 -*-
"""
Physical quantities relevant to other modules\packages
"""
import numpy as np

c = 299792458           #m/s
h = 6.63*10**(-34)      #J*s.
e = 1.602*10**(-19)     #in C
m0 = 9.109*10**(-31)    #in kg

def lorentz(kV):
    """
    Relativistic factor.

    Parameters
    ----------
    kV : array_like or float
        Electron gun voltage [kV].
    
    Returns
    -------
    out : array_like or float
    """
    return np.sqrt(1 + (e*kV*1e3)/(2*m0*c**2))

def electron_wavelength(kV):
    """ 
    Relativistic wavelength of an accelerated electron.
        
    Parameters
    ----------
    kV : float or ndarray
        Voltage in kilovolts of the instrument
    
    Returns
    -------
    out : float
        Electron wavelength [Angs]
    """
    return (h/np.sqrt(2*m0*e*kV*1e3))/lorentz(kV)*1e10

def interaction_parameter(kV):
    """
    Interaction parameter from relativistic electron wavelength.

    Parameters
    ----------
    kV : float
        Electron gun voltage [kV].
    
    Returns
    -------
    out : float
        Interaction parameter [rad/(V*Angs)]

    References
    ----------
    .. Kirkland 2010 Eq. 5.6
    """
    l = electron_wavelength(kV)
    V = kV * 1e3

    return (2*np.pi)/(electron_wavelength(kV)*V)*(m0*c**2 + e * V)/(2*m0*c**2 + e * V)