# -*- coding: utf-8 -*-
""" Electron properties """

import numpy as np

c = 299792458           # speed of light [m/s]
h = 6.63*10**(-34)      # Planck's constant [J*s]
e = 1.602*10**(-19)     # electron charge [C]
m0 = 9.109*10**(-31)    # electron mass [kg]

def lorentz(keV):
    """
    Relativistic factor.

    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].
    
    Returns
    -------
    out : array_like or float
    """
    return 1/np.sqrt(1 + (e*keV*1e3)/(2*m0*c**2))

def electron_wavelength(keV):
    """ 
    Relativistic wavelength of an accelerated electron.
        
    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].
    
    Returns
    -------
    out : float
        Electron wavelength [Angs]
    """
    return (h/np.sqrt(2*m0*e*keV*1e3))*lorentz(keV)*1e10

def interaction_parameter(keV):
    """
    Interaction parameter from relativistic electron wavelength.

    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].
    
    Returns
    -------
    out : float
        Interaction parameter [rad/(V*Angs)]

    References
    ----------
    .. Kirkland 2010 Eq. 5.6
    """
    l = electron_wavelength(keV)
    V = keV * 1e3

    return (2*np.pi)/(electron_wavelength(keV)*V)*(m0*c**2 + e * V)/(2*m0*c**2 + e * V)