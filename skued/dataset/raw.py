# -*- coding: utf-8 -*-
"""
Handling UED Datasets
=====================
"""
from abc import abstractmethod, abstractproperty

from . import ExperimentalParameter, MetaRawDataset

class AbstractRawDataset(metaclass = MetaRawDataset):
    """
    Abstract base class for ultrafast electron diffraction data set. 
    RawDatasetBase allows for enforced metadata types and values, 
    as well as a standard interface. For example, AbstractRawDataset
    implements the context manager interface.

    Minimally, the following method must be implemented in subclasses:

        * raw_data

    It is suggested to also implement the following magic methods:

        * __init__ 
        * __exit__
    
    The call signature should remain the same for all overwritten methods.
    """

    # List of valid metadata below
    # Using the ExperimentalParameter allows for automatic registering
    # of the parameters as valid.
    # These attributes can be accessed using the usual property access
    pump_wavelength = ExperimentalParameter('pump_wavelength', int, default = 800)
    fluence =         ExperimentalParameter('fluence', float, default = 0)
    time_zero_shift = ExperimentalParameter('time_zero_shift', float, default = 0)
    temperature =     ExperimentalParameter('temperature', float, default = 293)
    exposure =        ExperimentalParameter('exposure', float, default = 1)
    resolution =      ExperimentalParameter('resolution', tuple, default = (2048, 2048))
    time_points =     ExperimentalParameter('time_points', tuple, default = tuple())
    scans =           ExperimentalParameter('scans', tuple, default = (1,))
    notes =           ExperimentalParameter('notes', str, default = '')

    def __init__(self, source = None, metadata = dict()):
        """
        Parameters
        ----------
        source : object
            Data source, for example a directory or external file.
        metadata : dict or None, optional
            Metadata and experimental parameters. Dictionary keys that are
            not valid metadata, they are ignored. Metadata can also be
            set directly later.

        Raises
        ------
        TypeError : if an item from the metadata has an unexpected type.
        """
        self.source = source

        for k, v in metadata.items():
            if k in self.valid_metadata:
                setattr(self, k, v)
    
    def __repr__(self):
        string = '< RawDataset object. '
        for k, v in self.metadata.items():
            string.join('\n {key}: {value} '.format(key = k, value = v))
        string.join(' >')
        return string
    
    def __enter__(self):
        """ Return `self` upon entering the runtime context. """
        return self
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        """ Raise any exception triggered within the runtime context. """
        # Perform cleanup operations here
        pass

    @property
    def metadata(self):
        """ Experimental parameters and dataset metadata as a dictionary. """
        # This property could be generated from some metadata file
        return {k:getattr(self, k) for k in self.valid_metadata}
    
    @abstractmethod
    def raw_data(self, timedelay, scan = 1, **kwargs):
        """
        Returns an array of the image at a timedelay and scan.
        
        Parameters
        ----------
        timdelay : float
            Acquisition time-delay.
        scan : int, optional
            Scan number. Default is 1.
        
        Returns
        -------
        arr : `~numpy.ndarray`, ndim 2
        
        Raises
        ------
        IOError : Filename is not associated with an image/does not exist.
        """ 
        pass
