# -*- coding: utf-8 -*-
"""
Siwick Research Group RawDataset class as an example use
of RawDatasetBase
"""
from contextlib import suppress
from glob import iglob
from os import listdir
from os.path import isdir, join
from re import search

from cached_property import cached_property
from skimage.io import imread

from npstreams import imean, last

from . import AbstractRawDataset


class McGillRawDataset(AbstractRawDataset):
    """
    Raw dataset
    """

    @staticmethod
    def parse_tagfile(path):
        """ Parse a tagfile.txt from a raw dataset into a dictionary of values """
        metadata = dict()
        with open(path) as f:
            for line in f:
                key, value = re.sub('\s+', '', line).split('=')
                try:
                    value = float(value.strip('s'))    # exposure values have units
                except:
                    value = None    # value might be 'BLANK'
                metadata[key.lower()] = value
        
        return metadata

    def __init__(self, directory, **kwargs):
        super().__init__(*args, **kwargs)

        if not isdir(directory):
            raise ValueError('{} does not point to an existing directory'.format(directory))
        self.directory = directory

        # Populate experimental parameters
        # from a metadata file called 'tagfile.txt'
        _metadata = self.parse_tagfile(join(self.directory, 'tagfile.txt'))
        self.fluence = _metadata.get('fluence', 0)
        self.resolution = (2048, 2048)
        self.current = _metadata.get('current', 0)
        self.exposure = _metadata.get('exposure', 0)
        self.energy = _metadata.get('energy', 90)
        
        # Determine acquisition date
        # If directory name doesn't match the time pattern, the
        # acquisition date will be the default value
        with suppress(AttributeError):
            self.acquisition_date = search('(\d+[.])+', self.directory).group()[:-1]      #Last [:-1] removes a '.' at the end

        # To determine the scans and time-points, we need a list of all files
        image_list = [f for f in listdir(self.directory) 
                      if isfile(join(self.directory, f)) 
                      and f.endswith(('.tif', '.tiff'))]

        # Determine the number of scans
        # by listing all possible files
        scans = [search('[n][s][c][a][n][.](\d+)', f).group() for f in image_list if 'nscan' in f]
        self.scans = tuple({int(string.strip('nscan.')) for string in scans})

        # Determine the time-points by listing all possible files
        time_data = [search('[+-]\d+[.]\d+', f).group() for f in image_list if 'timedelay' in f]
        time_list =  list(set(time_data))     #Conversion to set then back to list to remove repeated values
        time_list.sort(key = float)
        self.time_points = tuple(map(float, time_list))
    
    # Only three methods/properties, beyond __init__, need to be
    # implemented to have a valid RawDataset:
    #   * pumpon_background
    #   * pumpoff_background
    #   * raw_data

    @cached_property
    def pumpon_background(self):
        backgrounds = map(imread, iglob(join(self.directory, 'background.*.pumpon.tif')))
        return last(imean(backgrounds))
    
    @cached_property
    def pumpoff_background(self):
        backgrounds = map(imread, iglob(join(self.directory, 'background.*.pumpoff.tif')))
        return last(imean(backgrounds))

    def raw_data(self, timedelay, scan = 1, bgr = True, **kwargs): 
        """
        Returns an array of the image at a timedelay and scan.
        
        Parameters
        ----------
        timedelay : float
            Time-delay in picoseconds.
        scan : int, optional
            Scan number. 
        bgr : bool, optional
            If True, pump-on background is removed from the pattern
            before being returned.
        
        Returns
        -------
        arr : ndarray, shape (N,M)
        
        Raises
        ------
        ImageNotFoundError
            Filename is not associated with an image/does not exist.
        """ 
        #Template filename looks like:
        #    'data.timedelay.+1.00.nscan.04.pumpon.tif'
        sign = '' if float(timedelay) < 0 else '+'
        str_time = sign + '{0:.2f}'.format(float(timedelay))
        filename = 'data.timedelay.' + str_time + '.nscan.' + str(int(scan)).zfill(2) + '.pumpon.tif'

        im = imread(join(self.directory, filename))
        if bgr:
            return uint_subtract_safe(im, self.pumpon_background)
        return im
