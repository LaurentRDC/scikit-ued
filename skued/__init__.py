# -*- coding: utf-8 -*-
__author__ = 'Laurent P. René de Cotret'
__email__ = 'laurent.renedecotret@mail.mcgill.ca'
__license__ = 'MIT'
__version__ = '0.4.1' # TODO: automatic versioning

from .array_utils import repeated_array
from .parallel import pmap, preduce
from .plot_utils import spectrum_colors
from .quantities import lorentz, electron_wavelength, interaction_parameter
from .transformations import (affine_map, transform, change_of_basis, is_basis, translation_matrix,
							  is_rotation_matrix, rotation_matrix, translation_rotation_matrix,
							  change_basis_mesh, minimum_image_distance)
from .voigt import gaussian, lorentzian, pseudo_voigt