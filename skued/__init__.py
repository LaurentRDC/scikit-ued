# -*- coding: utf-8 -*-
__author__ = 'Laurent P. René de Cotret'
__email__ = 'laurent.renedecotret@mail.mcgill.ca'
__license__ = 'MIT'
__version__ = '0.4.8' # TODO: automatic versioning?

from .affine import (affine_map, change_basis_mesh, change_of_basis, is_basis,
                     is_rotation_matrix, minimum_image_distance,
                     rotation_matrix, transform, translation_matrix,
                     translation_rotation_matrix)
from .array_utils import mirror, repeated_array, cart2polar, polar2cart
from .baseline import baseline_dt, baseline_dwt, dtcwt, idtcwt
from npstreams import chunked, last, linspace, multilinspace
from npstreams import pmap, preduce
from .plot_utils import spectrum_colors, rgb_sweep
from .quantities import electron_wavelength, interaction_parameter, lorentz
from .structure import Crystal, Atom, Lattice
from .voigt import gaussian, lorentzian, pseudo_voigt