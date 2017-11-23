# -*- coding: utf-8 -*-
__author__ = 'Laurent P. René de Cotret'
__email__ = 'laurent.renedecotret@mail.mcgill.ca'
__license__ = 'MIT'
__version__ = '0.4.10' # TODO: automatic versioning?

from .affine import (affine_map, change_basis_mesh, change_of_basis, is_basis,
                     is_rotation_matrix, minimum_image_distance,
                     rotation_matrix, transform, translation_matrix,
                     translation_rotation_matrix)
from .array_utils import (cart2polar, cart2spherical, mirror, plane_mesh,
                          polar2cart, repeated_array, spherical2cart)
from .baseline import baseline_dt, baseline_dwt, dtcwt, idtcwt
from .eproperties import (electron_velocity, electron_wavelength,
                          interaction_parameter, lorentz)
from .image import (align, azimuthal_average, combine_masks, diff_register,
                    ialign, isnr, mask_from_collection, mask_image, mnxc2,
                    nfold, powder_center, shift_image, snr_from_collection,
                    triml, trimr)
from .plot_utils import rgb_sweep, spectrum_colors
from .simulation import (affe, bounded_reflections, electrostatic,
                         pelectrostatic, powdersim, structure_factor)
from .structure import Atom, Crystal, Lattice
from .voigt import gaussian, lorentzian, pseudo_voigt
