# -*- coding: utf-8 -*-
__author__ = 'Laurent P. René de Cotret'
__email__ = 'laurent.renedecotret@mail.mcgill.ca'
__license__ = 'MIT'
__version__ = '1.0.1.1'

from .affine import (affine_map, change_basis_mesh, change_of_basis, is_basis,
                     is_rotation_matrix, minimum_image_distance,
                     rotation_matrix, transform, translation_matrix,
                     translation_rotation_matrix)
from .array_utils import (cart2polar, cart2spherical, complex_array, mirror,
                          plane_mesh, polar2cart, repeated_array,
                          spherical2cart)
from .baseline import (available_dt_filters, available_first_stage_filters,
                       baseline_dt, baseline_dwt, dt_max_level, dtcwt, idtcwt)
from .eproperties import (electron_velocity, electron_wavelength,
                          interaction_parameter, lorentz)
from .image import (align, azimuthal_average, combine_masks, diff_register,
                    ialign, isnr, itrack_peak, mask_from_collection,
                    mask_image, mnxc2, nfold, powder_calq, powder_center,
                    reflection, shift_image, snr_from_collection, triml, trimr,
                    xcorr)
from .io import diffread, diffshow, dmread, imibread, mibheader, mibread
from .plot_utils import rgb_sweep, spectrum_colors
from .potential_map import potential_map
from .simulation import (affe, bounded_reflections, electrostatic,
                         pelectrostatic, powdersim, structure_factor)
from .structure import (Atom, AtomicStructure, Crystal, Lattice, LatticeSystem,
                        lattice_system, symmetry_expansion)
from .thin_films import film_optical_coefficients
from .time_series import (biexponential_decay, exponential_decay, mad, nfft,
                          nfftfreq, register_time_shift, register_time_shifts,
                          time_shift, time_shifts)
from .voigt import gaussian, lorentzian, pseudo_voigt
