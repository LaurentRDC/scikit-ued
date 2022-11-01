# -*- coding: utf-8 -*-
__author__ = "Laurent P. René de Cotret"
__email__ = "laurent.renedecotret@mail.mcgill.ca"
__license__ = "GPLv3"
__version__ = "2.1.13"

from .affine import (
    affine_map,
    change_basis_mesh,
    change_of_basis,
    is_basis,
    is_rotation_matrix,
    minimum_image_distance,
    rotation_matrix,
    transform,
    translation_matrix,
    translation_rotation_matrix,
)
from .array_utils import (
    cart2polar,
    cart2spherical,
    complex_array,
    mirror,
    plane_mesh,
    polar2cart,
    repeated_array,
    spherical2cart,
)
from .baseline import (
    available_dt_filters,
    available_first_stage_filters,
    baseline_dt,
    baseline_dwt,
    dt_max_level,
    dtcwt,
    idtcwt,
)
from .eproperties import (
    electron_velocity,
    electron_wavelength,
    interaction_parameter,
    lorentz,
)
from .image import (
    align,
    autocenter,
    azimuthal_average,
    brillouin_zones,
    combine_masks,
    detector_scattvectors,
    ialign,
    isnr,
    itrack_peak,
    bragg_peaks,
    bragg_peaks_persistence,
    mask_from_collection,
    mask_image,
    nfold,
    powder_calq,
    reflection,
    snr_from_collection,
    triml,
    trimr,
)
from .io import diffread, diffshow, dmread, imibread, mibheader, mibread
from .patterson import patterson
from .plot_utils import rgb_sweep, spectrum_colors, spectrum_cmap, indices_to_text
from .potential_map import potential_map, potential_synthesis
from .simulation import (
    affe,
    electrostatic,
    pelectrostatic,
    powdersim,
    structure_factor,
    kinematicsim,
)
from .thin_films import film_optical_coefficients
from .time_series import (
    biexponential,
    exponential,
    with_irf,
    mad,
    nfft,
    nfftfreq,
    register_time_shift,
    register_time_shifts,
    Selection,
    ArbitrarySelection,
    RectSelection,
    DiskSelection,
    RingSelection,
    RingArcSelection,
)
from .voigt import gaussian, lorentzian, pseudo_voigt
