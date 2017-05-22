"""
The Partial Alignment (pa) module is use to align residual dipolar couplings
(RDCs) and residual anisotropic chemical shifts (RACSs) measured by a partially
aligned sample by NMR to one or more molecular structure(s).
"""
from .data_types import RDC, RACS
from .process_molecule import Process
from .data_readers import read_pa_file
from .svd import calc_pa_SVD
from .analysis import find_outliers, calc_summary
from .reports import report_tables
from . import settings

# Load plugins
from .plugin import PA
pa = PA()

from mollib.core import register_settings
register_settings(settings)

# TODO: Implement a format flag to select the '-a' format
# TODO: add fixer for methyl groups
# TODO: Round re-scaled RDC/RACS in the NHScaleFixer
# TODO: Add fixer to flatten or rename subunits
