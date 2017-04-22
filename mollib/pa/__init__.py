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


# TODO: add --set option to select an RDC dataset
# TODO: check the 'mr' file for 2OED
# TODO: write test for 2MJB '.mr' file
# TODO: add --summary option to only output summary (no data).
