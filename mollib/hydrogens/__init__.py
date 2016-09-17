from . import settings
from .hydrogenate import (add_hydrogens, add_one_sp2_h, add_two_sp2_h,
                          add_one_sp3_h, add_two_sp3_h, add_three_sp3_h)

# Load the plugins
from .plugin import Hydrogenate
hydrogenate = Hydrogenate()

from mollib.core import register_settings
register_settings(settings)