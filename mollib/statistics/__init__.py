from .statistics import Statistics, BuildData
from .ramachandran_statistics import RamachandranStatistics
from . import settings

# Load the plugins
from mollib.core import register_settings
register_settings(settings)