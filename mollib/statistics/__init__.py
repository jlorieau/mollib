"""Plugin to calculate statistic datasets from molecular structures."""

from .statistics import Statistics, BuildData
from .ramachandran_statistics import RamachandranStatistics
from .hbond_statistics import HbondStatistics
from . import settings

# Load the plugins
from mollib.core import register_settings
register_settings(settings)