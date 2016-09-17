from .markdown import MDTable
from . import settings

# Register the settings
from mollib.core import register_settings
register_settings(settings)