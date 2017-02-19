from .markdown import MDTable
from .formatted_str import FormattedStr
from .data_types import Datum
from . import settings

# Register the settings
from mollib.core import register_settings
register_settings(settings)