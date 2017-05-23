"""
Utility functions for reading and writing files.
"""
import os
import tempfile

from .formatted_str import FormattedStr
from .exceptions import ParameterError
from . import settings


def write_file(text, filepath=None, temporary=False, overwrite=None):
    """Write to an ASCII file.
    
    Parameters
    ----------  
    text: str
        The text to write to the file
    filepath: str (optional)
        The path and filename of the file to write to. If the path
        directory doesn't exist, it will be created. Either the filepath or the
        temporary argument should be specified.
    temporary: bool (optional)
        Write to a temporary file. Either the filepath or temporary should be
        specified
    overwrite: bool (optional)
        If True, the file will be overwritten if it exists.
        If False, a number for a non-existent filename will be appended after
        the base filename. ex: "text.txt" will become "text.1.txt".
    """
    # Check the pass parameters to make sure they're consisten
    if filepath is not None and not temporary:
        pass  # Use the filepath as is
    elif filepath is None and temporary:
        _, filepath = tempfile.mkstemp()
    else:
        msg = ("write_file: either a filepath must be specified of temporary "
               "must be True--not both.")
        raise ParameterError(msg)

    # If the overwrite option is not specified, get its value from the settings
    if overwrite is None:
        overwrite = settings.overwrite_files

    # Break the filepath into a path and filename
    basepath, filename = os.path.split(filepath)

    # Check to see whether the basepath exists, if not, create it
    if not os.path.isdir(basepath):
        os.makedirs(basepath)

    # Check to see if the filename exists, rename it if overwrite is False
    if os.path.exists(filepath) and not overwrite:
        # Find which filename does not exist
        file_pieces = os.path.splitext(filename)

        # Set the maximum number of file versions
        max_versions = settings.max_file_versions

        for i in range(1, max_versions + 1):
            # If the max_versions has been reached, raise an exception
            if i == max_versions:
                msg = ("The maximum number of version ({}) has been reached "
                       "for '{}'")
                raise OSError(msg.format(max_versions, filepath))

            # Try a new filename with the next version
            new_filename = '.' + str(i)
            new_filename = new_filename.join(file_pieces)

            # See if the the new_filepath exists, if it does, keep going in the
            # cycle
            new_filepath = os.path.join(basepath, new_filename)

            if not os.path.exists(new_filepath):
                filename = new_filename
                break

    # At this point, we are allowed to write to the basepath and filename
    new_filepath = os.path.join(basepath, filename)

    # Prepare the string to write. If it's a formatted string, then strip it
    # of its formatting.
    text = FormattedStr(text).stripped_str()

    with open(new_filepath, 'w') as f:
        f.write(text)
        return new_filepath
