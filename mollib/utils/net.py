"""
Utility functions to access files from the internet.
"""

import os
import tempfile
import shutil
import logging
import ssl
try:
    from urllib2 import urlopen, URLError
except ImportError:
    from urllib.request import urlopen, URLError

from . import settings

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


def get_or_fetch(identifier, extensions=None, urls=None, load_cached=True,
                 critical=True):
    """Retrieve a local file path, either from the specified path of from a
    local copy retrieved from the internet.
    
    Parameters
    ----------
    identifier: str
        The filename or identifier (without extension) for the file to 
        download.
    extensions: str or tuple of str (optional)
        The extension(s) to try for the files to download.
    urls: str or tuple of str (optional)
        The urls to try, in order, for the files to download.
    load_cached: bool (optional)
        If True, the path of a cached file will be used if it's available.
    critical: bool (optional)
        If True, then a file not found is considered a critical problem and
        will cause the program to exit.
    
    Returns
    -------
    path: str or None
        The path to the locally stored file. If the file couldn't be found
        locally or downloaded from the web, then None is returned.
    
    Examples
    --------
    >>> # Retrieve a file that exists    
    >>> temp_path = get_or_fetch('index', 'html', 'http://www.google.com/',
    ...                          critical=False)
    >>> temp_path is not None
    True
    >>> # Try to retrieve a file that doesn't exist. None is returned.
    >>> temp_path = get_or_fetch('21312', 'html', 'http://www.google.com/',
    ...                          critical=False)
    >>> temp_path is None
    True
    """
    # Setup the SSL context
    global ctx

    # Setup the message in case the file was not found
    msg = ("Could not find file or identifier '{}'. A suitable file must be "
           "specified to continue.")

    # See if the file exists at the path already. If it does, return this
    # path
    if os.path.isfile(identifier):
        # See if a local copy should be saved
        _save_locally(identifier)

        # Return the identifier filename
        return identifier

    # If the extensions and url aren't specified, then there's nothing else
    # that can be done: the file cannot be found. Return none
    if extensions is None or urls is None:
        if critical:
            print(msg.format(identifier))
            exit()
        return None

    # Convert urls and extensions to iterators
    if isinstance(extensions, str):
        extensions = (extensions, )
    if isinstance(urls, str):
        urls = (urls, )

    # At this point, the file couldn't be found at a local path. In this case,
    # see if it can be retrieved from a locally cached version or retrieved
    # from the web.

    # Setup the temporary directory
    temp_dir = os.path.join(tempfile.gettempdir(), 'mollib')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Go through the possible filenames
    for ext in extensions:
        # Get the filename
        filename = '.'.join((identifier, ext))

        # Get the path for the temporary file
        temp_path = os.path.join(temp_dir, filename)

        # See if a cached version exists and whether it's a valid file. If so,
        # return it.
        if load_cached and _is_valid(temp_path):
            # See if a local copy should be saved
            _save_locally(temp_path)

            # Return the cached copy
            return temp_path

        # The file isn't available, try retrieving it. Go through the possible
        # urls
        for url in urls:
            # Check the url and filename. An exception is raised if the
            # response is a 404
            try:
                response = urlopen('/'.join((url, filename)), context=ctx)
            except URLError:
                continue

            # At this point, the url was successful retrieved. Save it to a
            # temporary path.
            with open(temp_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)

            # See if a copy in the local/present should be saved
            _save_locally(temp_path)

            # Return its local locations
            return temp_path

    # If a temporary file was not produced, then the file couldn't be
    # found locally or retrieved from the web.

    if critical:
        print(msg.format(identifier))
        exit()

    # In this case, return None
    return None


def _save_locally(filepath):
    """Saves the file locally, if specified in the settings, and logs 
    information on the saved file.
    
    Parameters
    ----------
    filepath: str
        The path of the cached file to save locally.
    
    Returns
    -------
    bool
        True if the file was saved locally, False otherwise.
    """
    if settings.save_fetched_files_locally:
        filename = os.path.basename(filepath)
        msg = "Saving file '{}' to the current directory."
        logging.info(msg.format(filename))

        # copy the file locally
        shutil.copy2(filepath, '.')
        return True
    return False


def _is_valid(filepath):
    """Test whether the file at the given filepath is a valid file.
    
    Parameters
    ----------
    filepath: str
        The path of the cached file to test.
    
    Returns
    -------
    bool
        True if the file is valid, False otherwise.
    """
    try:
        # See if the file contains data.
        return os.stat(filepath).st_size > 0
    except OSError:
        # File not found.
        return False



def clear_cache():
    """Clears the temporary cache for mollib."""
    temp_path = os.path.join(tempfile.gettempdir(), 'mollib')
    for filename in os.listdir(temp_path):
        filepath = os.path.join(temp_path, filename)
        try:
            if os.path.isfile(filepath):
                os.unlink(filepath)
        except Exception as e:
            raise e
