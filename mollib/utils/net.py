"""
Utility functions to access files from the internet.
"""

import os
import tempfile
try:
    from urllib.request import URLopener
except ImportError:
    from urllib import URLopener


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
    >>> temp_path = get_or_fetch('index', 'html', 'http://www.google.com/')
    >>> temp_path is not None
    True
    >>> # Try to retrieve a file that doesn't exist. None is returned.
    >>> temp_path = get_or_fetch('21312', 'html', 'http://www.glg222.com/',
    ...                          critical=False)
    >>> temp_path is None
    True
    """
    # Setup the message in case the file was not found
    msg = ("Could not find file or identifier '{}'. A suitable file must be "
           "specified to continue.")

    # See if the file exists at the path already. If it does, return this
    # path
    if os.path.isfile(identifier):
        return identifier

    # If the extensions and url aren't specified, then there's nothing else
    # that can be done: the file cannot be found. Return none
    if extensions is None or urls is None:
        if critical:
            print(msg.format(identifier))
            exit()
        return None

    # Convert urls and extensions to iterators
    if not hasattr(extensions, '__iter___'):
        extensions = (extensions, )
    if not hasattr(urls, '__iter__'):
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

        # See if a cached version exists and whether it should be returned
        if load_cached and os.path.isfile(temp_path):
            return temp_path

        # The file isn't available, try retrieving it. Go through the possible
        # urls
        for url in urls:
            # Check the url and filename. An exception is raised if the
            # response is a 404
            opener = URLopener()
            try:
                opener.retrieve('/'.join((url, filename)), temp_path)
            except:
                continue

            # At this point, the url was successful retrieved. Return its
            # local locations
            return temp_path

    # If a temporary file was not produced, then the file couldn't be
    # found locally or retrieved from the web.

    if critical:
        print(msg.format(identifier))
        exit()

    # In this case, return None
    return None


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
