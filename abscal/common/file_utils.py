"""
This module includes data file utility functions.

Authors
-------
- Brian York

Use
---
Individual functions from this module are intended to be imported where
needed::

    from abscal.common.file_utils import get_data_file
"""

import os
from pathlib import Path


INTERNAL_PATH = Path(__file__).parent.parent


def get_base_data_dir():
    """
    Find the location of the ABSCAL data files.
    
    ABSCAL stores both default parameters and exposure-specific corrections and settings 
    in a number of data files. There are default copies stored internally, but they can 
    also be stored elsewhere. The ABSCAL_DATA environment variable points to that 
    location, although there is always a fallback to local data if a specified file does 
    not exist elsewhere.
    
    Returns
    -------
    data_dir : str
        The directory pointed to by ABSCAL_DATA (if both the environment variable and 
        the directory to which it points exist)
    """
    if "ABSCAL_DATA" in os.environ:
        # First look for all-capitals
        data_dir = Path(os.environ["ABSCAL_DATA"])
        if data_dir.is_dir():
            return data_dir
    elif "abscal_data" in os.environ:
        # Next look for all lower-case
        data_dir = Path(os.environ["abscal_data"])
        if data_dir.is_dir():
            return data_dir
    
    # Fall back to internal
    return INTERNAL_PATH


def get_data_path(module, subdir=None):
    """
    Get the relative path from the base data directory (whatever it may be) to a
    particular data directory.

    Parameters
    ----------
    module : str
        The module to search in, using standard dot separators (e.g.
        abscal.wfc3)
    subdir : str, default None
        Subdirectory (within the data directory)

    Returns
    -------
    data_path : pathlib.Path
        Relative path to the data directory.
    """
    # Create a relative path from 
    module = module.replace("abscal.", "")
    module_components = module.split(".")
    if len(module_components) == 0:
        data_path = Path("data")
    elif len(module_components) == 1:
        data_path = Path(module_components[0]) / "data"
    else:
        data_path = Path(module_components[0])
        for comp in module_components[1:]:
            data_path /= comp
        data_path /= "data"
    if subdir is not None:
        data_path /= subdir
    return data_path


def check_data_path(module, subdir=None, fname=None):
    """
    Take a relative path, and, in order:
    
    - If the relative path exists from CWD, return that
    - If the relative path exists from `base_data_dir()`, return that
    - If the relative path exists in internal data, return that
    - If none of the above, return None

    Parameters
    ----------
    module : str
        The module to search in, using standard dot separators (e.g.
        abscal.wfc3)
    subdir : str, default None
        Subdirectory (within the data directory)
    fname : str, default None
        File to add at the end

    Returns
    -------
    data_path : pathlib.Path or None
        The resulting (absolute) path, or None if not found
    """
    # Create a relative path from module and subdir
    module = module.replace("abscal.", "")
    module_components = module.split(".")
    if len(module_components) == 0:
        data_path = Path("data")
    elif len(module_components) == 1:
        data_path = Path(module_components[0]) / "data"
    else:
        data_path = Path(module_components[0])
        for comp in module_components[1:]:
            data_path /= comp
        data_path /= "data"
    if subdir is not None:
        data_path /= subdir
    if fname is not None:
        data_path /= fname

    # Compare the relative path against the canonical options
    if (Path.cwd() / data_path).exists():
        return (Path.cwd() / data_path)
    elif (get_base_data_dir() / data_path).exists():
        return (get_base_data_dir() / data_path)
    elif (INTERNAL_PATH / data_path).exists():
        return (INTERNAL_PATH / data_path)
    return None


def get_data_dir(module, subdir=None, make=False):
    """
    Find an internal data directory.
    
    Returns the path to an internal data directory based on module and (optional) 
    subdirectory. Keeps there from being a lot of repetitive code in order 
    to find data file paths.

    Parameters
    ----------
    module : str
        The module to search in, using standard dot separators (e.g.
        abscal.wfc3)
    subdir : str, default None
        Subdirectory (within the data directory)
    make : bool, default False
        If the data directory doesn't exist, should it be created?

    Returns
    -------
    data_path : str or None
        Full path to the data directory. If no directory is found at the generated path, 
        None will be returned. This is not necessarily a failure state, because (for 
        example) a function may call for a known-issues file even though there are no 
        current known issues (and thus no file to return). This way, the code doesn't need 
        to be changed when there *is* a file, and a user can add a file to their local 
        directory without needing to alter the code, because the code will just 
        transparently find the file.
    """
    data_path = check_data_path(module, subdir=subdir)

    if data_path is None and make:
        data_path = get_data_path(module, subdir)
        data_path.mkdir(parents=True, exist_ok=True)

    return data_path


def get_data_file(module, fname, subdir=None, optional=False):
    """
    Find an internal data file.
    
    Returns the path to a file named `fname` in the data directory of the
    (sub)module named `module`. Keeps there from being a lot of repetitive code in order 
    to find data file paths.

    Parameters
    ----------
    module : str
        The module to search in, using standard dot separators (e.g.
        abscal.wfc3)
    fname : str
        The file name of interest
    subdir : str, default None
        Subdirectory (within the data directory)
    optional : bool, default False
        Is the file optional? If so, no error on not found.

    Returns
    -------
    data_file : str or None
        Full path to the data file. If no file is found at the generated path, None will 
        be returned. This is not necessarily a failure state, because (for example) a 
        function may call for a known-issues file even though there are no current known 
        issues (and thus no file to return). This way, the code doesn't need to be 
        changed when there *is* a file, and a user can add a file to their local directory 
        without needing to alter the code, because the code will just transparently find 
        the file.
    """
    data_file = check_data_path(module, subdir=subdir, fname=fname)
    if data_file is not None:
        return data_file

    # This is a convenience. The exposure parameters file associated with a particular
    # python file is that file's name with ".py" replaced with ".yaml". For example, the 
    # util_grism_cross_correlate.py exposure parameters file would be 
    # util_grism_cross_correlate.yaml.
    #
    # As such, this basically says that, if you're asking for a python file (e.g. by just 
    # sending in __file__ to get_data_file) (recognized by the extension being ".py"), 
    # then search for a ".yaml" file with the same name if the previous search failed.
    if Path(fname).suffix == ".py":
        fname = fname.replace(".py", ".yaml")
        data_file = check_data_path(module, subdir=subdir, fname=fname)
        if data_file is not None:
            return data_file

    if not optional:
        msg = f"ERROR: File {module}.{fname} not found at {data_file} or {local_file}.\n"
        raise FileNotFoundError(msg)
    
    # If nothing was found, return None
    return None
