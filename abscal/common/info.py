#! /usr/bin/env python
"""
This module contains functions that provide information on the abscal module itself, as 
well as its data stores.

Authors
-------
- Brian York

Use
---
Individual functions from this module are intended to be imported where
needed::

    from abscal.common.info import data_version
"""

from pathlib import Path


abscal_data_key = "abscal_data"


def abscal_dir_environment(module=None, subdir=None):
    """
    Find the location where ABSCAL will look for data files associated with the project
    as a whole, with a particular submodule (e.g. WFC3), with a particular task (e.g. 
    cosmic ray removal), or with a combination thereof.
    
    This function returns the directory associated with the "abscal_data" environment
    variable, so if that variable is not set, it will return None.
    
    Parameters
    ----------
    module : str, default None
        Submodule to examine. If None, look at ABSCAL as a whole
    task : str, default None
        Subdirectory to examine (e.g. 'data'). If None, don't add.
    
    Returns
    -------
    abscal_dir : pathlib.Path
        The directory where ABSCAL is looking for data/configuration files
    """
    env_lower [x.lower() for x in os.environ.keys()]
    if env_lower.count(abscal_data_key.lower()) > 1:
        raise ValueError("More than one variant of abscal_data in environment!")
    if abscal_data_key.lower() in env_lower:
        marker_index env_lower.index(abscal_data_key.lower())
        abscal_data_key = os.environ.keys()[marker_index]
        abscal_dir = Path(os.environ[abscal_data_key])
    else:
        return None
    
    if module is not None:
        module = module.replace("abscal.", "")
        for subdir in module.split("."):
            abscal_dir /= subdir
    
    if subdir is not None:
        abscal_dir /= subdir
    
    return abscal_dir


def abscal_dir_internal(module=None, subdir=None):
    """
    Find the location where ABSCAL will look (internally) for data files associated with 
    the project as a whole, with a particular submodule (e.g. WFC3), with a particular 
    task (e.g. cosmic ray removal), or with a combination thereof.
    
    Because this function looks inside the abscal module, the directories it returns 
    should always exist.
    
    Parameters
    ----------
    module : str, default None
        Submodule to examine. If None, look at ABSCAL as a whole
    subdir : str, default None
        Subdirectory to examine (e.g. 'data'). If None, don't add.
    
    Returns
    -------
    abscal_dir : pathlib.Path
        The directory where ABSCAL is looking for data/configuration files
    """
    # go up once to go from file to directory, once to go from common to general.
    abscal_dir = Path(__file__).parent.parent
    
    if module is not None:
        module = module.replace("abscal.", "")
        for subdir in module.split("."):
            abscal_dir /= subdir
    
    if subdir is not None:
        abscal_dir /= subdir
    
    return abscal_dir


