#! /usr/bin/env python
"""
This module contains the STISDataTable class. The class is derived from the common 
AbscalDataTable class, with modifications to deal with the specific header keywords that
are relevant to STIS exposures.

Authors
-------
    - Brian York

Use
---

This module provides the STISDataTable class, which holds an astropy table
of exposure metadata, along with the ability to present that data either in
a form expected by IDL, or in a form more suitable for passing on to
later scripts in abscal::

    from abscal.stis.stis_data_table import STISDataTable
    t = STISDataTable()
    
    t.add_exposure(<exposure metadata dict>)
    ... # repeat as needed
    t.write_to_file(<file>, <idl_mode>, <other arguments>)
    
The goal is to create an astropy.table.Table subclass that has the ability to successfully 
read in IDL-formatted tables, and to write out IDL-formatted tables *as an additional 
option*, but that generally uses its own format.
"""

import argparse
import datetime
import glob
import io
import os

import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.time import Time
from copy import deepcopy
from datetime import datetime as dt
from pathlib import Path

from abscal.common.exposure_data_table import AbscalDataTable


def mglobal_formatter(mglobal):
    """
    Format the IDL table 'MGLOBAL' column.
    
    This column may be either an integer (possibly 0), or a string of asterisk (*)
    characters.
    
    Parameters
    ----------
    mglobal : string-or-float
        The MGLOBAL value to format
    
    Returns
    -------
    mglobal_str : str
        The MGLOBAL value, formatted as a string.
    """
    if '*' in str(mglobal):
        return "{:<7}"
    return "{:>7}"


def postarg_formatter(postarg):
    """
    Format the IDL table 'POSTARG' column.
    
    This column may be either an floating point or a string.
    
    Parameters
    ----------
    postarg : string-or-float
        The POSTARG value to format
    
    Returns
    -------
    postarg_str : str
        The MGLOBAL value, formatted as a string.
    """
    if 'pos' in str(postarg):
        return "{:<5}"
    return "{:.1f"


class STISDataTable(AbscalDataTable):
    """
    A class to represent STIS exposure metadata for ABSCAL
    """

    column_mappings = {
                        "ROOT": "root",
                        "MODE": "mode",
                        "APER": "aperture",
                        "CENWAV": "central_wavelength",
                        "DETECTOR": "detector",
                        "TARGET": "target",
                        "OBSMODE": "observation_mode",
                        "MGLOBAL": "mglobal",
                        # DATE handled specially
                        # TIME handled specially
                        "PROPID": "proposal",
                        "EXPTIME": "exptime",
                        "CR": "cr",
                        "MINW": "min_wavelength",
                        "MAXW": "max_wavelength",
                        "POSTARG": "postarg"
                      }
    
    default_format = 'ascii.ipac'
    
    default_search_str = 'o*raw.fits'
    idl_str = 'STISDIR'
    idl_columns = (0, 11, 19, 28, 33, 42, 55, 62, 70, 79, 88, 96, 104, 107, 112, 118)        
    idl_exts = ['raw']
    idl_column_formats = {
                            "ROOT": "<10",
                            "MODE": "<7",
                            "APER": "<15",
                            "CENWAV": ">4",
                            "DETECTOR": "<8",
                            "TARGET": "<12",
                            "OBSMODE": "<5",
                            "MGLOBAL": mglobal_formatter,
                            "DATE": "<8",
                            "TIME": "<8",
                            "PROPID": "<6",
                            "EXPTIME": "7.1f",
                            "CR": "<1",
                            "MINW": ">4",
                            "MAXW": ">5",
                            "POSTARG": postarg_formatter
                         }
    
    standard_columns = {
                            'root': {'dtype': 'O', 'idl': True},
                            'obset': {'dtype': 'S6', 'idl': True},
                            'mode': {'dtype': 'O', 'idl': True},
                            'aperture': {'dtype': 'O', 'idl': True},
                            'detector': {'dtype': 'O', 'idl': True},
                            'target': {'dtype': 'O', 'idl': True},
                            'obsmode': {'dtype': 'O', 'idl': True},
                            'mglobal': {'dtype': 'i4', 'idl': True},
                            'date': {'dtype': 'O', 'idl': True},
                            'proposal': {'dtype': 'i4', 'idl': True},
                            'exptime': {'dtype': 'f8', 'idl': True},
                            'cr': {'dtype': 'i4', 'idl': True},
                            'central_wavelength': {'dtype': 'f8', 'idl': True},
                            'min_wavelength': {'dtype': 'f8', 'idl': True},
                            'max_wavelength': {'dtype': 'f8', 'idl': True},
                            'raw_postarg': {'dtype': 'f8', 'idl': False},
                            'postarg': {'dtype': 'O', 'idl': True, 'default': '0.0'},
                            'subarray': {'dtype': '?', 'idl': False, 'default': False},
                            'instrument': {'dtype': 'O', 'idl': False},
                            'gain': {'dtype': 'f8', 'idl': False},

                            'path': {'dtype': 'O', 'idl': False, 'default': 'N/A'},
                            'use': {'dtype': '?', 'idl': False, 'default': True},
                            'filename': {'dtype': 'O', 'idl': False,'default': 'N/A'},
                            'ra_targ': {'dtype': 'f8', 'idl': False, 'default': 0.},
                            'dec_targ': {'dtype': 'f8', 'idl': False, 'default': 0.},
                            'wl_offset': {'dtype': 'f8', 'idl': False, 'default': -1.},
                            'flatfielded': {'dtype': 'O', 'idl': False, 'default': ''},
                            'extracted': {'dtype': 'O', 'idl': False, 'default': ''},
                            'notes': {'dtype': 'O', 'idl': False, 'default': 'N/A'},
                       }
