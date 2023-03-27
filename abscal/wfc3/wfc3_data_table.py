#! /usr/bin/env python
"""
This module contains the WFC3DataTable class. The class is derived from the common 
AbscalDataTable class, with modifications to deal with the specific header keywords that
are relevant to WFC3 exposures.

Authors
-------
    - Brian York

Use
---

This module provides the WFC3DataTable class, which holds an astropy table
of exposure metadata, along with the ability to present that data either in
a form expected by IDL, or in a form more suitable for passing on to
later scripts in abscal::

    from abscal.stis.stis_data_table import WFC3DataTable
    t = WFC3DataTable()
    
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
from distutils.util import strtobool
from pathlib import Path
from simpleeval import simple_eval

from abscal.common.utils import build_expr
from abscal.common.exposure_data_table import AbscalDataTable


def scan_rate_formatter(scan_rate):
    """
    Format the IDL table 'SCAN_RAT' column.
    
    In IDL at least, the scan rate column has a very particular format (0.0000 for cases 
    where there is no scan being done, left-aligned). There's no good one-step way of 
    doing that formatting in python, so here's what we're stuck with.
    
    Parameters
    ----------
    scan_rate : float
        The scan rate to format
    
    Returns
    -------
    scan_rate_str : str
        The scan rate, formatted as a four-decimal-place floating point, and left-aligned.
    """
    scan_rate_str = "{:6.4f}".format(scan_rate)
    return "{:<9}".format(scan_rate_str)


class WFC3DataTable(AbscalDataTable):
    """
    A class to represent HST WFC3 IR grism exposure metadata for ABSCAL
    """

    @classmethod
    def from_idl_table(cls, table_file):
        """
        Adds a call to the `set_filter_images()` function particular to WFC3 data, after
        importing the table in general. Parameters and return values are as for the 
        superclass.
        """
        table = super().from_idl_table(cls, table_file)

        # Now that we've made the table, set the filter images if possible.
        if cls.instrument == 'wfc3ir':
            table.set_filter_images()
        
        return table
    
    
    def set_filter_images(self):
        """
        For any grism exposures, look for the nearest associated filter
        exposure and, if one is found, associate it using the filter_root
        column. A filter image is considered associated if it:
        
        - is from the same program
        - is from the same visit
        - has the same POSTARG values
        
        If there is more than one appropriate exposure, take the one that's closest in 
        time.
        """
        filter_table = deepcopy(self)
        filter_list = [f[0] != 'F' for f in filter_table['filter']]
        filter_table.remove_rows(filter_list)
        
        grism_table = deepcopy(self)
        grism_list = [f[0] != 'G' for f in grism_table['filter']]
        grism_table.remove_rows(grism_list)
        
        for row in grism_table:
            root = row['root']
            obset = row['obset']
            postarg1 = row['postarg1']
            postarg2 = row['postarg2']
            base_mask = [r == root for r in self['root']]
            
            check_mask = [(o == obset) & (p1 == postarg1) & (p2 == postarg2) for \
                          o,p1,p2 in zip(filter_table['obset'], \
                          filter_table['postarg1'], filter_table['postarg2'])]
            check_table = filter_table[check_mask]
            
            if len(check_table) > 0:
                check_dates = [abs(row['date']-t) for t in check_table['date']]
                minimum_time_index = check_dates.index(min(check_dates))
                minimum_time_row = check_table[minimum_time_index]
                self["filter_root"][base_mask] = minimum_time_row["root"]
            else:
                msg = "No corresponding filter exposure found."
                self["notes"][base_mask] += " {}".format(msg)
                self["filter_root"][base_mask] = "NONE"
        return
       
    
    @staticmethod
    def _filter_table(table, filter):
        """
        Applies WFC3-IR-specific filter logic to the general filtering function. Parameters
        and return values are as per the superclass.
        """
        table = super()._filter_table(table, filter)
        
        if isinstance(filter, str):
            if filter == 'stare':
                table.remove_rows(table['scan_rate'] > 0.)
            elif filter == 'scan':
                table.remove_rows(table['scan_rate'] == 0.)
            elif filter == 'filter':
                filter_list = [f[0] != 'F' for f in table['filter']]
                table.remove_rows(filter_list)
            elif filter == 'grism':
                table = WFC3DataTable._grism_filter(table)
            else:
                print("Warning: Unknown filter {}".format(filter))
        else:
            print("Warning: Unknown filter {}".format(filter))
            table = filter(table)

        return table


    @staticmethod
    def _grism_filter(table):
        """
        Filter the table as follows:
        
        - keep all grism exposures
        - for each grism exposure,
        
          - if there is at least one filter exposure from the same program and visit, 
            keep the filter exposure closest in time to the grism exposure
          - else annotate the grism that no corresponding filter was found
        
        In order to do this, the function uses the above set_filter_images function.
        
        Parameters
        ----------
        table : astropy.table.Table
            The table to filter.
        
        Returns
        -------
        filtered_table : astropy.table.Table
            The filtered table.
        """
        grism_table = deepcopy(table)
        grism_table.set_filter_images()
        filter_images = [f[0] != 'G' for f in grism_table['filter']]
        grism_images = [f not in grism_table['filter_root'] for f in grism_table['root']]
        grism_list = np.array(filter_images) & np.array(grism_images)
        grism_table.remove_rows(grism_list)
        
        return grism_table
    
    
    column_mappings = {
                        "ROOT": "root",
                        "MODE": "filter",
                        "APER": "aperture",
                        "TYPE": "exposure_type",
                        "TARGET": "target",
                        # "IMG SIZE" handled specially
                        # DATE handled specially
                        # TIME handled specially
                        "PROPID": "proposal",
                        "EXPTIME": "exptime",
                        # "POSTARG X,Y" handled specially
                        "SCAN_RAT": "scan_rate"
                      }
    
    default_format = 'ascii.ipac'
    
    instrument = 'wfc3ir'
    default_search_str = 'i*flt.fits'
    idl_str = 'WFCDIR'
    idl_columns = (0, 11, 19, 35, 41, 54, 64, 73, 82, 89, 98, 112)
    idl_exts = ['flt', 'ima', 'raw']
    idl_column_formats = {
                            "ROOT": "<10",
                            "MODE": "<7",
                            "APER": "<15",
                            "TYPE": "<5",
                            "TARGET": "<12",
                            "IMG SIZE": "<9",
                            "DATE": "<8",
                            "TIME": "<8",
                            "PROPID": "<6",
                            "EXPTIME": "7.1f",
                            "POSTARG X,Y": ">14",
                            "SCAN_RAT": scan_rate_formatter
                         }
    
    standard_columns = {
                            'root': {'dtype': 'O', 'idl': True},
                            'obset': {'dtype': 'S6', 'idl': True},
                            'filter': {'dtype': 'O', 'idl': True},
                            'aperture': {'dtype': 'O', 'idl': True},
                            'exposure_type': {'dtype': 'O', 'idl': True},
                            'target': {'dtype': 'O', 'idl': True},
                            'xsize': {'dtype': 'i4', 'idl': True},
                            'ysize': {'dtype': 'i4', 'idl': True},
                            'date': {'dtype': 'O', 'idl': True},
                            'proposal': {'dtype': 'i4', 'idl': True},
                            'exptime': {'dtype': 'f8', 'idl': True},
                            'postarg1': {'dtype': 'f8', 'idl': True},
                            'postarg2': {'dtype': 'f8', 'idl': True},
                            'scan_rate': {'dtype': 'f8', 'idl': True},
                            'path': {'dtype': 'O', 'idl': False, 
                                     'default': 'N/A'},
                            'use': {'dtype': '?', 'idl': False, 
                                    'default': True},
                            'filter_root': {'dtype': 'S10', 'idl': False,
                                            'default': 'unknown'},
                            'filename': {'dtype': 'O', 'idl': False,
                                         'default': 'N/A'},
                            'xc': {'dtype': 'f8', 'idl': False, 'default': -1.},
                            'yc': {'dtype': 'f8', 'idl': False, 'default': -1.},
                            'xerr': {'dtype': 'f8', 'idl': False, 
                                     'default': -1.},
                            'yerr': {'dtype': 'f8', 'idl': False, 
                                     'default': -1.},
                            'crval1': {'dtype': 'f8', 'idl': False, 
                                       'default': 0.},
                            'crval2': {'dtype': 'f8', 'idl': False, 
                                       'default': 0.},
                            'ra_targ': {'dtype': 'f8', 'idl': False, 
                                        'default': 0.},
                            'dec_targ': {'dtype': 'f8', 'idl': False, 
                                         'default': 0.},
                            'extracted': {'dtype': 'O', 'idl': False, 
                                          'default': ''},
                            'coadded': {'dtype': 'O', 'idl': False, 
                                        'default': ''},
                            'wl_offset': {'dtype': 'f8', 'idl': False, 'default': -1.},
                            'planetary_nebula': {'dtype': '?', 'idl': False,
                                                 'default': False},
                            'notes': {'dtype': 'O', 'idl': False, 
                                      'default': 'N/A'},
                       }
