#! /usr/bin/env python
"""
This module includes MAST utility functions to retrieve data.

Authors
-------
- Brian York

Use
---
Individual functions from this module are intended to be imported where
needed::

    from abscal.common.mast import download_mast_program
"""

import argparse
import glob
import os
import re
import shutil
import sys
import yaml

import numpy as np

from copy import deepcopy
from datetime import datetime
from distutils.util import strtobool
from astroquery.mast import Observations
from ruamel.yaml import YAML
from simpleeval import simple_eval


mast_instrument_prefixes = {'wfc3': 'i', 'stis': 'o'}
exp_str = "[a-zA-Z0-9]{9}_[a-zA-Z0-9]{3}\\.fits"


def download_mast_program(proposal_ids, download_dir, instruments='all', exts='all', force=False, verbose=False):
    """
    Download the HST data related to a MAST program, with the option to choose files of 
    only a given instrument (or instruments), and to choose files of only a given 
    extension (or extensions). Optionally, only download files that are not already 
    present on disk.
    
    Parameters
    ----------
    proposal_ids : str or int
        The proposal or proposals to download. Single value or comma-separated list.
    download_dir : str
        The directory that the data should be downloaded to
    instruments: str, default 'all'
        The instruments to download, single value or comma-separated list
    exts: str, default 'all'
        The extensions to download, single value or comma-separated list. Because abscal
        uses support (spt) files for almost all reduction paths, these files will always
        be included.
    force : bool, default False
        Force the re-downloading of files already in "download_dir"?
    """
    if verbose:
        print("Retrieving program(s) {}".format(proposal_ids))
        print("Retrieving to {}".format(download_dir))
        print("Retrieving instrument(s) {} and extension(s) {}".format(instruments, exts))
    proposal_ids = str(proposal_ids).split(",")
    result_table = Observations.query_criteria(proposal_id=proposal_ids, project='hst')
#     if verbose:
#         print("Result Table:")
#         print(result_table)
#     obs_mask = [x[:4] != 'hst_' for x in flux_table['obs_id']]
#     obs_filter = [id for id in flux_table['obs_id'] if id[:4] != 'hst_']
#     flux_table = flux_table[obs_mask]
    obs_ids = list(result_table['obs_id'])
    print("obs_ids list: {}".format(obs_ids))
    if not force:
        i = 0
        while i < len(obs_ids):
            # This is an idiom for going through a list and potentially removing items 
            # from it. If you're going through a list in a for loop, the result of 
            # removing items from the list during the loop is not well-defined, so this 
            # method is used instead.
            if len(glob.glob(os.path.join(download_dir, obs_ids[i]+"*"))) > 0:
                obs_ids.remove(obs_ids[i])
            else:
                i += 1
    mask = [x in obs_ids for x in result_table['obs_id']]
    result_table = result_table[mask]
    if verbose:
        print("Masked Result Table")
        print(result_table)
    if len(result_table) == 0:
        if verbose:
            print("No records left to retrieve")
        return
    data_products = Observations.get_product_list(result_table[mask])
    mask = [x in obs_ids for x in data_products['obs_id']]    
    data_products = data_products[mask]
    mask = [x[-5:] == ".fits" for x in data_products['productFilename']]
    data_products = data_products[mask]
    if exts != 'all':
        mask = np.array([False for x in data_products['productFilename']])
        for ext in exts.split(",") + ["spt", "asn", "wav"]:
            mask = np.bitwise_or(mask, np.array(["{}.fits".format(ext) in x for x in data_products['productFilename']]))
        data_products = data_products[mask]
    if verbose:
        print("Data Products:")
        print(data_products)
    if len(data_products) == 0:
        if verbose:
            print("No files left to retrieve")
        return
    manifest = Observations.download_products(data_products, download_dir=download_dir)
    if verbose:
        print("Manifest:")
        print(manifest)
    for file_name in manifest['Local Path']:
        base_file = os.path.basename(file_name)
        if verbose:
            print("Copying {}".format(base_file))
        shutil.copy(os.path.join(download_dir, file_name), os.path.join(download_dir, base_file))
    shutil.rmtree(os.path.join(download_dir, "mastDownload"))
  

def main():
    """
    Main function that sets up arguments and runs the download function. Included so that
    the `download_mast` command can just directly invoke the interface setup here.
    """
    description = "Download MAST data for chosen program(s)."
    
    prop_help = "Comma-separated list of proposals to download"
    
    path_help = "Download directory (default $CWD)."
    
    ins_help = "Download files for specific instrument(s). Comma-separated list of "
    ins_help += "instruments to download (default 'all')."
    
    exts_help = "Download specific calibration extension(s). Comma-separated list of "
    exts_help += "extensions to download (default 'all')."

    force_help = "Force the re-downloading of files that are already in the download path. (default False)"

    verbose_help = "Print diagnostic information while running."

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("proposal_ids", help=prop_help)
    parser.add_argument("-p", "--path", dest='path', help=path_help, default=os.getcwd())
    parser.add_argument("-i", "--instruments", dest="instruments", help=ins_help, default='all')
    parser.add_argument("-e", "--exts", dest="exts", help=exts_help, default='all')
    parser.add_argument('-f', '--force', dest='force', action='store_true', 
                        default=False, help=force_help)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                        default=False, help=verbose_help)
    
    res = parser.parse_args()
    
    download_mast_program(res.proposal_ids, res.path, instruments=res.instruments, 
                          exts=res.exts, force=res.force, verbose=res.verbose)


if __name__ == "__main__":
    main()
