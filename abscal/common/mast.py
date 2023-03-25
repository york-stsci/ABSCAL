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

import glob
import os
import shutil
import sys
import yaml

import numpy as np

from astropy.time import Time
from copy import deepcopy
from datetime import datetime
from distutils.util import strtobool
from astroquery.mast import Observations
from ruamel.yaml import YAML
from simpleeval import simple_eval


def download_mast_program(proposal_id, download_dir, skip_existing_files=True):
    flux_table = Observations.query_criteria(proposal_id=proposal_id)
    obs_mask = [x[:4] != 'hst_' for x in flux_table['obs_id']]
    obs_filter = [id for id in flux_table['obs_id'] if id[:4] != 'hst_']
    flux_table = flux_table[obs_mask]
    obs_ids = flux_table['obsid']
    if skip_existing_files:
        i = 0
        while i < len(obs_filter):
            # This is an idiom for going through a list and potentially removing items from it. If you're going 
            # through a list in a for loop, the result of removing items from the list during the loop is not 
            # well-defined, so this method is used instead.
            if len(glob.glob(os.path.join(download_dir, obs_filter[i]+"*"))) > 0:
                obs_filter.remove(obs_filter[i])
            else:
                i += 1
    data_products = Observations.get_product_list(obs_ids)
    if len(data_products) > 0 and len(obs_filter) > 0:
        manifest = Observations.download_products(data_products, download_dir=download_dir, extension=["fits"], 
                                                  obs_id=obs_filter)

        for file_name in manifest['Local Path']:
            base_file = os.path.basename(file_name)
            print("Copying {}".format(base_file))
            shutil.copy(os.path.join(download_dir, file_name), os.path.join(download_dir, base_file))