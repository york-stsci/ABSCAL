#! /usr/bin/env python
"""
This submodule takes an input metadata table and extracts STIS spectra from the files in 
the table. For each exposure, it does the following:

- For CCD exposures,

    - Subtract the bias level
    - Subtract the bias image
    - Do cosmic ray rejection with standard paramemters
    
        - Iterate cosmic ray parameters until the user is satisfied that cosmic rays have
          been removed without removing any of the spectral trace

- For MAMA exposures

    - Lo-Res Conversion
    - Data Initialization

- Subtract dark
- Subtract flatfield

Authors
-------
- Brian York (all python code)
- Ralph Bohlin (original IDL code)
- STIS team (calstis code imported and used)

Use
---

This submodule is intended to be run as part of the ABSCAL STIS reduction suite. In order
to run it, you can use the form::

    from abscal.stis.reduce_extract import reduce
    
    output_table = reduce(input_table, command_line_arg_namespace, override_dict)

The override dict allows for many of the default input parameters to be overriden (as 
defaults -- individual per-exposure overrides defined in the data files will still take 
priority). Parameters that can be overriden in coadd are:

xc: default -1
    X centre of zeroth order image. If set to a negative value, the submodule will find 
    and fit the centre itself, either from a corresponding filter exposure (preferred) or 
    from the grism exposure directly.
"""

import datetime
import glob
import os
import shutil
import yaml

import matplotlib.pyplot as plt
import numpy as np

from astropy import wcs
from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.time import Time
from copy import deepcopy
from crds import assign_bestrefs
from matplotlib.colors import LogNorm
from matplotlib.widgets import TextBox, Button
from pathlib import Path
from stistools import basic2d
from stistools import ocrreject
from stistools import x1d

from abscal.common.args import parse
from abscal.common.utils import find_value
from abscal.common.utils import get_data_file
from abscal.common.utils import get_defaults
from abscal.common.utils import handle_parameter
from abscal.common.utils import initialize_value
from abscal.common.utils import set_params
from abscal.common.utils import set_image
from abscal.stis.stis_data_table import STISDataTable

def reduce_flatfield(input_table, **kwargs):
    """
    Perform the calstis "basic2d" data reduction, including iterating on cosmic ray
    rejection. This is done by repeatedly calling the calstis "basic2d" on the raw STIS
    data with various calibration flag settings, and by calling "ocrreject" separately in 
    order to have finer control over cosmic ray rejection parameters.

    Parameters
    ----------
    input_table : abscal.stis.stis_data_table.STISDataTable
        Single-row input

    Returns
    -------
    final_file : str
        Name (including path) of the CR-rejected 2D image.
    """
    task = "stis_reduce_flatfield"
    cr_flag_value = 8192
    root = input_table['root']
    mode = input_table['mode']
    target = input_table['target']
    detector = input_table['detector']
    path = input_table['path']
    fname = input_table['filename']
    raw_file = os.path.join(path, fname)
    final_file = raw_file.replace("_raw", "_{}_{}_2d".format(target, mode))
    preamble = "{}: {}: {} ({})".format(task, root, detector, mode)
    verbose = kwargs.get("verbose", False)

    issues = {}
    param_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if param_file is not None:
        with open(param_file, 'r') as inf:
            issues = yaml.safe_load(inf)

    # stistools.basic2d.basic2d(input, 
    #                           output='', 
    #                           outblev='', 
    #                           dqicorr='perform', 
    #                           atodcorr='omit', 
    #                           blevcorr='perform', 
    #                           doppcorr='perform', 
    #                           lorscorr='perform', 
    #                           glincorr='perform', 
    #                           lflgcorr='perform', 
    #                           biascorr='perform', 
    #                           darkcorr='perform', 
    #                           flatcorr='perform', 
    #                           shadcorr='omit', 
    #                           photcorr='perform', 
    #                           statflag=True, 
    #                           darkscale='', 
    #                           verbose=False, 
    #                           timestamps=False, 
    #                           trailer='', 
    #                           print_version=False, 
    #                           print_revision=False)

    if verbose:
        print("{}: starting".format(preamble))
    
    if "MAMA" in detector:
        # UV reduction. No cosmic rays issues. Just reduce
        if verbose:
            print("{}: UV reduction".format(preamble))
        basic2d.basic2d(raw_file, output=final_file, verbose=verbose)
    else:
        if verbose:
            print("{}: CCD reduction: first pass".format(preamble))

        initial_mulnoise, mulnoise = initialize_value("mulnoise", issues, input_table, 0., verbose=verbose)
        _, skysub = initialize_value("skysub", issues, input_table, "", verbose=verbose)
        _, wave = initialize_value("wavecal", issues, input_table, "file", verbose=verbose)
        
        if wave != "file":
            # We have an overridden wavecal value. Put it in the file.
            with fits.open(raw_file, mode="update") as raw:
                raw[0].header['WAVECAL'] = wave
        
        interim_file = os.path.join(path, root+"_interim.fits")
        if os.path.isfile(interim_file):
            os.remove(interim_file)
        basic2d.basic2d(raw_file, output=interim_file, verbose=verbose)
        with fits.open(interim_file) as interim_fits:
            initial_dat = interim_fits[1].data
        crj_file = os.path.join(path, root+"_interim_crj.fits")
        crj_params = {
                        "scalense": "{}".format(mulnoise*100), 
                        "initgues": "", 
                        "skysub": skysub, 
                        "crsigmas": "8,6,4",
                        "crradius": 0., 
                        "crthresh": 0.5, 
                        "crmask": ""}
        trial_params = {"accepted": False,
                        "custom": False}
        while not trial_params['accepted']:
            if os.path.isfile(crj_file):
                os.remove(crj_file)
            # Iterate cosmic ray rejection, potentially editing parameters each time,
            # until the result is accepted
            ocrreject.ocrreject(interim_file, crj_file, verbose=verbose, **crj_params)
            with fits.open(crj_file) as cosmic_ray_file:
                crj_dat = cosmic_ray_file['SCI'].data
                cr_flag_dat = cosmic_ray_file['DQ'].data & cr_flag_value
            # Display 
            #   - initial image
            #   - image after rejection
            #   - difference image
            # and collect user input on whether the CR rejection has hit the actual
            # spectrum.
            fig = plt.figure()
            ax = fig.add_subplot(1, 3, 1)
            ax.matshow(initial_dat, origin='lower', norm=LogNorm())
            ax.set_title("Before Cosmic Ray Rejection")
            ax = fig.add_subplot(1, 3, 2)
            ax.matshow(crj_dat, origin='lower', norm=LogNorm())
            ax.set_title("After Cosmic Ray Rejection")
            ax = fig.add_subplot(1, 3, 3)
            ax.imshow(cr_flag_dat, origin='lower')
            ax.set_title("Pixels Flagged as cosmic rays")
            fig.subplots_adjust(bottom=0.2)
            
            def accept(event):
                trial_params['accepted'] = True
                plt.close()
            accept_axis = fig.add_axes([0.1, 0.05, 0.15, 0.075])
            accept_button = Button(accept_axis, 'Accept')
            accept_button.on_clicked(accept)
            
            def tighten(event):
                mulnoise = float(crj_params['scalense'])/100.
                mulnoise -= 0.01
                crj_params['scalense'] = "{}".format(mulnoise*100)
                plt.close()
            more_axis = fig.add_axes([0.35, 0.05, 0.15, 0.075])
            more_button = Button(more_axis, 'Reject more CRs')
            more_button.on_clicked(tighten)
            
            def loosen(event):
                mulnoise = float(crj_params['scalense'])/100.
                mulnoise += 0.01
                crj_params['scalense'] = "{}".format(mulnoise*100)
                plt.close()
            less_axis = fig.add_axes([0.55, 0.05, 0.15, 0.075])
            less_button = Button(less_axis, 'Reject Fewer CRs')
            less_button.on_clicked(loosen)
            
            def custom(event):
                trial_params['accepted'] = True
                trial_params['custom'] = True
                plt.close()
            custom_axis = fig.add_axes([0.75, 0.05, 0.15, 0.075])
            custom_button = Button(custom_axis, 'Enter Custom Value')
            custom_button.on_clicked(custom)

            plt.show()
        # end iteration on cosmic ray rejection
        
        if trial_params['custom']:
            got_custom = False
            while not got_custom:
                try:
                    mulnoise = float(input("Enter custom mulnoise value: "))
                except ValueError as e:
                    print("Error: {}".format(e))
            if verbose:
                msg = "{}: Repeating cosmic ray rejection with mulnoise = {}"
                print(msg.format(preamble, mulnoise))
            crj_params['scalense'] = "{}".format(mulnoise*100)
            ocrreject.ocrreject(interim_file, final_file, verbose=verbose, **crj_params)            
        else:
            shutil.copy(crj_file, final_file)
        
        if verbose:
            print("{}: Finished Cosmic Ray Iteration".format(preamble))
        
        mulnoise = float(crj_params['scalense'])/100.
        handle_parameter('mulnoise', initial_mulnoise, mulnoise, param_file, input_table)

    # Finished 2d extraction preparation
    if verbose:
        print("{}: finished".format(preamble))
    
    return os.path.basename(final_file)


def reduce_extract(input_table, **kwargs):
    """
    Perform the calstis "x1d" spectral extraction, including all custom parameters 
    derived from data (or chosen by the user).

    Parameters
    ----------
    input_table : abscal.stis.stis_data_table.STISDataTable
        Single-row input

    Returns
    -------
    image : np.ndarray
        Flatfielded science image
    """
    task = "stis_reduce_extract"
    root = input_table['root']
    mode = input_table['mode']
    target = input_table['target']
    detector = input_table['detector']
    path = input_table['path']
    fname = input_table['filename']
    flt_file = os.path.join(path, input_table['flatfielded'])
    final_file = os.path.join(path, root+"_{}_{}_x1d.fits".format(target, mode))    
    verbose = kwargs.get('verbose', False)
    preamble = "{}: {}: {} ({})".format(task, root, detector, mode)
    
    issues = {}
    exposure_parameter_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if exposure_parameter_file is not None:
        with open(exposure_parameter_file, 'r') as inf:
            issues = yaml.safe_load(inf)

    # Vignetting info
    dist1 = "{}".format(100 + abs(float(input_table['raw_postarg']))/0.05)

    if verbose:
        print("{}: starting".format(preamble))
    
    x1d_params = {"output": final_file, 
                  "backcorr": 'perform', 
                  "ctecorr": 'perform', 
                  "dispcorr": 'perform', 
                  "helcorr": 'perform', 
                  "fluxcorr": 'omit', 
                  "sporder": None, 
                  "a2center": None, 
                  "maxsrch": None, 
                  "globalx": False, 
                  "extrsize": None, 
                  "bk1size": None, 
                  "bk2size": None, 
                  "bk1offst": None, 
                  "bk2offst": None, 
                  "bktilt": None, 
                  "backord": None, 
                  "bksmode": 'median', 
                  "bksorder": 3, 
                  "blazeshift": None, 
                  "algorithm": 'unweighted', 
                  "xoffset": None, 
                  "verbose": verbose, 
                  "timestamps": False, 
                  "trailer": '', 
                  "print_version": False, 
                  "print_revision": False}
    trial_params = {"accepted": False}
    
    for param in x1d_params:
        # Find any parameter values
        found, value = find_value(param, 
                                  issues, 
                                  input_table, 
                                  default=x1d_params[param], 
                                  verbose=verbose)
        if found:
            x1d_params[param] = value
    
    initial_params = deepcopy(x1d_params)
    
    x1d.x1d(flt_file, **x1d_params)

    return os.path.basename(final_file)


def reduce(input_table, **kwargs):
    """
    Reduces STIS spectroscopic data.
    
    Takes a table of data and extracts spectra for individual exposures, using the STIS 
    pipeline while overriding the default values in order to match Ralph's IDL pipeline
    defaults and exposure-specific values.
    
    Parameters
    ----------
    input_table : abscal.common.exposure_data_table.AbscalDataTable
        Table of exposures to be extracted.
    kwargs : dict
        Dictionary of overrides to the default reduction parameters, and command-line 
        options.
    
    Returns
    -------
    input_table : abscal.stis.stis_data_table.STISDataTable
        Updated table of exposures
    """
    task = "stis: reduce"
    
    default_values = get_defaults('abscal.common.args')
    base_defaults = default_values | get_defaults(kwargs.get('module_name', __name__))
    verbose = kwargs.get('verbose', base_defaults['verbose'])
    show_plots = kwargs.get('plots', base_defaults['plots'])
    force = kwargs.get('force', base_defaults['force'])

    if 'out_file' in kwargs:
        out_file = kwargs['out_file']
        out_dir, out_table = os.path.split(out_file)
        if out_dir == '':
            out_dir = os.getcwd()
    else:
        out_dir = os.getcwd()
    spec_name = kwargs.get('spec_dir', base_defaults['spec_dir'])
    spec_dir = os.path.join(out_dir, spec_name)

    if verbose:
        print("{}: Starting STIS data reduction for spectroscopic data.".format(task))

    issues = {}
    exposure_parameter_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if exposure_parameter_file is not None:
        with open(exposure_parameter_file, 'r') as inf:
            issues = yaml.safe_load(inf)
    
    # Make sure the files are set to use the appropriate CRDS references, and that the
    # references have been appropriately retrieved.
    task_verbosity = -1
    if verbose:
        print("{}: Checking Reference Data".format(task))
        task_verbosity = 10
    files = [os.path.join(p,f) for p,f in zip(input_table['path'], input_table['filename'])]
    assign_bestrefs(files, sync)_references=True, verbosity=task_verbosity

    if verbose:
        print("{}: Starting individual file reductions".format(task))
    for row in input_table:
        root = row['root']
        target = row['target']
        mode = row['mode']
        preamble = "{}: {}".format(task, root)

        # Don't extract if there's already an extracted version of
        #   the file present.
        if row['extracted'] != '':
            ext_file = os.path.join(out_dir, row['extracted'])
            if os.path.isfile(ext_file):
                if force:
                    if verbose:
                        msg = "{}: {}: extracted file exists. Re-extracting."
                        print(msg.format(task, root))
                else:
                    if verbose:
                        msg = "{}: {}: skipping extraction because file exists."
                        print(msg.format(task, root))
                    continue
        else:
            extracted_file_name = "{}_{}_x1d.fits".format(root, target)
            extracted_dest = os.path.join(spec_dir, extracted_file_name)

            # If there is already an extracted file for this input, skip.
            if os.path.isfile(extracted_dest):
                ext_str = os.path.join(spec_dir, extracted_file_name)
                if force:
                    if verbose:
                        msg = "{}: {}: extracted file exists. Re-extracting."
                        print(msg.format(task, root))
                else:
                    row['extracted'] = ext_str
                    if verbose:
                        msg = "{}: {}: skipping extraction because file exists."
                        print(msg.format(task, root))
                    continue

        # Only reduce grism data in the reduce function.
        if row['use']:
            if verbose:
                print("{}: Starting {}".format(task, root))

            if verbose:
                print("{}: Flatfielding {} ({})".format(task, root, row['mode']))
            reduced_2d = reduce_flatfield(row, **kwargs)
            row['flatfielded'] = reduced_2d
            
            if verbose:
                print("{}: Extracting {} ({})".format(task, root, row['mode']))
            reduced_extracted = reduce_extract(row, **kwargs)
            row['extracted'] = reduced_extracted
            
            if verbose:
                print("{}: Finished {}".format(task, root))

        elif verbose:
            msg = "{}: Skipping {} because it's been set to don't use (reason: {})."
            print(msg.format(task, root, row['notes']))

    return input_table


def additional_args(**kwargs):
    """
    Additional command-line arguments. 
    
    Provides additional command-line arguments that are unique to the extraction process.
    
    Returns
    -------
    additional_args : dict
        Dictionary of tuples in the form (fixed,keyword) that can be passed to an argument 
        parser to create a new command-line option
    """
    module_name = kwargs.get('module_name', __name__)
    base_defaults = get_defaults(module_name)

    additional_args = {}

    table_help = "The input metadata table to use."
    table_args = ['table']
    table_kwargs = {'help': table_help}
    additional_args['table'] = (table_args, table_kwargs)

    plots_help = "Include result plots while running (default False)."
    plots_args = ["-p", "--plots"]
    plots_kwargs = {'dest': 'plots', 'action': 'store_true', 
                    'default': base_defaults['plots'], 'help': plots_help}
    additional_args['plots'] = (plots_args, plots_kwargs)

    return additional_args


def parse_args(**kwargs):
    """
    Parse command-line arguments.
    
    Gets the custom arguments for extractions, and passes them to the common command-line 
    option function.
    
    Returns
    -------
    res : namespace
        parsed argument namespace
    """
    description_str = 'Process files from metadata table.'
    default_out_file = kwargs.get('default_output_file', 'stis_extracted.log')
    default_in_file = kwargs.get('default_input_file', 'dirtemp.log')

    args = additional_args(**kwargs)

    res = parse(description_str, default_out_file, args, **kwargs)

    if res.paths is not None:
        if "," in res.paths:
            res.paths = res.paths.split(",")
        else:
            res.paths = [res.paths]
    else:
        res.paths = []

    if res.table is None:
        res.table = "dirtemp.log"

    if len(res.paths) == 0:
        res.paths.append(os.getcwd())

    return res


def main(**kwargs):
    """
    Run the extract function.
    
    Runs the extract function if called from the command line, with command-line arguments 
    added in.
    
    Parameters
    ----------
    kwargs : dict
        Dictionary of parameters to override when running.
    """
    kwargs['default_output_file'] = 'stis_extracted.log'
    parsed = parse_args(**kwargs)

    for key in kwargs:
        if hasattr(parsed, key):
            setattr(parsed, key, kwargs[key])

    input_table = STISDataTable(table=parsed.table,
                                duplicates='both',
                                search_str='',
                                search_dirs=parsed.paths)

    output_table = reduce(input_table, **vars(parsed), **kwargs)

    table_fname = parsed.out_file
    output_table.write_to_file(table_fname, parsed.compat)


if __name__ == "__main__":
    main(module_name='abscal.stis.reduce_extract')
