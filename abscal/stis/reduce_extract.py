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

from astropy import wcs
from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.time import Time
from collections import defaultdict
from copy import deepcopy
from crds import assign_bestrefs
from distutils.util import strtobool
import datetime
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import PySimpleGUI as sg
import shutil
from stistools import basic2d
from stistools import ocrreject
from stistools import x1d as stis_extract
import traceback
import yaml

from abscal.common.args import parse
from abscal.common.plots import draw_figure
from abscal.common.plots import make_img_fig
from abscal.common.plots import make_figure_window
from abscal.common.plots import make_params_window
from abscal.common.plots import make_spectrum_window
from abscal.common.utils import get_value
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
    exp_info = "{} {} {}".format(root, target, mode)
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
        crj_params = {"scalense": "{}".format(mulnoise*100), 
                      "initgues": "",
                      "skysub": skysub,
                      "crsigmas": "8,6,4",
                      "crradius": 0.,
                      "crthresh": 0.5,
                      "crmask": ""}
        crj_types = defaultdict(lambda: str)
        crj_types["crradius"] = float
        crj_types["crthresh"] = float
        initial_crj_params = deepcopy(crj_params)

        pre_window = make_figure_window("{}: Before Cosmic Ray Removal".format(exp_info))
        pre_fig = make_img_fig(np.log10(np.where(initial_dat>=0.1, initial_dat, 0.1)))
        draw_figure(pre_window, pre_fig)

        post_window = make_figure_window("{}: After Cosmic Ray Removal".format(exp_info))
        flag_window = make_figure_window("{}: Cosmic Ray Flags".format(exp_info))

        # Create Parameter Window
        cr_buttons = [sg.Button('Flag More CRs'), sg.Push(), sg.Button('Flag Fewer CRs')]
        param_window = make_params_window(exp_info, "ocrreject", crj_params, [cr_buttons])
                
        changed = True
        while True:
            if changed:
                # If a value was changed, re-do the cosmic ray rejection.
                if os.path.isfile(crj_file):
                    os.remove(crj_file)
                ocrreject.ocrreject(interim_file, crj_file, verbose=verbose, **crj_params)
                with fits.open(crj_file) as cosmic_ray_file:
                    crj_dat = cosmic_ray_file['SCI'].data
                    cr_flag_dat = cosmic_ray_file['DQ'].data & cr_flag_value
                
                # Create the new figures
                post_fig = make_img_fig(np.log10(np.where(crj_dat>=0.1, crj_dat, 0.1)))
                draw_figure(post_window, post_fig)
                flag_fig = make_img_fig(cr_flag_dat)
                draw_figure(flag_window, flag_fig)
                changed = False
            window, event, values = sg.read_all_windows()
            if window is None and event != sg.TIMEOUT_EVENT:
                print('exiting because no windows are left')
                break
            elif event == sg.WIN_CLOSED or event == 'Exit':
                # Closing the parameter window is the equivalent of "Accept"
                if window == param_window:
                    break
                window.close()
            elif event == 'Accept':
                break
            elif event == "Reset":
                changed = True
                crj_params = deepcopy(initial_crj_params)
            elif event == "Re-plot":
                for item in crj_params:
                    if crj_params[item] != crj_types[item](values[item]):
                        changed = True
                    crj_params[item] = crj_types[item](values[item])
            elif event == "Flag More CRs":
                crj_params["scalense"] -= 1.
                changed = True
            elif event == "Flag Fewer CRs":
                crj_params["scalense"] += 1.
                changed = True
            else:
                print("{}: Got event {} for window {} with values {}".format(preamble, event, window, values))
        pre_window.close()
        post_window.close()
        flag_window.close()
        param_window.close()

        if changed:
            # At least one CRJ parameter was changed without re-running ocrreject
            ocrreject.ocrreject(interim_file, final_file, verbose=verbose, **crj_params)            
        else:
            shutil.copy(crj_file, final_file)
        
        if verbose:
            print("{}: Finished Cosmic Ray Iteration".format(preamble))
        
        if os.path.isfile(interim_file):
            os.remove(interim_file)
        if os.path.isfile(crj_file):
            os.remove(crj_file)
        
        for item in crj_params:
            if item == "scalense":
                initial_mulnoise = float(initial_crj_params["scalense"])/100.
                mulnoise = float(crj_params["scalense"])/100.
                handle_parameter("mulnoise", initial_mulnoise, mulnoise, param_file, 
                                 input_table)
            else:
                handle_parameter(item, initial_crj_params[item], crj_params[item], 
                                 param_file, input_table)

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
    exp_info = "{} {} {}".format(root, target, mode)
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
                  "fluxcorr": 'perform', 
                  "sporder": None, 
                  "a2center": None, 
                  "maxsrch": None, 
                  "globalx": False,  # Whether to use global cross-correlation offset for all orders
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
                  "verbose": verbose}
    x1d_types = {'sporder': int,
                 'a2center': float,
                 'maxsrch': float,
                 'extrsize': float,
                 'bk1size': float,
                 'bk2size': float,
                 'bk1offset': float,
                 'bk2offset': float,
                 'bktilt': float,
                 'backord': int,
                 'bksorder': int,
                 'blazeshift': float,
                 'xoffset': float}
    
    for param in x1d_params:
        # Find any parameter values
        found, value = get_value(param, issues, input_table, default=x1d_params[param], 
                                 verbose=verbose)
        if found:
            x1d_params[param] = value
    
    initial_params = deepcopy(x1d_params)
    
    stis_extract.x1d(flt_file, **x1d_params)
    
    # For plotting what's going on in the 2D file:
    #   - "backcorr":   f[0].header["BACKCORR"]
    #   - "ctecorr":    f[0].header["CTECORR"]
    #   - "dispcorr":   f[0].header["DISPCORR"]
    #   - "helcorr":    f[0].header["HELCORR"]
    #   - "fluxcorr":   f[0].header["FLUXCORR"]
    #   - "sporder":    f[1].data["SPORDER"]
    #   - "a2center":   f[1].data["A2CENTER"] 
    #   - "maxsrch":    f[1].data["MAXSRCH"]
    #   - "globalx":    [*****probably not needed?*****]
    #   - "extrsize":   f[1].data["EXTRSIZE"]
    #   - "bk1size":    f[1].data["BK1SIZE"]
    #   - "bk2size":    f[1].data["BK2SIZE"]
    #   - "bk1offst":   f[1].data["BK1OFFST"]
    #   - "bk2offst":   f[1].data["BK2OFFST"]
    #   - "bktilt": 
    #                   - look at f[0].header["XTRACTAB"]
    #                   - get aperture from f[0].header["PROPAPER"]
    #                   - get opt_elem from f[0].header["OPT_ELEM"]
    #                   - get cenwave from int(f[0].header["CENWAVE"])
    #                   - get the BKTCOEFF values matching these
    #                   - somehow figure out what each of these 7 coefficients means?
    #   - "backord":    as "bktilt" but for "BACKORD" column. It's either 0 or 1. It 
    #                   defines whether to use the average of lower/upper background 
    #                   region (backord=0, default), or do a linear interpolation
    #                   (backord=1)
    #   - "bksmode":    how to smooth background (mode or median). Default is median.
    #   - "bksorder":   polynomial to fit to smoothed background. Integer. Defaults to 3.
    #   - "blazeshift": echelle blazeshift, based on SHIFTA1, SHIFTA2, and MJD of exposure
    #   - "algorithm":  as "bktilt" but XTRACALG
    #   - "xoffset":    f[0].header["SHIFTA1"]
    
    if os.path.isfile(final_file):
        with fits.open(flt_file) as fltf, fits.open(final_file) as x1df:
            x_arr = np.arange(1024, dtype=np.int32)
            spec_loc = x1df[1].data['EXTRLOCY'][0]
            spec_low = spec_loc - x1df[1].data["EXTRSIZE"][0]//2
            spec_high = spec_loc + x1df[1].data["EXTRSIZE"][0]//2
        
            b1_loc = spec_loc + x1df[1].data["BK1OFFST"][0]
            b1_low = b1_loc - x1df[1].data["BK1SIZE"][0]//2
            b1_high = b1_loc + x1df[1].data["BK1SIZE"][0]//2
        
            b2_loc = spec_loc + x1df[1].data["BK2OFFST"][0]
            b2_low = b2_loc - x1df[1].data["BK2SIZE"][0]//2
            b2_high = b2_loc + x1df[1].data["BK2SIZE"][0]//2
            
            ext_data = fltf[1].data

        extr_window = make_figure_window("{}: Extraction Region".format(exp_info))
        ext_fig = make_img_fig(np.log10(np.where(ext_data>=0.1,ext_data,0.1)))
        ax = ext_fig.axes[0]
        ax.plot(x_arr, spec_low, color='white', label='Extraction Region')
        ax.plot(x_arr, spec_high, color='white')
        ax.plot(x_arr, b1_low, color='red', label='Background 1 Region')
        ax.plot(x_arr, b1_high, color='red')
        ax.plot(x_arr, b2_low, color='blue', label='Background 2 Region')
        ax.plot(x_arr, b2_high, color='blue')
        ax.legend()
        draw_figure(extr_window, ext_fig)
        
        spec_window, spec_fig, lines = make_spectrum_window("{}: Extracted".format(exp_info),
                                                            final_file)
        draw_figure(spec_window, spec_fig)

        # Create Parameter Window
        param_window = make_params_window(exp_info, "x1d", x1d_params)
        
        changed = False
        while True:
            if changed:
                # If a value was changed, re-do the cosmic ray rejection.
                if os.path.isfile(final_file):
                    os.remove(final_file)
                stis_extract.x1d(flt_file, **x1d_params)
                with fits.open(final_file) as x1df:
                    x_arr = np.arange(1024, dtype=np.int32)
                    spec_loc = x1df[1].data['EXTRLOCY'][0]
                    spec_low = spec_loc - x1df[1].data["EXTRSIZE"][0]//2
                    spec_high = spec_loc + x1df[1].data["EXTRSIZE"][0]//2
        
                    b1_loc = spec_loc + x1df[1].data["BK1OFFST"][0]
                    b1_low = b1_loc - x1df[1].data["BK1SIZE"][0]//2
                    b1_high = b1_loc + x1df[1].data["BK1SIZE"][0]//2
        
                    b2_loc = spec_loc + x1df[1].data["BK2OFFST"][0]
                    b2_low = b2_loc - x1df[1].data["BK2SIZE"][0]//2
                    b2_high = b2_loc + x1df[1].data["BK2SIZE"][0]//2
                
                # Create the new figures
                ext_fig = make_img_fig(np.log10(np.where(ext_data>=0.1,ext_data,0.1)))
                ax = ext_fig.axes[0]
                ax.plot(x_arr, spec_low, color='white', label='Extraction Region')
                ax.plot(x_arr, spec_high, color='white')
                ax.plot(x_arr, b1_low, color='red', label='Background 1 Region')
                ax.plot(x_arr, b1_high, color='red')
                ax.plot(x_arr, b2_low, color='blue', label='Background 2 Region')
                ax.plot(x_arr, b2_high, color='blue')
                ax.legend()
                draw_figure(extr_window, ext_fig)
                
                for line in lines:
                    lines[line].set_visible(spec_window[line].get())
                # Toggle line visibility on the extracted figure
                draw_figure(spec_window, spec_fig)
                changed = False
            window, event, values = sg.read_all_windows()
            if window is None and event != sg.TIMEOUT_EVENT:
                print('exiting because no windows are left')
                break
            elif event == sg.WIN_CLOSED or event == 'Exit':
                # Closing the parameter window is the equivalent of "Accept"
                if window == param_window:
                    break
                window.close()
            elif event == 'Accept':
                break
            elif event == "Reset":
                changed = True
                x1d_params = deepcopy(initial_x1d_params)
            elif event == "Re-plot":
                for item in x1d_params:
                    if x1d_params[item] != x1d_types[item](values[item]):
                        if not (x1d_params[item] is None and values[item] == 'None'):
                            changed = True
                    if values[item] == 'None':
                        x1d_params[item] = None
                    else:
                        x1d_params[item] = x1d_types[item](values[item])
            elif event in lines.keys():
                lines[event].set_visible(not lines[event].get_visible())
                spec_fig.canvas.draw()
                changed = True
            else:
                print("{}: Got event {} for window {} with values {}".format(preamble, event, window, values))
        extr_window.close()
        spec_window.close()
        param_window.close()
    else:
        print("{}: ERROR: Spectral extraction failed.")

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
    update_refs = kwargs.get('ref_update', base_defaults['ref_update'])

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
    
    if update_refs:
        # Make sure the files are set to use the appropriate CRDS references, and that the
        # references have been appropriately retrieved.
        task_verbosity = -1
        if verbose:
            print("{}: Checking Reference Data".format(task))
            task_verbosity = 10
        files = [os.path.join(p,f) for p,f in zip(input_table['path'], input_table['filename'])]
        assign_bestrefs(files, sync_references=True, verbosity=task_verbosity)
    else:
        if verbose:
            print("{}: Skipping reference file update".format(task))

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

    plots_default = base_defaults['plots']
    if isinstance(plots_default, str):
        plots_default = strtobool(plots_default)
    plots_help = "Toggle interactive mode (default {}).".format(plots_default)
    plots_args = ["-p", "--plots"]
    plots_kwargs = {'dest': 'plots', 'default': plots_default, 'help': plots_help}
    if plots_default:
        plots_kwargs['action'] = 'store_false'
    else:
        plots_kwargs['action'] = 'store_true'
    additional_args['plots'] = (plots_args, plots_kwargs)
    
    ref_default = base_defaults['ref_update']
    if isinstance(ref_default, str):
        ref_default = strtobool(ref_default)
    ref_help = "Toggle whether to update reference files while running "
    ref_help += "(default {})".format(ref_default)
    ref_args = ["-r", "--ref_update"]
    ref_kwargs = {'dest': 'ref_update', 'default': ref_default, 'help': ref_help}
    if ref_default:
        ref_kwargs['action'] = 'store_false'
    else:
        ref_kwargs['action'] = 'store_true'
    additional_args['ref_update'] = (ref_args, ref_kwargs)

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
