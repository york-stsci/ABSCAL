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
from datetime import datetime
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import PySimpleGUI as sg
from scipy.signal import find_peaks_cwt
import shutil
from stistools import basic2d
from stistools import ocrreject
from stistools import x1d as stis_extract
import traceback
import yaml

from abscal.common.args import parse
from abscal.common.ui import handle_parameter_window
from abscal.common.ui import ImageWindow
from abscal.common.ui import run_task
from abscal.common.ui import SpectrumWindow
from abscal.common.ui import TaskWindow
from abscal.common.ui import TwoColumnWindow
from abscal.common.utils import check_params
from abscal.common.utils import get_value
from abscal.common.utils import get_data_dir
from abscal.common.utils import get_data_file
from abscal.common.utils import get_defaults
from abscal.common.utils import get_param_types
from abscal.common.utils import handle_parameter
from abscal.common.utils import initialize_value
from abscal.common.utils import set_params
from abscal.common.utils import set_image
from abscal.common.utils import setup_params
from abscal.stis.stis_data_table import STISDataTable


def make_hot_pixel_list(dark_file, threshold):
    """
    Create and return a list of hot pixels, given a dark reference file and a count rate
    threshold.
    
    Dark reference files contain a science (SCI) extension which shows the dark count rate
    at every detector pixel. In order to turn this into a list of hot pixels (pixels whose 
    count rates exceed the provided threshold), this function filters the data array to an
    array of pixels whose values exceed the provided threshold, and uses the pixel index 
    values to produce the pixel co-ordinates.
    
    Parameters
    ----------
    dark_file : str
        Name of the dark reference file, potentially including the STIS reference keyword
        location ("oref") as a path specifier
    threshold : float
        The count rate threshold for considering a pixel as a hot pixel
    
    Returns
    -------
    hot_list : list of tuple
        List of tuples in the form (x, y, rate) including hot pixels
    """
    if "oref$" in dark_file:
        dark_file = os.path.join(os.environ["oref"], dark_file.replace("oref$", ""))
    
    print("Creating hot pixel list")
    with fits.open(dark_file) as dark:
        sci_ext = dark['SCI'].data
        hot_pixels = np.where(sci_ext>=threshold)
        x_list = hot_pixels[1]
        y_list = hot_pixels[0]
        rate_list = sci_ext[hot_pixels]
        hot_list = [(x, y, r) for x, y, r in zip(x_list, y_list, rate_list)]
        print("Found {} hot pixels".format(len(hot_list)))

    return hot_list


class ExtractionInfoWindow(ImageWindow):
    """
    This is a version of the ImageWindow that also plots the extraction region and the 
    background regions.
    """
    def set_up_config(self, kwargs):
        self.extraction_file = self.handle_kwarg('extraction_file', None, kwargs)
        return super().set_up_config(kwargs)
    
    def get_line_data(self):
        line_data = {}
        with fits.open(self.extraction_file) as x1df:
            spec_loc = x1df[1].data['EXTRLOCY'][0]
            line_data['spec_low'] = spec_loc - x1df[1].data["EXTRSIZE"][0]//2
            line_data['spec_high'] = spec_loc + x1df[1].data["EXTRSIZE"][0]//2
    
            b1_loc = spec_loc + x1df[1].data["BK1OFFST"][0]
            line_data['b1_low'] = b1_loc - x1df[1].data["BK1SIZE"][0]//2
            line_data['b1_high'] = b1_loc + x1df[1].data["BK1SIZE"][0]//2
    
            b2_loc = spec_loc + x1df[1].data["BK2OFFST"][0]
            line_data['b2_low'] = b2_loc - x1df[1].data["BK2SIZE"][0]//2
            line_data['b2_high'] = b2_loc + x1df[1].data["BK2SIZE"][0]//2
        return line_data

    def create_figure(self):
        figure = super().create_figure()
        ax = figure.axes[0]
        line_data = self.get_line_data()
        extr_lines = {}
        x_arr = np.arange(1024, dtype=np.int32)
        extr_lines["spec_low"] = ax.plot(x_arr, line_data['spec_low'], color='green', 
            label='Extraction Region')[0]
        extr_lines["spec_high"] = ax.plot(x_arr, line_data['spec_high'], color='green')[0]
        extr_lines["b1_low"] = ax.plot(x_arr, line_data['b1_low'], color='red', 
            label='Background 1 Region')[0]
        extr_lines["b1_high"] = ax.plot(x_arr, line_data['b1_high'], color='red')[0]
        extr_lines["b2_low"] = ax.plot(x_arr, line_data['b2_low'], color='blue', 
            label='Background 2 Region')[0]
        extr_lines["b2_high"] = ax.plot(x_arr, line_data['b2_high'], color='blue')[0]
        ax.legend()
        self.lines = extr_lines
        return figure

    def update_figure(self, **kwargs):
        """
        Updates the figure image, using the file and log settings from the initial 
        creation. The keyword arguments aren't used here, but are provided for subclasses
        where it's possible that only a partial update might be needed (e.g. only updating
        plot lines, only updating data image, only updating visibility, etc.)
        """
        line_data = self.get_line_data()
        for key in line_data:
            self.lines[key].set_ydata(line_data[key])
        super().update_figure(**kwargs)


class CosmicRayWindow(TwoColumnWindow):
    """
    This is a version of the ImageWindow that includes a function to get the approximate
    spectrum location (based on doing a collapse of the image along the spectral direction),
    and offering the option of showing the (approximate) extraction region on the image of
    cosmic ray flags, to offer the option of a guide for seeing whether the cosmic ray 
    flagging has accidentally flagged some part of the spectral trace.
    """
    def set_up_config(self, kwargs):
        self.extrsize = self.handle_kwarg('extrsize', 11., kwargs)
        self.trace_file = self.handle_kwarg('trace_file', None, kwargs)
        return super().set_up_config(kwargs)
    
    def make_ui_column(self):
        trace_check = sg.Checkbox("Show Approximate Spectral Trace", default=True,
                                  text_color="red", key="trace", enable_events=True,
                                  background_color="white")
        trace_frame = sg.Frame("Trace:", [[trace_check]], background_color="white",
                               title_color="black", expand_y=True)
        trace_col = sg.Column([[trace_frame]], expand_y=True)
        return trace_col
    
    def get_trace_data(self):
        peaks = self.get_approximate_peak()
        y_vals = []
        if len(peaks) == 0:
            print("{}: WARNING: no spectral trace found".format(self.data_file))
            return y_vals
        if len(peaks) > 1:
            print("{}: WARNING: multiple possible traces found".format(self.data_file))
        for peak in peaks:
            y_vals.append(peak-self.extrsize/2)
            y_vals.append(peak+self.extrsize/2)
        return y_vals
    
    def get_image_data(self):
        cr_flag_value = 8192
        with fits.open(self.data_file) as dat:
            cr_flag_dat = np.zeros_like(dat["DQ"].data, dtype=np.int32)
            for ext in dat:
                if dat[ext].name == 'DQ':
                    cr_flag_dat += np.where(dat[ext].data & cr_flag_value != 0, 1, 0)
        return cr_flag_dat
    
    def create_figure(self):
        figure = super().create_figure()
        ax = figure.axes[0]
        self.lines = []
        for value in self.get_trace_data():
            self.lines.append(ax.axhline(value, color="red", linestyle='dotted'))
        return figure

    def update_figure(self, **kwargs):
        update_data = kwargs.get("update_data", True)
        for line in self.lines:
            line.set_visible(self['trace'].Get())
        if update_data:
            y_vals = self.get_trace_data()
            for line, val in zip(self.lines, y_vals):
                line.set_ydata(val)
            super().update_figure()
        else:
            self.figure.canvas.draw()
    
    def get_approximate_peak(self):
        """
        Gets an approximate spectrum location by 
        - summing the spectrum along the X axis
        - subtracting the medium
        - setting any value < 0.1 of the maximum value to 0.
        - running a peak-finding routine
        """
        with fits.open(self.trace_file) as data_file:
            dat = np.median(data_file[1].data, axis=1)
            dat -= np.median(dat)
            dat = np.where(dat>=0.1*np.max(dat), dat, 0.)
            peaks = find_peaks_cwt(dat, [10])
        return peaks


class CosmicRayTaskWindow(TaskWindow):
    """
    Modified version of the TaskWindow that adds in buttons for tweaking up and down
    the scalense parameter.
    """
    def add_to_layout(self, layout):
        layout.append([sg.Button('Flag More CRs'), sg.Push(), sg.Button('Flag Fewer CRs')])
        return layout
    
    def handle_ui_event(self, event, values):
        if event == "Flag More CRs":
            self.params["scalense"] = "{}".format(float(self.params["scalense"]) - 1.)
            self["scalense"].Update(self.params["scalense"])
            return True
        elif event == "Flag Fewer CRs":
            self.params["scalense"] = "{}".format(float(self.params["scalense"]) + 1.)
            self["scalense"].Update(self.params["scalense"])
            return True
        return super().handle_ui_event(event, values)


def reduce_flatfield(input_row, **kwargs):
    """
    Perform the calstis "basic2d" data reduction, including iterating on cosmic ray
    rejection. This is done by repeatedly calling the calstis "basic2d" on the raw STIS
    data with various calibration flag settings, and by calling "ocrreject" separately in 
    order to have finer control over cosmic ray rejection parameters.

    Parameters
    ----------
    input_row : abscal.stis.stis_data_table.STISDataTable
        Single-row input

    Returns
    -------
    final_file : str
        Name (including path) of the CR-rejected 2D image.
    """
    task = "stis_reduce_flatfield"
    cr_flag_value = 8192
    hot_pixel_flag_value = 16
    root = input_row['root']
    mode = input_row['mode']
    target = input_row['target']
    detector = input_row['detector']
    path = input_row['path']
    fname = input_row['filename']
    raw_file = os.path.join(path, fname)
    final_file = raw_file.replace("_raw", "_{}_{}_2d".format(target, mode))
    preamble = "{}: {}: {} ({})".format(task, root, detector, mode)
    exp_info = "{} {} {}".format(root, target, mode)
    verbose = kwargs.get("verbose", False)

    settings = {}
    setting_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if setting_file is not None and os.path.isfile(setting_file):
        with open(setting_file, 'r') as inf:
            settings = yaml.safe_load(inf)
    
    # Do we want to make basic2d interactive? My thought is currently "no"

    if verbose:
        print("{}: starting".format(preamble))
    
    if "MAMA" in detector:
        # UV reduction. No cosmic rays issues. Just reduce
        if verbose:
            print("{}: MAMA reduction".format(preamble))
        basic2d.basic2d(raw_file, output=final_file, verbose=verbose)
    else:
        if verbose:
            print("{}: CCD reduction: first pass".format(preamble))
        
        with fits.open(raw_file, mode="update") as exposure:
            exposure[0].header['HISTORY'] = "ABSCAL: Starting 2d reduction"
        
        # Do the 2D reduction
        interim_file = os.path.join(path, root+"_interim.fits")
        if os.path.isfile(interim_file):
            os.remove(interim_file)
        basic2d.basic2d(raw_file, output=interim_file, verbose=verbose)
        
        with fits.open(interim_file, mode="update") as exposure:
            exposure[0].header['HISTORY'] = "ABSCAL: Finished running basic2d"
        
        # Set up and run the cosmic ray rejection task
        figure_windows = []
        standard_kwargs = {"do_log": True, "draw_toolbar": True}
        figure_windows.append({"title": "{}: Before Cosmic Ray Removal".format(exp_info),
                               "static": True,
                               "data_file": interim_file,
                               "kwargs": standard_kwargs})
        figure_windows.append({"title": "{}: After Cosmic Ray Removal".format(exp_info),
                               "data_file": "working",
                               "kwargs": standard_kwargs})
        x1d_params = setup_params("x1d", "stis", settings, input_row, verbose)
        cr_kwargs = {"draw_toolbar": True, "extrsize": x1d_params['extrsize'], "trace_file": "working"}
        figure_windows.append({"title": "{}: Cosmic Ray Flags".format(exp_info),
                               "window_class": CosmicRayWindow,
                               "data_file": interim_file,
                               "kwargs": cr_kwargs})

        try:
            run_task("stis", "ocrreject", setting_file, input_row, ocrreject.ocrreject, 
                     interim_file, final_file, verbose=verbose, figure_windows=figure_windows, 
                     task_window_class=CosmicRayTaskWindow, output_is_param=False)
        except Exception as e:
            print("{}: ERROR: {}".format(preamble, e))
            return None
        
        with fits.open(final_file, mode="update") as exposure:
            exposure[0].header["HISTORY"] = "ABSCAL: Finished running OCRREJECT"
        
        # Interpolate across hot pixels
        if verbose:
            print("{}: Averaging over hot pixels".format(preamble))
        hpix_params = setup_params("hpix", "stis", settings, input_row, verbose)
        hot_pixel_mapping = kwargs['hot_pixel_mapping']
        with fits.open(interim_file) as exposure:
            dark_file = exposure[0].header["DARKFILE"]
        if dark_file not in hot_pixel_mapping:
            if verbose:
                print("{}: Creating mapping for {}".format(preamble, dark_file))
            hot_pixel_mapping[dark_file] = make_hot_pixel_list(dark_file, hpix_params['threshold'])
        if verbose:
            print("{}: Hot pixel interpolation".format(preamble))
        # 
        # The DQ "magic value" of 125 is from the original hot pixel interpolation routine
        # in IDL calstis.
        # 
        with fits.open(final_file, mode='update') as exposure:
            for (x, y, rate) in hot_pixel_mapping[dark_file]:
                if x != 0 and x != exposure['SCI'].data.shape[1]-1:
                    if (exposure['DQ'].data[y, x-1] < 125) and (exposure['DQ'].data[y, x+1] < 125):
                        exposure['SCI'].data[y, x] = (exposure['SCI'].data[y, x-1] + exposure['SCI'].data[y, x+1])/2.
                    elif exposure['DQ'].data[y, x-1] < 125:
                        exposure['SCI'].data[y, x] = exposure['SCI'].data[y, x-1]
                    elif exposure['DQ'].data[y, x+1] < 125:
                        exposure['SCI'].data[y, x] = exposure['SCI'].data[y, x+1]
                    else:
                        exposure['SCI'].data[y, x] = 0.
                        exposure['DQ'].data[y, x] = 252
                exposure['DQ'].data[y, x] |= hot_pixel_flag_value
            exposure[0].header["HISTORY"] = "ABSCAL: Finished hot pixel interpolation"
        
    # Finished 2d extraction preparation

    if verbose:
        print("{}: finished".format(preamble))
    
    return os.path.basename(final_file)


def reduce_extract(input_row, **kwargs):
    """
    Perform the calstis "x1d" spectral extraction, including all custom parameters 
    derived from data (or chosen by the user).

    Parameters
    ----------
    input_row : abscal.stis.stis_data_table.STISDataTable
        Single-row input

    Returns
    -------
    image : np.ndarray
        Flatfielded science image
    """
    task = "stis_reduce_extract"
    root = input_row['root']
    mode = input_row['mode']
    target = input_row['target']
    detector = input_row['detector']
    path = input_row['path']
    fname = input_row['filename']
    flt_file = os.path.join(path, input_row['flatfielded'])
    final_file = os.path.join(path, root+"_{}_{}_x1d.fits".format(target, mode))    
    exp_info = "{} {} {}".format(root, target, mode)
    verbose = kwargs.get('verbose', False)
    preamble = "{}: {}: {} ({})".format(task, root, detector, mode)
    
    settings = {}
    setting_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if setting_file is not None and os.path.isfile(setting_file):
        with open(setting_file, 'r') as inf:
            settings = yaml.safe_load(inf)

    # Vignetting info
    dist1 = "{}".format(100 + abs(float(input_row['raw_postarg']))/0.05)

    if verbose:
        print("{}: starting".format(preamble))

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

    figure_windows = []
    extr_kwargs = {"extraction_file": "working", "do_log": True, "draw_toolbar": True}
    figure_windows.append({"title": "{}: Extraction Region".format(exp_info),
                           "window_class": ExtractionInfoWindow,
                           "data_file": flt_file,
                           "kwargs": extr_kwargs})
    spec_kwargs = {"draw_toolbar": True, "plot_x": "wave"}
    figure_windows.append({"title": "{}: Extracted".format(exp_info),
                           "window_class": SpectrumWindow,
                           "data_file": "working",
                           "kwargs": spec_kwargs})
    
    try:
        run_task("stis", "x1d", setting_file, input_row, stis_extract.x1d, flt_file,
                 final_file, verbose=verbose, figure_windows=figure_windows, 
                 output_is_param=True)
    except Exception as e:
        print("{}: ERROR: {}".format(preamble, e))
        return None

    if verbose:
        print("{}: finished".format(preamble))

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

    settings = {}
    setting_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if setting_file is not None:
        with open(setting_file, 'r') as inf:
            settings = yaml.safe_load(inf)
    
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
    kwargs['hot_pixel_mapping'] = {}

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

        # Only reduce spectroscopic data in the reduce function.
        if row['use']:
            if verbose:
                print("{}: Starting {}".format(task, root))

            if verbose:
                print("{}: Flatfielding {} ({})".format(task, root, row['mode']))
            reduced_2d = reduce_flatfield(row, **kwargs)
            if reduced_2d is not None:
                row['flatfielded'] = reduced_2d
            
                if verbose:
                    print("{}: Extracting {} ({})".format(task, root, row['mode']))
                reduced_extracted = reduce_extract(row, **kwargs)
                if reduced_extracted is not None:
                    row['extracted'] = reduced_extracted
            else:
                if verbose:
                    print("{}: Flatfielding failed. Skipping extraction.".format(task))
            
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
