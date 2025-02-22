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
from functools import partial
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from scipy.signal import find_peaks_cwt
import shutil
import tkinter as tk
from tkinter import ttk
import traceback
import yaml

from stistools import basic2d
from stistools import ocrreject
from stistools import x1d as stis_extract
from stistools.wavecal import wavecal
from stistools.defringe import defringe
from stistools.defringe import mkfringeflat
from stistools.defringe import normspflat

from abscal.common.args import parse
from abscal.common.file_utils import get_data_file
from abscal.common.logging import DEFAULT_LOGGER as logger
from abscal.common.ui import AbscalFrame
from abscal.common.ui import AbscalTask
from abscal.common.ui import ImageFrame
from abscal.common.ui import SpectrumFrame
from abscal.common.ui import SpectrumWindow
from abscal.common.ui import TaskFrame
from abscal.common.ui import TwoColumnWindow
from abscal.common.utils import check_params
from abscal.common.utils import get_value
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
    
    logger.info("Creating hot pixel list")
    with fits.open(dark_file) as dark:
        sci_ext = dark['SCI'].data
        hot_pixels = np.where(sci_ext>=threshold)
        x_list = hot_pixels[1]
        y_list = hot_pixels[0]
        rate_list = sci_ext[hot_pixels]
        hot_list = [(x, y, r) for x, y, r in zip(x_list, y_list, rate_list)]
        logger.info(f"Found {len(hot_list)} hot pixels")

    return hot_list


class ExtractionInfoFrame(ImageFrame):
    """
    This is a version of the ImageFrame that also plots the extraction region and the 
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


class CosmicRayTraceFrame(AbscalFrame):
    def __init__(self, parent, visibility_callback, starting_value, *args, **kwargs):
        self.visibility_callback = visibility_callback
        self.starting_value = starting_value
        super().__init__(parent, "Spectral Trace", *args, **kwargs)

    def layout_frame(self):
        self.trace_var = tk.BooleanVar(master=self, value=self.starting_value)
        b = ttk.Checkbutton(
            self,
            text="Display Approximate Trace",
            variable=self.trace_var,
            onvalue=True,
            offvalue=False,
            command=partial(self.handle_ui_event, "visibility")
        )
        b.grid(row=0, column=0, sticky="news")
        self.grid_rowconfigure(0, weight=0)
        self.grid_columnconfigure(0, weight=0)

    def handle_ui_event(self, event):
        if event == "visibility":
            self.visibility_callback(self.trace_var.get())


class CosmicRayFrame(ImageFrame):
    """
    This is a version of the ImageFrame that includes a function to get the approximate
    spectrum location (based on doing a collapse of the image along the spectral direction),
    and offering the option of showing the (approximate) extraction region on the image of
    cosmic ray flags, to offer the option of a guide for seeing whether the cosmic ray 
    flagging has accidentally flagged some part of the spectral trace.
    """
    def set_up_config(self, kwargs):
        self.extrsize = self.handle_kwarg('extrsize', 11., kwargs)
        self.trace_file = self.handle_kwarg('trace_file', None, kwargs)
        self.extraction_file = self.handle_kwarg('extraction_file', None, kwargs)
        self.show_trace = self.handle_kwarg('trace_visible', True, kwargs)
        return super().set_up_config(kwargs)
    
    def get_trace_data(self):
        y_vals = []
        if (self.extraction_file is not None) and (self.extraction_file.is_file()):
            with fits.open(self.extraction_file) as x1df:
                spec_loc = x1df[1].data['EXTRLOCY'][0]
                y_vals.append(spec_loc - x1df[1].data["EXTRSIZE"][0]//2)
                y_vals.append(spec_loc + x1df[1].data["EXTRSIZE"][0]//2)
        else:
            peaks = self.get_approximate_peak()
            if len(peaks) == 0:
                logger.warning(f"{self.data_file}: WARNING: no spectral trace found")
                return y_vals
            if len(peaks) > 1:
                logger.warning(f"{self.data_file}: WARNING: multiple possible traces found")
            for peak in peaks:
                y_vals.append(np.array([peak-self.extrsize//2 for x in range(self.x_size)]))
                y_vals.append(np.array([peak+self.extrsize//2 for x in range(self.x_size)]))
        return y_vals
    
    def get_figure_data(self):
        cr_flag_value = 8192
        with fits.open(self.data_file) as dat:
            cr_flag_dat = np.zeros_like(dat["DQ"].data, dtype=np.int32)
            self.x_size = cr_flag_dat.shape[1]
            for ext in dat:
                if dat[ext].name == 'DQ':
                    cr_flag_dat += np.where(dat[ext].data & cr_flag_value != 0, 1, 0)
        return cr_flag_dat
    
    def create_figure(self):
        figure = super().create_figure()
        xax = np.arange(self.x_size)
        ax = figure.axes[0]
        self.lines = []
        for value in self.get_trace_data():
            line, = ax.plot(xax, value, color="red", linestyle='dotted')
            self.lines.append(line)
            self.lines[-1].set_visible(self.show_trace)
        return figure

    def update_figure(self, **kwargs):
        update_data = kwargs.get("update_data", True)
        for line in self.lines:
            line.set_visible(self.show_trace)
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

    def toggle_trace(self, show_trace):
        self.show_trace = show_trace
        self.update_figure(update_data=False)


class CosmicRayWindow(TwoColumnWindow):
    def setup_ui(self, frame_class, frame_args, frame_kwargs):
        figure_frame = frame_class(self, *frame_args, **frame_kwargs)
        trace_frame = CosmicRayTraceFrame(
            self,
            figure_frame.toggle_trace,
            figure_frame.show_trace
        )
        super().setup_ui(trace_frame, figure_frame)
        self.grid_columnconfigure(0, weight=0)


class CosmicRayTaskFrame(TaskFrame):
    """
    Modified version of the TaskFrame that adds in buttons for tweaking up and down
    the scalense parameter.
    """
    def add_to_layout(self, current_row):
        b = ttk.Button(
            self,
            text="Flag Fewer CRs",
            command=partial(self.handle_ui_event, "flag_fewer")
        )
        b.grid(row=current_row, column=0, sticky="news")
        b = ttk.Button(
            self,
            text="Flag More CRs",
            command=partial(self.handle_ui_event, "flag_more")
        )
        b.grid(row=current_row, column=2, sticky="news")
        self.grid_rowconfigure(current_row, weight=0)
        current_row += 1
        return current_row
    
    def handle_ui_event(self, event):
        if event == "flag_more":
            scalense = float(self.params["scalense"])
            scalense -= 1
            self.params["scalense"] = f"{scalense}"
            self.param_variables["scalense"].set(self.params["scalense"])
            return super().handle_ui_event("run")
        elif event == "flag_fewer":
            scalense = float(self.params["scalense"])
            scalense += 1
            self.params["scalense"] = f"{scalense}"
            self.param_variables["scalense"].set(self.params["scalense"])
            return super().handle_ui_event("run")
        return super().handle_ui_event(event)


def reduce_flatfield(input_row, **kwargs):
    """
    Perform the calstis "basic2d" data reduction. For MAMA exposures, this can be done 
    with a `stistools.wavecal.wavecal` call to derive the individual SHIFTA1 and SHIFTA2 
    values, followed by a single call to `stistools.basic2d`. For CCD exposures, the 
    procedure is considerably more involved.
    
    In order to obtain the maximum possible signal from the exposures, ABSCAL requires 
    interactive iteration on cosmic ray rejection. This is done by repeatedly calling the
    `stistools.ocrreject` function with different parameters until it can be confirmed 
    that the cosmic ray rejection is not rejecting parts of the spectrum as cosmic rays.
    
    In addition, because the STIS absolute flux calibration exposures use CR-SPLIT, all
    exposures in a single visit are taken at the same pointing (i.e. without any 
    intervening dithers), so hot pixels will be at the same place in each exposure. 
    Because the STIS dark correction does not necessarily perfectly remove hot pixels, 
    and because the spectra of the targets is well-known, the ABSCAL calibration 
    introduces the additional step of interpolating over hot pixels in the spectral 
    direction, something that must be done after the general dark correction step.
    
    All of the above is done by repeatedly calling the calstis "basic2d" on the raw STIS
    data with various calibration flag settings, in order to add the additional operations
    in to the appropriate place. Because the STIS calibration program is re-entrant (i.e. 
    you can run part of it on a given input file, do some other work, and then run the rest 
    on that same file without any theoretical loss of correctness, unless your own changes 
    introduced errors), this function calls the "basic2d" function multiple times, with 
    different sets of flags enabled. In particular, the standard CCD basic2d flow is:
    
    - Initialize DQI
    - Bias correction and subtraction
    - Cosmic ray rejection
    - Dark subtraction
    - Flatfield division
    - Populate photometric keywords
    
    The revised flow (including the additional steps added for ABSCAL, and the steps that 
    are usually done as part of calstis, but are not done by basic2d itself) is:
    
    - Initialize DQI, bias correction, bias subtraction
    - Interactive cosmic ray rejection
    - Dark subtraction
    - Hot pixel interpolation
    - Flatfield division
    - Populate photometric keywords
    - Populate SHIFT keywords by calibrating wavecal exposure

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
    crj_file = raw_file.replace("_raw", "_crj")
    final_file = raw_file.replace("_raw", f"_{target}_{mode}_2d")
    wave_file = raw_file.replace("_raw", "_wav")
    preamble = f"{task}: {root}: {detector} ({mode})"
    exp_info = f"{root} {target} {mode}"
    verbose = kwargs.get("verbose", False)
    interp_hot = kwargs.get("interp_hot", True)

    settings = {}
    setting_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if setting_file is not None and os.path.isfile(setting_file):
        with open(setting_file, 'r') as inf:
            settings = yaml.safe_load(inf)
    
    logger.info(f"{preamble}: starting")
    
    if "MAMA" in detector:
        # UV reduction. No cosmic rays issues. Just reduce
        logger.info(f"{preamble}: MAMA reduction")
        basic2d.basic2d(raw_file, output=final_file, verbose=verbose)
    else:
        logger.info(f"{preamble}: CCD reduction.")
        
        # Default calibration flag values
        # dqicorr='perform', 
        # blevcorr='perform', 
        # biascorr='perform', 
        # darkcorr='perform', 
        # flatcorr='perform', 
        # photcorr='perform', 

        # First pass: DQI and Bias
        with fits.open(raw_file, mode="update") as exposure:
            exposure[0].header['HISTORY'] = "ABSCAL: Starting 2d reduction"
            exposure[0].header['DARKCORR'] = "omit"
            exposure[0].header['FLATCORR'] = "omit"
            exposure[0].header['PHOTCORR'] = "omit"
        
        # Do the 2D reduction
        interim_file = os.path.join(path, root+"_firstpass.fits")
        if os.path.isfile(interim_file):
            os.remove(interim_file)
        basic2d.basic2d(raw_file, output=interim_file, verbose=verbose)
        
        with fits.open(interim_file, mode="update") as exposure:
            exposure[0].header['HISTORY'] = "ABSCAL: Finished running basic2d first pass"

        # Set up and run the cosmic ray rejection task
        task_info = {
            "name": "ocrreject",
            "function": ocrreject.ocrreject,
            "frame_class": CosmicRayTaskFrame,
            "metadata": setting_file,
            "start_file": interim_file,
            "end_file": crj_file,
            "output_is_param": False
        }

        standard_kwargs = {"do_log": True, "draw_toolbar": True}
        x1d_params = setup_params("x1d", "stis", settings, input_row, verbose)
        cr_kwargs = {"draw_toolbar": True, "extrsize": x1d_params['extrsize'], "trace_file": "working"}

        figure_windows = {}
        figure_windows["before"] = {
            "title": f"{exp_info}: Before Cosmic Ray Removal",
            "frame_class": ImageFrame,
            "data_file": interim_file,
            "frame_args": [],
            "frame_kwargs": standard_kwargs
        }
        figure_windows["after"] = {
            "title": f"{exp_info}: After Cosmic Ray Removal",
            "frame_class": ImageFrame,
            "data_file": "working",
            "frame_args": [],
            "frame_kwargs": standard_kwargs
        }
        figure_windows["flags"] = {
            "title": f"{exp_info}: Cosmic Ray Flags",
            "window_class": CosmicRayWindow,
            "frame_class": CosmicRayFrame,
            "data_file": interim_file,
            "frame_args": [],
            "frame_kwargs": cr_kwargs
        }

        cr_app = AbscalTask(
            "stis",
            task_info,
            input_row,
            figure_windows,
            verbose
        )
        cr_app.mainloop()
        
        # Intermediate pass for G750L: defringing
        if mode == "G750L" and "fringe" in kwargs and kwargs["fringe"] is not None:
            # Do de-fringing
            with fits.open(crj_file, mode="update") as exposure:
                msg = f"ABSCAL: Running de-fringing with {kwargs['fringe']}"
                exposure[0].header["HISTORY"] = msg
            fringe_file = os.path.join(path, kwargs["fringe"])
            interim_fringe = fringe_file.replace("raw", "nsp")
            if Path(interim_fringe).is_file():
                Path(interim_fringe).unlink()
            normspflat(fringe_file, do_cal=True, wavecal=wave_file)
            fringe_flat = fringe_file.replace("raw", "frr")
            if Path(fringe_flat).is_file():
                Path(fringe_flat).unlink()
            mkfringeflat(crj_file, interim_fringe, fringe_flat)
            defringe(crj_file, fringe_flat)
            crj_file = crj_file.replace("crj", "drj")
            with fits.open(crj_file, mode="update") as exposure:
                exposure[0].header["HISTORY"] = "ABSCAL: Finished de-fringing"
        
        # Second pass: dark correction
        with fits.open(crj_file, mode="update") as exposure:
            exposure[0].header["HISTORY"] = "ABSCAL: Finished running OCRREJECT"
            exposure[0].header["DARKCORR"] = "perform"
            exposure[0].header["HISTORY"] = "ABSCAL: Running dark correction"
        
        interim_file = os.path.join(path, root+"_secondpass.fits")
        if os.path.isfile(interim_file):
            os.remove(interim_file)
        basic2d.basic2d(crj_file, output=interim_file, verbose=verbose)
        
        # Don't do interpolation on flats
        if target.lower() == 'tungsten':
            with fits.open(interim_file, mode='update') as exposure:
                exposure[0].header["HISTORY"] = "ABSCAL: Skipping hot pixel interpolation for TUNGSTEN"
                exposure[0].header["FLATCORR"] = "perform"
                exposure[0].header["PHOTCORR"] = "perform"
                exposure[0].header["HISTORY"] = "ABSCAL: Third pass"
        elif not interp_hot:
            with fits.open(interim_file, mode='update') as exposure:
                exposure[0].header["HISTORY"] = "ABSCAL: Skipping hot pixel interpolation (flagged OMIT)"
                exposure[0].header["FLATCORR"] = "perform"
                exposure[0].header["PHOTCORR"] = "perform"
                exposure[0].header["HISTORY"] = "ABSCAL: Third pass"
        else:
            # Interpolate across hot pixels
            logger.debug(f"{preamble}: Averaging over hot pixels")
            hpix_params = setup_params("hpix", "stis", settings, input_row, verbose)
            hot_pixel_mapping = kwargs['hot_pixel_mapping']
            with fits.open(interim_file) as exposure:
                dark_file = exposure[0].header["DARKFILE"]
            if dark_file not in hot_pixel_mapping:
                logger.debug(f"{preamble}: Creating mapping for {dark_file}")
                hot_pixel_mapping[dark_file] = make_hot_pixel_list(dark_file, hpix_params['threshold'])
            logger.debug(f"{preamble}: Hot pixel interpolation")
            # 
            # The DQ "magic value" of 125 is from the original hot pixel interpolation routine
            # in IDL calstis.
            # 
            with fits.open(interim_file, mode='update') as exposure:
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
                exposure[0].header["FLATCORR"] = "perform"
                exposure[0].header["PHOTCORR"] = "perform"
                exposure[0].header["HISTORY"] = "ABSCAL: Third pass"
        
        # Third pass: flatfield and photometric keywords
        basic2d.basic2d(interim_file, output=final_file, verbose=verbose)
        
    # Finished 2d extraction preparation
    
    # Assign wavelength offsets
    if os.path.isfile(wave_file):
        logger.debug("{}: Running wavecal to assign SHIFT keywords")
        wavecal(final_file, wave_file, verbose=verbose)

    logger.info(f"{preamble}: finished")
    
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
    final_file = os.path.join(path, root+f"_{target}_{mode}_x1d.fits")    
    exp_info = f"{root} {target} {mode}"
    verbose = kwargs.get('verbose', False)
    preamble = f"{task}: {root}: {detector} ({mode})"
    
    settings = {}
    setting_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if setting_file is not None and os.path.isfile(setting_file):
        with open(setting_file, 'r') as inf:
            settings = yaml.safe_load(inf)

    # Vignetting info
    vignetting = 100 + abs(float(input_row['raw_postarg'])) / 0.05
    dist1 = f"{vignetting}"

    logger.info(f"{preamble}: starting")

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

    # Set up and run the extraction task
    task_info = {
        "name": "x1d",
        "function": stis_extract.x1d,
        "frame_class": TaskFrame,
        "metadata": setting_file,
        "start_file": flt_file,
        "end_file": final_file,
        "output_is_param": True
    }


    figure_windows = {}
    extr_kwargs = {"extraction_file": "working", "do_log": True, "draw_toolbar": True}
    figure_windows["region"] = {
        "title": f"{exp_info}: Extraction Region",
        "frame_class": ExtractionInfoFrame,
        "data_file": flt_file,
        "frame_args": [],
        "frame_kwargs": extr_kwargs
    }
    spec_kwargs = {"draw_toolbar": True, "plot_x": "wave"}
    figure_windows["extracted"] = {
        "title": f"{exp_info}: Extracted Spectrum",
        "window_class": SpectrumWindow,
        "frame_class": SpectrumFrame,
        "data_file": "working",
        "frame_args": [],
        "frame_kwargs": spec_kwargs
    }

    ext_app = AbscalTask(
        "stis",
        task_info,
        input_row,
        figure_windows,
        verbose
    )
    ext_app.mainloop()

    logger.info("{preamble}: finished")

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

    logger.info(f"{task}: Starting STIS data reduction for spectroscopic data.")

    settings = {}
    setting_file = get_data_file("abscal.stis", os.path.basename(__file__))
    if setting_file is not None:
        with open(setting_file, 'r') as inf:
            settings = yaml.safe_load(inf)
    
    if update_refs:
        # Make sure the files are set to use the appropriate CRDS references, and that the
        # references have been appropriately retrieved.
        task_verbosity = -1
        logger.debug(f"{task}: Checking Reference Data")
        if verbose:
            task_verbosity = 10
        files = [os.path.join(p,f) for p,f in zip(input_table['path'], input_table['filename'])]
        assign_bestrefs(files, sync_references=True, verbosity=task_verbosity)
    else:
        logger.debug(f"{task}: Skipping reference file update")
    kwargs['hot_pixel_mapping'] = {}

    logger.info(f"{task}: Starting individual file reductions")
    for row in input_table:
        root = row['root']
        target = row['target']
        mode = row['mode']
        preamble = f"{task}: {root}"
        # lamp exposure, for de-fringing
        if target == 'TUNGSTEN':
            continue
        #exposure that needs a fringe flat
        elif mode == 'G750L': 
            kwargs['fringe'] = None
            mask = (input_table['obset'] == row['obset']) & (input_table['target'] == 'TUNGSTEN')
            if len(input_table[mask]) > 0:
                found_flat = False
                for row in input_table[mask]:
                    if row['aperture'] == "0.3X0.09":
                        kwargs['fringe'] = row['filename']
                        found_flat = True
                if not found_flat:
                    kwargs['fringe'] = input_table[mask][0]['filename']

        # Don't extract if there's already an extracted version of
        #   the file present.
        if row['extracted'] != '':
            ext_file = os.path.join(out_dir, row['extracted'])
            if os.path.isfile(ext_file):
                if force:
                    logger.debug(f"{task}: {root}: extracted file exists. Re-extracting.")
                else:
                    logger.debug(f"{task}: {root}: skipping extraction because file exists.")
                    continue
        else:
            extracted_file_name = f"{root}_{target}_x1d.fits"
            extracted_dest = os.path.join(spec_dir, extracted_file_name)

            # If there is already an extracted file for this input, skip.
            if os.path.isfile(extracted_dest):
                ext_str = os.path.join(spec_dir, extracted_file_name)
                if force:
                    logger.debug(f"{task}: {root}: extracted file exists. Re-extracting.")
                else:
                    row['extracted'] = ext_str
                    logger.debug(f"{task}: {root}: skipping extraction because file exists.")
                    continue

        # Only reduce spectroscopic data in the reduce function.
        if row['use']:
            logger.info(f"{task}: Starting {root}")
            logger.debug("{task}: Flatfielding {root} ({row['mode']})")
            reduced_2d = reduce_flatfield(row, **kwargs)
            if reduced_2d is not None:
                row['flatfielded'] = reduced_2d
            
                logger.info(f"{task}: Extracting {root} ({row['mode']})")
                reduced_extracted = reduce_extract(row, **kwargs)
                if reduced_extracted is not None:
                    row['extracted'] = reduced_extracted
            else:
                logger.warning(f"{task}: Flatfielding failed. Skipping extraction.")
            
            logger.info(f"{task}: Finished {root}")
        else:
            msg = f"{task}: Skipping {root} because it's been set to don't use (reason:"
            msg += f" {row['notes']})."
            logger.info(msg)

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
    plots_help = f"Toggle interactive mode (default {plots_default})."
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
    ref_help += f"(default {ref_default})"
    ref_args = ["-r", "--ref_update"]
    ref_kwargs = {'dest': 'ref_update', 'default': ref_default, 'help': ref_help}
    if ref_default:
        ref_kwargs['action'] = 'store_false'
    else:
        ref_kwargs['action'] = 'store_true'
    additional_args['ref_update'] = (ref_args, ref_kwargs)
    
    interp_default = base_defaults['interpolate_hot_pixels']
    if isinstance(interp_default, str):
        interp_default = strtobool(interp_default)
    interp_help = f"Toggle whether to interpolate over hot pixels (default"
    interp_help += f" {interp_default})"
    interp_args = ["--interpolate_hot_pixels"]
    interp_kwargs = {'dest': 'interp_hot', 'default': interp_default, 'help': interp_help}
    if interp_default:
        interp_kwargs['action'] = 'store_false'
    else:
        interp_kwargs['action'] = 'store_true'
    additional_args['interp_hot'] = (interp_args, interp_kwargs)

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
