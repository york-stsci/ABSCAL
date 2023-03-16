#! /usr/bin/env python
"""
This module includes utility plotting functions, 

Authors
-------
- Brian York

Use
---
Individual classes from this module are intended to be imported where
needed::

    from abscal.common.plots import ImageWindow
"""

from astropy.io import ascii, fits
from collections import defaultdict
from copy import deepcopy
import datetime
from distutils.util import strtobool
import glob
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.colors import to_hex
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import PySimpleGUI as sg
from ruamel.yaml import YAML
import shutil
import yaml
matplotlib.use('TkAgg')

from .utils import check_params
from .utils import get_param_types
from .utils import set_override
from .utils import setup_params


class Toolbar(NavigationToolbar2Tk):
    """
    Subclass of matplotlib/TK toolbar to let you add it to a window canvas, while 
    linking it to a figure.
    """
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)


class AbscalWindow(sg.Window):
    """
    Abstract base class setting up functions for smart windows. Contains functions that
    would otherwise need to be replicated across matplotlib windows, task windows, and 
    save-changes windows
    """
    def handle_kwarg(self, key, default, kwargs, remove=True, add=False):
        """
        If the key is in the dictionary, return the associated value (and optionally 
        remove it from the dictionary so that the super() init doesn't have an "unknown 
        keyword" issue). Otherwise take the default value. Either way, return the value.
        """
        value = default
        if key in kwargs:
            value = kwargs[key]
            if remove:
                del kwargs[key]
        elif add:
            kwargs[key] = value
        return value

    def handle_ui_event(self, event, values):
        """
        Stub function for handling events targeting the window (presented as a stub so 
        that for subclasses with (e.g.) checkboxes that toggle line visibility the event
        can be handled by the window (because it knows what's going on)).
        
        Always returns boolean false because that's used to figure out if anything has 
        changed WRT the event.
        """
        return False


class ImageWindow(AbscalWindow):
    """
    Basic class for creating a PySimpleGUI window that contains a figure, along with an
    optional matplotlib toolbar, and an optional "Done" button that will close the window
    when pressed (well, it will do so if it's set to do so in the event loop).
    """
    def __init__(self, window_title, data_file, *args, **kwargs):
        # Set up configuration
        kwargs = self.set_up_config(kwargs)
        # Set data file
        self.data_file = data_file
        # Set up matplotlib figure
        self.figure = self.create_figure()
        # Set up window layout
        layout = self.create_layout()
        # Create window and draw figure
        super().__init__(window_title, layout, *args, **kwargs)
        self.draw_figure()

    def set_up_config(self, kwargs):
        self.draw_toolbar = self.handle_kwarg("draw_toolbar", True, kwargs)
        self.add_done = self.handle_kwarg("add_done", False, kwargs)
        self.file_ext = self.handle_kwarg("file_ext", "SCI", kwargs)
        self.do_log = self.handle_kwarg("do_log", False, kwargs)
        self.handle_kwarg("font", "Helvetica 16", kwargs, remove=False, add=True)
        self.handle_kwarg("finalize", True, kwargs, remove=False, add=True)
        return kwargs

    def create_figure(self):
        data = self.get_image_data()
        figure = matplotlib.figure.Figure(figsize=(11, 11), dpi=100)
        ax = figure.add_axes((0.04, 0.04, 0.95, 0.95))
        img = ax.imshow(data, origin='lower', cmap='Greys', resample=False)
        return figure

    def create_layout(self):
        layout = [[sg.Canvas(key='canvas')]]
        if self.draw_toolbar:
            layout.append([sg.Canvas(key='canvas_controls')])
        if self.add_done:
            layout.append([sg.Button('Done')])
        return layout

    def get_image_data(self):
        with fits.open(self.data_file) as inf:
            dat = inf[self.file_ext].data
        if self.do_log:
            dat = np.log10(np.where(dat>0.0001, dat, 0.0001))
        return dat

    def draw_figure(self):
        """
        Draw the matplotlib figure stored in the ImageWindow class on the TK canvas 
        present in the window. If the window is set up for a matplotlib toolbar, draw 
        the toolbar.
        """
        canvas = self["canvas"].TKCanvas
        control_key = "canvas_controls"
        if control_key in self.AllKeysDict:
            canvas_toolbar = self[control_key].TKCanvas
            if canvas.children:
                for child in canvas.winfo_children():
                    child.destroy()
            if canvas_toolbar.children:
                for child in canvas_toolbar.winfo_children():
                     child.destroy()
        self.figure_canvas_agg = FigureCanvasTkAgg(self.figure, canvas)
        self.figure_canvas_agg.draw()
        if control_key in self.AllKeysDict:
            toolbar = Toolbar(self.figure_canvas_agg, canvas_toolbar)
            toolbar.update()
        self.figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)

    def update_figure(self, **kwargs):
        """
        Updates the figure image, using the file and log settings from the initial 
        creation. The keyword arguments aren't used here, but are provided for subclasses
        where it's possible that only a partial update might be needed (e.g. only updating
        plot lines, only updating data image, only updating visibility, etc.)
        """
        ax = self.figure.axes[0]
        data = self.get_image_data()
        img = ax.get_images()[0]
        img.set_array(data)
        self.figure.canvas.draw()


class TwoColumnWindow(ImageWindow):
    """
    Window that has a plot in the right column, and something else in the left column.
    """
    def create_layout(self):
        ui_column = self.make_ui_column()
        plot_frame = [[sg.Canvas(key='canvas')]]
        if self.draw_toolbar:
            plot_frame.append([sg.Canvas(key='canvas_controls')])
        plot_col = sg.Column([[sg.Frame('Plot:', plot_frame)]])
        layout = [[ui_column, plot_col]]
        if self.add_done:
            layout.append([sg.Button('Done')])
        return layout
    
    def make_ui_column(self):
        """
        Stub function
        """
        pass


class SpectrumWindow(TwoColumnWindow):
    """
    Window that puts spectral plot visibility selectors in the left column.
    """
    def set_up_config(self, kwargs):
        self.plot_x = self.handle_kwarg('plot_x', 'wave', kwargs)
        return super().set_up_config(kwargs)
    
    def make_ui_column(self):
        check_list = []
        for line in self.spec_lines:
            color = self.spec_lines[line].get_color()
            val = True
            if line == "DQ":
                val = False
            checkbox = sg.Checkbox(line, default=val, text_color=color, key=line, 
                                   enable_events=True, background_color='white')
            check_list.append([checkbox])
        check_col = sg.Column([[sg.Frame('Lines:', check_list, background_color='white', 
                                         title_color='black')]])
        return check_col
    
    def get_spec_data(self):
        spec_data = {}
        with fits.open(self.data_file) as x1d:
            if self.plot_x == 'wave':
                spec_data['x_axis'] = x1d[1].data['WAVELENGTH'][0]
            else:
                spec_data['x_axis'] = np.arange(x1d[1].data['SPORDER'][0])
            for dataset in ['GROSS', 'BACKGROUND', 'NET', 'ERROR', 'NET_ERROR', 'DQ']:
                if dataset in x1d[1].data.columns.names:
                    spec_data[dataset] = x1d[1].data[dataset][0]
        return spec_data
    
    def create_figure(self):
        spec_data = self.get_spec_data()
        self.spec_lines = {}
        figure = matplotlib.figure.Figure(figsize=(12,8), dpi=100)
        ax = figure.add_subplot(1, 1, 1)
        for item in spec_data:
            if item != "x_axis":
                line = ax.plot(spec_data['x_axis'], spec_data[item], label=item)
                if item == 'DQ':
                    line[0].set_visible(False)
                self.spec_lines[item] = line[0]
        if self.plot_x == 'wave':
            ax.set_xlabel("Wavelength (Angstroms)")
        else:
            ax.set_xlabel("Pixel")
        ax.relim(visible_only=True)
        ax.autoscale_view()
        return figure

    def update_figure(self, **kwargs):
        """
        Updates figure, either by changing line visibility or by changing line data and,
        potentially, visibility.
        """
        update_data = kwargs.get("update_data", True)
        ax = self.figure.axes[0]
        if kwargs["update_data"]:
            spec_data = self.get_spec_data()
        for dataset in ['GROSS', 'BACKGROUND', 'NET', 'ERROR', 'NET_ERROR', 'DQ']:
            if kwargs["update_data"]:
                self.spec_lines[dataset].set_ydata(spec_data[dataset])
            self.spec_lines[dataset].set_visible(self.find_element(dataset).Get())
        ax.relim(visible_only=True)
        ax.autoscale_view()
        self.figure.canvas.draw()


class TaskWindow(AbscalWindow):
    """
    Create a window for displaying and editing task parameters
    """
    def __init__(self, task, module, metadata, row, *args, **kwargs):
        # Set up configuration
        kwargs = self.set_up_config(kwargs)
        # Set task info
        self.task = task
        self.module = module
        self.params = setup_params(task, module, metadata, row, self.verbose)
        self.starting_params = deepcopy(self.params)
        self.param_types = get_param_types(task, module, self.verbose)
        # Set up window layout
        layout = self.create_layout()
        # Create window and draw figure
        super().__init__("{} parameters".format(task), layout, *args, **kwargs)

    def set_up_config(self, kwargs):
        self.verbose = self.handle_kwarg('verbose', False, kwargs)
        self.handle_kwarg("font", "Helvetica 16", kwargs, remove=False, add=True)
        self.handle_kwarg("finalize", True, kwargs, remove=False, add=True)
        return kwargs
    
    def create_layout(self):
        layout = [[sg.Text(self.task)]]
        for item in self.params:
            if item in self.param_types:
                item_type_str = self.param_types[item].replace("_", " ").replace("none", "None")
                item_label = "{} ({})".format(item, item_type_str)
            else:
                item_label = "{}".format(item)
            label = sg.Text(item_label)
            field = sg.InputText(default_text="{}".format(self.params[item]), 
                                 key='{}'.format(item))
            layout.append([label, sg.Push(), field])
        layout = self.add_to_layout(layout)
        layout.append([sg.Button('Plot'), sg.Button('Reset'), sg.Push(), sg.Button('Accept')])
        return layout
    
    def add_to_layout(self, layout):
        """
        Stub function to add additional buttons to the layout (but have them above the
        standard buttons)
        """
        return layout

    def handle_ui_event(self, event, values):
        """
        Handle a UI event directed to the window by comparing all of the values found in 
        the window's fields to the values stored in the internal parameters.
        """
        if event == "Reset":
            self.params = deepcopy(self.starting_params)
            return True
        return check_params(self.task, self.module, self.params, values, self.verbose)


def handle_parameter_window(name, start_value, current_value, file, row):
    """
    Check whether the user has changed a parameter. If they have, ask whether they want to
    update the reference file. If they do, update the file.
    """
    if start_value != current_value:
        msg = "You have changed {} from {} to {}.".format(name, start_value, current_value)
        layout = [[sg.Text("Change {}?".format(name), font='Helvetica 18')],
                  [sg.Text(msg, font='Helvetica 16')],
                  [sg.Text("Name:"), sg.Push(), sg.InputText(key='name')],
                  [sg.Text("Reason:"), sg.Push(), sg.InputText(key='reason')],
                  [sg.Push()],
                  [sg.Button('Save Changes'), sg.Push(), sg.Button('Discard Changes')]]
        window = sg.Window("Save Changes to {}?".format(name), layout, finalize=True)
        while True:
            event, values = window.read()
            if event == "Save Changes":
                if values["name"] == "":
                    sg.popup_error("Name field must not be blank.")
                    continue
                elif values["reason"] == "":
                    sg.popup_error("Reason field must not be blank.")
                    continue
                yaml = YAML()
                with open(file) as inf:
                    config = yaml.load(inf)
                set_override(name, config, row, values["name"], current_value, values["reason"])
                with open(file, mode="w") as outf:
                    yaml.dump(config, outf)
                break
            elif event == "Discard Changes":
                break
        window.close()


def run_task(module, task, settings_file, row, task_fn, start_file, end_file, **kwargs):
    """
    Run a calibration task, including displaying output data and allowing the user to 
    adjust parameters. Then check the final parameters against the original parameters, 
    and offer to save any changes.
    """
    verbose = kwargs.get("verbose", False)
    figure_windows = kwargs.get("figure_windows", [])
    task_info = kwargs.get("task_info", row["root"])
    task_window_class = kwargs.get("task_window_class", TaskWindow)
    output_is_param = kwargs.get("output_is_param", True)
    
    if verbose:
        print("{}.{}: Starting".format(module, task))
    
    # Remove final output file (if present)
    if os.path.isfile(end_file):
        os.remove(end_file)

    # Set up working and final files
    working_dir, end_file_name = os.path.split(end_file)
    interim_file_name = "interim_"+end_file_name
    interim_file = os.path.join(working_dir, interim_file_name)
    
    # Set up metadata dictionary
    with open(settings_file) as input_file:
        settings = yaml.safe_load(input_file)

    # Set up lists for plot windows and figures
    static_windows = []
    dynamic_windows = []
    
    task_window = task_window_class(task, module, settings, row, finalize=False)
    
    # Do first task run
    if os.path.isfile(interim_file):
        os.remove(interim_file)
    if output_is_param:
        task_fn(start_file, output=interim_file,  **task_window.params)
    else:
        task_fn(start_file, interim_file, **task_window.params)
    
    # Check for the output file actually being created
    if not os.path.isfile(interim_file):
        raise ValueError("Task {}.{} did not produce an output file".format(module, task))
    
    # Set up figure windows
    for window_dict in figure_windows:
        window_title = window_dict["title"]
        static = window_dict.get("static", False)
        window_class = window_dict.get("window_class", ImageWindow)
        window_kwargs = window_dict.get("kwargs", {})
        for item in window_kwargs:
            if window_kwargs[item] == "working":
                window_kwargs[item] = interim_file
        if window_dict["data_file"] == "working":
            data_file = interim_file
        else:
            data_file = window_dict["data_file"]
        window = window_class(window_title, data_file, **window_kwargs)
        if static:
            static_windows.append(window)
        else:
            dynamic_windows.append(window)
    
    # Set up the parameter window
    task_window.finalize()
    
    changed = False
    while True:
        if changed:
            if verbose:
                print("\t{}: Re-running task with new parameter values".format(task))
            # Re-run task and dynamic plots
            if os.path.isfile(interim_file):
                os.remove(interim_file)
            if output_is_param:
                task_fn(start_file, output=interim_file,  **task_window.params)
            else:
                task_fn(start_file, interim_file, **task_window.params)
            # Check for the output file actually being created
            if not os.path.isfile(interim_file):
                raise ValueError("Task {}.{} did not produce an output file".format(module, task))
            for window in dynamic_windows:
                window.update_figure(update_data=True)

        window, event, values = sg.read_all_windows()
        if window is None and event != sg.TIMEOUT_EVENT:
            print('\t{}: Exiting because no windows are left'.format(task))
            break
        elif event == sg.WIN_CLOSED or event == 'Exit':
            # Closing the parameter window is the equivalent of "Accept"
            if window == task_window:
                break
            window.close()
        elif window in dynamic_windows:
            window.update_figure(update_data=False)
        elif event == 'Accept':
            break
        elif window == task_window:
            changed = task_window.handle_ui_event(event, values)
        else:
            msg = "\t{}: Got unknown event {} for window {} with values {}"
            print(msg.format(preamble, event, window, values))


    # Get the final parameters
    final_params = task_window.params
    initial_params = task_window.starting_params
    
    task_window.close()
    for window in static_windows:
        window.close()
    for window in dynamic_windows:
        window.close()
    
    if changed:
        if output_is_param:
            task_fn(start_file, output=end_file,  **task_window.params)
        else:
            task_fn(start_file, end_file, **task_window.params)
    else:
        shutil.copy(interim_file, end_file)
        
    if os.path.isfile(interim_file):
        os.remove(interim_file)

    for item in final_params:
        handle_parameter_window(item, initial_params[item], final_params[item], settings_file, row)

    if verbose:
        print("{}.{}: Finished".format(module, task))
    # Done.
