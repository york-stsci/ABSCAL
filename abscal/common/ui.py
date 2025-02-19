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
from functools import partial
import glob
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.colors import to_hex
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from ruamel.yaml import YAML
import shutil
import tkinter as tk
from tkinter import ttk
import yaml
matplotlib.use('TkAgg')

from .utils import check_params
from .utils import get_param_types
from .utils import get_data_file
from .utils import get_default_params
from .utils import return_type
from .utils import set_override
from .utils import setup_params
from .utils import update_param


class AbscalRoot(tk.Tk):
    """
    Essentially a base class for tasks (organized sets of calibration that can be done in
    a single live GUI instance)
    """

    def run_task(self, task_info, callback=None):
        """
        Very basic implementation
        """
        task_fn = task_info["function"]
        task_name = task_info["name"]
        task_args = task_info["args"]
        task_kwargs = task_info["kwargs"]
        loading_window = LoadingWindow(self, task_name, f"Running {task_name}")
        return_value = task_fn(*task_args, **task_kwargs)
        loading_window.withdraw()
        loading_window.destroy()
        if callback is not None:
            callback(return_value)

    def invoke_event(self, event, values=None):
        """
        Very basic implementation
        """
        if event in self.known_events:
            return self.known_events[event](values)
        return None

    def clear_outputs(self, output_dir, output_prefix):
        """
        Very basic implementation. Search for files in the provided directory with the
        provided prefix and delete them.
        """
        file_list = Path(output_dir).rglob(f"{output_prefix}*")
        for file in file_list:
            file.unlink()

    def finish_up(self, *args, **kwargs):
        """
        Just finish and quit
        """
        self.destroy()


class AbscalWindow(tk.Toplevel):
    def __init__(self, parent, cnf={}, **kwargs):
        self.parent = parent
        self.known_events = {}
        super().__init__(parent, cnf, **kwargs)

    def setup_ui(self):
        """
        Stub
        """
        pass

    def invoke_event(self, event, values=None):
        """
        If the event can be handled locally, handle it locally. Otherwise, pass it to
        parent.
        """
        if event in self.known_events:
            return self.known_events[event](values)
        return self.parent.invoke_event(event, values)

    def run_task(self, task_info, callback=None):
        """
        Run the provided task and, if a callback is provided, invoke the callback
        """
        self.parent.run_task(task_info, callback=None)

    def trigger_update(self, event, values=None):
        """
        Stub
        """
        pass


class LoadingWindow(tk.Toplevel):
    def __init__(self, parent, title, message, cnf={}, **kwargs):
        super().__init__(parent, cnf, **kwargs)
        self.title(title)
        self.geometry("400x200")
        self.resizable(False, False)
        ttk.Label(self, text=message).grid(row=0, column=0, sticky="news")

        # Grab focus and keep topmost while open
        self.wait_visibility()
        self.grab_set()
        self.transient(parent)
        self.attributes('-topmost', True)


class OneColumnWindow(AbscalWindow):
    """
    Window that has a single content frame.
    """
    def setup_ui(self, frame_class, frame_args, frame_kwargs):
        self.content = frame_class(self, *frame_args, **frame_kwargs)
        self.content.grid(row=0, column=0, sticky="news")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

    def trigger_update(self, event, values=None):
        """
        Pass the update to our content
        """
        self.content.trigger_update(event, values)


class TwoColumnWindow(AbscalWindow):
    """
    Window that has a frame each in two columns.
    """
    def setup_ui(self, left_col, right_col):
        self.left_col = left_col
        self.left_col.grid(row=0, column=0, sticky="news")
        self.right_col = right_col
        self.right_col.grid(row=0, column=1, sticky="news")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

    def trigger_update(self, event, values=None):
        """
        Pass the update to our content
        """
        self.left_col.trigger_update(event, values)
        self.right_col.trigger_update(event, values)


class AbscalFrame(ttk.LabelFrame):
    """
    Abstract base class setting up functions for smart windows. Contains functions that
    would otherwise need to be replicated across matplotlib windows, task windows, and 
    save-changes windows
    """
    def __init__(self, parent, name, *args, **kwargs):
        """
        Set up the frame
        """
        self.parent = parent

        # Set up configuration
        kwargs = self.set_up_config(kwargs)

        # Superclass init
        super().__init__(parent, text=name, *args, **kwargs)

        # Set up window layout
        self.layout_frame()

        # Draw contents
        self.draw_frame()

    def set_up_config(self, kwargs):
        """
        Configure Keywords and do basic setup. Stub.
        """
        return kwargs

    def layout_frame(self):
        """
        Create widget contents, and set up their position and weighting. Stub.
        """
        pass

    def draw_frame(self):
        """
        Any frame-specific content drawing goes here. Stub.
        """
        pass

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

    def handle_ui_event(self, event):
        """
        Stub function for handling events targeting the window (presented as a stub so 
        that for subclasses with (e.g.) checkboxes that toggle line visibility the event
        can be handled by the window (because it knows what's going on)).
        
        Always returns boolean false because that's used to figure out if anything has 
        changed WRT the event.
        """
        return False

    def trigger_update(self, event, values=None):
        """
        Handle an update event sent from the parent window. Stub.
        """
        pass


class ImageFrame(AbscalFrame):
    """
    Basic class for creating a window that contains a figure, along with an
    optional matplotlib toolbar, and an optional "Done" button that will close the window
    when pressed (well, it will do so if it's set to do so in the event loop).
    """
    def __init__(self, parent, name, data_file, *args, **kwargs):
        # Set up frame data
        self.data_file = data_file

        # Superclass
        super().__init__(parent, name, *args, **kwargs)

        # Draw figure
        self.canvas.draw()

    def set_up_config(self, kwargs):
        self.draw_toolbar = self.handle_kwarg("draw_toolbar", True, kwargs)
        self.add_done = self.handle_kwarg("add_done", False, kwargs)
        self.file_ext = self.handle_kwarg("file_ext", "SCI", kwargs)
        self.do_log = self.handle_kwarg("do_log", False, kwargs)
        return kwargs

    def create_figure(self):
        data = self.get_figure_data()
        figure = matplotlib.figure.Figure()
        ax = figure.add_subplot(1, 1, 1)
        img = ax.imshow(data, origin='lower', cmap='Greys', resample=False)
        return figure

    def layout_frame(self):
        """
        Set up the frame layout
        """
        self.figure = self.create_figure()
        self.canvas = FigureCanvasTkAgg(self.figure.figure, master=self)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky='news')
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        current_row = 1
        if self.draw_toolbar:
            self.toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
            self.toolbar.grid(row=current_row, column=0, sticky="news")
            self.toolbar.update()
            self.grid_rowconfigure(current_row, weight=0)
            current_row += 1
        if self.add_done:
            b = ttk.Button(self, text="Done", command=partial(self.handle_ui_event, "done"))
            b.grid(row=current_row, column=0, sticky="news")
            self.grid_rowconfigure(current_row, weight=0)

    def draw_frame(self):
        """
        Draw the frame
        """
        self.canvas.draw()

    def get_figure_data(self):
        with fits.open(self.data_file) as inf:
            dat = inf[self.file_ext].data
        if self.do_log:
            dat = np.log10(np.where(dat>0.0001, dat, 0.0001))
        return dat

    def update_figure(self, **kwargs):
        """
        Updates the figure image, using the file and log settings from the initial 
        creation. The keyword arguments aren't used here, but are provided for subclasses
        where it's possible that only a partial update might be needed (e.g. only updating
        plot lines, only updating data image, only updating visibility, etc.)
        """
        ax = self.figure.axes[0]
        data = self.get_figure_data()
        img = ax.get_images()[0]
        img.set_array(data)
        self.figure.canvas.draw()

    def handle_ui_event(self, event):
        if event == "done":
            self.withdraw()

    def trigger_update(self, event, values=None):
        if event == "new_data":
            self.update_figure()


class LegendFrame(AbscalFrame):
    """
    Frame that takes an input figure and creates a legend for that figure.
    """
    def __init__(self, parent, figure, *args, **kwargs):
        self.figure = figure
        super().__init__(parent, "Legend", *args, **kwargs)

    def layout_frame(self):
        """
        Create a text label for each item in the figure.
        """
        legend_style = ttk.Style()
        legend_style.configure('LegendFrame.TLabelFrame', background='white')
        self.configure(style="LegendFrame.TLabelFrame")
        current_row = 0
        for item in self.figure.get_children():
            if hasattr(item, get_label()):
                artist_label = item.get_label()
                if (len(artist_label) > 0) and (artist_label[0] != "_"):
                    color = item.get_color()
                    label = ttk.Label(self, text=artist_label, foreground=color)
                    label.grid(row=current_row, column=0)
                    self.grid_rowconfigure(current_row, weight=0)
                    current_row += 1
        self.grid_columnconfigure(0, weight=0)


class PlotLegendWindow(TwoColumnWindow):
    def setup_ui(self, figure_class, figure_args, figure_kwargs):
        figure_frame = figure_class(self, *figure_args, **figure_kwargs)
        legend_frame = LegendFrame(self, figure_frame.figure)
        super().setup_ui(legend_frame, figure_frame)
        self.grid_columnconfigure(0, weight=0)


class SpectrumFrame(ImageFrame):
    """
    Frame that plots a spectrum with optional element visibility.
    """
    def set_up_config(self, kwargs):
        self.plot_x = self.handle_kwarg('plot_x', 'wave', kwargs)
        self.line_visibility = {
            'GROSS': True,
            'BACKGROUND': True,
            'NET': True,
            'ERROR': True,
            'NET_ERROR': True,
            'DQ': False
        }
        return super().set_up_config(kwargs)

    def get_figure_data(self):
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
        spec_data = self.get_figure_data()
        self.spec_lines = {}
        figure = matplotlib.figure.Figure()
        ax = figure.add_subplot(1, 1, 1)
        for item in spec_data:
            if item != "x_axis":
                line, = ax.plot(spec_data['x_axis'], spec_data[item], label=item)
                self.spec_lines[item] = line
        self._set_line_visibility()
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
        if update_data:
            spec_data = self.get_spec_data()
            for dataset in ['GROSS', 'BACKGROUND', 'NET', 'ERROR', 'NET_ERROR', 'DQ']:
                self.spec_lines[dataset].set_ydata(spec_data[dataset])
        self._set_line_visibility()
        ax.relim(visible_only=True)
        ax.autoscale_view()
        self.figure.canvas.draw()

    def toggle_line(self, line, value):
        self.line_visibility[line] = value
        self.update_figure(update_data=False)    

    def _set_line_visibility(self):
        for line_key in self.spec_lines:
            if line_key in self.line_visibility:
                self.spec_lines[line_key].set_visible(self.line_visibility[line_key])
            else:
                self.spec_lines[line_key].set_visibility(True)

    def trigger_update(self, event, values=None):
        if event == "new_data":
            self.update_figure(update_data=True)


class VisibilityFrame(AbscalFrame):
    def __init__(self, parent, lines, visibility_callback, *args, **kwargs):
        self.visibility_callback = visibility_callback
        self.lines = lines
        self.line_list = {}
        super().__init__(parent, "Legend", *args, **kwargs)

    def layout_frame(self):
        current_row = 0
        for line in self.lines:
            v = tk.BooleanVar(master=self, value=self.lines[line])
            b = ttk.Checkbutton(
                self, 
                text=line, 
                variable=v, 
                onvalue=True, 
                offvalue=False,
                command=partial(self.handle_ui_event, line)
            )
            b.grid(row=current_row, column=0, sticky="news")
            self.grid_rowconfigure(current_row, weight=0)
            current_row += 1
            self.line_list[line] = v
        self.grid_columnconfigure(0, weight=0)

    def handle_ui_event(self, event):
        self.visibility_callback(event, self.line_list[event].get())


class SpectrumWindow(TwoColumnWindow):
    def setup_ui(self, figure_class, figure_args, figure_kwargs):
        figure_frame = figure_class(self, *figure_args, **figure_kwargs)
        visibility_frame = VisibilityFrame(
            self,
            figure_frame.line_visibility,
            figure_frame.toggle_line
        )
        super().setup_ui(visibility_frame, figure_frame)
        self.grid_columnconfigure(0, weight=0)


class TaskFrame(AbscalFrame):
    """
    Create a frame for displaying and editing task parameters
    """
    def __init__(self, parent, task, module, metadata, row, verbose, *args, **kwargs):
        self.verbose = verbose
        self.task = task
        self.module = module
        self.params = setup_params(task, module, metadata, row, self.verbose)
        self.starting_params = deepcopy(self.params)
        self.param_types = get_param_types(task, module, self.verbose)
        self.param_variables = {}
        super().__init__(parent, f"{module}.{task}", *args, **kwargs)

    def layout_frame(self):
        current_row = 0
        for param in self.params:
            if param in self.param_types:
                type_str = self.param_types[param].replace("_", " ")
                type_str = type_str.replace("none", "None")
                label = ttk.Label(self, text=f"{param} ({type_str})")
                label.grid(row=current_row, column=0, sticky="news")
                if self.param_types[param] == "int":
                    v = tk.IntVar(self, value=self.params[param])
                elif self.param_types[param] == "float":
                    v = tk.DoubleVar(self, value=self.params[param])
                else:
                    v = tk.StringVar(self, value=self.params[param])
            else:
                label = ttk.Label(self, text=param)
                label.grid(row=current_row, column=0, sticky="news")
                v = tk.StringVar(self, value=self.params[param])
            entry = ttk.Entry(self, textvariable=v)
            entry.grid(row=current_row, column=1, columnspan=2, sticky="news")
            self.grid_rowconfigure(current_row, weight=1)
            self.param_variables[param] = v
            current_row += 1
        current_row = self.add_to_layout(current_row)
        b = ttk.Button(
            self,
            text="Run with Current Parameters",
            command=partial(self.handle_ui_event, "run")
        )
        b.grid(row=current_row, column=0, sticky="news")
        b = ttk.Button(
            self,
            text="Reset Parameters to Default",
            command=partial(self.handle_ui_event, "reset")
        )
        b.grid(row=current_row, column=1, sticky="news")
        b = ttk.Button(
            self,
            text="Accept Current Parameters",
            command=partial(self.handle_ui_event, "accept")
        )
        b.grid(row=current_row, column=2, sticky="news")
        self.grid_rowconfigure(current_row, weight=0)
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

    def add_to_layout(self, current_row):
        """
        Stub for anything additional that gets added after the parameter list.
        """
        return current_row

    def handle_ui_event(self, event):
        """
        Handle a UI event directed to the window by comparing all of the values found in 
        the window's fields to the values stored in the internal parameters.
        """
        if event == "reset":
            self.params = deepcopy(self.starting_params)
            for param in self.params:
                self.param_variables[param].set(self.params[param])
            return
        elif event == "run":
            self.parent.invoke_event("run", self.params)
        elif event == "accept":
            self.parent.invoke_event("accept", self.params)

    def check_params(self):
        """
        Update the parameter dictionary to make sure that the dictionary contains values
        that are the same as the UI, but of the type that the dictionary expects.
        """
        changed = False
        for param in params:
            item_type = self.param_types[param].split("_")[0]
            current_value = self.param_variables[param].get()
            if "none" in self.param_types[param]:
                gen_type = self.param_types[param].split("_")[0]
                if self.params[param] is None and current_value == "None":
                    continue
                elif self.params[param] is None:
                    changed = True
                    self.params[param] = return_type(item_type, current_value)
                elif current_value == 'None':
                    changed = True
                    self.params[param] = None
                else:
                    typed_current_value = return_type(item_type, current_value)
                    if self.params[param] != typed_current_value:
                        changed = True
                        self.params[param] = typed_current_value
            else:
                typed_current_value = return_type(item_type, current_value)
                if self.params[param] != typed_current_value:
                    changed = True
                    self.params[param] = typed_current_value

        return changed


class ConfirmChangeFrame(AbscalFrame):
    """
    Frame that shows how a parameter has been changed, and asks how the change should be
    applied.
    """
    def __init__(self, parent, params, starting_params, *args, **kwargs):
        self.params = params
        self.starting_params = starting_params
        super().__init__(parent, f"Changed Parameters", *args, **kwargs)

    def layout_frame(self):
        # File Information
        file_frame = ttk.Frame(self)
        file_frame.grid(row=0, column=0, sticky="news")

        self.grid_rowconfigure(0, weight=0)

        # Layout name and notes
        name_frame = ttk.Frame(self)
        name_frame.grid(row=1, column=0, sticky="news")

        l = ttk.Label(name_frame, text="Name:")
        l.grid(row=0, column=0, sticky="news")
        self.name_var = tk.StringVar(self, value="")
        w = ttk.Entry(name_frame, textvariable=self.name_var)
        w.grid(row=0, column=1, sticky="news")
        name_frame.grid_rowconfigure(0, weight=1)

        l = ttk.Label(name_frame, text="Comments:")
        l.grid(row=1, column=0, sticky="news")
        self.comment_var = tk.StringVar(self, value="")
        w = ttk.Entry(name_frame, textvariable=self.comment_var)
        w.grid(row=1, column=1, sticky="news")
        name_frame.grid_rowconfigure(1, weight=1)

        name_frame.grid_columnconfigure(0, weight=0)
        name_frame.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

        # Layout changed parameters
        param_frame = ttk.Frame(self)
        param_frame.grid(row=2, column=0, sticky="news")
        l = ttk.Label(param_frame, text="Parameter")
        l.grid(row=0, column=0, sticky="news")
        l = ttk.Label(param_frame, text="Initial")
        l.grid(row=0, column=1, sticky="news")
        l = ttk.Label(param_frame, text="Current")
        l.grid(row=0, column=2, sticky="news")
        l = ttk.Label(param_frame, text="Action")
        l.grid(row=0, column=3, columnspan=3, sticky="news")
        param_frame.grid_rowconfigure(0, weight=0)
        current_row = 1
        self.param_selections = {}
        for param in self.params:
            if self.params[param] != self.starting_params[param]:
                l = ttk.Label(param_frame, text=param)
                l.grid(row=current_row, column=0, sticky="news")
                l = ttk.Label(param_frame, text=f"{self.starting_params[param]}")
                l.grid(row=current_row, column=1, sticky="news")
                l = ttk.Label(param_frame, text=f"{self.params[param]}")
                l.grid(row=current_row, column=2, sticky="news")
                v = tk.StringVar(self, value="current")
                b = ttk.Radiobutton(
                    self,
                    text="Set for current file",
                    variable=v,
                    value="current"
                )
                b.grid(row=current_row, column=3, sticky="news")
                b = ttk.Radiobutton(
                    self,
                    text="Set as default",
                    variable=v,
                    value="default"
                )
                b.grid(row=current_row, column=4, sticky="news")
                b = ttk.Radiobutton(
                    self,
                    text="Revert to initial value",
                    variable=v,
                    value="revert"
                )
                b.grid(row=current_row, column=5, sticky="news")
                self.param_selections[param] = v
                param_frame.grid_rowconfigure(current_row, weight=0)
        param_frame.grid_columnconfigure(0, weight=0)
        param_frame.grid_columnconfigure(1, weight=0)
        param_frame.grid_columnconfigure(2, weight=0)
        param_frame.grid_columnconfigure(3, weight=0)
        param_frame.grid_columnconfigure(4, weight=0)
        param_frame.grid_columnconfigure(5, weight=0)
        self.grid_rowconfigure(2, weight=0)

        # Layout Buttons
        button_frame = ttk.Frame(self)
        button_frame.grid(row=3, column=0, sticky="news")
        b = ttk.Button(
            button_frame,
            text="Cancel",
            command=partial(self.handle_ui_event, "cancel")
        )
        b.grid(row=0, column=0, sticky="news")
        b = ttk.Button(
            button_frame,
            text="Confirm",
            command=partial(self.handle_ui_event, "confirm")
        )
        b.grid(row=0, column=1, sticky="news")
        button_frame.grid_rowconfigure(0, weight=0)
        button_frame.grid_columnconfigure(0, weight=0)
        button_frame.grid_columnconfigure(1, weight=0)

        self.grid_rowconfigure(3, weight=0)

        self.grid_columnconfigure(0, weight=1)

    def set_up_config(self, kwargs):        
        self.module = self.handle_kwarg('module', '', kwargs)
        self.file = self.handle_kwarg('file', '', kwargs)
        self.row = self.handle_kwarg('row', '', kwargs)

        self.default_file = get_data_file(self.module, "tasks.yaml")
        return kwargs

    def handle_ui_event(self, event):
        """
        Handle a UI event directed to the window by comparing all of the values found in 
        the window's fields to the values stored in the internal parameters.
        """
        if event == "confirm":
            self.param_dict["name"] = self.name_var.get()
            self.param_dict["comments"] = self.comment_var.get()
            self.parent.invoke_event("confirm_changes", self.param_dict)
        elif event == "cancel":
            self.parent.invoke_event("cancel_changes")

    @property
    def param_dict(self):
        """
        Turn the dictionary of tk variables into a dictionary of values
        """
        param_values = {}
        for param in self.param_selections:
            self.param_values[param] = self.param_selections[param].get()
        return param_values


class AbscalTask(AbscalRoot):
    """
    GUI App instance to run a calibration task, display output and parameters, allow the
    user to change the parameters and re-run the task, and when the output is accepted,
    prompt the user about how to 
    """
    def __init__(self, module, task_info, row, figure_windows, verbose, *args, **kwargs):
        self.verbose = verbose
        self.module = module
        self.task_name = task_info["name"]
        self.task_function = task_info["function"]
        self.task_frame_class = task_info["frame_class"]
        self.metadata_file = task_info["metadata"]
        self.task_meta = self.read_meta(self.metadata_file)
        self.start_file = Path(task_info["start_file"])
        self.end_file = Path(task_info["end_file"])
        self.interim_prefix = f"interim_{self.task_name}_"
        self.interim_file_name = f"{self.interim_prefix}{self.end_file.name}"
        self.interim_file = self.end_file.parent / self.interim_file_name
        self.output_is_param = task_info["output_is_param"]
        self.row = row
        self.figure_windows = figure_windows
        self.confirm_changes_window = None

        super().__init__(*args, **kwargs)
        self.withdraw()
        self.display_windows = {}
        self.clear_outputs()
        self.setup_task_window()

        # So my plan is this:
        #   - create the task window
        #   - "Accept" means
        #       - Ask for confirmation of parameter changes
        #       - Run the task
        #       - Quit
        #   - "Try" means run the task, display the output
        #   - "Cancel" means quit

    def setup_task_window(self):
        self.task_window = OneColumnWindow(self)
        self.task_window.setup_ui(
            self.task_frame_class,
            (self.task_name, self.module, self.task_meta, self.row, self.verbose),
            {}
        )
        self.starting_parameters = deepcopy(self.task_window.content.starting_params)
        self.accepted_parameters = {}

    def read_meta(self, meta_file):
        """
        Open and return a metadata dictionary from a metadata YAML file
        """
        with open(meta_file, "rt") as input_file:
            metadata = yaml.safe_load(input_file)
        return metadata

    def invoke_event(self, event, values=None):
        """
        Receive an event passed up by a window.
        """
        if event == "run":
            # Run the task, with the supplied values as the task parameters
            self.run_task(values)
        elif event == "accept":
            # Run the task, with the supplied values as the task parameters, and then
            # run the "confirm changes" interface.
            self.accepted_parameters = deepcopy(values)
            self.task_window.destroy()
            self.task_window = None
            self.run_task(values, callback=self.confirm_changes)
        elif event == "confirm_changes":
            for param in values:
                update_param(
                    self.task,
                    self.module,
                    self.row,
                    self.metadata_file,
                    param,
                    self.accepted_parameters[param],
                    values["name"],
                    values["comment"],
                    values[param]
                )
            self.finish_up()
            self.destroy()
        elif event == "cancel_changes":
            self.accepted_parameters = deepcopy(self.starting_parameters)
            self.run_task(self.starting_parameters, callback=self.finish_up)
            self.destroy()
        else:
            title = "Error: Unknown Event"
            message = f"Error: Unknown event {event} with values {values}"
            tk.messagebox.showerror(title=title, message=message)

    def clear_outputs(self):
        output_path = self.end_file.parent
        prefix = f"interim_{self.task_name}"
        if self.end_file.is_file():
            self.end_file.unlink()
        super().clear_outputs(output_path, prefix)

    def run_task(self, task_info, callback=None):
        self.clear_outputs()
        interim_prefix = f"interim_{self.task_name}_"
        interim_file = self.end_file.parent / f"{interim_prefix}{self.end_file.name}"
        task_dict = {
            "name": self.task_name,
            "function": self.task_function,
            "args": [self.start_file.as_posix()],
            "kwargs": task_info
        }
        if self.output_is_param:
            task_dict["kwargs"]["output"] = self.interim_file.as_posix()
        else:
            task_dict["args"].append(self.interim_file.as_posix())
        super().run_task(task_dict, callback)
        self.check_task_completion()

    def check_task_completion(self):
        if not self.interim_file.is_file():
            title = "Error: No output file created"
            message = f"{self.module}.{self.task} failed to create an output file."
            tk.messagebox.showerror(title=title, message=message)
            return
        self.refresh_figure_windows()

    def refresh_figure_windows(self):
        for window_name in self.figure_windows:
            if window_name in self.display_windows:
                self.display_windows[window_name].trigger_update("new_data")
            else:
                window_info = self.figure_windows[window_name]
                window_title = window_info["title"]
                window_class = window_info.get("window_class", OneColumnWindow)
                frame_class = window_info["frame_class"]
                frame_kwargs = window_info["frame_kwargs"]
                for item in frame_kwargs:
                    if frame_kwargs[item] == "working":
                        frame_kwargs[item] = self.interim_file
                if window_info["data_file"] == "working":
                    data_file = self.interim_file
                else:
                    data_file = window_info["data_file"]
                frame_args = [window_name, data_file] + window_info["frame_args"]
                window = window_class(self)
                window.setup_ui(frame_class, frame_args, frame_kwargs)
                self.display_windows[window_name] = window

    def finish_up(self):
        self.end_file.write_bytes(self.interim_file.read_bytes())
        if self.task_window is not None:
            self.task_window.destroy()
        if self.confirm_changes_window is not None:
            self.confirm_changes_window.destroy()
        for item in self.display_windows:
            self.display_windows[item].destroy()

    def confirm_changes(self, return_values):
        """
        Run the interface to confirm parameter changes.
        """
        any_changes = False
        for item in self.starting_parameters:
            if self.accepted_parameters[item] != self.starting_parameters[item]:
                any_changes = True
        if any_changes:
            self.confirm_changes_window = OneColumnWindow(self)
            cc_args = [self.accepted_parameters, self.starting_parameters]
            cc_kwargs = {}
            self.confirm_changes_window.setup_ui(ConfirmChangeFrame, cc_args, cc_kwargs)
            return
        self.finish_up()
        self.destroy()
