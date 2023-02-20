#! /usr/bin/env python
"""
This module includes utility plotting functions, 

Authors
-------
- Brian York

Use
---
Individual functions from this module are intended to be imported where
needed::

    from abscal.common.plots import make_figure_window
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
import shutil
import yaml
matplotlib.use('TkAgg')


class Toolbar(NavigationToolbar2Tk):
    """
    Subclass of matplotlib/TK toolbar to let you add it to a window canvas, while 
    linking it to a figure.
    """
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)


def draw_figure(window, figure, canvas_key='canvas', draw_toolbar=True):
    """
    Given a window with a TK canvas, draw a matplotlib figure in that window (or update
    an already existing figure in the window). Optionally, if the window is set up for 
    a matplotlib toolbar to be drawn, draw said toolbar.
    
    Parameters
    ----------
    window : PySimpleGUI.Window
        The window that has an appropriate canvas
    figure : matplotlib.figure.Figure
        The figure to draw on the canvas
    canvas_key : str, default 'canvas'
        The key to select the canvas from the window
    draw_toolbar : bool, default False
        Whether to draw a matplotlib toolbar for the figure
    
    Returns
    -------
    figure_canvas_agg : matplotlib.backends.backend_tkagg.FigureCanvasTkAgg
        A reference to the canvas as a TK widget
    """
    canvas = window[canvas_key].TKCanvas
    if draw_toolbar:
        control_key = canvas_key+'_controls'
        print("Looking for canvas controls {}".format(control_key))
        canvas_toolbar = window[control_key].TKCanvas
        if canvas.children:
            for child in canvas.winfo_children():
                child.destroy()
        if canvas_toolbar.children:
            for child in canvas_toolbar.winfo_children():
                 child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    if draw_toolbar:
        toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
        toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def make_figure_window(window_title, draw_toolbar=True, add_done=False):
    """
    Creates a window that has a canvas element (for drawing a figure) and, optionally,
    a matplotlib toolbar to allow for matplotlib interactions with the figure.
    
    Parameters
    ----------
    window_title : str
        The title to put on the window
    draw_toolbar : bool, default True
        Whether to include a matplotlib toolbar widget
    add_done : bool, default False
        Whether to put a "Done" button on the figure that will close the window.
    
    Returns
    -------
    window : PySimpleGUI.Window
        The created window
    """
    layout = [[sg.Canvas(key='canvas')]]
    if draw_toolbar:
        layout.append([sg.Canvas(key='canvas_controls')])
    if add_done:
        layout.append([sg.Button('Done')])
    window = sg.Window(window_title, layout, finalize=True, element_justification='center', 
                       font='Helvetica 18')
    return window


def make_spectrum_window(window_title, fits_file, plot_x='wave', draw_toolbar=True, add_done=False):
    """
    Creates a window that has two columns:
    
    - A lefthand column with a list of spectral elements (e.g. gross, net, flux, etc.)
      that can be plotted, provided as checkboxes
    - A righthand column with a matplotlib spectrum, showing the selected elements of the
      spectrum.
    
    The window also optionally has two bottom frames:
    
    - A "Done" button to close the window
    - A matplotlib toolbar
    
    The function returns the window object itself, as well as a list of matplotlib lines,
    so that their visibilities can be toggled.
    
    Parameters
    ----------
    window_title : str
        The title to put on the window
    fits_file : str
        The file name and path to a FITS file that has a spectrum data table in its 'SCI'
        extension
    plot_x : str, default 'wave'
        The x axis to use on the plot, either pixel or wavelength
    draw_toolbar : bool, default True
        Whether to include a matplotlib toolbar widget
    add_done : bool, default False
        Whether to put a "Done" button on the figure that will close the window.
    
    Returns
    -------
    window : PySimpleGUI.Window
        The created window
    fig : matplotlib.figure.Figure
        The created figure, to be used for things like draw_figure
    line_list : dict
        A list linking checkbox titles/line names to their associated matplotlib line
        objects.
    """
    fig, line_list = make_spectrum_figure(fits_file, plot_x=plot_x)
    check_list = []
    for item in line_list:
        color = line_list[item].get_color()
        default_val = True
        if item == 'DQ':
            default_val = False
        check_list.append([sg.Checkbox(item, default=default_val, text_color=color, 
            key=item, enable_events=True, background_color='white')])
    check_col = sg.Column([[sg.Frame('Lines:', check_list, background_color='white',
        title_color='black')]])
    plot_frame = [[sg.Canvas(key='canvas')]]
    if draw_toolbar:
        plot_frame.append([sg.Canvas(key='canvas_controls')])
    plot_col = sg.Column([[sg.Frame('Plot:', plot_frame)]])
    layout = [[check_col, plot_col]]
    if add_done:
        layout.append([sg.Button('Done')])
    window = sg.Window(window_title, layout, finalize=True)
    return window, fig, line_list


def make_spectrum_figure(fits_file, plot_x='wave'):
    """
    Create a matplotlib figure from a FITS file containing a table spectrum in its first
    extension, and return a set of lines and colors to be used to make the spectrum
    selectable.
    
    Parameters
    ----------
    fits_file : str
        The file name and path to a FITS file that has a spectrum data table in its 'SCI'
        extension
    plot_x : str, default 'wave'
        The x axis to use on the plot, either pixel or wavelength
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure, to be used for things like draw_figure
    line_list : dict
        A list linking checkbox titles/line names to their associated matplotlib line
        objects
    """
    spec_data = {}
    spec_lines = {}
    with fits.open(fits_file) as x1d:
        if plot_x == 'wave':
            x_axis = x1d[1].data['WAVELENGTH'][0]
        else:
            x_axis = np.arange(x1d[1].data['SPORDER'][0])
        for dataset in ['GROSS', 'BACKGROUND', 'NET', 'ERROR', 'NET_ERROR', 'DQ']:
            if dataset in x1d[1].data.columns.names:
                spec_data[dataset] = x1d[1].data[dataset][0]
    fig = matplotlib.figure.Figure(figsize=(12,8), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    for item in spec_data:
        line = ax.plot(x_axis, spec_data[item], label=item)
        if item == 'DQ':
            line[0].set_visible(False)
        spec_lines[item] = line[0]
    if plot_x == 'wave':
        ax.set_xlabel("Wavelength (Angstroms)")
    else:
        ax.set_xlabel("Pixel")
    ax.relim(visible_only=True)
    ax.autoscale_view()
    return fig, spec_lines


def make_img_fig(data):
    """
    Create a matplotlib figure that is based around an imshow of a 2D image (likely a
    FITS file of an observation). The data should be pre-processed however is required
    before being sent to this function (e.g. if you want to plot log data, take the log
    before sending in the data).
    
    Note that additional commands can be sent to the return figure to add additional 
    material to the plot.
    
    Note that the figure size is optimized around showing a 1024x1024 image.
    
    Parameters
    ----------
    data : numpy.array
        2D numpy array of data (or anything else that can be sent to pyplot.imshow)
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig = matplotlib.figure.Figure(figsize=(11, 11), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    img = ax.imshow(data, origin='lower')
    return fig    


def make_params_window(window_title, param_title, param_dict, extra_layout=[]):
    """
    Create a window for editing task parameters.
    
    Parameters
    ----------
    window_title : str
        The title for the window
    param_title : str
        The text to put above the parameter list
    param_dict : dict
        The dictionary of parameters (and their current values)
    extra_layout : list, default []
        Extra items to add below the parameter list and above the bottom buttons
    
    Returns
    -------
    window : PySimpleGUI.Window
        The created window
    """
    # Edit items UI
    layout = [[sg.Text(param_title)]]
    for item in param_dict:
        label = sg.Text("{}:".format(item))
        field = sg.InputText(default_text="{}".format(param_dict[item]), 
                             key='{}'.format(item))
        layout.append([label, sg.Push(), field])
    for item in extra_layout:
        layout.append(item)
    layout.append([sg.Button('Plot'), sg.Button('Reset'), sg.Push(), sg.Button('Accept')])
    window = sg.Window(window_title, layout, finalize=True)
    return window
