#! /usr/bin/env python
"""
This module includes general utility functions.

Authors
-------
- Brian York

Use
---
Individual functions from this module are intended to be imported where
needed::

    from abscal.common.utils import absdate
"""

import glob, os, sys, yaml

import numpy as np

from astropy.time import Time
from astropy.io import fits
from copy import deepcopy
from datetime import datetime
from dateutil.parser import parse as date_parse
from distutils.util import strtobool
from ruamel.yaml import YAML
import shutil
from simpleeval import simple_eval

from .file_utils import get_data_dir, get_data_file
from .logging import DEFAULT_LOGGER as logger


def absdate(pstrtime):
    """
    Get the date in decimal years. 
    
    This is used to figure out a target position, given that we have the target 
    co-ordinates, the co-ordinate epoch, the annual proper motion, and the observation 
    time.
    
    Parameters
    ----------
    pstrtime : str or Time or datetime
        The observation start time
    
    Returns
    -------
    dt : float
        The observation year + fractional (decimal) year
    """
    # '2013.057:04:24:48'
    # pstrtime is in the format yyyy.ddd:hh:mm:ss where 'ddd' is the decimal
    #   date (i.e. January 1 is 001, January 31 is 031, February 1 is 032,
    #   December 31 is 365 or 366 depending on leap year status)
    if isinstance(pstrtime, datetime):
        dt = pstrtime
    elif isinstance(pstrtime, Time):
        dt = pstrtime.datetime
    else:
        dt = datetime.strptime(pstrtime, "%Y.%j:%H:%M:%S")
    next_year = datetime(year=dt.year+1, month=1, day=1)
    this_year = datetime(year=dt.year, month=1, day=1)
    year_part = dt - this_year
    year_length = next_year - this_year
    return dt.year + year_part/year_length


def _extract_dict(input_dict, output_dict, input_keys):
    """
    Recursively extract values from a defaults dictionary.
    
    A defaults dictionary consists of:
    
    - an optional "all" key
    - zero or more other keys, each of whose values is a defaults dictionary
    
    The goal is to add any matching values to an output dictionary, with more specific 
    matching values overriding less specific matching values. As such, given an input 
    dictionary and a list of keywords,
    
    - Add all key/value pairs from the "all" dictionary (if present) to the output 
      dictionary.
    - For each keyword in the list, if that keyword is in the dictionary, call this 
      function recursively on the value of that key, which is (see above) a dictionary.
    - Don't check on whether a value already exists in the output dictionary, because 
      more-specific overrides less-specific (if you need a default for a specific value to 
      definitely override a more general default, nest that value as a keyword inside the 
      more general dictionary).
    
    Parameters
    ----------
    input_dict : dict
        The dictionary to search
    output_dict : dict
        The dictionary to build from
    input_keys : list
        A list of keys to search for
        
    Returns
    -------
    output_dict : dict
        The edited output dictionary
    """
    if "all" in input_dict:
        for keyword in input_dict["all"].keys():
            output_dict[keyword] = input_dict["all"][keyword]
    
    for keyword in input_keys:
        if keyword in input_dict:
            output_dict = _extract_dict(input_dict[keyword], output_dict, input_keys)
    
    return output_dict


def get_defaults(module, *args):
    """
    Find an internal defaults data file, load it using YAML, and return the resulting 
    dictionary.
    
    Takes the dot-separated module path (e.g. "abscal.wfc3.reduce_grism_extract"), splits 
    off the last item (e.g. ["abscal.wfc3", "reduce_grism_extract"]), adds ".yaml" to the 
    end of the second item (e.g. ["abscal.wfc3", "reduce_grism_extract.yaml"]), adds 
    ".defaults" to the first item 
    (e.g. ["abscal.wfc3.defaults", "reduce_grism_extract.yaml"]), and feeds the result 
    into :code:`get_data_file()`. Then loads the resulting file as a dictionary, and 
    builds a new dictionary consisting of:
    
    - All key/value pairs in the "all" dictionary
    - All key/value pairs in any dictionary matching any of the keyword arguments
    - The above two items from any dictionary matching any of the keyword arguments, 
      extending recursively into the depths of the dictionary.
    
    The result will be a flat (i.e. single-level) dictionary.

    Parameters
    ----------
    module : str
        The module to search in, using standard dot separators (e.g. abscal.wfc3)
    args : list
        A list of specific keyword arguments, provided to ensure the inclusion of 
        specific sub-values or sub-dictionaries.

    Returns
    -------
    defaults : dict
        Dictionary of default parameters.
    """
    items = module.split(".")
    module = ".".join(items[:-1])
    file_name = items[-1]+".yaml"
    
    defaults_file = get_data_file(module, file_name, subdir="defaults")
    
    with open(defaults_file, "r") as inf:
        defaults_dict = yaml.safe_load(inf)
    
    defaults = _extract_dict(defaults_dict, {}, args)
    
    return defaults


def build_expr(value):
    """
    Build an expression for ``simple_eval()`` out of a dictionary that's basically in the
    form of a tree boolean expression.
    """
    wildcard_token = "*"
    source = value.get("source", "")
    reason = value.get("reason", "")
    
    if value["type"].upper() == "NODE":
        # Leaf node. Make the check
        search_column = value["column"]
        search_key = value["key"]
        if "special" in value:
            # We're doing something odd. Currently set up for "length:l" which means 
            # we only care about the first l characters, and "contains", which means
            # we're checking if the value contains a supplied string.
            if "length" in value["special"]:
                l = int(value["special"].split(":")[1])
                expr = f"{search_column}[:{l}] == '{search_key}'"
            elif value["special"] == "contains":
                expr = f"'{search_key}' in {search_column}"
            else:
                logger.error(f"ERROR: Unknown special column type {inputs['special']}")
        else:
            if wildcard_token in search_key:
                expr = f"{search_column} in '{search_key.replace(wildcard_token, "")}'"
            else:
                expr = f"{search_column} {search_key}"
    elif value["type"].upper() == "NOT":
        # Special case -- there will only be a single entry.
        e, s, r = build_expr(value["entries"][0])
        expr = f"not ({e})"
        if source != "":
            source = f"{source}. {s}"
        else:
            source = s
        if reason != "":
            reason = f"{reason}. {r}"
        else:
            reason = r
    elif value["type"].upper() == "AND":
        # AND -- combine all of the sub-elements
        e, s, r = build_expr(value["entries"][0])
        expr = f"({e})"
        if source != "":
            source = f"{source}. {s}"
        else:
            source = s
        if reason != "":
            reason = f"{reason}. {r}"
        else:
            reason = r
        for entry in value["entries"][1:]:
            e, s, r = build_expr(entry)
            expr = f"{expr} and ({e})"
            if source != "":
                source = f"{source}. {s}"
            else:
                source = s
            if reason != "":
                reason = f"{reason}. {r}"
            else:
                reason = r
    elif value["type"].upper() == "OR":
        # OR -- combine all of the sub-elements
        e, s, r = build_expr(value["entries"][0])
        expr = f"({e})"
        if source != "":
            source = f"{source}. {s}"
        else:
            source = s
        if reason != "":
            reason = f"{reason}. {r}"
        else:
            reason = r
        for entry in value["entries"][1:]:
            e, s, r = build_expr(entry)
            expr = f"{expr} or ({e})"
            if source != "":
                source = f"{source}. {s}"
            else:
                source = s
            if reason != "":
                reason = f"{reason. {r}"
            else:
                reason = r
    
    return expr, source, reason


def value_eval(var, expr, pre="", verbose=False):
    """
    Return the output when evaluating the variable ``var`` as part of an arithmetic 
    string. This uses ``simple_eval()``, which is intended to allow only arithmetic 
    expressions that won't cause performance bottlenecks or system crashes, as a 
    replacement for the (much more general and dangerous) ``eval()`` function.
    
    In general, var will be either a string, a number (integer or floating-point), or a 
    boolean value, and numbers and boolean values can be treated the same way (strings 
    must have quotation marks added around ``var``'s value)
    """
    wildcard_char = '*'
    negation_char = '!'
    
    if isinstance(var, str):
        if negation_char in expr:
            expr = expr.replace(negation_char, "")
            if wildcard_char in expr:
                expr = expr.replace(wildcard_char, "")
                final = f"not ('{expr}' in '{var}')"
                logger.debug(f"{pre}: evaluating '{final}'")
                return simple_eval(final)
            else:
                final = f"not ('{var}' {expr)"
                logger.debug(f"{pre}: evaluating '{final}'")
                return simple_eval(final)
        elif wildcard_char in expr:
            expr = expr.replace(wildcard_char, "")
            final = f"'{expr}' in '{var}'"
            logger.debug(f"{pre}: evaluating '{final}'")
            return simple_eval(final)
        final = f"'{var}' {expr}"
        logger.debug(f"{pre}: evaluating '{final}'")
        return simple_eval(final)
    final = f"{var} {expr}"
    logger.debug(f"{pre}: evaluating '{final}'")
    return simple_eval(final)


def initialize_value(name, inputs, row, default, verbose):
    """
    Initialize a settable parameter, by
    - Getting an initial value based on a supplied default and a check of the data files
    - Setting the to-be-used value to the initial value
    - Returning both
    """
    initial_value = default
    found, value = get_value(name, inputs, row, default=default, verbose=verbose)
    if found:
        return value
    return initial_value


def get_value(name, inputs, row, default=None, pre="", verbose=False):
    """
    Stub function because the base-level needs to check for the presence of the keyword,
    whereas the recursive function does not necessarily need to.
    """
    if name in inputs:
        return find_value(name, inputs, row, default, pre, verbose)
    return False, default


def find_value(name, inputs, row, default=None, pre="", verbose=False):
    """
    Find a parameter from a nested dictionary, given an optional default value, and
    return the most applicable value for that parameter.
    """
    value = default
    match_found = False
    
    if name in inputs:
        logger.debug(f"{pre}Found {name}")
        inputs = inputs[name]
        if 'parameters' in inputs:
            # Trunk node
            logger.debug(f"{pre}Parameter Trunk")
            if 'value' in inputs:
                if (name in row.columns) and (value_eval(row[name], inputs['value'], pre, verbose)):
                    logger.debug(f"{pre}Matched {name} ({row[name]})")
                    match_found = True
            elif 'values' in inputs:
                if (name in row.columns) and (row[name] in inputs['values']):
                    logger.debug(f"{pre}Matched {name} ({row[name]})")
                    match_found = True
            else:
                match_found = True
            if match_found:
                if 'default' in inputs:
                    value = inputs['default']
                    logger.debug(f"{pre}: Setting value to default {value}")
                ordering = inputs['eval_order']
                inputs = inputs['parameters']
                for key in ordering:
                    logger.debug(f"{pre}Checking {key}")
                    local_match_found, value_found = find_value(key, inputs, row, default=value, pre=f"\t{pre}", verbose=verbose)
                    if local_match_found:
                        match_found = True
                        value = value_found
                        break
        else:
            if isinstance(inputs, list):
                for item in inputs:
                    logger.debug(f"{pre}Checking list item")
                    m, v = find_value(name, item, row, default=value, pre=f"\t{pre}", verbose=verbose)
                    if m:
                        return m, v
            elif "default" in inputs:
                source = inputs.get("source", "[NONE PROVIDED]")
                reason = inputs.get("reason", "[NONE PROVIDED]")
                msg = f"{pre}: Using default {inputs['default']} from {source} because "
                msg += f"{reason}"
                logger.debug(msg)
                return True, inputs["default"]
            else:
                logger.warning(f"{pre}: ERROR: Invalid node {inputs}")
    elif 'parameters' in inputs:
        # Trunk node
        logger.debug(f"{pre}Parameter Trunk")
        if 'value' in inputs:
            if (name in row.columns) and (value_eval(row[name], inputs['value'], pre, verbose)):
                logger.debug(f"{pre}Matched {name} ({row[name]})")
                match_found = True
        elif 'values' in inputs:
            if (name in row.columns) and (row[name] in inputs['values']):
                logger.debug(f"{pre}Matched {name} ({row[name]})")
                match_found = True
        else:
            match_found = True
        if match_found:
            if 'default' in inputs:
                value = inputs['default']
                logger.debug(f"{pre}: Setting value to default {value}")
            ordering = inputs['eval_order']
            inputs = inputs['parameters']
            for key in ordering:
                logger.debug(f"{pre}Checking {key}")
                local_match_found, value_found = find_value(key, inputs, row, default=value, pre=f"\t{pre}", verbose=verbose)
                if local_match_found:
                    match_found = True
                    value = value_found
                    break
    elif 'value' in inputs:
        # Leaf node, single-expression
        logger.debug(f"{pre}Leaf node {inputs}")
        if (name in row.columns) and (value_eval(row[name], inputs['value'], pre, verbose)):
            match_found = True
            if "user" in inputs:
                # user override
                name = inputs['user']['name']
                date = inputs['user']['date']
                result = inputs['user']['param_value']
                reason = inputs['user']['reason']
                msg = f"{pre}User Override: {name} on {date} set value to {result}"
                msg += f" because {reason}"
                logger.debug(msg)
                value = inputs['user']['param_value']
            else:
                value = inputs['param_value']
                msg = f"{pre}Setting {name} to {value} based on {inputs['source']}"
                msg += f" because {inputs['reason']}"
                logger.debug(msg)
    elif 'values' in inputs:
        # Leaf node, multi-value list
        logger.debug(f"{pre}Leaf node {inputs}")
        if (name in row.columns) and (row[name] in inputs['values']):
            match_found = True
            if "user" in inputs:
                # user override
                name = inputs['user']['name']
                date = inputs['user']['date']
                result = inputs['user']['param_value']
                reason = inputs['user']['reason']
                msg = f"{pre}User Override: {name} on {date} set value to {result}"
                msg += f" because {reason}"
                logger.debug(msg)
                value = inputs['user']['param_value']
            else:
                value = inputs['param_value']
                msg = f"{pre}Setting {name} to {value} based on {inputs['source']}"
                msg += f" because {inputs['reason']}"
                logger.debug(msg)
    elif "default" in inputs:
        value = default
        match_found = True
        source = inputs.get("source", "[NONE PROVIDED]")
        reason = inputs.get("reason", "[NONE PROVIDED]")
        logger.debug(f"{pre}: Using default {value} from {source} because {reason}")
    else:
        msg = f"{pre}: ERROR: Invalid node {inputs} with no values or parameters."
        msg += f" Returning {value}"
        logger.error(msg)
    
    logger.debug(f"{pre}Returning {value} ({match_found})")
    return match_found, value


def set_override(name, inputs, row, uname, value, reason):
    """
    Set a user override on a parameter
    """
    
    user_dict = {'name': uname,
                 'date': datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                 'param_value': value,
                 'reason': reason}
    
    if name not in inputs:
        inputs[name] = {'eval_order': ['root'], 'parameters': {'root': []}}
    p = inputs[name]['parameters']
    if 'root' not in p:
        p['root'] = []
    leaf_found = False
    for item in p['root']:
        if 'value' in item:
            # Leaf node, single-expression
            if value_eval(row['root'], item['value'], verbose=True):
                item['user'] = user_dict
                return
        elif 'values' in item:
            # Leaf node, multi-value list. Since we're overriding this-one-only, REMOVE
            # from list.
            if row['root'] in item['values']:
                item['values'].remove(row['root'])
    item_entry = {f'value': "== '{row[\'root\']}'", 'user': user_dict}
    p['root'].append(item_entry)
    
    return inputs


def handle_parameter(name, start_value, current_value, file, row):
    """
    Check whether the user has changed a parameter. If they have, ask whether they want to
    update the reference file. If they do, update the file.
    """
    if start_value != current_value:
        done = False;
        while not done:
            logger.info(f"You changed {name} from {start_value} => {current_value}.")
            response = input("Do you want to save this change to the data file? (y/N)")
            if response.lower() in ["y", "yes"]:
                rname = input("Please enter a name for the change: ")
                reason = input("Please enter a reason for the change: ")
                yaml = YAML()
                with open(file) as inf:
                    config = yaml.load(inf)
                set_override(name, config, row, rname, current_value, reason)
                with open(file, mode="w") as outf:
                    yaml.dump(config, outf)
                logger.info(f"Changed value for {row['root']}")
                done = True
            elif response.lower() in ["n", "no"]:
                logger.info("Value not saved.")
                done = True
            else:
                logger.info(f"Unknown response {response}")


def set_param(param, default, row, issues, pre, overrides={}, verbose=False):
    """
    Set a parameter value
    
    Given a parameter name, that parameter's default value, a data table
    row, and a JSON dictionary which may have an entry for the current row that
    will override the parameter, return the parameter value that should be used.

    Parameters
    ----------
    param : str
        The parameter to check and return
    default : object
        The default value for that parameter
    row : abscal.common.exposure_data_table.AbscalDataTable
        A single-row table containing the data of interest.
    issues : dict
        A dictionary containing a set of parameters (one of which may be
        param), along with information to identify files whose parameters
        should be adjusted .
    overrides : dict
        A dictionary containing any parameters whose value is being overridden
        by the user.
    verbose : bool
        Whether or not informational output should be printed.

    Returns
    -------
    value : object
        The appropriate value for the parameter given
    """
    value = default
    try:
        if param in issues:
            issue_list = issues[param]
            for item in issue_list:
                comp_value = row[item["column"]]
                if isinstance(item["value"], str):
                    comp_value = comp_value[:len(item["value"])]
                if comp_value == item["value"]:
                    reason = item["reason"]
                    source = item["source"]
                    msg = f"{pre}: changed {param} {value}=>{item['param_value]} because"
                    msg += f" {reason} from {source}"
                    logger.debug(msg)
                    value = item["param_value"]
        if param in overrides:
            msg = f"{pre}: changed {param} {value}=>{overrides[param]} by user override"
            logger.debug(msg)
            value = overrides[param]
    except Exception as e:
        logger.error("ERROR in set_param(). Arguments are:")
        logger.error(f"\tparam={param}")
        logger.error(f"\tdefault={default}")
        logger.error("\trow:\n")
        logger.error(row)
        logger.error("\n")
        logger.error("\tissues:\n")
        logger.error(issues)
        logger.error(f"\tpre={pre}")
        logger.error("\toverrides:\n")
        logger.error(overrides)
        logger.error(f"\tverbose={verbose}")
        raise e

    return value


def set_params(defaults, row, issues, pre, overrides={}, verbose=False):
    """
    Set multiple parameter values
    
    Given a dictionary of default values, a metadata row, a dictionary of
    known issues and overrides, a dictionary of user-supplied overrides,
    and a verbose flag, produce a dictionary of parameters (all with the
    appropriate value) which also contains a 'set' key (an array of all
    the parameters that have been overridden from their default values).

    Parameters
    ----------
    defaults : dict
        A dictionary of default values (also names the parameters)
    row : abscal.common.exposure_data_table.AbscalDataTable
        A single-row table containing the data of interest.
    issues : dict
        A dictionary containing a set of parameters, along with information to
        identify files whose parameters should be adjusted.
    overrides : dict
        A dictionary containing any parameters whose value is being overridden
        by the user.
    verbose : bool
        Whether or not informational output should be printed.

    Returns
    -------
    params : dict
        The supplied parameters, each with its value.
    """
    params = {'set': []}

    for param in defaults:
        default = defaults[param]
        value = set_param(param, default, row, issues, pre, overrides, verbose)
        params[param] = value
        if value != default:
            params['set'].append(param)

    return params

def set_image(images, row, issues, pre, overrides={}, verbose=False):
    """
    Update an image based on known issues.
    
    Given an image, image metadata, and a set of known issues, determine if any
    of the known issues apply to the image in question and, if they do, make
    the appropriate edits to the image.

    Parameters
    ----------
    images : dict
        Dict of (SCI, ERR, DQ) np.ndarray images
    row : abscal.common.exposure_data_table.AbscalDataTable
        A single-row table containing metadata on the image
    issues : dict
        A dictionary containing a set of parameters, along with information to
        identify files whose parameters should be adjusted.
    overrides : dict
        A dictionary containing any parameters whose value is being overridden
        by the user.
    verbose : bool
        Whether or not informational output should be printed.

    Returns
    -------
    image : tuple
        Tuple of (SCI, ERR, DQ) np.ndarray images, as edited.
    """

    for issue in issues:
        found = False
        if issue["column"] in row:
            if isinstance(issue["column"], str):
                issue_len = len(issue["value"])
                if issue["value"] == row[issue["column"]][:issue_len]:
                    found = True
            else:
                if issue["value"] == row[issue["column"]]:
                    found = True
        if found:
            if len(issue["x"]) > 1:
                x1, x2 = issue["x"][0], issue["x"][1]
            else:
                x1, x2 = issue["x"][0], issue["x"][0]+1
            if len(issue["y"]) > 1:
                y1, y2 = issue["y"][0], issue["y"][1]
            else:
                y1, y2 = issue["y"][0], issue["y"][0]+1
            images[issue["ext"]][y1:y2,x1:x2] = issue["value"]
            reason = issue["reason"]
            source = issue["source"]
            value = issue["value"]
            msg = f"{pre}: changed ({y1}:{y2},{x1}:{x2}) to {value} because {reason}"
            msg += f" from {source}"
            logger.debug(msg)


    return images


def air2vac(air):
    """
    Convert a set of wavelengths from air to vacuum.

    Parameters
    ----------
    air : array-like
        Air wavelengths

    Returns
    -------
    vac : array-like
        Vacuum wavelengths
    """
    vac = []
    c0 = 0.00008336624212083
    c1 = 0.02408926869968
    c2 = 0.0001599740894897

    for wl_air in air:
        s = 1.e4/wl_air
        n = 1 + c0 + c1 / (130.1065924522 - s*s) + c2 / (38.92568793293 - s*s)
        wl_vac = wl_air * n
        vac.append(wl_vac)

    return np.array(vac)


def smooth_model(wave, flux, fwhm):
    """
    Smooth a model spectrum with a non-uniform sampling interval. 
    
    Based on Ralph Bohlin's "smomod.pro", which itself references "tin.pro"

    Parameters
    ----------
    wave : array-like
        Wavelength array

    flux : array-like
        Flux array
    fwhm : float
        FWHM of delta function.

    Returns
    -------
    smoothed : array-like
        Smoothed flux array.
    """
    wmin = wave - fwhm/2.
    wmax = wave + fwhm/2.

    good = np.where((wmin > wave[0]) & (wmax < wave[-1]))
    smoothed = deepcopy(flux)
    smoothed[good] = trapezoidal(wave, flux, wmin[good], wmax[good])
    smoothed[good] = trapezoidal(wave, smoothed, wmin[good], wmax[good])

    return smoothed

def trapezoidal(wave, flux, wmin, wmax):
    """
    Make a trapezoidal integral
    
    Trapezoidal 'integral' (really an average) from Ralph Bohlin's 'tin.pro'
    and 'integral.pro'. Uses wmin and wmax to set limits

    Parameters
    ----------
    wave : array-like
        Wavelength array
    flux : array-like
        Flux array
    wmin : array-like
        Wavelength array shifted bluewards by FWHM/2
    wmax : array-like
        Wavelength array shifted redwards by FWHM/2

    Returns
    -------
    trapint : array-like
        Flux array after trapezoidal integral
    """
    trapint = integral(wave, flux, wmin, wmax)/(wmax - wmin)
    return trapint

def integral(x, y, xmin, xmax):
    """
    Return the approximate integral of y over x for the range (xmin, xmax)
    
    Parameters
    ----------
    x : array-like
        X value array
    y : array-like
        Y value array. Must be the same length as x
    xmin : float
        minimum x value for integral
    xmax : float
        maximum x value for integral
        
    Returns
    -------
    int : float
        integral value
    """
    rmin = tabinv(x, xmin)
    rmax = tabinv(x, xmax)
    n = len(x)
    dx = np.roll(x, -1) - x
    if (np.max(xmin) > np.max(xmax)) or (np.min(xmax) < np.min(xmin)) or (np.min(xmax - xmin) < 0.):
        raise ValueError("Invalid integral limits")
    dx = np.roll(x, -1) - x
    dy = np.roll(y, -1) - y

    dint = (np.roll(y, -1) + y)/(2.*dx)
    imin = np.floor(rmin).astype(np.int32)
    imax = np.floor(rmax).astype(np.int32)

    dxmin = xmin - x[imin]
    ymin = y[imin] + dxmin*(y[imin+1]-y[imin])/(dx[imin] + np.where(dx[imin]==0, 1, 0))
    dxmax = xmax - x[imax]
    coeff0 = y[np.where(imax+1<len(y)-1, imax+1, len(y)-1)] - y[imax]
    coeff1 = dx[imax] + np.where(dx[imax]==0, 1, 0)
    ymax = y[imax] + dxmax*coeff0/coeff1

    int = np.zeros_like(xmin)
    for i in range(len(xmin)):
        if imax[i] != imin[i]:
            int[i] = np.sum(dint[imin[i]:imax[i]]-1)

    int -= (y[imin] + ymin)/(2.*dxmin)
    int += (y[imax] + ymax)/(2.*dxmax)

    return int


def tabinv(xarr, x):
    """
    Find the effective index in xarr of each element in x.
    
    The effective index for each element j in x is the value i such that 
    :math:`xarr[i] <= x[j] <= xarr[i+1]`, to which is added an interpolation fraction 
    based on the size of the intervals in xarr.
    
    Parameters
    ----------
    x_arr : array-like
        The array of values to search
    x : float or array-like
        Value (or list of values) to look for in x_arr
    
    Returns
    -------
    ieff : float
        Effective index
    """
    npoints, npt = len(xarr), len(xarr) - 1
    if npoints <= 1:
        raise ValueError("Search array must contain at least 2 elements")

    if not (np.all(np.diff(xarr) >= 0) or (np.all(np.diff(xarr) <= 0))):
        raise ValueError("Search array must be monotonic")
    
    if not isinstance(x, (list, tuple, np.ndarray)):
        x = np.array([x])
    
    # ieff contains values j1, ..., jn such that
    #   ji = x where xarr[x-1] <= ji < xarr[x]
    #   If no position is found, ji = len(xarr)
    ieff = np.searchsorted(xarr, x, side='right').astype(np.float64)
    g = np.where((ieff >= 0) & (ieff < (len(xarr) - 1)))
    if len(g) > 0 and len(g[0] > 0):
        neff = ieff[g].astype(np.int32)
        x0 = xarr[neff].astype(np.float64)
        diff = x[g] - x0
        ieff[g] = neff + diff / (xarr[neff+1] - x0)
    ieff = np.where(ieff>0., ieff, 0.)
    return ieff


def linecen(wave, spec, cont):
    """
    Find the centre of an emission line.
    
    Computes the centroid of an emission line over the range of
    
    :math:
        xapprox \pm fwhm/2

    after subtracting any continuum and half value at the remaining peak. After
    clipping at zero, the weights of the remaining spectral wings approach zero,
    so any marginally missed or included point matters little.

    Parameters
    ----------
    wave : np.ndarray
        1-d array of x values
    spec : np.ndarray
        1-d array of y values
    cont : float
        Approximate continuum value

    Returns
    -------
    centroid : float
        The x value of the centroid
    badflag : bool
        False for good data, true for bad data
    """
    badflag = False
    profile = spec - cont
    clip = (profile - np.max(profile)*0.5)
    n_points = len(clip)
    midpoint = (n_points-1)//2
    low_bad = 0
    while low_bad < len(clip) and clip[low_bad] < 0:
        low_bad += 1
    if low_bad > 0:
        clip[:low_bad+1] = 0.
    high_bad = len(clip) - 1
    while high_bad >= 0 and clip[high_bad] < 0:
        high_bad -= 1
    if high_bad < len(clip) - 1:
        clip[high_bad:] = 0.
    clip = np.where(clip>0., clip, 0.)
    good = np.where(clip > 0.)
    if len(good) > 0:
        n_good = len(good[0])
    if n_good <= 1:
        logger.warning("WARNING: LINECEN: Bad profile, centroid set to midpoint.")
        return wave[midpoint], "bad"
    centroid = np.sum(wave*clip)/np.sum(clip)
    return centroid, "good"


def return_type(var_type, var_value):
    """
    Handles updating a variable given a type string. Currently this depends on the type
    being one of "float", "int", "bool", or "string".
    
    Parameters
    ----------
    var_type : str
        Variable type indicator
    var_value : str
        Variable value, formatted as a string
    
    Returns
    -------
    var_value : object
        The variable value, in the specified type
    """
    if var_type == "float":
        return float(var_value)
    elif var_type == "int":
        return int(var_value)
    elif var_type == "bool":
        if isinstance(var_value, str):
            return strtobool(var_value)
        elif not isinstance(var_value, bool):
            return bool(var_value)
    return var_value


def check_params(task, module, params, values, verbose):
    """
    Given a task and module, check the current values (in values) against the existing
    values (in params), update params as needed, and return a value indicating whether 
    any have changed.
    
    Parameters
    ----------
    task : str
        Task name
    module : str
        Abscal module
    params : dict
        Parameter values at the start of the iteration
    values : dict
        Current parameter values
    verbose : bool
        Whether to print verbose output
    
    Returns
    -------
    changed : bool
        Whether any settings where changed
    """
    param_types = get_param_types(task, module, verbose)

    changed = False
    for item in params:
        item_type = param_types[item].split("_")[0]
        if "none" in param_types[item]:
            gen_type = param_types[item].split("_")[0]
            if params[item] is None and values[item] == "None":
                continue
            elif params[item] is None:
                changed = True
                params[item] = return_type(item_type, values[item])
            elif values[item] == 'None':
                changed = True
                params[item] = None
            else:
                current_value = return_type(item_type, values[item])
                if params[item] != current_value:
                    changed = True
                    params[item] = current_value
        else:
            current_value = return_type(item_type, values[item])
            if params[item] != current_value:
                changed = True
                params[item] = current_value
        
    return changed


def get_param_types(task, module, verbose):
    """
    Given a task name, create a dictionary that has type info for every parameter in that
    task.
    
    Parameters
    ----------
    task : str
        Task name
    module : str
        Abscal module
    verbose : bool
        Whether to print verbose output
    
    Returns
    -------
    types : dict
        The task types dictionary
    """
    param_types = {}
    task_params = {}
    task_file = get_data_file(module, "tasks.yaml")
    if task_file is not None:
        with open(task_file) as in_file:
            task_dict = yaml.safe_load(in_file)
        if task in task_dict:
            task_params = task_dict[task]
    
    for item in task_params:
        param_types[item] = task_params[item]["type"]

    return param_types


def get_task_defaults(task, module, verbose=False):
    """
    Given a task name, create a dictionary that holds the default parameter set for
    running that task, as well as type and history information.
    
    Parameters
    ----------
    task : str
        Task name
    module : str
        Abscal module
    verbose : bool, default False
        Whether to print verbose output
    
    Returns
    -------
    task_params : dict
        The set parameter dictionary
    """
    task_params = {}
    task_file = get_data_file(module, "tasks.yaml")
    if task_file is not None:
        with open(task_file) as in_file:
            task_dict = yaml.safe_load(in_file)
        if task in task_dict:
            task_params = task_dict[task]

    return task_params


def get_default_params(task, module, verbose=False):
    """
    Given a task name, create a dictionary that holds the default parameter values for
    that task.

    Parameters
    ----------
    task : str
        Task name
    module : str
        Abscal module
    verbose : bool, default False
        Whether to print verbose output
    
    Returns
    -------
    params : dict
        The set parameter dictionary
    """
    params = {}
    param_types = get_param_types(task, module, verbose)
    task_params = get_task_defaults(task, module, verbose=verbose)
    
    for item in task_params:
        default = task_params[item]["value"]
        if "none" in task_params[item]["type"] and default == "None":
            default = None
        if "string" in param_types[item]:
            default = f"{default}"
        params[item] = default

    return params    


def setup_params(task, module, metadata, row, verbose=False):
    """
    Given a task name, create a dictionary that uses the current exposure of interest to 
    populate a dictionary of parameter values, then return that dictionary
    
    Parameters
    ----------
    task : str
        Task name
    module : str
        Abscal module
    metadata : dict-like
        Data file to search for exposure-specific settings
    row : ABSCAL.ExposureDataTable row
        Row of information on the current exposure
    verbose : bool, default False
        Whether to print verbose output
    
    Returns
    -------
    params : dict
        The set parameter dictionary
    """
    params = {}
    param_types = get_param_types(task, module, verbose)
    task_params = get_task_defaults(task, module, verbose=verbose)
    
    for item in task_params:
        default = task_params[item]["value"]
        if "none" in task_params[item]["type"] and default == "None":
            default = None
        params[item] = initialize_value(item, metadata, row, default, verbose)
        if "string" in param_types[item]:
            params[item] = f"{params[item]}"

    return params

def update_param(task, module, row, settings_file, param, value, name, notes, update_type):
    """
    Update a parameter for a task. The update_type can be "default", to update the overall
    default, "current" to update the task version for this particular file, or "revert" to
    do nothing.
    """
    if update_type == "current":
        yaml = YAML()
        with open(settings_file) as inf:
            config = yaml.load(inf)
        set_override(param, config, row, name, value, notes)
        with open(file, mode="w") as outf:
            yaml.dump(config, outf)
    elif update_type == "default":
        task_file = get_data_file(module, "tasks.yaml")
        if task_file is not None:
            yaml = YAML()
            with open(task_file) as inf:
                config = yaml.load(inf)
            update_dict = {
                "value": deepcopy(config[param]["value"]),
                "date": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                "name": name,
                "reason": notes
            }
            update_history = [update_dict] + deepcopy(config[param]["history"])
            config["value"] = value
            config["name"] = name
            config["reason"] = notes
            with open(file, mode="w") as outf:
                yaml.dump(config, outf)
