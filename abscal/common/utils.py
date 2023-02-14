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
from copy import deepcopy
from datetime import datetime
from ruamel.yaml import YAML
from simpleeval import simple_eval


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


def get_base_data_dir():
    """
    Find the location of the ABSCAL data files.
    
    ABSCAL stores both default parameters and exposure-specific corrections and settings 
    in a number of data files. There are default copies stored internally, but they can 
    also be stored elsewhere. The ABSCAL_DATA environment variable points to that 
    location, although there is always a fallback to local data if a specified file does 
    not exist elsewhere.
    
    Returns
    -------
    data_dir : str
        The directory pointed to by ABSCAL_DATA (if both the environment variable and 
        the directory to which it points exist)
    """
    if "ABSCAL_DATA" in os.environ:
        # First look for all-capitals
        if os.path.isdir(os.environ["ABSCAL_DATA"]):
            return os.environ["ABSCAL_DATA"]
    elif "abscal_data" in os.environ:
        # Next look for all lower-case
        if os.path.isdir(os.environ["abscal_data"]):
            return os.environ["abscal_data"]
    
    # Fall back to internal
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get_data_file(module, fname, defaults=False):
    """
    Find an internal data file.
    
    Returns the path to a file named `fname` in the data directory of the
    (sub)module named `module`. Keeps there from being a lot of repetitive code in order 
    to find data file paths.

    Parameters
    ----------
    module : str
        The module to search in, using standard dot separators (e.g.
        abscal.wfc3)
    fname : str
        The file name of interest
    defaults : bool, default False
        Whether to append a "defaults" directory to the final path

    Returns
    -------
    data_file : str or None
        Full path to the data file. If no file is found at the generated path, None will 
        be returned. This is not necessarily a failure state, because (for example) a 
        function may call for a known-issues file even though there are no current known 
        issues (and thus no file to return). This way, the code doesn't need to be 
        changed when there *is* a file, and a user can add a file to their local directory 
        without needing to alter the code, because the code will just transparently find 
        the file.
    """
    # /path/to/abscal (with /common/utils.py stripped off)
    local_loc = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    current_loc = get_base_data_dir()
        
    module = module.replace("abscal.", "")
    
    # Replace '.' with path separator
    module_path = module.replace(".", "/")

    data_path = os.path.join(current_loc, module_path, "data")
    if defaults:
        data_path = os.path.join(data_path, "defaults")
    data_file = os.path.join(data_path, fname)
    
    if os.path.isfile(data_file):
        # Try for the data file (with potential user-supplied path)
        return data_file
    elif os.path.isfile(data_file.replace(current_loc, local_loc)):
        # Fall back to the local version
        return data_file.replace(current_loc, local_loc)

    # This is a convenience. The exposure parameters file associated with a particular
    # python file is that file's name with ".py" replaced with ".yaml". For example, the 
    # util_grism_cross_correlate.py exposure parameters file would be 
    # util_grism_cross_correlate.yaml.
    #
    # As such, this basically says that, if you're asking for a python file (e.g. by just 
    # sending in __file__ to get_data_file) (recognized by the extension being ".py"), 
    # then search for a ".yaml" file with the same name if the previous search failed.
    if os.path.splitext(fname)[1] == ".py":
        fname = os.path.splitext(fname)[0] + ".yaml"
        data_file = os.path.join(data_path, fname)
    
        if os.path.isfile(data_file):
            # Try for the data file (with potential user-supplied path)
            return data_file
        elif os.path.isfile(data_file.replace(current_loc, local_loc)):
            # Fall back to the local version
            return data_file.replace(current_loc, local_loc)
    
    msg = "ERROR: File {}.{} not found at {} or {}.\n"
    msg = msg.format(module, fname, data_file, data_file.replace(current_loc, local_loc))
    sys.stderr.write(msg)
#     for i in range(5):
#         data_file = os.path.dirname(data_file)
#         search_str = os.path.join(data_file, "*")
#         try:
#             names = [os.path.basename(x) for x in glob.glob(search_str)]
#             msg += "Directory {} contains: {}\n".format(data_file, names)
#         except Exception as e:
#             break
#     raise FileNotFoundError(msg)
    
    # If nothing was found, return None
    return None
    
    
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
    
    defaults_file = get_data_file(module, file_name, defaults=True)
    
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
            if "length" in inputs["special"]:
                l = int(inputs["special"].split(":")[1])
                expr = "{}[:{}] == '{}'".format(search_column, l, search_key)
            elif inputs["special"] == "contains":
                expr = "'{}' in {}".format(search_key, search_column)
            else:
                print("ERROR: Unknown special column type {}".format(inputs["special"]))
        else:
            if wildcard_token in search_key:
                expr = "{} in '{}'".format(search_column, search_key.replace(wildcard_token, ""))
            else:
                expr = "{} {}".format(search_column, search_key)
    elif value["type"].upper() == "NOT":
        # Special case -- there will only be a single entry.
        e, s, r = build_expr(value["entries"][0])
        expr = "not ({})".format(e)
        if source != "":
            source = "{}. {}".format(source, s)
        else:
            source = s
        if reason != "":
            reason = "{}. {}".format(reason, r)
        else:
            reason = r
    elif value["type"].upper() == "AND":
        # AND -- combine all of the sub-elements
        e, s, r = build_expr(value["entries"][0])
        expr = "({})".format(e)
        if source != "":
            source = "{}. {}".format(source, s)
        else:
            source = s
        if reason != "":
            reason = "{}. {}".format(reason, r)
        else:
            reason = r
        for entry in value["entries"][1:]:
            e, s, r = build_expr(entry)
            expr = "{} and ({})".format(expr, e)
            if source != "":
                source = "{}. {}".format(source, s)
            else:
                source = s
            if reason != "":
                reason = "{}. {}".format(reason, r)
            else:
                reason = r
    elif value["type"].upper() == "OR":
        # OR -- combine all of the sub-elements
        e, s, r = build_expr(value["entries"][0])
        expr = "({})".format(e)
        if source != "":
            source = "{}. {}".format(source, s)
        else:
            source = s
        if reason != "":
            reason = "{}. {}".format(reason, r)
        else:
            reason = r
        for entry in value["entries"][1:]:
            e, s, r = build_expr(entry)
            expr = "{} or ({})".format(expr, e)
            if source != "":
                source = "{}. {}".format(source, s)
            else:
                source = s
            if reason != "":
                reason = "{}. {}".format(reason, r)
            else:
                reason = r
    
    return expr, source, reason


def value_eval(var, expr):
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
                return simple_eval("not ({} in '{}')".format(expr, var))
            else:
                return simple_eval("not ('{}' {})".format(var, expr))
        elif wildcard_char in expr:
            expr = expr.replace(wildcard_char, "")
            return simple_eval("{} in '{}'".format(expr, var))
        return simple_eval("'{}' {}".format(var, expr))
    return simple_eval("{} {}".format(var, expr))


def initialize_value(name, inputs, row, default, verbose):
    """
    Initialize a settable parameter, by
    - Getting an initial value based on a supplied default and a check of the data files
    - Setting the to-be-used value to the initial value
    - Returning both
    """
    initial_value = default
    found, value = find_value(name, inputs, row, default=default, verbose=verbose)
    if found:
        initial_value = value
    value = initial_value
    return initial_value, value


def find_value(name, inputs, row, default=None, pre="", verbose=False):
    """
    Find a parameter from a nested dictionary, given an optional default value, and
    return the most applicable value for that parameter.
    """
    value = default
    match_found = False
    
    if name in inputs:
        if verbose:
            print("{}Found {}".format(pre, name))
        inputs = inputs[name]
        if 'parameters' in inputs:
            # Trunk node
            if 'value' in inputs:
                if (name in row) and (value_eval(row[name], inputs['value'])):
                    if verbose:
                        print("{}Matched {} ({})".format(pre, name, row[name]))
                    match_found = True
            elif 'values' in inputs:
                if (name in row) and (row[name] in inputs['values']):
                    if verbose:
                        print("{}Matched {} ({})".format(pre, name, row[name]))
                    match_found = True
            else:
                match_found = True
            if match_found:
                if 'default' in inputs:
                    value = inputs['default']
                ordering = inputs['eval_order']
                inputs = inputs['parameters']
                for key in ordering:
                    if verbose:
                        print("{}Checking {}".format(pre, key))
                    local_match_found, value_found = find_value(key, inputs, row, default=value, pre="\t{}".format(pre), verbose=verbose)
                    if local_match_found:
                        match_found = True
                        value = value_found
                        break
        elif 'user' in inputs:
            # leaf note, user override
            if verbose:
                name = inputs['user']['name']
                date = inputs['user']['date']
                result = inputs['user']['param_value']
                reason = inputs['user']['reason']
                msg = "{}User Override: {} on {} set value to {} because {}"
                print(msg.format(pre, name, date, result, reason))
            value = inputs['user']['param_value']
            match_found = True
        elif 'value' in inputs:
            # Leaf node, single-expression
            if (name in row) and (value_eval(row[name], inputs['value'])):
                value = inputs['param_value']
                match_found = True
                if verbose:
                    msg = "{}Setting {} to {} based on {} because {}"
                    print(msg.format(pre, name, value, inputs['source'], inputs['reason']))
        elif 'values' in inputs:
            # Leaf node, multi-value list
            if (name in row) and (row[name] in inputs['values']):
                value = inputs['param_value']
                match_found = True
                if verbose:
                    msg = "{}Setting {} to {} based on {} because {}"
                    print(msg.format(pre, name, value, inputs['source'], inputs['reason']))
        else:
            if verbose:
                msg = "{}: ERROR: Invalid node {} with no values or parameters. Returning {}"
                print(msg.format(pre, inputs, value))
    
    if verbose:
        print("{}Returning {} ({})".format(pre, value, match_found))
    return match_found, value


def set_override(name, inputs, row, uname, value, reason):
    """
    Set a user override on a parameter
    """
    user_dict = {'name': uname,
                 'date': datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                 'param_value': value,
                 'reason': reason}
    
    if name not in inputs:
        inputs[name] = {'eval_order': ['root'], 'parameters': {}}
    p = inputs[name]['parameters']
    if root not in p:
        p['root'] = []
    leaf_found = False
    for item in p['root']:
        if 'value' in item:
            # Leaf node, single-expression
            if value_eval(row['root'], item['value']):
                item['user'] = user_dict
                return
        elif 'values' in item:
            # Leaf node, multi-value list. Since we're overriding this-one-only, REMOVE
            # from list.
            if row['root'] in item['values']:
                item['values'].remove(row['root'])
    item_entry = {'value': "== '{}'".format(row['root']),
                  'user': user_dict}
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
            print("You changed {} from {} => {}.".format(name, start_value, current_value))
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
                print("Changed value for {}".format(row['root']))
            elif response.lower() in ["n", "no"]:
                print("Value not saved.")
            else:
                print("Unknown response {}".format(response))


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
                    if verbose:
                        reason = item["reason"]
                        source = item["source"]
                        msg = "{}: changed {} {}=>{} because {} from {}"
                        print(msg.format(pre, param, value, item["param_value"], reason, source))
                    value = item["param_value"]
        if param in overrides:
            if verbose:
                msg = "{}: changed {} {}=> by user override"
                print(msg.format(pre, param, value, overrides[param]))
            value = overrides[param]
    except Exception as e:
        print("ERROR in set_param(). Arguments are:")
        print("\tparam={}".format(param))
        print("\tdefault={}".format(default))
        print("\trow:\n")
        print(row)
        print("\n")
        print("\tissues:\n")
        print(issues)
        print("\n")
        print("\tpre={}".format(pre))
        print("\toverrides:\n")
        print(overrides)
        print("\n")
        print("\tverbose={}".format(verbose))
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
#     print("set_image with {}, {}, {}, {}".format(images, row, issues, overrides))

    for issue in issues:
#         print(issue)
#         print(issue["column"], type(issue["column"]))
#         print(row)
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
            if verbose:
                reason = issue["reason"]
                source = issue["source"]
                value = issue["value"]
                msg = "{}: changed ({}:{},{}:{}) to {} because {} from {}"
                print(msg.format(pre, y1, y2, x1, x2, value, reason, source))


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
#     print("\tlinecen called with:")
#     print("\t\twave={}".format(wave))
#     print("\t\tspec={}".format(spec))
#     print("\t\tcont={}".format(cont))
    badflag = False
    profile = spec - cont
    clip = (profile - np.max(profile)*0.5)
#     print("\tClip starts at {}".format(clip))
    n_points = len(clip)
    midpoint = (n_points-1)//2
#     print("\tProfile has {} points, midpoint at {}.".format(n_points, midpoint))
    low_bad = 0
    while low_bad < len(clip) and clip[low_bad] < 0:
        low_bad += 1
    if low_bad > 0:
        clip[:low_bad+1] = 0.
#         print("\t\tClipping points [:{}]".format(low_bad+1))
    high_bad = len(clip) - 1
    while high_bad >= 0 and clip[high_bad] < 0:
        high_bad -= 1
    if high_bad < len(clip) - 1:
        clip[high_bad:] = 0.
#         print("\t\tClipping points [{}:]".format(high_bad))
    clip = np.where(clip>0., clip, 0.)
    good = np.where(clip > 0.)
#     print("\tClip is now {}, with good points {}".format(clip, good))
    if len(good) > 0:
        n_good = len(good[0])
    if n_good <= 1:
        print("WARNING: LINECEN: Bad profile, centroid set to midpoint.")
        return wave[midpoint], "bad"
    centroid = np.sum(wave*clip)/np.sum(clip)
    return centroid, "good"
