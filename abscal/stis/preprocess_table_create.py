#! /usr/bin/env python
"""
This module takes a file path and file type specification that matches to one
or more STIS exposures. In general, these files should be observations of
STIS standard stars. The module will open these files, retrieve information from
their headers, and produce an output table of these files. 

For compatibility with the original program, the module can be run from the 
command line directly (and given a file specification as a positional argument). 
In that case, it will work exactly the same as 'stisdir.pro' as created by Ralph
Bohlin. It also includes additional options that allow the output file name to
be specified, allow the stdout content (and verbosity) to be specified, and 
allow multiple input directories to be specified with the positional argument
acting as a template for files to search for rather than as the full search
path.

Authors
-------
- Brian York (all python code)
- Ralph Bohlin (original IDL code)

Use
---
This module is intended to be either run from the command line or used by
other module code as the first step (creating an annotated list of files 
for flux calibration)::

    python preprocess_table_create.py

Alternately there are a set of python scripts placed in your path by setup.py that will
import and call this module. The general command is::

    stis_setup

whilst the IDL version is::

    stisdir

If you use this module from within python, it is recommended to import the 
`populate_table()` function (if you need any custom table values, you can either pass in 
an AbscalDataTable directly or pass in keyword parameters that will, in turn, be passed 
along to the table creation)::

    from abscal.stis import populate_table
    
    output_table = populate_table(some_arguments=some_values)
"""

import datetime
import glob
import os
import yaml

import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.time import Time
from copy import deepcopy

from abscal.common.args import parse
from abscal.common.standard_stars import find_star_by_name, find_closest_star
from abscal.common.utils import absdate, get_data_file, get_defaults
from abscal.stis.stis_data_table import STISDataTable

def get_target_name(header):
    """
    Find the canonical target name.
    
    Return a standardized target name despite any inconsistencies in target
    naming by different PIs. For now, use the existing stisdir.pro checks to
    make a standard target name. In the future, potentially use RA and DEC to
    do a lookup to figure out the target and fill in appropriately.
    
    Parameters
    ----------
    header : astropy.io.fits header
        The header containing target information. In this header, the keys
        
        - targname (target name)
        - ra_targ (target RA)
        - dec_targ (target DEC)
        -  sclamp (active lamp)
        
        are (or may be) used by the function.
    
    Returns
    -------
    name : str
        Target name, standardized.
    """
    targname = header['TARGNAME']
    ra = header['RA_TARG']
    dec = header['DEC_TARG']
    lamp = header.get('SCLAMP', None)

    star = find_star_by_name(targname)
    if star is not None:
        return star['name']
    star = find_closest_star(ra, dec, max_distance=0.01)
    if star is not None:
        return star['name']
    if "NONE" in targname:
        if lamp is not None:
            return lamp
        return "LAMP"
    return targname


def populate_table(data_table=None, **kwargs):
    """
    Search a directory and produce a table of exposures.
    
    Uses glob to search all directories in the table's `search_dirs` array for
    files matching the table's `search_str` template, and adds rows to the table
    containing metadata based on the files that were found.
    
    Parameters
    ----------
    data_table : abscal.stis.stis_data_table.STISDataTable, default None
        The table (which may contain existing data) to which the new data should
        be added.
    kwargs : dict
        A dictionary of optional keywords. Currently checked keywords are:
        
        verbose : bool
            Flag to indicate whether or not the program should print out 
            informative text whilst running.
        compat : bool
            Whether to operate in strict IDL compatibility mode
        
        In addition to these keywords, any default parameters may be set through passing 
        a keyword argument.
        
    If data_table is None, a new table will be created in the function. In that case, the 
    kwargs dict will be passed to that table, so any table-creation keywords will be sent 
    through.
    
    Returns
    -------
    data_table : abscal.stis.stis_data_table.STISDataTable
        A table containing an entry for each input file and necessary metadata
        obtained from the FITS header of that file.
    """
    default_values = get_defaults('abscal.common.args')
    base_defaults = default_values | get_defaults(kwargs.get('module_name', __name__))
    
    if data_table is None:
        data_table = STISDataTable(**kwargs)
    
    paths = data_table.search_dirs
    file_template = data_table.search_str
    idl_strict = kwargs.get('compat', base_defaults['compat'])
    verbose = kwargs.get('verbose', base_defaults['verbose'])
    task = "create_table"

    for path in paths:
        if verbose:
            print("{}: searching {}...".format(task, path))
        all_files = glob.glob(os.path.join(path, file_template))
        
        for file_name in all_files:
            if verbose:
                print("{}: adding {}".format(task, file_name))
            loc = "START"
            file_metadata = {}

            file_path, base_name = os.path.split(file_name)
            file_ext = base_name[-8:-5]
            
            root = base_name[:9]
            file_metadata['root'] = root
            file_metadata['obset'] = base_name[:6]
            file_metadata['path'] = file_path
            file_metadata['filename'] = base_name
            file_metadata['use'] = True
            
            try:
                with fits.open(file_name) as fits_file:
                    phdr = fits_file[0].header
                    
                    cr_count = 0
                    for ext in fits_file:
                        hdr = fits_file[ext].header
                        if "EXTNAME" in hdr and hdr["EXTNAME"].strip() == 'SCI':
                            cr_count += 1
                    file_metadata['cr'] = cr_count
                    
                    file_metadata['notes'] = ''
                    
                    file_metadata['mode'] = phdr['OPT_ELEM'].upper()
                    if file_metadata['mode'].strip() == '0':
                        file_metadata['mode'] = phdr['optmode']

                    file_metadata['aperture'] = phdr['PROPAPER']
                    if phdr['PROPAPER'].strip() == '':
                        file_metadata['aperture'] = phdr['APERTURE']
                    if file_metadata['aperture'].strip() == '0':
                        file_metadata['aperture'] = phdr['SLITSIZE']

                    file_metadata['subarray'] = False
                    if phdr['subarray'] == 1:
                        file_metadata['subarray'] = True
                    file_metadata['target'] = get_target_name(phdr)
                    file_metadata['central_wavelength'] = phdr.get('CENWAVE', -1.)
                    file_metadata['min_wavelength'] = int(round(phdr['MINWAVE']))
                    file_metadata['max_wavelength'] = int(round(phdr['MAXWAVE']))
                    file_metadata['obsmode'] = phdr['OBSMODE']
                    file_metadata['instrument'] = phdr['INSTRUME']
                    file_metadata['proposal'] = phdr["proposid"]
                    file_metadata['detector'] = phdr["DETECTOR"]
                    if "TEXPTIME" not in phdr or phdr["TEXPTIME"] == 0:
                        file_metadata['exptime'] = phdr['exptime']
                    else:
                        file_metadata['exptime'] = phdr["texptime"]
                    file_metadata['gain'] = phdr.get('CCDGAIN', 0.)
                    if file_metadata['gain'] == 0.:
                        file_metadata['gain'] = phdr.get('CCDGAIN4', 0.)
                    if 'CCD' in file_metadata['detector'].upper():
                        file_metadata['detector'] = 'CCDgain{}'.format(file_metadata['gain'])
                        if phdr['CCDAMP'][:1].upper() != 'D':
                            file_metadata['root'] += phdr['CCDAMP'][:1].upper()
                    file_metadata['raw_postarg'] = phdr['POSTARG2']
                    if file_metadata['raw_postarg'] == 0.:
                        m1 = phdr.get("MOFFSET1", 0.)
                        m2 = phdr.get("MOFFSET2", 0.)
                        if m1 != 0. and m2 != 0.:
                            file_metadata['postarg'] = "{},{}px".format(m1, m2)
                        # FUV MAMA POSTARG
                        if "G140L" in file_metadata['mode'].upper():
                            spt_file_name = base_name.replace(file_ext, 'spt')
                            spt_file = os.path.join(file_path, spt_file_name)
                            if os.path.isfile(spt_file):
                                with fits.open(spt_file) as spt_inf:
                                    msmpos = spt_inf[1].header['OMSCYL1P']
                                    if msmpos < 800:
                                        file_metadata['postarg'] = '-3pos'
                                    else:
                                        file_metadata['postarg'] = '+3pos'
                            else:
                                msg = "{}: {}: WARNING: spt file not found"
                                print(msg.format(task, root))
                                msg = "G140L x1d file with no postarg information."
                                file_metadata['notes'] += msg
                    
                    mglobal = 0.
                    if "MAMA" in file_metadata["detector"].upper():
                        mglobal = np.sum(fits_file['SCI'].data)/file_metadata['exptime']
                    file_metadata['mglobal'] = mglobal

                    loc = "PARSING DATE"
                    date = fits_file[1].header['date-obs']
                    time = fits_file[1].header['time-obs']
                    date_str = "{}T{}".format(date, time)
                    file_metadata['date'] = Time(date_str)

                    ra = phdr['RA_TARG']
                    dec = phdr['DEC_TARG']
                
                loc = "ASSEMBLING WRITE DATA"
                new_target = (file_metadata['target'], 'Updated by ABSCAL')
                standard_star = find_star_by_name(new_target[0])
                if standard_star is None:
                    # Try by distance
                    standard_star = find_closest_star(ra, dec, max_distance=1.)
                if standard_star is not None:
                    delta_from_epoch = absdate(file_metadata['date']) - 2000.
                    epoch_ra = standard_star['ra']
                    epoch_dec = standard_star['dec']
                    pm_ra = standard_star['pm_ra']
                    pm_dec = standard_star['pm_dec']
                    delta_ra = pm_ra * delta_from_epoch/1000.
                    delta_dec = pm_dec * delta_from_epoch/1000.
                    corrected_ra = epoch_ra + delta_ra/3600.
                    corrected_dec = epoch_dec + delta_dec/3600.
                    new_ra = (corrected_ra, 'Updated for PM by ABSCAL')
                    new_dec = (corrected_dec, 'Updated for PM by ABSCAL')
                    if verbose:
                        msg = "{}: {}: Target Star: {}"
                        print(msg.format(task, root, file_metadata['target']))
                        print("\tEpoch RA,DEC = {},{}".format(epoch_ra, epoch_dec))
                        print("\tTime Since Epoch = {}".format(delta_from_epoch))
                        print("\tDelta RA,DEC = {},{}".format(delta_ra, delta_dec))
                        print("\tFinal RA,DEC = {},{}".format(corrected_ra, corrected_dec))
                else:
                    msg = file_metadata['target'] + " not a STIS standard star" 
                    new_target = (file_metadata['target'], msg)
                    new_ra = (phdr["RA_TARG"], msg)
                    new_dec = (phdr["DEC_TARG"], msg)
                file_metadata['ra_targ'] = new_ra[0]
                file_metadata['dec_targ'] = new_dec[0]
                loc = "WRITING NEW FILE"
                with fits.open(file_name, mode="update") as fits_file:
                    fits_file[0].header["TARGNAME"] = new_target
                    fits_file[0].header["RA_TARG"] = new_ra
                    fits_file[0].header["DEC_TARG"] = new_dec
                loc = "DONE"
                    
            except Exception as e:
                print("{}: {}: ERROR: {} {}".format(task, file_name, e, loc))
                for key in data_table.columns:
                    if key not in file_metadata:
                        print("\t{} missing".format(key))
                msg = "ERROR: Exception {} while processing.".format(str(e))
                file_metadata['notes'] += " {}".format(msg)
            
            data_table.add_exposure(file_metadata)
    
    # Adjust the values in the table based on the known metadata edits
    adjustments_file = get_data_file("abscal.stis", "metadata.yaml")
    if adjustments_file is not None:
        with open(adjustments_file, 'r') as inf:
            adjustments = yaml.safe_load(inf)
        data_table.adjust(adjustments)

    if data_table.n_exposures == 0:
        error_str = "Error: no files found for filespec {}"
        raise ValueError(error_str.format(file_template))
    
    data_table.sort(['root'])
    return data_table


def additional_args(**kwargs):
    """
    Adds process-specific command-line arguments.
    
    This function generates arguments (in a form understandable by 
    abscal.common.args.parse) to handle items unique to table creation.

    - How duplicate entries should be handled (important because this is process is the 
      one that adds new entries to a table)
    - The search template
    
    Returns
    -------
    args : dict
        Dictionary of tuples of arguments for building a module command-line argument 
        list.
    """
    module_name = kwargs.get('module_name', __name__)
    base_defaults = get_defaults(module_name)

    args = {}
    
    dup_help = "How to handle duplicate entries (entries defined as "
    dup_help += "duplicates if they have the same ipppssoot). Valid values are "
    dup_help += "'both' (keep both), 'preserve' (keep first), 'replace' (keep "
    dup_help += "second), and 'neither' (delete both). Duplicates should only "
    dup_help += "be an issue if an input table is specified. Default: {}"
    dup_help = dup_help.format(base_defaults['duplicates'])
    dup_args = ['--duplicates']
    dup_kwargs = {'dest': 'duplicates', 'help': dup_help, 
                  'default': base_defaults['duplicates']}
    args['duplicates'] = (dup_args, dup_kwargs)
    
    template_help = "The file template to match, or the path to search "
    template_help += "and the file template to match within that directory. "
    template_help += "If the '--paths' option is used to specify one or more input "
    template_help += "paths, then file file template will be joined to "
    template_help += "each input path."
    template_args = ['template']
    template_kwargs = {'help': template_help}
    args['template'] = (template_args, template_kwargs)
    
    return args


def parse_args(**kwargs):
    """
    Parse command-line arguments.
        
    Returns
    -------
    res : namespace
        A namespace populated by the command-line arguments.
    """    
    description_str = "Build metadata table from input files."
    default_output_file = kwargs.get('default_output_file', 'dirtemp.log')
    default_input_file = kwargs.get('default_input_file', 'dirirstare.log')

    args = additional_args(**kwargs)

    res = parse(description_str, default_output_file, args, **kwargs)
    
    if res.paths is not None: 
        if "," in res.paths:
            res.paths = res.paths.split(",")
        else:
            res.paths = [res.paths]
    else:
        res.paths = []
    
    if res.template is None:
        res.template = "o*raw.fits"
    
    if os.path.sep in res.template and len(res.paths) == 0:
        (template_path, template_value) = os.path.split(res.template)
        res.paths.append(template_path)
        res.template = template_value
    
    if len(res.paths) == 0:
        res.paths.append(os.getcwd())
    
    return res


def main(**kwargs):
    """
    Run the process.
    
    This function is called if the script is run directly (i.e. __name__ == "__main__"), 
    and is also imported by the binary command scripts as a way to run this process as a 
    standalone application.
    
    Parameters
    ----------
    kwargs : dict
        Contains keys named after keyword parameters (whether command-line arguments or 
        parameters used by table creation) that will override whatever value is set there.
        Note that specific exposure-specific values from data files will still override
        values specified here.
    """
    dt_str = datetime.datetime.now().strftime("%Y-%m-%d")
    kwargs['default_output_file'] = 'stis_spec_{}.log'.format(dt_str)
    
    res = parse_args(**kwargs)
    
    for key in kwargs:
        if hasattr(res, key):
            setattr(res, key, kwargs[key])
    
    table = STISDataTable(search_str=res.template,
                          search_dirs=res.paths,
                          table=res.in_file,
                          idl_mode=res.compat,
                          duplicates=res.duplicates)
    
    table = populate_table(data_table=table, verbose=res.verbose, compat=res.compat, **kwargs)
    table.write_to_file(res.out_file, res.compat)


if __name__ == "__main__":
    main(module_name='abscal.stis.preprocess_table_create')
