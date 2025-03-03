#! /usr/bin/env python
"""
This module contains the ExposureTable class. The class contains a table of
exposure metadata, as well as the ability to export that table in formats
suitable for IDL compatibility, or for further use by later scripts.

Authors
-------
    - Brian York

Use
---

This module provides the AbscalDataTable class, which holds an astropy table
of exposure metadata, along with the ability to present that data either in
a form expected by IDL, or in a form more suitable for passing on to
later scripts in abscal::

    from abscal.common.exposure_data_table import AbscalDataTable
    t = AbscalDataTable()
    
    t.add_exposure(<exposure metadata>)
    ... # repeat as needed
    t.write_to_file(<file>,<idl_mode>,<other arguments>)
    
The goal is to create an astropy.table.Table subclass that has the ability to successfully 
read in IDL-formatted tables, and to write out IDL-formatted tables *as an additional 
option*, but that generally uses its own format.
"""

import argparse
import datetime
import glob
import io
import os

import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.time import Time
from copy import deepcopy
from datetime import datetime as dt
from distutils.util import strtobool
from pathlib import Path
from simpleeval import simple_eval

from .logging import DEFAULT_LOGGER as logger
from .utils import build_expr


class AbscalDataTable(Table):
    """
    A class to represent HST exposure metadata for ABSCAL
    """

    def __init__(self, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, rows=None, copy_indices=True, **kwargs):
        """
        Initialize the table, look for (and extract) abscal-specific keywords/values, 
        and shove everything else in to the superclass.
        
        Parameters
        ----------
        kwargs : dict
            Dictionary of optional keyword arguments. Known kwargs include
                search_str : str
                    String used to search for files.
                search_dirs : list of str
                    Directories in which to search for search_str
                table : string or None
                    Initial table to read in on initialization
                duplicates : string
                    How to handle duplicate entries
                create_time : datetime.datetime
                    When the table was created (if creating based on an existing table)
        
        Parameters Passed to astropy.table.Table
        ----------------------------------------
        data : object, default None
        masked : object, default None
        names : object, default None
        dtype : object, default None
        meta : object, default None
        copy : bool, default True
        rows : object, default None
        copy_indices : bool, default True
        """
        # Persistent Metadata
        self.create_time, kwargs = self._get_kwarg('create_date', datetime.datetime.now(),
                                                   kwargs)
        self.search_str, kwargs = self._get_kwarg('search_str', self.default_search_str, 
                                                  kwargs)
        search_dirs, kwargs = self._get_kwarg('search_dirs', os.getcwd(), kwargs)
        if isinstance(search_dirs, str):
            search_dirs = [search_dirs]
        self.search_dirs = search_dirs
        self.duplicates, kwargs = self._get_kwarg('duplicates', 'both', kwargs)
        
        # Creation Flags
        initial_table, kwargs = self._get_kwarg('table', None, kwargs)
        
        if initial_table is not None:
            if data is not None:
                msg = f"ERROR: Supplied both initial data {data} and input file "
                msg += f"{initial_table}"
                raise ValueError(msg)
            data = self._read_file_to_table(initial_table, False, **kwargs)
        else:
            if data is None:
                kwargs = {}
                names = self.standard_columns.keys()
                dtype = [self.standard_columns[x]['dtype'] for x in names]

        # Initialize superclass
        super().__init__(data=data, masked=masked, names=names, dtype=dtype,
                         meta=meta, copy=copy, rows=rows,
                         copy_indices=copy_indices, **kwargs)
        self._fix_columns()
    
    
    @classmethod
    def from_idl_table(cls, table_file):
        """
        Reads in an IDL-formatted text table, and turns it into a standard table. Because 
        it uses the "add_exposure()" function (defined below), rows that just aren't 
        present in the IDL table will have default values.
        
        Parameters
        ----------
        cls : class
            The class type to return
        table_file : str
            The file to read.
        
        Returns
        -------
        table : instantiated class cls
            The populated table
        """
        parse_str = "%y-%m-%d %H:%M:%S"
        parse_creation_date = f"{cls.idl_str} %d-%b-%Y %H:%M:%S"

        # Read in the IDL-format table as undifferentiated ASCII
        idl_table = ascii.read(table_file, 
                               format="fixed_width",
                               col_starts=cls.idl_columns)
        
        # Read in the available metadata
        with open(table_file, 'r') as inf:
            date_line = inf.readline().strip()
            try:
                initial_date = dt.strptime(date_line, parse_creation_date)
            except ValueError:
                # try again with potentially commented line
                initial_date = dt.strptime(date_line, '# '+parse_creation_date)
            # metatext is "SEARCH FOR [path]/[search string]"
            search_path_items = inf.readline().strip().split()[-1]
            search_path, search_pattern = os.path.split(search_path_items)

        table = cls(search_str=search_pattern, search_dirs=search_path, 
                    duplicates='replace', create_time=initial_date)
        
        for row in idl_table:
            row_data = {}
            for col_name in cls.column_mappings.keys():
                row_data[cls.column_mappings[col_name]] = row[col_name]
            row_data['obset'] = row_data['root'][:6]
            file_name = '{row_data["root"]}_{search_pattern[-8:]}'
            if len(glob.glob(os.path.join(search_path, row_data['root']+'*.fits'))) > 0:
                row_data['path'] = search_path
                for ext in cls.idl_exts:
                    file_name = f'{row_data["root"]}_{ext}.fits'
                    if os.path.isfile(os.path.join(search_path, file_name)):
                        row_data['filename'] = file_name
                        break
            if 'IMG SIZE' in row:
                row_data['xsize'] = int(row['IMG SIZE'][:4].strip())
                row_data['ysize'] = int(row['IMG SIZE'][5:].strip())
            row_data['date'] = dt.strptime(f"{row['DATE']} {row['TIME']}", parse_str) 
            if 'POSTARG X,Y' in row:
                row_data['postarg1'] = float(row["POSTARG X,Y"][:5].strip())
                row_data['postarg2'] = float(row["POSTARG X,Y"][8:].strip())
            elif 'POSTARG' in row:
                row_data['postarg'] = float(row["POSTARG"].strip())
            row_data['notes'] = "Imported from IDL-style Table"
            table.add_exposure(row_data)
        
        return table
    
    
    def read(*args, **kwargs):
        """
        Override the read method to make sure that strings have the Object type. That way 
        they don't have a fixed character count, and we can append to them as needed.
        
        Parameters
        ----------
        args : list
            List of positional arguments. Passed directly to super.
        kwargs : dictionary
            Dictionary of keyword arguments. Passed directly to super.
        
        Returns
        -------
        out : astropy.table.Table
            The table that was read in.
        """
        parse_str = "%Y-%m-%dT%H:%M:%S"
        
        out = super().read(*args, **kwargs)
        for col in out.itercols():
            if col.dtype.kind in 'SU':
                out.replace_column(col.name, col.astype('object'))
            if col.name == 'date':
                date_data = [dt.strptime(x, parse_str) for x in col]
                date_col = Column(data=date_data)
                out.replace_column(col.name, date_col)
        return out


    def add_exposure(self, metadata_dict):
        """
        Add a new row to the table from a dictionary built up from an exposure.
        The dictionary should contain all of the standard abscal columns. If the
        value of the root column is already in the table, then deal with it as
        provided for in the internal duplicates metadata.
        
        Handling duplicates can currently be any of:
        
        - both: keep the existing entry, add the new entry.
        - preserve: keep the existing entry, discard the new entry.
        - replace: remove the existing entry, add the new entry.
        - neither: remove the existing entry, discard the new entry.
        
        Parameters
        ----------
        metadata_dict : dict
            A dictionary containing information on a single exposure.
        
        Returns
        -------
        success : bool
            True if the column was added, False if not. Will generally only be
            False in the case of duplicate columns.
        """
        for column in self.standard_columns:
            if column not in metadata_dict:
                # At this point, there is a column from the abscal table that doesn't have
                # a provided value. So there are a few options:
                #   - If the column should be present in either the IDL table or the full
                #     table, do nothing.
                #   - If the column is present in the full table, but not in the IDL 
                #     table, *and* if the column has a default value that's not just 
                #     'N/A', then add an entry to the metadata dict with its value set to 
                #     the default value.
                idl = self.standard_columns[column]['idl']
                default = 'N/A'
                if 'default' in self.standard_columns[column]:
                    default = self.standard_columns[column]['default']
                if (not idl) and (default != 'N/A'):
                    metadata_dict[column] = default
        
        for column in metadata_dict:
            if column not in self.standard_columns:
                logger.error(f"ERROR: Unknown Column {column}")
        
        # If this is a duplicate entry (same root exists in table already), treat it 
        # according to the duplicate policy. Otherwise, pass along to add_row().
        if metadata_dict['root'] in self['root']:
            if self.duplicates == 'both':
                self.add_row(metadata_dict)
            elif self.duplicates == 'preserve':
                return False
            elif self.duplicates == 'replace':
                remove_mask = [self['root'] == metadata_dict['root']]
                self.remove_rows(remove_mask)
                self.add_row(metadata_dict)
            elif self.duplicates == 'neither':
                remove_mask = [self['root'] == metadata_dict['root']]
                self.remove_rows(remove_mask)
                return False
            else:
                raise ValueError(f"ERROR: Unknown duplicate policy '{self.duplicates}'")
        else:
            self.add_row(vals=metadata_dict)
        
        # After adding anything to the table, make sure that string columns are 
        # variable-length.
        self._fix_columns()
        return True
    
    
    def adjust(self, adjustments):
        """
        Adjust the table based on an adjustment dictionary. The dictionary may contain 
        entries that edit columns or delete rows. This is provided because there are 
        several observations in the calibration set which are apparently just *wrong* in 
        some fundamental way (visit number incorrect, proposal incorrect, etc.) These 
        exposures need to be put into the correct place if they're found.
        
        This is a fairly generic function because we don't know when there might be 
        another exposure that needs to have something done to it.
        
        The adjustments dictionary consists of two sub-dictionaries, one of items to be 
        deleted, and the other of items to be edited. The dictionaries are structured as 
        follows:
        
        column: str
            column to look at
        key: str
            value to search for in the column
        operation: str
            operation to perform. ONLY PRESENT FOR EDIT
        value: str
            value to use when applying the operation. ONLY PRESENT FOR EDIT
        source: str
            Where the change came from
        reason: str
            Why the change needed to be made
        
        Parameters
        ----------
        adjustments : dict
            Dictionary of adjustments
        """
        # There may be exposures that are just bad, and need to not be used. Rather than
        # remove those rows, we set the "use" flag to False. We also add the reason that
        # the exposure needed to be removed, as found in the adjustment dictionary.
        for item in adjustments['delete']:
            expression, source, reason = build_expr(item)
            explanation = f"Removed for {source} ({reason})."
            for row in self:
                if simple_eval(expression, names=row):
                    row["use"] = False
                    if row["notes"] != "":
                        row["notes"] = f"{row['notes']}. {explanation}"
                    else:
                        row["notes"] = explanation
        
        # Currently the only supported edits are replacement or appending.
        for item in adjustments['edit']:
            value = item["key"]
            
            # The key provides just enough information to specify the right place to edit.
            # So we basically need to look for items whose first n characters match the 
            # key value, where n is the length of the key. Doing a search like that from
            # the full column is difficult, so we make a temporary column to search 
            # instead.
            prefix = [c[:len(value)] for c in self[item["column"]]]
            if value in prefix:
                mask = [x == value for x in prefix]
                masked = self[mask]
                if operation == "replace":
                    # Value has the form [old]->[new]
                    column = masked[item["column"]]
                    original, revised = item["value"].split("->")
                    new_col = [x.replace(original, revised) for x in column]
                elif operation == "append":
                    # Value is just the text to be appended
                    column = masked[item["column"]]
                    new_col = ["{x}{item['value']}" for x in column]
                elif operation == "add":
                    # Value is interpreted as a float
                    column = masked[item["column"]]
                    new_col = [x + float(item["value"]) for x in column]
                masked[item["column"]][:] = new_col[:]
                reason = f"Edited {item['column']} by {item['operation']}"
                reason += f" {item['value']} because {item['reason']} by {item['source]}."
                new_notes = [f"{x} {reason}" for x in masked["notes"]]
                masked["notes"][:] = new_notes[:]
        
        return
    
    
    def filtered_copy(self, filter_or_filters):
        """
        Get a copy of the table, filtered by the specified filters. The is a convenience 
        function that implicitly creates a deep copy, specifically so that edits can be 
        made to the now-filtered table without altering the original (or even needing to 
        know that there *is* and original).
        
        Parameters
        ----------
        filter_or_filters : list of filters or single filter
            The filters to apply
        
        Returns
        -------
        table : AbscalDataTable
            The filtered table
        """
        if not isinstance(filter_or_filters, list):
            filter_or_filters = [filter_or_filters]
        
        table = deepcopy(self)
        for filter in filter_or_filters:
            table = self._filter_table(table, filter)
        
        return table
    
    
    def write_to_file(self, file_name, idl_mode, **kwargs):
        """
        Write the table to a file, optionally using IDL strict compatibility.
        
        Parameters
        ----------
        file_name : str
            The file to write the table to.
        idl_mode : bool
            Whether to write the table in IDL compatibility mode.
        kwargs : dict
            Optional keyword arguments. Known arguments are:
                filters : list
                    List of filters to apply to the table before writing it.
        """
        date_str = self.create_time.strftime("%d-%b-%Y %H:%M:%S")
        format = kwargs.get('format', self.default_format)
        filters = kwargs.get('filters', [])
        
        table = self.filtered_copy(filters)
        if len(table) == 0:
            return

        table.comments = []
        table.comments.append(f"Start time: {date_str}")
        for dir in self.search_dirs:
            search_str = os.path.join(dir, self.search_str)
            table.comments.append(f"Searched {search_str}")
        for filter in filters:
            table.comments.append(f"Filtered by: {filter}")

        # Write the table out
        table.metadata = table.comments        
        table.write(file_name, format=format, overwrite=True)
        
        if idl_mode:
            # Also write an IDL-compatible version of the table
            file_ext = Path(file_name).suffix
            file_name = file_name.replace(file_ext, f'_idl{file_ext}')
            self._write_to_idl(file_name, table, create_time=self.create_time,
                               search_dirs=self.search_dirs,
                               search_str = self.search_str)
    

    @property
    def n_exposures(self):
        """
        Return the number of rows in the table.
        """
        return len(self)


    @staticmethod
    def _filter_table(table, filter):
        """
        Filter a table by means of a particular filter. The filter should be
        either a known string value, or a callable which takes an astropy table
        and returns a filtered version of the table.
        
        Parameters
        ----------
        table : astropy.table.Table
            The table to filter. Should probably in fact be an AbscalDataTable.
        filter : str or callable
            Either a keyword string or a callable. Known strings are:
            
            stare: 
                remove scan_rate > 0
            scan: 
                remove scan_rate == 0
            filter: 
                remove grism data
            grism: 
                keep all grism data and the closest-in-time filter exposure from the 
                same visit as the grism exposure.
            obset:<value>: 
                keep all data where the root column begins with <value>
        
        Returns
        -------
        filtered_table : astropy.table.Table
            Filtered table.
        """
        if isinstance(filter, str):
            if "obset" in filter:
                val = filter[filter.find(":")+1:]
                len_val = len(val)
                mask = [r[:len_val] == val for r in table['root']]
                table = table[mask]
        else:
            table = filter(table)

        return table


    def _read_file_to_table(self, file_name, idl_mode, **kwargs):
        """
        Read an abscal table from a file, optionally using IDL compatibility mode. Note 
        that, even if idl_mode is set to False, if reading the table fails with the 
        standard method, then the program will attempt to read the table with the IDL 
        reader.

        Parameters
        ----------
        file_name : str
            The table file to read.
        idl_mode : bool
            Whether to operate in IDL compatibility mode.
        kwargs : dict
            Optional keyword arguments. Known arguments are:
            
            format: str
                Table format to attempt to read.
        
        Returns
        -------
        t : astropy.table.Table
            The table that was read in.
        """
        if idl_mode:
            return self.from_idl_table(file_name)
        format = kwargs.get('format', self.default_format)
        
        t = Table()
        try:
            t = t.read(file_name, format=format)
        except Exception as e:
            logger.error(f"ERROR while reading: {e}")
            logger.error("Falling back to IDL-format reader")
            t = AbscalDataTable.from_idl_table(file_name)
        
        return t

        
    def _write_to_idl(self, file_name, table, **kwargs):
        """
        Write the table in strict IDL mode. This mode uses different columns (which 
        are sometimes combinations of several variables that need to be separated for 
        use, and are thus stored in separate columns in standard mode).
        
        Parameters
        ----------
        file_name : str
            The file to write to.
        table : astropy.table.Table
            The (filtered) table to write
        kwargs : dict
            Dictionary of optional keywords. Known keywords include:
            
            create_time: datetime.DateTime
                Table creation date
        """
        if len(table) == 0:
            return
        parse_str = "%Y-%m-%dT%H:%M:%S.000"
        
        column_formats = self.idl_column_formats
        columns = column_formats.keys()

        if isinstance(table["date"][0], str):
            dates = [dt.strptime(d, parse_str) for d in table["date"]]
        else:
            dates = table["date"]
        
        idl_table = Table()
        for column in columns:
            if column in table.column_mappings:
                data = table[table.column_mappings[column]]
            elif column == "IMG SIZE":
                xsize = table["xsize"]
                ysize = table["ysize"]
                data = [f"{x:4d}x{y:4d}" for x,y in zip(xsize,ysize)]
            elif column == "DATE":
                data = [d.strftime("%y-%m-%d") for d in dates]
            elif column == "TIME":
                data = [d.strftime("%H:%M:%S") for d in dates]
            elif column == "POSTARG X,Y":
                p1 = table["postarg1"]
                p2 = table["postarg2"]
                data = [f"{x:6.1f}, {y:6.1f}" for x,y in zip(p1,p2)]
            else:
                raise ValueError(f"Trying to create nonexistent column {column}")
            format = column_formats[column]
            idl_table[column] = Column(data, name=column, format=format)
        
        # Write the table itself
        idl_table.write(file_name,
                        format='ascii.fixed_width',
                        delimiter=' ',
                        delimiter_pad=None,
                        bookend=False,
                        overwrite=True)
        
        # Re-open the table and write the preamble lines
        create_time = table.create_time
        date_str = create_time.strftime("%d-%b-%Y %H:%M:%S")
        if isinstance(table.search_dirs, str):
            search_dirs = [table.search_dirs]
        else:
            search_dirs = table.search_dirs
        search_str = table.search_str
        with open(file_name, 'r+') as table_file:
            content = table_file.read()
            table_file.seek(0, 0)
            table_file.write(f"# {self.idl_str} {date_str}\n")
            for search_dir in search_dirs:
                full_search_str = os.path.join(search_dir, search_str)
                table_file.write(f"# SEARCH FOR {full_search_str}\n")
            table_file.write(content)

        return

    
    @staticmethod
    def _get_kwarg(key, default, kwargs):
        """
        Get an abscal-specific keyword argument from the kwargs dictionary (if present), 
        and pop that keyword out of the kwargs dictionary. If not present, return a 
        (supplied) default value.
        
        Parameters
        ----------
        key : str
            The kwarg key to look for
        default : obj
            The default value of key
        kwargs : dict
            The keyword dictionary to search.
        
        Returns
        -------
        value : obj
            The value found (or the default value)
        kwargs : dict
            The keyword dictionary with 'key' removed from it (if it was present).
        """
        if key in kwargs:
            return kwargs.pop(key), kwargs
        return default, kwargs
    
    
    def _fix_columns(self):
        """
        Adjust all columns to replace fixed-length strings with objects.
        """
        for col in self.itercols():
            if col.dtype.kind in 'SU':
                self.replace_column(col.name, col.astype('object'))
            if self.standard_columns[col.name]['dtype'] == '?':
                if col.dtype.kind in 'fib':
                    new_values = [bool(x) for x in col]
                else:
                    new_values = [bool(strtobool(x)) for x in col]
                self.replace_column(col.name, new_values)
    

    column_mappings = {}
    
    default_format = 'ascii.ipac'
    
    instrument = ''
    default_search_str = '*.fits'
    idl_str = ''
    idl_columns = (0, 11, 19, 35, 41, 54, 64, 73, 82, 89, 98, 112)
    idl_exts = ['flt', 'ima', 'raw']
    idl_column_formats = {}
    
    standard_columns = {}
