#! /usr/bin/env python
"""
This module takes the name of an input metadata table, groups the exposures in 
that table by program and visit, and then:
    - calibrates each exposue
    - coadds together all exposures that have the same program/visit/star

Authors
-------
    - Brian York (all python code)
    - Ralph Bohlin (original IDL code)

Use
---
    This module is intended to be either run from the command line or used by
    other module code as the first step (creating an annotated list of files 
    for flux calibration).
    ::
        python coadd_grism.py <input_file>

Dependencies
------------
    - ``astropy``
"""

__all__ = ['coadd']

import datetime
import glob
import json
import os

import matplotlib.pyplot as plt
import numpy as np

from astropy import wcs
from astropy.convolution import convolve, Box1DKernel
from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.time import Time
from copy import deepcopy
from pathlib import Path
from photutils.centroids import centroid_sources
from photutils.centroids import centroid_1dg as g1
from photutils.centroids import centroid_2dg as g2
from photutils.centroids import centroid_com as com
from photutils.detection import DAOStarFinder
from rebin import rebin
from scipy import ndimage
from scipy.interpolate import interp2d

from abscal.common.args import parse
from abscal.common.utils import get_data_file, set_params, set_image
from abscal.common.exposure_data_table import AbscalDataTable
from abscal.wfc3.util_filter_locate_image import locate_image


def make_monotonic(wave, indx):
    """
    Make the supplied wavelength array (and the index of the wavelength values
    that were just altered) monotonically increasing.
    
    Parameters
    ----------
    wave : np.ndarray
        Array of wavelength values
    indx : np.ndarray
        Index of wave entries that were just replaced by a different order.
    """
    ibreak = np.max(indx)
    bad = np.where(wave[ibreak+1:1013] < wave[ibreak])
    npts = len(bad[0])
    if npts > 0:
        print("Non-monotonic: {}".format(bad))
        delwl = np.where((wave[ibreak+1+npts] - wave[ibreak])>1.,
                         (wave[ibreak+1+npts] - wave[ibreak]), 1.)
        delwl /= npts
        wave[ibreak:ibreak+npts] = wave[ibreak] + delwl*np.arange(1, npts+1)
    return wave


def set_hdr(hdr, params):
    """
    Set header parameters based on metadata input into or found from extraction.
    
    Parameters
    ----------
    h : astropy.io.fits header
        The header to edit. Must be writable.
    params : dict
        Dictionary of parameters to write.
    
    Returns
    -------
    h : astropy.io.fits header
        Edited header.
    """
    hdr["XACTUAL"] = (params['xc'], "Actual Found Direct Image X Position")
    hdr["YACTUAL"] = (params['yc'], "Actual Found Direct Image Y Position")
    hdr["XCAMERR"] = (params['xerr'], "Pointing error XACTUAL-Predicted (px)")
    hdr["YCAMERR"] = (params['yerr'], "Pointing error YACTUAL-Predicted (px)")
    if params['axeflg']:
        hdr['HISTORY'] = 'WAVELENGTH solution per AXE coef. No Z-order.'
    else:
        hdr['HISTORY'] = 'WAVELENGTH solution per 0-order position ISR'
        hdr["XZORDER"] = (params['xc_f'], "Zero order found X-Position")
        hdr["YZORDER"] = (params['yc_f'], "Zero order found Y-Position")
        xzerr = params['xastr'] - params['xc_f']
        hdr["XZERR"] = (xzerr, "Zero order PREDICT-FOUND (px)")
        yzerr = params['yastr'] - params['yc_f']
        hdr["YZERR"] = (yzerr, "Zero order PREDICT-FOUND (px)")
    hdr["AVGBKGR"] = (params['avgbkg'], "Average Background from wfc_flatscl")
    hdr["FLATFILE"] = (params['flatfile'], "Flatfield File Used.")
    hdr["EXTC0"] = (params['coef_0'], "Extraction Position: Constant Term")
    hdr["EXTC1"] = (params['coef_1'], "Extraction Position: Linear Term")
    hdr["DIRIMAGE"] = (params['ref'], "Direct Image")
    hdr["GWIDTH"] = (params['gwidth'], "Ext. Width")
    hdr["BWIDTH"] = (params['bwidth'], "Bkg. Ext. Width")
    hdr["UBDIST"] = (params['ubdist'], "Up. Bkg. Dist")
    hdr["LBDIST"] = (params['lbdist'], "Lo. Bkg. Dist")
    hdr["BMEDIAN"] = (params['bmedian'], "Bkg. Median filter width")
    hdr["BMEAN1"] = (params['bmean1'], "Bkg. Boxcar Smooth width 1")
    hdr["BMEAN2"] = (params['bmean2'], "Bkg. Boxcar Smooth width 2")
    hdr["WLOFFSET"] = (params['wl_offset'], "Wavelength correction.")
    hdr["ANGLE"] = (params['angle'], "Found angle of spectrum WRT x-axis")
    now = datetime.datetime.now()
    now_str = now.strftime("%Y-%m-%d %H:%M:%S")
    hdr["HISTORY"] = "Written by calwfc_spec.pro {}".format(now_str)
    
    return hdr


def reduce_flatfield(input_table, params):
    """
    Scale the flatfield image based on the coefficients in the flatfield cube.
    
    Parameters
    ----------
    input_table : abscal.common.exposure_data_table.AbscalDataTable
        Single-row input
    params : dict
        Parameter dictionary. Contains
            image : np.ndarray
                Science image
            hdr : astropy.io.fits.header
                FITS header for the science image
            flat : np.ndarray
                image file of the WFC3 flatfield
            wave : np.ndarray
                wavelength mapping for the science image
            verbose : bool
                Whether to print diagnostic output
    
    Returns
    -------
    image : np.ndarray
        Flatfielded science image
    flatfile : str
        File name of the flatfield used.
    """
    task = "wfc3_grism_reduce_flatfield"
    
    root = params['root']
    image = params['image']
    hdr = params['hdr']
    flat = params['flat']
    wave = params['wave']
    ltv1 = params['ltv1']
    ltv2 = params['ltv2']
    filter = params['filter']
    verbose = params['verbose']
    
    grating = input_table['filter']

    np_formatter = {'float_kind':lambda x: "{:8.4f}".format(x)}
    np_opt = {'max_line_width': 175, 'threshold': 2000000,#}#,
               'formatter': np_formatter}
    
    preamble = "{}: {}: {}".format(task, root, grating)
    
    if verbose:
        print("{}: starting".format(preamble))
    
    # Create an array for quadratic fits of SED
    coef3 = np.zeros((1014, 1014), dtype='float64')

    # Get flatfield data cube
    cal_data = get_data_file("abscal.wfc3", "calibration_files.json")
    with open(cal_data, 'r') as inf:
        cal_files = json.load(inf)
    flat_file_name = cal_files["flatfield_cube"]["wfc3"][filter]
    flat_file = get_data_file("abscal.wfc3", flat_file_name)

    with fits.open(flat_file) as in_flat:
        # Figure out extension of first coefficient extension of flatfield
        coeff_offset = 0
        if "sed" in flat_file_name:
            # Header in ext 0, coefficients start in ext 1
            coeff_offset = 1
        coef0 = in_flat[0+coeff_offset].data
        coef1 = in_flat[1+coeff_offset].data
        coef2 = in_flat[2+coeff_offset].data
        
        flat_hdr = in_flat[0].header

    wmin, wmax = float(flat_hdr["WMIN"]), float(flat_hdr["WMAX"])
    abswl = abs(wave)
    x = (abswl - wmin)/(wmax - wmin)
    if grating == 'G102':
        good = np.where((abswl >= 7000) & (abswl <= 12000))
    elif grating == 'G141':
        good = np.where((abswl >= 9000) & (abswl <= 18000))
    else:
        raise ValueError("Unknown Grating {}".format(grating))
    xgood = x[good]
    ff = np.ones_like(image)
    ns, nl = hdr['naxis1'], hdr['naxis2']
    xpx = ((good[0]+ltv1).astype('int32'),)
    
    for i in range(nl):
        if abswl.ndim == 2:
            if grating == 'G102':
                good = np.where((abswl[:,i] >= 7000) & (abswl[:,i] <= 12000))
            elif grating == 'G141':
                good = np.where((abswl[:,i] >= 9000) & (abswl[:,i] <= 18000))
            else:
                raise ValueError("Unknown Grating {}".format(grating))
            xgood = x[:,i][good]
        ypx = i + int(ltv2)
        coef0i = coef0[ypx,:][xpx]
        coef1i = coef1[ypx,:][xpx]
        coef2i = coef2[ypx,:][xpx]
        coef3i = coef3[ypx,:][xpx]
        ff[i,:][good] = coef0i + coef1i*xgood + coef2i*xgood**2 + coef3i*xgood**3
    
    if ns < 1014:
        if verbose:
            msg = "{}: Subarray w/ ltv1,ltv2 = ({},{})"
            print(msg.format(preamble, ltv1, ltv2))
    ff = np.where(ff<=0.5, 1., ff)
    image = image/ff
    if verbose:
        print("{}: flatfield for number of wave points {}".format(preamble, ns))

    if verbose:
        print("{}: finished".format(preamble))
    
    return image, flat_file


def calculate_order(params, xc, yc):
    """
    Calculate a spectral order coefficient given the corresponding coefficients.
    
    Parameters
    ----------
    params : dict
        Dictionary of coefficients. Contains 6 coefficients labeled as
        a, b, c, d, e, f, which correspond to a constant (a), linear terms in 
        x (b) and y (c), and quadratic terms in x (d), xy (e), and y (f)
    xc : float
        Central x pixel of spectrum.
    yc : float
        Central y pixel of spectrum.
    
    Returns
    -------
    result : float
        The result of the calculation.
    """
    result = params['a'] + params['b']*xc + params['c']*yc
    result += params['d']*xc*xc + params['e']*xc*yc + params['f']*yc*yc
    return result


def reduce_wave_axe(row, params):
    """
    Calculate the wavelength vector for a WFC3 grism image. Per the IDL code,
    this should only be used for the AXE solution. The Zero-order (Z-ORD)
    solution should use a different function.
    
    Parameters
    ----------
    input_table : abscal.common.exposure_data_table.AbscalDataTable
        Table containing the input data
    params : dict
        Dictionary of parameters, including:
            xc, yc : float
                predicted central location
            verbose : bool
                whether to print diagnostic output
    
    Returns
    -------
    x_arr : np.ndarray
        X pixel co-ordinate array (just effectively np.arange(1014))
    wave : np.ndarray
        Wavelength vector, customized for each order
    angle : float
        Average slope of spectrum
    wav1st : np.ndarray
        First-order wavelengths for flatfielding with FF cube data.
    """
    task = "wfc3_grism_reduce_wave_axe"
    file = os.path.join(row["path"], row["filename"])
    xcin = params['xin']
    ycin = params['yin']
    root = row['root']
    preamble = '{}: {}'.format(task, root)
    
    if verbose:
        print("{}: starting.".format(preamble))

    with fits.open(file) as inf:
        image = inf['SCI'].data
        filter = row['filter']
        x_arr = np.arange(1014, dtype='i16')
        xc, yc = xcin, ycin
        ns = inf[0].header['naxis1']
        if ns < 1014:
            ltv1, ltv2 = inf[1].header['ltv1'], inf[1].header['ltv2']
            xc, yc = xcin + ltv1, ycin + ltv2
        
        if row['filter'] == 'G102':
            # Order information taken from WFC3 ISR 2016-15 - AXE.
        
            # First Order
            params = { 'a': 6344.081, 'b': 0.20143085, 'c': 0.080213136,
                       'd': -0.00019613, 'e': 0.0000301396, 'f': -0.0000843157 }
            a0 = calculate_order(params, xc, yc)
            params = { 'a': 24.00123, 'b': -0.00071606, 'c': 0.00084115,
                       'd': 8.9775481e-7, 'e': -3.160441e-7, 'f': 7.1404362e-7 }
            a1 = calculate_order(params, xc, yc)
            a0p1, a1p1 = a0, a1
            wave = a0 + a1*(x - xc)
            
            # Zeroth Order: 
            #   'Susana says the ISR is crap and seems true to me... 
            #   Ignore 0-order.'
            
            # -1st order:
            params = { 'a': -6376.843, 'b': 0., 'c': -1.11124775,
                       'd': 0., 'e': 0., 'f': 0.00095901 }
            a0 = calculate_order(params, xc, yc)
            params = { 'a': -24.27561, 'b': 0., 'c': -0.003251784,
                       'd': 0., 'e': 0., 'f': 1.4988e-6 }
            a1 = calculate_order(params, xc, yc)
            a0m1, a1m1 = a0, a1
            
            # Paraphrased from calwfc_spec.pro, 'the -1 formula seems better 
            #   than the +1 formula for finding 0-order pixels.
            # As such, find probably 0-order pixels (wavelength<300), and 
            #   replace their values with values calculated with these new
            #   coefficients. Then test to ensure the wavelength array is still
            #   monotonically increasing.
            
            indx = np.where(wave < 300)
            if len(indx[0]) > 0:
                wave[indx] = -(a0 + a1*(indx[0]-xc))
                wave = make_monotonic(wave, indx)
            wav1st = wave
            
            # 2nd order
            params = { 'a': 3189.9195, 'b': 0.291324446, 'c': 0.039748254,
                       'd': 0.000405844, 'e': 6.5079365e-6, 'f': -0.00003221 }
            a0 = calculate_order(params, xc, yc)
            params = { 'a': 12.08004, 'b': -0.00046352, 'c': 0.000670315,
                       'd': 5.7894508e-7, 'e': 0., 'f': 5.667e-8 }
            a1 = calculate_order(params, xc, yc)
            a0p2, a1p2 = a0, a1
            
            # Replace points with wave>14000 with 2nd order
            indx = np.where(wave > 14000)
            if len(indx[0]) > 0:
                wave[indx] = 2*(a0 + a1*(indx[0]-xc))
                wav1st[indx] = a0 + a1*(indx[0]-xc)
                wave = make_monotonic(wave, indx)
                wav1st = make_monotonic(wav1st, indx)
                        
            # 3rd order is commented out in IDL code, provided below:
#             params = { 'a': 2.17651e+03, 'b': 0., 'c': 5.01084e-02,
#                        'd': 0., 'e': 0., 'f': 0. }
#             a0 = calculate_order(params, xc, yc)
#             params = { 'a': 8.00453, 'b': 0., 'c': 4.28339e-04,
#                        'd': 0., 'e': 0., 'f': 0. }
#             a1 = calculate_order(params, xc, yc)
#             a0p3, a1p3 = a0, a1
# 
#             index = np.where(wave > 23500)
#             if len(indx[0]) > 0:
#                 wave[indx] = 2*(a0 + a1*(indx[0]-xc))
#                 wav1st[indx] = a0 + a1*(indx[0]-xc)

            angle = np.radians(0.66)
            
        elif row['filter'] == 'G141':
            # WFC3 ISR 2009-17

            # First Order
            params = { 'a': 8951.386, 'b': 0.08044033, 'c': -0.00927970,
                       'd': 0.000021857, 'e': -0.000011048, 'f': 0.000033527 }
            a0 = calculate_order(params, xc, yc)
            params = { 'a': 44.972279, 'b': 0.000492789, 'c': 0.00357824,
                       'd': -9.175233345e-7, 'e': 2.235506e-7, 'f': -9.25869e-7 }
            a1 = calculate_order(params, xc, yc)
            a0p1, a1p1 = a0, a1
            wave = a0 + a1*(x - xc)
            
            #Negative First Order
            params = { 'a': -46.4855, 'b': 0., 'c': -0.8732184,
                       'd': 0., 'e': 0., 'f': 0.0009233797 }
            a0 = calculate_order(params, xc, yc)
            params = { 'a': 44.972279, 'b': 0., 'c': -0.004813895,
                       'd': 0., 'e': 0., 'f': 2.0768286663e-6 }
            a1 = calculate_order(params, xc, yc)
            a0p1, a1p1 = a0, a1

            indx = np.where(wave < 300)
            if len(indx[0]) > 0:
                wave[indx] = -(a0 + a1*(indx[0]-xc))
                wave = make_monotonic(wave, indx)
            wav1st = wave

            # 2nd order
            params = { 'a': 4474.5297, 'b': 0.17615670, 'c': 0.046354019,
                       'd': -0.00012965, 'e': 0.00001513, 'f': -0.00002961 }
            a0 = calculate_order(params, xc, yc)
            params = { 'a': 22.8791467, 'b': -0.0002159637, 'c': 0.00133454,
                       'd': 4.277729e-8, 'e': -8.522518e-8, 'f': 6.08125e-8 }
            a1 = calculate_order(params, xc, yc)
            a0p2, a1p2 = a0, a1

            indx = np.where(wave > 19000)
            if len(indx[0]) > 0:
                wave[indx] = 2*(a0 + a1*(indx[0]-xc))
                wav1st[indx] = a0 + a1*(indx[0]-xc)
                wave = make_monotonic(wave, indx)
                wav1st = make_monotonic(wav1st, indx)

            # 3rd order is commented out in IDL code, provided below:
#             params = { 'a': 3.00187e+03, 'b': 1.04205e-01, 'c': -1.18134e-03,
#                        'd': 0., 'e': 0., 'f': 0. }
#             a0 = calculate_order(params, xc, yc)
#             params = { 'a': 1.52552e+01, 'b': -2.08555e-04, 'c': 9.55645e-04,
#                        'd': 0., 'e': 0., 'f': 0. }
#             a1 = calculate_order(params, xc, yc)
#             a0p3, a1p3 = a0, a1
# 
#             index = np.where(wave > 31000)
#             if len(indx[0]) > 0:
#                 wave[indx] = 2*(a0 + a1*(indx[0]-xc))
#                 wav1st[indx] = a0 + a1*(indx[0]-xc)

            angle = np.radians(0.44)

        # Check for monotonic
        dx = np.diff(wave)
        if not (np.all(dx<=0) or (np.all(dx>=0))):
            if verbose:
                for i in range(1, len(wave)):
                    if wave[i] < wave[i-1]:
                        msg = "{}: Not monotonic at {}: {}-{}"
                        print(msg.format(preamble, i, wave[i-1], wave[i]))
            # For now just die on an exception if not monotonic.
            msg = "Wavelength array {} not monotonic".format(wave)
            raise ValueError(msg)
        
        wave += params['wl_offset']
        wav1st += params['wl_offset']
        
        if ns < 1014:
            ibeg = ltv1
            iend = ltv1 + ns - 1
            wave = wave[ibeg:iend]
            wav1st = wav1st[ibeg:iend]
            if verbose:
                msg = "{}: wave minmax = {},{} Ref. (1014)px at ({},{})"
                print(msg.format(preamble, min(wave), max(wave), xc, yc))
       
    with fits.open(file, mode='update') as f:
        f[0].header['XC'] = (xcin, 'Dir img ref X position used for AXE WLs')
        f[0].header['YC'] = (ycin, 'Dir img ref Y position used for AXE WLs')
        f[0].header['A0+1ST'] = (a0p1, 'Constant Term of the +1st order disp.')
        f[0].header['A1+1ST'] = (a1p1, 'Linear Term of the +1st order disp.')
        f[0].header['A0-1ST'] = (a0m1, 'Constant Term of the -1st order disp.')
        f[0].header['A1-1ST'] = (a1m1, 'Linear Term of the -1st order disp.')
        f[0].header['A0+2ND'] = (a0p2, 'Constant Term of the +2nd order disp.')
        f[0].header['A1+2ND'] = (a1p2, 'Linear Term of the +2nd order disp.')
#         f[0].header['A0+3RD'] = (a0p3, 'Constant Term of the +3rd order disp.')
#         f[0].header['A1+3RD'] = (a1p3, 'Linear Term of the +3rd order disp.')
    
    if verbose:
        print("{}: finishing.".format(preamble))

    return x_arr, wave, angle, wav1st


def reduce_wave_zord(row, params):
    """
    Calculate the wavelength vector for a WFC3 grism image. Per the IDL code,
    this should only be used for the Zeroth Order (ZORD) solution. The AXE
    solution should use a different function.
    
    Parameters
    ----------
    row : abscal.common.exposure_data_table.AbscalDataTable.Row
        Row of astropy Table containing the input data
    params : dict
        Dictionary of parameters, including:
            xc, yc : float
                predicted central location
            verbose : bool
                whether to print diagnostic output
    
    Returns
    -------
    x_arr : np.ndarray
        X pixel co-ordinate array (just effectively np.arange(1014))
    wave : np.ndarray
        Wavelength vector, customized for each order
    angle : float
        Average slope of spectrum
    wav1st : np.ndarray
        First-order wavelengths for flatfielding with FF cube data.
    """
    task = "wfc3_grism_reduce_wave_zord"
    root = row['root']
    preamble = "{}: {}".format(task, root)
    file = os.path.join(row["path"], row["filename"])
    zxposin = params['xin']
    zyposin = params['yin']
    verbose = params['verbose']

    np_formatter = {'float_kind':lambda x: "{:7.1f}".format(x)}
    np_opt = {'max_line_width': 175, 'formatter': np_formatter,
              'threshold': 2000000}
    
    if verbose:
        print("{}: starting.".format(preamble))
    
    with fits.open(file) as inf:
        image = inf['SCI'].data
        filter = row['filter']
        x_arr = np.arange(1014, dtype='int32')
        zxpos, zypos = zxposin, zyposin

        if image.shape[1] < 1014:
            ltv1, ltv2 = inf[1].header['ltv1'], inf[1].header['ltv2']
            zxpos, zypos = zxposin + ltv1, zyposin + ltv2
            if verbose:
                msg = "{}: subarray shifts Z-ord ref px by ({},{})"
                print(msg.format(preamble, ltv1, ltv2))

        if filter == 'G102':

            # First Order        
            b_coeff = [148.538, 0.145605, -0.008558]
            m_coeff = [23.8796, -0.000332, 0.001489]
            b = b_coeff[0] + b_coeff[1]*zxpos + b_coeff[2]*zypos
            m = m_coeff[0] + m_coeff[1]*zxpos + m_coeff[2]*zypos
            bp1, mp1 = b, m
            wave = b + m*(x_arr - zxpos)
#             print("Wave after first order")
#             print(np.array2string(wave, **np_opt))
            
#             print("Wave: ", end='')
#             for i in range(len(wave)):
#                 print("{:.1f}".format(wave[i]), end=', ')
#             print("done wave.")
            
            # Negative First Order
            b_coeff = [205.229, -0.015426, -0.019207]
            m_coeff = [24.7007, 0.000047, 0.001478]
            b = b_coeff[0] + b_coeff[1]*zxpos + b_coeff[2]*zypos
            m = m_coeff[0] + m_coeff[1]*zxpos + m_coeff[2]*zypos
            bm1, mm1 = b, m

            indx = np.where(wave < -7000)
            if len(indx[0]) > 0:
                wave[indx] = b + m*(indx[0] - zxpos)
#                 print("Wave after negative first order")
#                 print(np.array2string(wave, **np_opt))
                wave = make_monotonic(wave, indx)
            wav1st = deepcopy(wave)
#             print("Wave: ", end='')
#             for i in range(len(wave)):
#                 print("{:.1f}".format(wave[i]), end=', ')
#             print("done wave.")
            
            # Second Order
            b_coeff = [213.571, 0.561877, -0.040419]
            m_coeff = [23.9983, -0.000797, 0.001532]
            b = b_coeff[0] + b_coeff[1]*zxpos + b_coeff[2]*zypos
            m = m_coeff[0] + m_coeff[1]*zxpos + m_coeff[2]*zypos
            bp2, mp2 = b, m

            indx = np.where(wave > 14000)
            if len(indx[0]) > 0:
                wave[indx] = b + m*(indx[0] - zxpos)
#                 print("Wave after second order")
#                 print(np.array2string(wave, **np_opt))
                wav1st[indx] = (b + m*(indx[0] - zxpos))//2
                wave = make_monotonic(wave, indx)
                wav1st = make_monotonic(wav1st, indx)
#             print("Wave: ", end='')
#             for i in range(len(wave)):
#                 print("{:.1f}".format(wave[i]), end=', ')
#             print("done wave.")
#             print("Wave after monotonic Correction")
#             print(np.array2string(wave, **np_opt))

            angle = np.radians(0.61)
        
        elif filter == 'G141':

            # First Order        
            b_coeff = [156.339, 0.111342, -0.010926]
            m_coeff = [45.3203, -0.000408, 0.002818]
            b = b_coeff[0] + b_coeff[1]*zxpos + b_coeff[2]*zypos
            m = m_coeff[0] + m_coeff[1]*zxpos + m_coeff[2]*zypos
            bp1, mp1 = b, m
            wave = b + m*(x_arr - zxpos)
            
            # Negative First Order
            b_coeff = [165.764, 0.055688, 0.016568]
            m_coeff = [46.4521, 0.000382, 0.002960]
            b = b_coeff[0] + b_coeff[1]*zxpos + b_coeff[2]*zypos
            m = m_coeff[0] + m_coeff[1]*zxpos + m_coeff[2]*zypos
            bm1, mm1 = b, m

            indx = np.where(wave < -8000)
            if len(indx[0]) > 0:
                wave[indx] = b + m*(indx[0] - zxpos)
                wave = make_monotonic(wave, indx)
            wav1st = deepcopy(wave)
            
            # Second Order
            b_coeff = [193.093, 0.164508, -0.017031]
            m_coeff = [45.6062, -0.000499, 0.002840]
            b = b_coeff[0] + b_coeff[1]*zxpos + b_coeff[2]*zypos
            m = m_coeff[0] + m_coeff[1]*zxpos + m_coeff[2]*zypos
            bp2, mp2 = b, m

            indx = np.where(wave > 18000)
            if len(indx[0]) > 0:
                wave[indx] = b + m*(indx[0] - zxpos)
                wav1st[indx] = (b + m*(indx[0] - zxpos))//2
                wave = make_monotonic(wave, indx)
                wav1st = make_monotonic(wav1st, indx)

            angle = np.radians(0.42)

        # Check for monotonic
        dx = np.diff(wave)
        if not (np.all(dx<=0) or (np.all(dx>=0))):
            if verbose:
                for i in range(1, len(wave)):
                    if wave[i] < wave[i-1]:
                        msg = "{}: Not monotonic at {}: {}-{}"
                        print(msg.format(preamble, i, wave[i-1], wave[i]))
            # For now just die on an exception if not monotonic.
            msg = "Wavelength array {} not monotonic".format(wave)
            raise ValueError(msg)
        
        if verbose:
            print("{}: Adding offset {}".format(preamble, params['wl_offset']))
        wave += params['wl_offset']
        wav1st += params['wl_offset']

        if image.shape[1] < 1014:
            ibeg = ltv1
            iend = ltv1 + ns - 1
            wave = wave[ibeg:iend]
            wav1st = wav1st[ibeg:iend]
            if verbose:
                msg = "{}: Wave minmax = {},{} Ref. (1014)px at ({},{})"
                print(msg.format(preamble, min(wave), max(wave), xc, yc))

    with fits.open(file, mode='update') as f:
        f[0].header['B+1ST'] = (bp1, 'Constant Term of the +1st order disp.')
        f[0].header['M+1ST'] = (mp1, 'Linear Term of the +1st order disp.')
        f[0].header['B-1ST'] = (bm1, 'Constant Term of the -1st order disp.')
        f[0].header['M-1ST'] = (mm1, 'Linear Term of the -1st order disp.')
        f[0].header['B+2ND'] = (bp2, 'Constant Term of the +2nd order disp.')
        f[0].header['M+2ND'] = (mp2, 'Linear Term of the +2nd order disp.')
    
    if verbose:
        print("{}: finishing.".format(preamble))
   
    return x_arr, wave, angle, wav1st


def reduce_scan(input_table, params, arg_list):
    """
    Reduces scan-mode grism data
    """
    verbose = arg_list.verbose
    interactive = arg_list.trace
    bkg_flat_order = arg_list.bkg_flat_order
    
    file = os.path.join(row["path"], row["filename"])

    with fits.open(file) as inf:
        image = inf['SCI'].data
        filter = row['filter']
        xsize, ysize = image.shape[1], image.shape[0]
        err = inf['ERR'].data
        time = inf['TIME'].data
        dq = inf['DQ'].data

    return input_table

# 	if scnrat gt 0 then strput,file,'ima',pos-3

# if scnrat gt 0 then begin		; trim crapola for scanned data
# 	image=image(5:1018,5:1018)
# 	err=err(5:1018,5:1018)
# 	time=time(5:1018,5:1018)
# 	dq=dq(5:1018,5:1018)
# ; Clean false dq for 2048-signal in 0-read & 8192-CR detected:
# 	dq=dq-(dq and (2048+8192))
# 	endif

# if strtrim(sxpar(h,'imagetyp'),2) eq 'FLAT' then begin	; 05nov15 rcb
# 	good=where(time gt 0)
# 	image(good)=image(good)/time(good)	; flats not div by exptm
# 	err(good)=err(good)/time(good)		; 05nov28
# 	endif

# sxaddpar,h,'xactual',xco,'Actual Found direct image X Position'
# sxaddpar,h,'yactual',yco,'Actual Found direct image Y Position'
# sxaddpar,h,'xcamerr',xerr,'Pointing error XACTUAL-PREDICT (px)'
# sxaddpar,h,'ycamerr',yerr,'Pointing error YACTUAL-PREDICT (px)'

# if n_elements(crval1) gt 0 then begin
# 	extast,h,astr,status	;get astrometry information
# 	if status eq -1 then begin
# 		print,'CALWFC_SPEC: Error - Invalid astrometry ' + $
# 		    		'info in the image '+file
# 		retall
# 		endif
# 	ad2xy,crval1,crval2,astr,x1,y1   ;center position of targ acq
# 						 ;image in the current image

# 	refpx1=sxpar(h,'crpix1')-1	; minus 1 for IDL
# 	refpx2=sxpar(h,'crpix2')-1
# 	xdither = x1[0]-refpx1		;  pos from direct image
# 	ydither = y1[0]-refpx2		; 506,506 is usual ref px for x1,y1
# 	print,'calwfc_spec grsm x,y-dither from dir image '+	$
# 		'at'+string([refpx1,refpx2],'(2i4)')+' in px=',xdither,ydither
#      end else begin
# 	xdither = 0.0
# 	ydither = 0.0
# 	endelse
# fit_slope = 1				;default is to fit the slope
# if keyword_set(slope) then begin			;06jun28-rcb
# 	fit_slope = 0		; do not fit slope
# 	if slope ne 1 then angle = slope	; specified slope=slope
# 	endif

# ; ################################SCANNED################################
# ; 
# ;SCANNED SPECTRA have no direct image. Read SPT file, get scan rate & est.xc,yc
# ;
# if scnrat gt 0 then begin		; start trailed processing
# 	sxaddpar,h,'scan_rat',scnrat,'Trail scan rate'
# 	xpostarg=sxpar(h,'postarg1')
# ;find zero order @ y=559. SCANNED ONLY, as locates Z-order in 1-D (X) only.
# 	ycmean=559  &  ytop=926		;G141 means from scnoffset.pro
# 	xzappr=303.6+ 7.336*xpostarg  ;G141 scnoffset 0-ord pos @y=558.8
# 	if filter eq 'G102' then zxappr=000
# 	bn=rebin(image(*,559-20:559+20),1014,1)	; binned spectrum @y=559
# 	plot,bn,xr=[xzappr-60,xzappr+60]
# 	oplot,[xzappr,xzappr],[0,1e4],line=2
# 	xind=dindgen(1014)
# 	mx=where(bn(xzappr-50:xzappr+50) eq 			$
# 			max(bn(xzappr-50:xzappr+50),mxpos))
# 	mxpos=fix(mxpos+xzappr-50+.5)
# 	good=where(bn(mxpos-50:mxpos+50) gt 0.5*bn(mxpos))+mxpos-50
# 	xzfound=total(xind(good)*bn(good))/total(bn(good))
# 	oplot,[xzfound,xzfound],[0,1e4]
# 	xyouts,.15,.2,'Xfound,at y-posn='+			$
# 			string([xzfound,ycmean],'(f6.1,i4)'),/norm
# ;	read,st
# ;find zero order @ y=926 (G141)
# 	bn=rebin(image(*,ytop-20:ytop+20),1014,1) ;binned spectrum @ytop
# 	xzfix=fix(xzfound+.5)
# 	plot,bn,xr=[xzfix-60,xzfix+60]
# 	oplot,[xzfix,xzfix],[0,1e4],line=2
# 	mx=where(bn(xzfix-50:xzfix+50) eq 			$
# 			max(bn(xzfix-50:xzfix+50),mxpos))
# 	mxpos=mxpos+xzfix-50
# 	good=where(bn(xzfix-50:xzfix+50) gt 0.5*bn(mxpos))+xzfix-50
# 	xztop=total(xind(good)*bn(good))/total(bn(good))
# 	oplot,[xztop,xztop],[0,1e4]
# ;help,xzfix,mxpos,xztop  &	stop
# 	xyouts,.15,.2,'Xztop,at y-posn='+			$
# 			string([xztop,ytop],'(f6.1,i4)'),/norm
# 	z0coef0=xzfound			   ;coef-0 for 0-ord xz vs yz
# 	z0coef1=(xztop-xzfound)/(ytop-559) ;coef-1 for del-y, ie (y-559)
# 	print,'coef for 0-order vs (Y-yzcent)',z0coef0,z0coef1
# ;	read,st
# ;amount to add to found x of 0-order for G141 at ycmean=559 from scnoffset.pro
# 	delx=187.63 + 0.0121233*xpostarg
# 	xc=xzfound+delx
# 	yc=ycmean				;avg y for center
# 	if filter eq 'G102' then begin
# 		xc=0	; run scnoffset.pro
# 		stop
# 		endif
# ; x is range of spec
# 	calwfc_spec_wave,h,xc,yc,x,wave,angle,wav1st	;master WLs wfc_wavecal?
# 		
# ; Differences btwn 0-order X znd xc posit of ref * from scnoffset.pro:
# 	delxbot=191.91+0.0096832*xpostarg	;G141 @ y=184
# 	delxtop=183.71+0.0134596*xpostarg	;G141 @ y=926
# 	yzord=[184.,559,926]
# 	xdelx=[delxbot,delx,delxtop]
# ; make wlimg for use in flat fielding
# 	wlimg=image*0.
# 	for i=0,1013 do begin
# 		xz=z0coef0+z0coef1*(i-ycmean)		; 0-ord x-posn
# 		del=interpol(xdelx,yzord,i)		; do all rows
# 		xref=xz+del				; ref * posn
# ; WLs @ y=i
# 		calwfc_spec_wave,h,xref,i,x,wli,dum,wl1st,/noprnt
# 		wlimg(*,i)=wl1st			;1-ord WL for FF
# 		endfor
# 	if strupcase(flatfile) ne 'NONE' then $
# 		calwfc_spec_flat,h,image,err,dq,x,wlimg,flatfile      ; FF image
# 	disp=24.5					; G102
# 	if filter eq 'G141' then disp=46.5
# ; make new image where each row is resampled to WL scale of avg Y position
# 	rectim=image*0.		; image rectified for disp & angle
# 	recter=rectim  &  rectdq=rectim
# 	for i=0,1013 do begin
# 		xz=z0coef0+z0coef1*(i-ycmean)		; 0-ord x-posn
# 		yfit=i+(x-xz)*sin(angle/!radeg)		;trace @ row=i
# 		yfit=round(yfit)		       ;nearest neighbor
# 		yfit=yfit>0<1013
# 		del=interpol(xdelx,yzord,i)		; do all rows
# 		xref=xz+del				; ref * posn
# 		calwfc_spec_wave,h,xref,i,x,wli,/noprnt	; full WLs @ y=i
# ; use disp of -1st order for corr; but the +1st change w/ Y is the same to <0.1%
# 		indx=where(wli ge -12000,npts)
# 		indx=indx(0)
# 		dsprow=wli(indx+1)-wli(indx)
# ; Accomodate slope of spectra to get a few more good rows at lowest good Y
# 		rowspc=image(x,yfit)*disp/dsprow      ;corr to +1st disp
# ; row on master WL scale:
# 		linterp,wli,rowspc,wave,spcterp,missing=0
# 		rectim(*,i)=spcterp
# 	if i eq 511 then stop
# 		rowerr=err(x,yfit)*disp/dsprow
# 		linterp,wli,rowerr,wave,errterp,missing=0
# 		recter(*,i)=errterp
# 		rowdq=dq(x,yfit)
# ; Expand dq to near neighbors & use dq=2 for missing, ie fill data:
# 		linterp,wli,rowdq,wave,dqnotrp,missing=2,/nointerp
# 		linterp,wli,rowdq,wave,dqterp,missing=2
# 		bad=where(dqterp ne 0,nbad)
# 		if nbad gt 0 then for ibad=0,nbad-1 do begin
# 			indx=bad(ibad)
# 			dqterp(indx)=max(dqnotrp((indx-1)>0:(indx+1)<1013))
# 			endfor
# 		rectdq(*,i)=dqterp
# 		endfor
# ; 2018may -see above now	fdecomp,file,disk,dir,root
# ;	root = gettok(root,'_')
# 
# 	sxaddhist,'Written by IDL calwfc_spec.pro '+!stime,h
# 	exptim=sxpar(h,'exptime')
# 	rectim=rectim*exptim*scnrat/0.13	;Texp=plate-scale/scnrat
# 	recter=recter*exptim*scnrat/0.13
# 	ext='fits'
# 	if strupcase(flatfile) eq 'NONE' then begin
# 		ext='_noff.fits
# 		sxaddhist,'z.scimage NOT flat-fielded.',h
# 		endif
# 	mwrfits,{wave:wave,scimage:rectim,err:recter,dq:rectdq},  $
# 					subdir+'/imag_'+root+ext,h,/create
# ;	read,st
# 	return						; END SCAN SECT
# 	endif


def reduce_stare(row, params, arg_list):
    """
    Reduces stare-mode grism data
    """
    verbose = arg_list.verbose
    interactive = arg_list.trace
    bkg_flat_order = arg_list.bkg_flat_order
    task = "wfc3_grism_reduce_stare"
    file = os.path.join(row["path"], row["filename"])
    root = row['root']
    target = row['target']
    ref = row['filter_root']
    filter = row['filter']
    preamble = '{}: {}'.format(task, root)

    np_formatter = {'float_kind':lambda x: "{:8.2f}".format(x)}
    np_opt = {'max_line_width': 175, 'formatter': np_formatter,
              'threshold': 2000000}
    
    if verbose:
        print("{}: starting {}.".format(preamble, row['filename']))

    with fits.open(file) as inf:
        image = inf['SCI'].data
        sci_hdr = inf['sci'].header
        filter = row['filter']
        xsize, ysize = image.shape[1], image.shape[0]
        err = inf['ERR'].data
        time = inf['TIME'].data
        dq = inf['DQ'].data
        ltv1, ltv2 = inf[1].header['LTV1'], inf[1].header['LTV2']
        
        if inf[0].header["IMAGETYP"] == "FLAT":
            if verbose:
                print("{}: IMAGETYP=FLAT".format(preamble))
            good = np.where(time>0.)
            image[good] = image[good] / time[good]
            err[good] = err[good] / time[good]
        
        img_wcs = wcs.WCS(inf['SCI'].header)
#         hdr = img_wcs.to_header()
#         for item in hdr.tostring(sep='\\n').split('\\n'):
#             print(item)
        ra, dec = float(row['crval1']), float(row['crval2'])
        targ = img_wcs.wcs_world2pix([ra], [dec], 0, ra_dec_order=True)
        x1, y1 = targ[1][0], targ[0][0]
        refpx1, refpx2 = inf[1].header['crpix1']-1, inf[1].header['crpix2']-1
        xdither, ydither = x1-refpx1, y1-refpx2
        if verbose:
            msg = "{}: x,y-dither from dir image at ({},{}) is ({},{}) px."
            print(msg.format(preamble, refpx1, refpx2, xdither, ydither))
        if (params['xc'] < 0 and 'xc' in params['set']) or \
           (params['yc'] < 0 and 'yc' in params['set']):
            msg = " No/Bad result in locating target from filter image."
            row['notes'] += msg
            for item in ['xc', 'yc']:
                if item in params['set']:
                    params['set'].remove(item)
        if 'xc' not in params['set'] or 'yc' not in params['set']:
            if verbose:
                msg = "{}: {}: manually searching for target image."
                print(msg.format(task, row['root']))
            tra = row['ra_targ']
            tdec = row['dec_targ']
            targ = row['target']
            targ = wcs.wcs_world2pix([tra, tdec], 0)
            ix, iy = targ[1], targ[0]
            xerr, yerr = 0., 0.
        else:
            ix, iy = params['xc'], params['yc']
            xerr, yerr = params['xerr'], params['yerr']
            
        if verbose:
            msg = "{}: {}: Predicted target image position ({},{})"
            print(msg.format(task, row['root'], ix, iy))
        ix += xerr
        iy += yerr
        if verbose:
            msg = "{}: {}: Including error, predicted position is ({},{})"
            print(msg.format(task, row['root'], ix, iy))
        xc, yc = ix - params['ix_shift'], iy - params['iy_shift']
        xastr, yastr = xc, yc
        
        if (xc < 4 or yc < 4 or xc > 998 or yc > 998):
            if verbose:
                msg = "{}: Predicted Z-ord position ({},{}) is off grism."
                print(msg.format(preamble, xc, yc))
                msg = "\tDistortion-sensitive measured position ({},{})"
                print(msg.format(xdither, ydither))
                msg = "\tUse Predicted position ({},{}) per astrometry"
                print(msg.format(ix, iy))
            wave_params = {
                            'xin': ix,
                            'yin': iy,
                            'wl_offset': params['wl_offset'],
                            'verbose': verbose
                          }
            x_arr, wave, angle, wav1st = reduce_wave_axe(row, wave_params)
            if verbose:
                msg = '{}: WAVELENGTH solution per AXE coef. No Z-order.'
                print(msg.format(preamble))
            axeflg = True
        else:
            if verbose:
                msg = '{}: 1st Z-ord centroid from astrom+ Petro starts at '
                msg += '({},{})'
                print(msg.format(preamble, xc, yc))

            xbeg, ybeg = int(round(xc)), int(round(yc))
            xb, xe = max(xbeg-20, 0), min(xbeg+20, image.shape[1])
            yb, ye = max(ybeg-20, 0), min(ybeg+20, image.shape[0])
            sbimg = image[yb:ye,xb:xe]
            # Threshold of 10 counts, FWHM of 2. Only return brightest result.
            star_finder = DAOStarFinder(10., 2., brightest=1)
            star_table = star_finder.find_stars(sbimg)
            star_x = star_table['xcentroid'][0] + xb
            star_y = star_table['ycentroid'][0] + yb
            if verbose:
                print("Testing centroiding methods...")
                print("\tDAOFind: xc={}, yc={}".format(star_x, star_y))

            np_formatter = {'float_kind':lambda x: "{:10.6f}".format(x)}
            star_opt = {'max_line_width': 375, 'formatter': np_formatter,
                        'threshold': 2000000}

            ywidth = params['ywidth']
            ns = 31
            ibeg = max(xbeg-ns//2, 0)
            iend = xbeg+ns//2
            xp, yp = centroid_sources(image, xc, yc, box_size=ns, centroid_func=g1)
            if verbose:
                print("\t1d Gaussian: xc,yc={},{}".format(xp, yp))
            xp, yp = centroid_sources(image, xc, yc, box_size=ns, centroid_func=g2)
            if verbose:
                print("\t2d Gaussian: xc,yc={},{}".format(xp, yp))
            xp, yp = centroid_sources(image, xc, yc, box_size=ns, centroid_func=com)
            if verbose:
                print("\t2d Moments: xc,yc={},{}".format(xp, yp))

            xc, yc = star_x, star_y
#             xc, yc = xpos, ypos
            if xc < ibeg or xc > iend:
                if verbose:
                    msg = "{}: Peak too close to edge. Use approx. pos ({},{})"
                    print(msg.format(preamble, xbeg, ybeg))
                xc, yc = xbeg, ybeg


            star_x, star_y = int(round(star_x)), int(round(star_y))
            
            wave_params = {
                            'xin': xc,
                            'yin': yc,
                            'wl_offset': params['wl_offset'],
                            'verbose': verbose
                          }
            x_arr, wave, angle, wav1st = reduce_wave_zord(row, wave_params)
            if verbose:
                msg = "{}: Zero-ord centroid at ({},{})"
                print(msg.format(preamble, xc, yc))
                print("\tSearch area {}X{} pix".format(ns, ns))
                print("\tX,Y astrom error=({},{})".format(xastr-xc, yastr-yc))
            axeflg = False

        fit_slope = True
        if 'slope' in params['set']:
            fit_slope = False
            angle = np.radians(params['slope'])
        
        if interactive:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            targ_str = "{} - {} ({})".format(root, target, filter)
            ax.set_title('{} Direct Image with Detected Source'.format(targ_str))
            plt.imshow(np.log10(np.where(image>=0.1,image,0.1)))
            plt.plot([xc, xc], [yc-8, yc-18], color='white')
            plt.plot([xc, xc], [yc+8, yc+18], color='white')
            plt.plot([xc-8, xc-18], [yc, yc], color='white')
            plt.plot([xc+8, xc+18], [yc, yc], color='white')
            plt.show()
        
        raw_image = deepcopy(image)
        
        # ***** Start of wfc_flatscal.pro
        #   Move to its own function?
        # *****

        # Get flatfield data cube
        if verbose:
            print("{}: Beginning flatfield.".format(preamble))
        cal_data = get_data_file("abscal.wfc3", "calibration_files.json")
        with open(cal_data, 'r') as inf:
            cal_files = json.load(inf)
        flat_file_name = cal_files["flatfield"]["wfc3"][row["filter"]]
        flat_file = get_data_file("abscal.wfc3", flat_file_name)
        if verbose:
            msg = "{}: wfc_flatscl: Using flatfield {}"
            print(msg.format(preamble, flat_file_name))

        with fits.open(flat_file) as in_flat:
            flat = in_flat[0].data
        ns, nl = image.shape[0], image.shape[1]
        flat = flat[int(ltv1):int(ltv1)+ns,int(ltv2):int(ltv2)+nl]
        xpos = np.arange(ns, dtype='int32')
        sclfac = nl/1014.
        sc_low, sc_high = int(round(100*sclfac)), int(round(900*sclfac))+1
        xmin, xmax = min(xpos), max(xpos)+1

        # Collapse and average along x axis
        imav = np.mean(image[:,xmin:xmax], axis=1)
        bkav = np.mean(flat[:,xmin:xmax], axis=1)
        # Remove spikes
        smimav = ndimage.median_filter(imav, size=61, mode='nearest')
        smimav[:30] = imav[:30]
        smimav[-30:] = imav[-30:]
        smbkav = ndimage.median_filter(bkav, size=21, mode='nearest')
        smbkav[:20] = bkav[:20]
        smbkav[-20:] = bkav[-20:]
        rat = smimav/smbkav
        avrat = np.mean(rat[sc_low:sc_high])
        sigma = np.std(rat[sc_low:sc_high])
        avgbkg = np.mean(smimav[sc_low:sc_high])*params['gwidth']
        if verbose:
            msg = "{}: wfc_flatscl: Avg bkg = {} +/- {} for gwidth={}"
            print(msg.format(preamble, avgbkg, sigma, params['gwidth']))
            msg = "{}: wfc_flatscl: Flatfield Scale Ratio = {}"
            print(msg.format(preamble, avrat))
        flat *= avrat
        if verbose:
            print("{}: Finished flatfield".format(preamble))
        # ***** End of wfc_flatscal.pro
        
        image -= flat

        if len(flat) > 1:
            sky_image = raw_image/(np.where(flat>0.001, flat, 0.001))
        else:
            sky_image = np.zeros_like(raw_image)
        
        if bkg_flat_order == 'flat_first':
            ff_params = {
                            "root": row['root'][0],
                            "image": image,
                            "hdr": sci_hdr,
                            "flat": flat,
                            "wave": wav1st,
                            "ltv1": ltv1,
                            "ltv2": ltv2,
                            "filter": filter,
                            "verbose": verbose
                        }
            image, flatfile = reduce_flatfield(row, ff_params)            
         
        # Extract
        y_offset = params['y_offset']
        y_shift = params['yshift']
        yapprox = yc + y_offset + (x_arr-xc)*np.sin(angle) + y_shift

        image_edit_file = get_data_file("abscal.wfc3", "image_edits.json")
        with open(image_edit_file, 'r') as inf:
            image_edits = json.load(inf)
        issues = []
        if "reduce" in image_edits:
            issues = image_edits["reduce"]
        overrides = params.get('overrides', {})
        img = {"sci": image, "err": err, "dq": dq}
        img = set_image(img, row, issues, preamble, overrides, verbose)
        image, err, dq = img["sci"], img["err"], img["dq"]
        
        fimage = deepcopy(image)

        for i in range(fimage.shape[0]):
            fimage[i,:] = ndimage.median_filter(image[i,:], size=7, mode='nearest')
            fimage[i,:4] = image[i,:4]
            fimage[i,-4:] = image[i,-4:]

        if fit_slope:
            x1bin, x2bin = {}, {}
            xfound, yfound = {0: xc}, {0: yc}
            if axeflg:
                yfound[0] = -yc
            
            for iord in [-1, 1, 2]:
                if verbose:
                    print("{}: Starting iord={}".format(preamble, iord))
                wlrangs = {
                            -1: ['wlrang_m1_low', 'wlrang_m1_high'],
                             1: ['wlrang_p1_low', 'wlrang_p1_high'],
                             2: ['wlrang_p2_low', 'wlrang_p2_high'],
                          }
                wlrang = wlrangs[iord]
                wl_low, wl_high = iord*params[wlrang[0]], iord*params[wlrang[1]]
                
                if iord > 0:
                    xrang = np.where((wave >= wl_low) & (wave <= wl_high))[0]
                else:
                    xrang = np.where((wave >= wl_high) & (wave <= wl_low))[0]
                npts = len(xrang)
                
                # if no points or minimum too high, order is too close to the
                #   edge, so set npts to zero.
                if (npts == 0) or (xrang[0] >= 985):
                    npts = 0
                
                if npts >= 50:
                    x1bin[iord] = xrang[0]
                    x2bin[iord] = min(max(xrang), 990)
                    xfound[iord] = np.mean(xrang)
                    
                    x1, x2 = x1bin[iord], x2bin[iord]
                    y1 = max(int(round(yapprox[(x1+x2)//2]))-ywidth//2, 0)
                    y2 = min((y1 + params['ywidth']-1), (ysize-1))
                    if verbose:
                        msg = "{}: iord={}; x,y search range=({},{}), ({},{})"
                        print(msg.format(preamble, iord, x1, x2, y1, y2))
                    y1pred = 0
                    if iord == 1:
                        y1pred = yapprox[(x1+x2)//2]
                    elif iord == 2:
                        ydel = yfound[iord-1]
                        ydel -= yapprox[(x1bin[iord-1]+x2bin[iord-1])//2]
                        y1 = int(round(y1+ydel))
                        y2 = int(round(y2+ydel))
                        if verbose:
                            msg = "\t2nd ord yposit tweak yrange={},{} by {}"
                            print(msg.format(y1, y2, ydel))
                    ny = y2 - y1 + 1

                    profile = np.zeros((ny,), dtype='float64')
                    iter = 0
                    iterate = True
                    while iterate:
                        iter += 1
                        for j in range(ny):
                            profile[j] = np.sum(fimage[y1+j,x1:x2+1])
                        
                        profile -= np.median(profile)
                        profile = np.where(profile>0., profile, 0.)
                        
#           ; patch for PN-G045.4-02 w/ low contin. & another spec just above:
# 			if keyword_set(star) then begin
#               if star eq 'PN' then begin
# 				    profile(0)=0  &  profile(ny-1)=0  
#               endif
#           endif
                        yprofile = np.arange(ny, dtype='int16') + y1
                        pmax, maxpos = np.max(profile), np.argmax(profile)
                        if pmax <= 0:
                            err = "{}: ERROR in Profile ".format(preamble)
                            err += "Iteration at {}. Maximum ".format(iter)
                            err += "value of profile is negative."
                            err += "Profile={}".format(profile)
                            raise ValueError(err)
                        if verbose:
                            msg = "{}: Profile maximum {} at {} ({})"
                            print(msg.format(preamble, pmax, maxpos, yprofile[maxpos]))
                        if maxpos <= 0 or maxpos > ny or pmax <= 2:
                            yfound[iord] = -999
                        else:
                            yp = yprofile[maxpos-1:maxpos+2]
                            pf = profile.astype('float64')[maxpos-1:maxpos+2]
                            yfound[iord] = np.sum(yp*pf)/np.sum(pf)
                        if verbose or interactive:
                            xf, yf = xfound[iord], yfound[iord]
                            msg = "{}: order {}, ".format(preamble, iord)
                            msg += "contin. x,y position=({},{})".format(xf, yf)
                            print(msg)
                            
                        if iord == 1 and yfound[iord] == -999 and iter == 1:
                            if maxpos == 0:
                                y1 = max((y1 - (params['ywidth']-3)), 0)
                            else:
                                guess = (y1 + params['ywidth']-3)
                                ceiling = (ysize - params['ywidth'] - 1)
                                y1 = min(guess, ceiling)
                            y2 = min((y1 + params['ywidth'] - 1), (ysize-1))
                            if verbose:
                                msg = "{}: Iterating with search ({},{})"
                                print(msg.format(preamble, y1, y2))
                        else:
                            iterate = False

#           ; std above continuum technique fails for PN em line and faint -1 order, SO,
#           ;		try the brightest px in the range:
# 			if iord eq -1 and star eq 'PN' then begin
# 				imax=where(image(x1:x2,y1:y2) eq max(image(x1:x2,y1:y2)))
# 				nsamp=fix(x2-x1+1)
# 				xfound(nbins)=x1+imax mod nsamp
# 				yfound(nbins)=y1+imax/nsamp
# 				print,'PN em line at ',xfound(nbins), yfound(nbins),' for ordr=-1'
# 			endif
                    
                    # Replaces the separate trace flag
                    if interactive:
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        targ_str = "{} - {} ({})".format(root, target, filter)
                        ax.set_title('{}: Trace Image'.format(targ_str))
                        plt.xlabel("Y (column) Pixel (px)")
                        plt.ylabel("Signal")
                        ax.scatter(yprofile, profile)
                        ax.text(.15, .88, 'order={}'.format(iord), 
                                transform=ax.transAxes)
                        plt.show()
                # end if npts > 50
            # end iord loop
            
            # Initially don't allow negative found positions, and don't
            #   count a result as good if aXe was used instead of zeroth-order
            good = []
            for iord in [-1, 0, 1, 2]:
                if iord in xfound and iord in yfound:
                    if xfound[iord] >= 4 and yfound[iord] >= 0:
                        good.append(iord)
            if verbose:
                print("{}: Initial good orders are {}".format(preamble, good))

            if len(good) == 1:
                if verbose:
                    msg = "{}: Only one order. ".format(preamble)
                    msg += "Try using predicted Z-order position. Adjust "
                    msg += "Z-order for Y error in first order. "
                    msg += "Found={}, pred={}, found-pred={}."
                    print(msg.format(yfound[1], y1pred, yfound[1]-y1pred))
                    msg = "\tWLs differ by <1A for corr. dir img at {},{}"
                    print(msg.format(ix, yfound[0]))
                yfound[0] = abs(yfound[0]) + (yfound[1]-y1pred)
            
            # This time allow negative found positions, and allow the zeroth
            #   order even via aXe.
            good = []
            for iord in [-1, 0, 1, 2]:
                if iord in xfound and iord in yfound:
                    if xfound[iord] >= -500 and yfound[iord] != 0 and yfound[iord] != -999:
                        good.append(iord)
            if verbose:
                print("{}: Final good orders are {}".format(preamble, good))

            if len(good) <= 1:
                msg = "{}: Only one order. ".format(preamble)
                msg += "Cannot measure spectral angle for subarray of "
                msg += "({},{}). Aborting.".format(xsize, ysize)
                raise ValueError(msg)
            
            if (len(good) >= 3) and axeflg:
                good = []
                for iord in [-1, 0, 1, 2]:
                    if iord in xfound and iord in yfound:
                        if yfound[iord] > 0:
                            good.append(iord)
                if verbose:
                    msg = "{}: {}: PROFILE yfound ".format(task, row['root'])
                    msg += "for AXE case={}".format(yfound)
                    print(msg)
            else:
                for iord in [-1, 0, 1, 2]:
                    if iord in yfound:
                        yfound[iord] = abs(yfound[iord])
            
            poly_x = np.array([xfound[i] for i in good])
            poly_y = np.array([yfound[i] for i in good])
            msg = "{}: X and Y found values: {}, {}"
            print(msg.format(preamble, poly_x, poly_y))
            coef = np.polyfit(poly_x, poly_y, 1)
            fit_angle = np.arctan(coef[0])*180./np.pi
            print("{}: Found coefficients {} angle {}".format(preamble, coef, fit_angle))
            if interactive:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                t_axes = ax.transAxes
                targ_str = "{} - {} ({})".format(root, target, filter)
                ax.set_title('{}: Fit Plot'.format(targ_str))
                plt.xlabel("X-pixel (px)")
                plt.ylabel("Y-pixel (px)")
                ax.scatter(poly_x, poly_y, label="Found Orders")
                ax.plot(x_arr, yapprox, label="Initial Trace Approximation")
                y_fit = coef[1] + coef[0]*x_arr
                ax.plot(x_arr, y_fit, label="Fitted Trace")
                ax.text(.7, .15, 'angle={:.3f}'.format(fit_angle), transform=t_axes)
                ax.legend()
                plt.show()
                if verbose:
                    msg = "{}: {}:, angle={}"
                    print(msg.format(preamble, filter, fit_angle))
        else:
            # Don't fit slope
            profile = np.zeros((ywidth,), dtype='float64')
            yoff = np.arange(ywidth, dtype='float64') - ywidth/2
            for i in range(ywidth):
                x_axis = np.arange(fimage.shape[0], dtype='int16')
                y_axis = np.arange(fimage.shape[1], dtype='int16')
                interp = interp2d(x_axis, y_axis, fimage, kind='linear')
                profile = np.sum(interp(x_arr, yapprox+yoff[i]))
            profile -= np.median(profile)
            profile -= np.max(profile)/4
            profile = np.where(profile>0, profile, 0.)
            ycent = np.sum(profile*yoff)/np.sum(profile)
            ypos = yapprox + ycent
            coef = np.polyfit(x, ypos, 1)
        
        if verbose:
            print("{}: Coef of trace fit={}".format(preamble, coef))
        yfit = coef[1] + coef[0]*x_arr
        
        # Extract Background Spectra
        if verbose:
            print("{}: Extract background spectra".format(preamble))
        ns = len(x_arr)
        nsb = ns
        x_b = deepcopy(x_arr)
        b_lower = np.zeros((nsb,), dtype='float64')
        b_upper = np.zeros((nsb,), dtype='float64')
        sky_b_lower = np.zeros((nsb,), dtype='float64')
        sky_b_upper = np.zeros((nsb,), dtype='float64')
        yfit_b = coef[1] + coef[0]*x_b
        
        half_bwidth = params['bwidth']//2
        for i in range(nsb):
            ypos = int(round(yfit_b[i] - params['lbdist']))
            y1 = max((ypos - half_bwidth), 0)
            y2 = min((ypos + half_bwidth + 1), image.shape[0])
            b_lower[i] = np.median(image[y1:y2,x_b[i]])
            sky_b_lower[i] = np.median(sky_image[y1:y2,x_b[i]])
            
            ypos = int(round(yfit_b[i] + params['ubdist']))
            y1 = max((ypos - half_bwidth), 0)
            y2 = min((ypos + half_bwidth + 1), image.shape[0])
            b_upper[i] = np.median(image[y1:y2,x_b[i]])
            sky_b_upper[i] = np.median(sky_image[y1:y2,x_b[i]])
        
        # Average upper and lower background, and smooth.
        raw_sky_back = (sky_b_lower + sky_b_upper)/2.

        bmedian = params['bmedian']
        half_bmedian = bmedian//2
        if bmedian > 1:
            sky_back = ndimage.median_filter(raw_sky_back, size=bmedian,
                                             mode='nearest')
            sky_back[:half_bmedian] = raw_sky_back[:half_bmedian]
            sky_back[-half_bmedian:] = raw_sky_back[-half_bmedian:]
        else:
            sky_back = raw_sky_back        

        if params['bmean1'] > 1:
            sky_back = convolve(sky_back, Box1DKernel(params['bmean1']))
        if params['bmean2'] > 1:
            sky_back = convolve(sky_back, Box1DKernel(params['bmean2']))
        
        back = (b_lower + b_upper)/2
        
        s_back_lo = ndimage.median_filter(b_lower, size=bmedian, mode='nearest')
        s_back_lo[:half_bmedian] = b_lower[:half_bmedian]
        s_back_lo[-half_bmedian:] = b_lower[-half_bmedian:]

        s_back_up = ndimage.median_filter(b_upper, size=bmedian, mode='nearest')
        s_back_up[:half_bmedian] = b_upper[:half_bmedian]
        s_back_up[-half_bmedian:] = b_upper[-half_bmedian:]
        
        x_b_f = x_b.astype('float64')

        # Do a two-pass iteration
        
        # First pass
        lo_res = np.polyfit(x_b_f, s_back_lo, 3, full=True)
        lo_coef = lo_res[0]
        lo_error = np.sqrt(lo_res[1]/(len(x_b_f)-2))[0]
        p = np.poly1d(lo_coef)
        lo_fit = p(x_b_f)

        up_res = np.polyfit(x_b_f, s_back_up, 3, full=True)
        up_coef = up_res[0]
        up_error = np.sqrt(up_res[1]/(len(x_b_f)-2))[0]
        p = np.poly1d(up_coef)
        up_fit = p(x_b_f)
        
        sigless = min(lo_error, up_error)
        if verbose:
            msg = "{}: loerr={}, uperr={}, sigless={}"
            print(msg.format(preamble, lo_error, up_error, sigless))
        
        # Second Pass
        good_lo = np.where(np.abs(s_back_lo-lo_fit) < sigless)
        n_good_lo = len(good_lo[0])
        if verbose:
            msg = "{}: {} of {} low bkg points have residual < standard error."
            print(msg.format(preamble, n_good_lo, len(lo_fit)))
        lo_res = np.polyfit(x_b_f[good_lo], s_back_lo[good_lo], 3, full=True)
        lo_coef = lo_res[0]
        lo_error = np.sqrt(lo_res[1]/(n_good_lo-2))[0]
        if verbose:
            msg = "{}: Lower fit has error {} on second pass."
            print(msg.format(preamble, lo_error))
        p = np.poly1d(lo_coef)
        lo_fit = p(x_b_f)

        good_up = np.where(np.abs(s_back_up-up_fit) < sigless)
        n_good_up = len(good_up[0])
        if verbose:
            msg = "{}: {} of {} high bkg points have residual < standard error."
            print(msg.format(preamble, n_good_up, len(up_fit)))
        up_res = np.polyfit(x_b_f[good_up], s_back_up[good_up], 3, full=True)
        up_coef = up_res[0]
        up_error = np.sqrt(up_res[1]/(n_good_up-2))[0]
        if verbose:
            msg = "{}: Upper fit has error {} on second pass."
            print(msg.format(preamble, up_error))
        p = np.poly1d(up_coef)
        up_fit = p(x_b_f)
        
        lo_slp = lo_coef[2] + 2*lo_coef[1]*x_b_f + 3*lo_coef[0]*x_b_f*x_b_f
        up_slp = up_coef[2] + 2*up_coef[1]*x_b_f + 3*up_coef[0]*x_b_f*x_b_f

        if verbose:
            msg = "{}: {} points, {} good low points, {} good high points"
            print(msg.format(preamble, nsb, n_good_lo, n_good_up))
        
        if np.mean(up_fit) < np.mean(lo_fit):
            if verbose:
                print("{}: Setting initial fit to upper fit".format(preamble))
            s_back = up_fit
        else:
            if verbose:
                print("{}: Setting initial fit to lower fit".format(preamble))
            s_back = lo_fit

        bettr = np.where((abs(up_fit)-abs(s_back) > 3*sigless) & (up_fit < lo_fit))
        if len(bettr[0]) > 0:
            if verbose:
                msg = "{}: Setting {} points to upper fit"
                print(msg.format(preamble, len(bettr[0])))
            s_back[bettr] = up_fit[bettr]
        bettr = np.where((abs(lo_fit)-abs(s_back) > 3*sigless) & (lo_fit < up_fit))
        if len(bettr[0]) > 0:
            if verbose:
                msg = "{}: Setting {} points to lower fit"
                print(msg.format(preamble, len(bettr[0])))
            s_back[bettr] = lo_fit[bettr]
        
        if len(good_lo[0]) < nsb/5:
            if verbose:
                print("{}: Setting fit to upper fit.".format(preamble))
            s_back = up_fit
        if len(good_up[0]) < nsb/5:
            if verbose:
                print("{}: Setting fit to lower fit.".format(preamble))
            s_back = lo_fit
        
        if interactive:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            targ_str = "{} - {} ({})".format(root, target, filter)
            ax.set_title("{} Background Fit".format(targ_str))
            t_axes = ax.transAxes
            ax.scatter(x_b_f[good_lo], s_back_lo[good_lo], label="Low Points")
            ax.plot(x_b_f, lo_fit, label="Low Fit")
            ax.scatter(x_b_f[good_up], s_back_up[good_up], label="High Points")
            ax.plot(x_b_f, up_fit, label='High Fit')
            ax.plot(x_b_f, s_back, label='Final Fit')
            ax.legend()
            plt.show()
        
        if (len(good_lo[0]) < nsb/5) and (len(good_up[0]) < nsb/5):
            raise ValueError("Not enough good points to fit.")
        if abs(np.mean(s_back)) > 3:
            raise ValueError("Average absolute fit is greater than 3(?)")
        
        # Full image subtraction
        for i in range(nsb):
            image[:,x_b[i]] = image[:,x_b[i]] - s_back[i]
        
        s_back = s_back * params['gwidth']
        sky_back = sky_back * params['gwidth']
        back = back * params['gwidth']
        raw_sky_back = raw_sky_back * params['gwidth']
        
        # Make Background results for optional output
        br_gross = rebin(raw_image[:,x_b[0]:x_b[0]+nsb-1], (5, 256), np.sum)
        br_net = rebin(image[:,x_b[0]:x_b[0]+nsb-1], (5, 256), np.sum)
        back_results = {
                        'sback': rebin(s_back, 5, func=np.sum)/params['gwidth'],
                        'gross': np.transpose(br_gross),
                        'net': np.transpose(br_net),
                        'wave': rebin(wave, 5),
                        'bupper': rebin(b_upper, 5, np.sum),
                        'blower': rebin(b_lower, 5, np.sum),
                        'x': rebin(x_arr, 5),
                        'y': rebin(yfit_b, 5),
                        'ubdist': params['ubdist'],
                        'lbdist': params['lbdist']
                       }
        
        # Create trace image
        if interactive:
            imagt = image
            ytmp = yfit
            if np.max(yfit) > 700:
                imagt = image[500:1013,:]
                ytmp -= 500
            fig = plt.figure()
            ax = fig.add_subplot(111)
            targ_str = "{} - {} ({})".format(root, target, filter)
            ax.set_title('{} Trace Fit on Source'.format(targ_str))
            plt.imshow(np.log10(np.where(imagt>=0.1,imagt,0.1)))
            ax.plot(x_arr, ytmp+params['gwidth']//2)
            ax.plot(x_arr, ytmp-params['gwidth']//2)
            ax.plot(x_arr, ytmp+params['ubdist']+params['bwidth']//2)
            ax.plot(x_arr, ytmp+params['ubdist']-params['bwidth']//2)
            ax.plot(x_arr, ytmp-params['lbdist']+params['bwidth']//2)
            ax.plot(x_arr, ytmp-params['lbdist']-params['bwidth']//2)
            plt.show()

        if bkg_flat_order == 'bkg_first':
            ff_params = {
                            "root": row['root'][0],
                            "image": image,
                            "hdr": sci_hdr,
                            "flat": flat,
                            "wave": wav1st,
                            "ltv1": ltv1,
                            "ltv2": ltv2,
                            "filter": filter,
                            "verbose": verbose
                        }
            image, flatfile = reduce_flatfield(row, ff_params)
        
        # Extract spectra
        flux = np.zeros((ns,))
        epsf = np.zeros((ns,), dtype='int16')
        errf = np.zeros((ns,))
        spec_time = np.zeros((ns,))

        nn_img = np.where(image>0, image, 0)        
        hwidth = params['gwidth']//2
        for i in range(ns):
            y1 = max((yfit[i] - hwidth), 0)
            y2 = max(min((yfit[i] + hwidth), image.shape[0]), y1+1)
            iy1 = int(round(y1))
            iy2 = int(round(y2))
            
            for irow in range(iy1, iy2+1):
                if (dq[i,irow] & 24) and (x_arr[i] != 0) and (x_arr[i] != ns-1):
                    image[i,irow] = (image[i-1,irow] + image[i+1,irow])/2
            
            frac1 = 0.5 + iy1-y1 # frac of pixel i1y
            frac2 = 0.5 + y2-iy2 # frac of pixel iy2
            ifull1 = iy1 + 1     # range of full pixels to extract - low
            ifull2 = iy2         # range of full pixels to extract - high
            
            if ifull2 >= ifull1:
                tot = np.sum(image[ifull1:ifull2,x_arr[i]])
                var = np.sum(err[ifull1:ifull2,x_arr[i]])**2
                tt = time[ifull1:ifull2,x_arr[i]]
                nn = nn_img[ifull1:ifull2,x_arr[i]]
                tot_time = np.sum(tt*nn)
                time_weight = np.sum(nn_img[ifull1:ifull2,x_arr[i]])
            else:
                tot = 0.
                var = 0.
                time_weight = 0.
                tot_time = 0.
            tot += frac1*image[iy1,x_arr[i]] + frac2*image[iy2,x_arr[i]]
            var += frac1*err[iy1,x_arr[i]]**2 + frac2*err[iy2,x_arr[i]]**2
            tot_time += frac1*time[iy1,x_arr[i]]*nn_img[iy1,x_arr[i]] 
            tot_time += frac2*time[iy2,x_arr[i]]*nn_img[iy2,x_arr[i]]
            time_weight += frac1*nn_img[iy1,x_arr[i]] 
            time_weight += frac2*nn_img[iy2,x_arr[i]]
            if time_weight > 0:
                ave_time = tot_time/time_weight
            else:
                ave_time = np.max(time[iy1:iy2+1,x_arr[i]])
            e = 0
            for j in range(iy1, iy2+1):
                e = e | dq[j,x_arr[i]]
            
            flux[i] = tot
            errf[i] = np.sqrt(var)
            epsf[i] = e
            spec_time[i] = ave_time
        
        # Correct response to +1 order dispersion at y-centre
        # This code exists in the IDL, but is followed by a direct assignment of
        #   dcorr which, as far as I can see, makes the rest of it effectively
        #   obsolete.
#         indx = np.where(wave >= 10900)
#         npts = len(indx[0])
#         indx = indx[0][0]
#         if npts <= 10:
#             indx = np.where(wave <= -12000)
#             npts = len(indx[0])
#             indx = indx[0][-1]
#         disp = wave[indx+1] - wave[indx]
#         if row['filter'] == 'G102':
#             dcorr = 24.5/disp
#         else:
#             dcorr = 46.5/disp
        dcorr = 1. # apparently big improvement as of 2018-06-19?
        if verbose:
            print("{}: extraction finished.".format(preamble))
        
        flux *= dcorr
        s_back *= dcorr
        sky_back *= dcorr
        errf *= dcorr
        b_lower *= dcorr
        b_upper *= dcorr
        sky_b_lower *= dcorr
        sky_b_upper *= dcorr
        gross = flux + s_back + sky_back

# if strpos(file,'icqw01') ge 0 then begin ;GD71 G102 contam by another star.
# 	bad=where(wave gt 11450 and wave lt 15400)
# 	spec_time(bad)=0.
# 	endif
# if strpos(file,'icqw02') ge 0 then begin ;GD71 G141 contam by another star.
# 	bad=where(wave gt 17000)
# 	spec_time(bad)=0.
# 	endif

        # Create trace image
        if interactive:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            targ_str = "{} - {} ({})".format(root, target, filter)
            ax.set_title('{}: Extracted Spectrum'.format(targ_str))
            ax.plot(wave, flux, label='flux')
            ax.plot(wave, gross, label='gross')
            ax.plot(wave, s_back, label='background')
            ax.plot(wave, sky_back, label='sky background')
            ax.legend()
            plt.show()
        
        # Write results
        extracted_file_name = "{}_{}_x1d.fits".format(row['root'], row['target'])
        
        x_col = fits.Column(name='x', format='E', array=x_arr)
        fit_col = fits.Column(name='y_fit', format='E', array=yfit)
        wave_col = fits.Column(name='wavelength', format='E', array=wave)
        net_col = fits.Column(name='net', format='E', array=flux)
        gross_col = fits.Column(name='gross', format='E', array=gross)
        back_col = fits.Column(name='background', format='E', array=s_back+sky_back)
        eps_col = fits.Column(name='eps', format='E', array=epsf)
        err_col = fits.Column(name='err', format='E', array=errf)
        x_back_col = fits.Column(name='x_background', format='E', array=x_b)
        b_lower_col = fits.Column(name='background_lower', format='E', array=b_lower)
        b_upper_col = fits.Column(name='background_upper', format='E', array=b_upper)
        time_col = fits.Column(name='time', format='E', array=spec_time)
        cols = fits.ColDefs([x_col, fit_col, wave_col, net_col, gross_col,
                             back_col, eps_col, err_col, x_back_col, 
                             b_lower_col, b_upper_col, time_col])
        table_hdu = fits.BinTableHDU.from_columns(cols)
        
        params['xc_f'] = xc
        params['yc_f'] = yc
        params['xastr'] = xastr
        params['yastr'] = yastr
        params['avgbkg'] = avgbkg
        params['flatfile'] = flatfile
        params['coef_0'] = coef[0]
        params['coef_1'] = coef[1]
        params['ref'] = ref
        params['angle'] = angle
        params['axeflg'] = axeflg
        
        hdr = fits.Header()
        hdr = set_hdr(hdr, params)
        primary_hdu = fits.PrimaryHDU(header=hdr)
        
        extracted_file = fits.HDUList([primary_hdu, table_hdu])
        
        out_dir, out_table = os.path.split(arg_list.out_file)
        if out_dir == '':
            out_dir = os.getcwd()
        spec_dir = os.path.join(out_dir, arg_list.spec_dir)
        Path(spec_dir).mkdir(parents=True, exist_ok=True)
        extracted_dest = os.path.join(spec_dir, extracted_file_name)
        extracted_file.writeto(extracted_dest, overwrite=True)
        row['extracted'] = os.path.join(arg_list.spec_dir, extracted_file_name)
        
    with fits.open(file, mode='update') as f:
        set_hdr(f[0].header, params)

    return row


def reduce(input_table, overrides, arg_list):
    """
    Reduces grism data
    """
    verbose = arg_list.verbose
    interactive = arg_list.trace
    task = "grism_reduce"
    
    if verbose:
        print("{}: Starting WFC3 data reduction for GRISM data.".format(task))
    
    known_issues_file = get_data_file("abscal.wfc3", "known_issues.json")
    with open(known_issues_file, 'r') as inf:
        known_issues = json.load(inf)
    input_table.adjust(known_issues['metadata'])
    issues = []
    if "reduction" in known_issues:
        if "reduce_grism" in known_issues["reduction"]:
            issues = known_issues["reduction"]["reduce_grism"]
    
    for row in input_table:
        root = row['root']
        target = row['target']
        ref = row['filter_root']
        preamble = "{}: {}".format(task, root)
        
        # Don't extract if there's already an extracted version of
        #   the file present.
        if row['extracted'] != '':
            out_dir, out_table = os.path.split(arg_list.out_file)
            if out_dir == '':
                out_dir = os.getcwd()
            ext_file = os.path.join(out_dir, row['extracted'])
            if os.path.isfile(ext_file):
                if arg_list.force:
                    if verbose:
                        msg = "{}: {}: extracted file exists. Re-extracting."
                        print(msg.format(task, root))
                else:
                    if verbose:
                        msg = "{}: {}: skipping extraction because file exists."
                        print(msg.format(task, root))
                    continue
        else:
            out_dir, out_table = os.path.split(arg_list.out_file)
            if out_dir == '':
                out_dir = os.getcwd()
            spec_dir = os.path.join(out_dir, arg_list.spec_dir)
            extracted_file_name = "{}_{}_x1d.fits".format(root, target)
            extracted_dest = os.path.join(spec_dir, extracted_file_name)
            
            # If there is already an extracted file for this input, skip.
            if os.path.isfile(extracted_dest):
                ext_str = os.path.join(arg_list.spec_dir, extracted_file_name)
                if arg_list.force:
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
        if row['use'] and row['filter'][0] == 'G':
            defaults = {
                        'xc': -1., 
                        'xerr': -1.,
                        'yc': -1., 
                        'yerr': -1., 
                        'ywidth': 11, 
                        'y_offset': 0,
                        'gwidth': 6.,
                        'bwidth': 13, 
                        'bmedian': 7,
                        'bmean1': 7, 
                        'bmean2': 7,
                        'slope': 1,
                        'yshift': 0,
                        'ix_shift': 0,
                        'iy_shift': 0,
                        'wl_offset': 0.,
                        'wlrang_m1_low': -1.,
                        'wlrang_m1_high': -1.,
                        'wlrang_p1_low': -1.,
                        'wlrang_p1_high': -1.,
                        'wlrang_p2_low': -1.,
                        'wlrang_p2_high': -1.
                       }
            defaults['bdist'] = 25 + defaults['bwidth']/2
            defaults['ubdist'] = defaults['bdist']
            defaults['lbdist'] = defaults['bdist']
            if row['filter'] == 'G102':
                defaults['ix_shift'] = 252
                defaults['iy_shift'] = 4
                defaults['wlrang_m1_low'] = 8000.
                defaults['wlrang_m1_high'] = 10000.
                defaults['wlrang_p1_low'] = 8800.
                defaults['wlrang_p1_high'] = 11000.
                defaults['wlrang_p2_low'] = 8000.
                defaults['wlrang_p2_high'] = 10000.
            elif row['filter'] == 'G141':
                defaults['ix_shift'] = 188
                defaults['iy_shift'] = 1
                defaults['wlrang_m1_low'] = 10800.
                defaults['wlrang_m1_high'] = 16000.
                defaults['wlrang_p1_low'] = 10800.
                defaults['wlrang_p1_high'] = 16000.
                defaults['wlrang_p2_low'] = 10800.
                defaults['wlrang_p2_high'] = 13000.
            params = set_params(defaults, row, issues, preamble, overrides, 
                                verbose)
            if 'bwidth' in params['set']:
                if 'bdist' not in params['set']:
                    params['bdist'] = 25 + params['bwidth']/2
                if 'ubdist' not in params['set']:
                    params['ubdist'] = 25 + params['bwidth']/2
                if 'lbdist' not in params['set']:
                    params['lbdist'] = 25 + params['bwidth']/2
                        
            if verbose:
                print("{}: Reducing {} ({})".format(task, root, row['filter']))
            
            print("{}: Reference Image is {}".format(task, ref))
            # Only get the position if there
            #   - is a known filter reference (not NONE or unknown)
            #   - that reference is in the table (so we can find it)
            if ref != 'NONE' and ref != 'unknown' and ref in input_table['root']:
                ref_row = input_table[input_table['root'] == ref]
                # Process the filter image only if its XC and YC haven't been
                #   set to an actual value yet.
                if float(ref_row['xc'][0]) < 0 or float(ref_row['yc'][0]) < 0:
                    processed = locate_image(ref_row, verbose, interactive)
                    ref_row['xc'] = float(processed['xc'][0])
                    ref_row['yc'] = float(processed['yc'][0])
                    ref_row['xerr'] = float(processed['xerr'][0])
                    ref_row['yerr'] = float(processed['yerr'][0])
                for item in ['xc', 'yc', 'xerr', 'yerr']:
                    params[item] = float(ref_row[item])
                    if item not in params['set']:
                        params['set'].append(item)

            # Check scan rate to determine extraction type
            scan_rate = row['scan_rate']

            if scan_rate > 0:
                print("{}: Reducing SCAN MODE data.".format(task))
                reduced = reduce_scan(row, params, arg_list)
            else:
                print("{}: Reducing STARE MODE data.".format(task))
                reduced = reduce_stare(row, params, arg_list)
                
            row['extracted'] = reduced['extracted']
                
        elif verbose:
            msg = "{}: Skipping {} because it's been set to don't use "
            msg += "(reason: {})."
            if row['filter'][0] != 'G':
                reason = 'Imaging exposure ({})'.format(row['filter'])
            else:
                reason = row['notes']
            print(msg.format(task, root, reason))
    
    return input_table


def additional_args():
    """
    Additional command-line arguments. Used when a single command may run
    another command, and need to add arguments from it.
    """
    
    additional_args = {}
    
    table_help = "The input metadata table to use."
    table_args = ['table']
    table_kwargs = {'help': table_help}
    additional_args['table'] = (table_args, table_kwargs)

    bkg_help = "Whether to subtract background before or after applying "
    bkg_help += "flatfield. Default is 'flat_first'. Available options are "
    bkg_help += "'flat_first', 'bkg_first' and 'bkg_only'."
    bkg_args = ['-b', '--bkg_flat_order']
    bkg_kwargs = {'dest': 'bkg_flat_order', 'default': 'flat_first',
                  'help': bkg_help}
    additional_args['bkg_flat_order'] = (bkg_args, bkg_kwargs)
    
    trace_help = "Include result plots while running."
    trace_args = ["-t", "--trace"]
    trace_kwargs = {'dest': 'trace', 'action': 'store_true', 'default': False,
                    'help': trace_help}
    additional_args['trace'] = (trace_args, trace_kwargs)

    return additional_args


def parse_args():
    """
    Parse command-line arguments.
    """
    description_str = 'Process files from metadata table.'
    default_output_file = 'ir_grism_stare_reduced.log'
    
    additional_args = additional_args()
    
    res = parse(description_str, default_output_file, additional_args)
    
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


def main(overrides={}):
    parsed = parse_args()
    
    for key in overrides:
        if hasattr(parsed, key):
            setattr(parsed, key, overrides[key])
    
    input_table = AbscalDataTable(table=parsed.table,
                                  duplicates=parsed.duplicates,
                                  search_str='',
                                  search_dirs=parsed.paths,
                                  idl=parsed.compat)

    output_table = reduce(input_table, overrides, parsed)
    
    table_fname = res.out_file
    output_table.write_to_file(table_fname, res.compat)


if __name__ == "__main__":
    main()
