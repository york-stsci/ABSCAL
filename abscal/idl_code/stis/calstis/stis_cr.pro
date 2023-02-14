pro stis_cr,lst,h,data,err,eps,nused,outfile=outfile, $
            MULT_NOISE = mult_noise, NSIG=nsig,  $
            MEDIAN_LOOP=median_loop, MEAN_LOOP=mean_loop, $
            MINIMUM_LOOP=minimum_loop, INIT_MED=init_med, INIT_MIN=init_min, $
            INIT_MEAN=init_mean, BIAS=bias, $
            DILATION=dilation,DFACTOR=dfactor,VERBOSE=verbose, $
            DISPLAY=display, NOSKYADJUST=noskyadjust, $
	    MASK_CUBE=mask_cube, NOCLEARMASK=noclearmask, $
	    XMEDSKY=xmedsky, NULL_VALUE=null_value, $
	    WEIGHTING=weighting, INPUT_MASK=input_mask, $
            SKYBOX=skybox, SELECT = select, scale_darks = scale_darks, $
	    DARK_SCALES = DARK_SCALES
;+
;  NAME:              
;    STIS_CR
;
; Routine to combine STIS CCD images with cosmic ray removal using
; routine CR_REJECT.
;
; CALLING SEQUENCE:
;   stis_cr,lst,h,data,err,eps,nused
;
; INPUTS:
;   lst - vector of observation id numbers of file names
;
; OUTPUTS:
;   h - header for combined image
;   data - combined image
;   err - statistical error image
;   nused - number of pixels of component images used in each
;               pixel of "data" and "err"
;   eps - pixel-by-pixel data flags
;            61 - 1 pixel rejected from stack
;            62 - 2 pixels rej.
;            63 - 3 pixels rej.
;            64 - 4 pixels rej.
;            65 - 5 or more pixels rej.
;           239 - all pixels rej
;
; KEYWORD INPUTS:
;   outfile    - name of output file to which to write results.
;                The default is to not write an output file.
;   display    - Flag.  If set, the routine will generate an
;                x-window display of the results
;   mult_noise - Coefficient for multiplicative noise term -- helps
;                account for differing PSFs or subpixel image shifts.
;                default = (0.03)
;   verbose    - If set, lots of output.
;   nsig       - Rejection limit in units of pixel-to-pixel noise
;                (sigma) on each input image.  Can be specified as
;                an array, in which case the dimension gives the
;                maximum number of iterations to run.  
;                (Default = [8,6,4])
;   bias       - Set if combining biases (divides through by number
;                of images at end, since exposure time is 0).
;   dilation   - Set to expand the cosmic rays around
;                the initially detected CR pixels.  The value gives
;                the width into which each CR pixel is expanded.
;                (Default = 0, no dilation)
;   dfactor    - Set to coefficient to be applied to error threshold
;                in region surrounding CRs if doing dilation.
;                (Default = 0.5 if dilation requested)
;   verbose    - Set to get messages.
;   noskyadjust- Set to disable adjusting the sky between the 
;                readouts for the CR rejection.  A non-sky-subtracted 
;                image is returned in either case.
;   xmedsky    - Flag.  If set, the sky is computed as a 1-d array
;                which is a column-by-column median.  This is intended
;                for STIS slitless spectra.  If sky adjustment is
;                disabled, this keyword has no effect.
;   skybox     - Bounds of image area to be used to compute sky:
;                [X0,X1,Y0,Y1].  If not supplied, default is returned.
;   noclearmask- By default, the mask of CR flags is reset before
;                every iteration, and a pixel that has been
;                rejected has a chance to get back in the game
;                if the average migrates toward its value.  If this
;                keyword is set, then any rejected pixel stays 
;                rejected in subsequent iterations.  Note that what 
;                stsdas.hst_calib.wfpc.crrej does is the same
;                as the default.  For this procedure, the default
;                was NOT to clear the flags, until 20 Oct. 1997.
;   null_value - Value to use in resultant image for pixels where
;                all input values rejected.  Default=0
;   select     - Vector of pixel indicies (e.g. from WHERE output) to be 
;                processed.   Useful if only a small number of pixels in the
;                array need to be processed.     If SELECT is used then the
;                output 'data' array will be the smallest square array large 
;                enough to  include all the selected pixels.    Use    
;                data(lindgen(N_elements(select))  to match with select vector.
;                
;   weighting    - Selects weighting scheme in final image
;                   combination:
;                    0 (default) - Poissonian weighting - co-add
;                        detected DN from non-CR pixels.  (Pixel-by-
;                        pixel scaling up to total exposure time,
;                        for pixels in stack where some rejected.)
;                        Equivalent to exptime weighting of rates.
;                    1 or more - Sky and read noise weighting of rates.
;                        Computed as weighted average of DN rates,
;                        with total exp time multiplied back in
;                        afterward.
;
;                   In all cases, the image is returned as a sum in
;                   DN with the total exposure time of the image
;                   stack, and with total sky added back in.
;
;     The following keywords control how the current guess at a CR-free
;     "check image" is recomputed on each iteration:
;
;       median_loop  - If set, the check image for each iteration is
;                      the pixel-by-pixel median.  (Default is mean.)
;       mean_loop    - If set, the check image for each iteration is
;                      the pixel-by-pixel mean.  (Same as the default.)
;       minimum_loop - If set, the check image for each iteration is
;                      the pixel-by-pixel minimum.  (Trivial case, mostly
;                      for testing.)
;
;     The following keywords allow the initial guess at a CR-free "check
;     image" to be of a different kind from the iterative guesses:
;
;       init_med  - If set, the initial check image is
;                   the pixel-by-pixel median.  (Default is minimum.)
;       init_mean - If set, the initial check image is
;                   the pixel-by-pixel mean.  (Default is minimum.)
;       init_min  - If set, the initial check image is
;                   the pixel-by-pixel minimum.  (Same as the default.)
;
; 	/scale_darks - scale each STIS CCD dark to 18 Deg. Housing Temperature
;		assuming 7% change per degree
;
; KEYWORD OUTPUTS:
;    	INPUT_MASK - Mask of saturated pixels used when combining the CR-SPLITS.
;		If then mask equals 0 then the pixel was saturated.
;       MASK_CUBE -  Mask of CR-rejected pixels; 1 means a good pixel, 0 means
;               a CR-rejected pixel 
;	DARK_SCALES - vector of scale factors used to scale each dark if
;		scale_darks is specified
;
; METHOD:
;   see CR_REJECT
;
; HISTORY:
;   version 1.0 D. Lindler  Mar. 13, 1997
;       Incorporated new features of cr_iterg.  First 4 args same
;       as before.  RSH, 17 Mar. 1997
;       26-mar-97 DJL added keyword NCOMBINE to output header
;       4 Apr. 1997 - Name of routine cr_iterg changed to cr_reject.  RSH
;   version 2.0 D. Lindler May 22, 1997, modified to accept a list
;       of file names, added display option
;   v2.1 Jul 2, 1997, Lindler, modified to use ccd parameter table to look
;       up value of readnoise
;       11 July, 1997, Plait, modified to display list of images in window
;   V2.2 Jul 14, 1997,  RSH.   Sky adjustment option, eps array.
;       18 July, 1997, NRC, if lst refers to one image (no cr_split), 
;                            return err and eps as arrays filled with
;                            zeros
;       12 Sep 1997, RSH.  Documentation corrected. 
;   V2.3 17 Sep 1997, RSH.  Underlying cr_reject changes:  /clearmask;
;                           intermediate means now exposure-time
;                           weighted like final sum.  
;                           Err for 1 image gets real error estimate.
;   V2.4 20 Oct 1997, RSH.  Underlying cr_reject changes:  /noclearmask
;                           instead of /clearmask; /xmedsky.
;   V2.5  4 Feb 1998, RSH.  Err and eps initialized with calstis_stat,
;                           instead of in this code.
;         6 Feb 1998, RSH.  Err wasn't getting initialized.  This is
;                           now fixed.
;         9 Feb 1998, RSH.  Eps array allocation for each individual image
;                           moved. 
;   V2.6 19 Mar 1998, RSH.  NULL_VALUE option added.  MAXITER option
;                           deleted because it was confusing.  Number
;                           of iterations is n_elements(nsig).
;        10 Apr 1998, RSH.  Returns mask_cube if requested.
;   V3.0 22 Sep 1998, RSH.  Weighting by sky and RON can be substituted
;                           for the usual summing of DN, in computing
;                           the final combined image. 
;   V3.1 22 Jan 1999, RCB/DJL  Added flagging of saturated pixels in input_mask.
;			    No change, if all or none of images are saturated.
;			    Add input_mask output keyword to allow saturated
;			    pixel mask to be returned to user.
;	 		    Fix bug of nused not set to 1 for single image.
;   V3.2 13 Aug 1999, RCB/DJL modified not to include epsilon values from
;			data points masked as saturated in output image
;			data quality.
;   V3.3 22 Sept 1999, DJL, no longer writes the _spt file (stisread
;			no longer requires it)
;   V3.4 22 Mar 2000, RSH.  Skybox keyword added.
;   V3.5 15 Sep 2001, WBL   Added Select keyword
;   V3.6 28 Jan 2002, DJL, added scale_darks and dark_scales keyword parameters
;   V3.7 06 Aug 2009, rcb-do an average of occdhtav from all extents, instead of
;	just from the last extent
;-
;------------------------------------------------------------------------------

VERSION = 3.7
;
; 
if n_params(0) lt 1 then begin
    print,'CALLING SEQUENCE: stis_cr,lst,h,data,err,eps,nused'
    print,'KEYWORD INPUTS: outfile, mult_noise, nsig, median_loop, '
    print,'                mean_loop, minimum_loop, init_med, init_min,'
    print,'                init_mean, bias, dilation, dfactor,'
    print,'                display, noskyadjust, noclearmask, xmedsky,'
    print,'                null_value, verbose, mask_cube, weighting, skybox'
    print,'                scale_darks, dark_scales
    return
endif

if n_elements(mult_noise) eq 0 then mult_noise = 0.03
if n_elements(dilation) eq 0 then begin
    dilation = 0
    dfactor = 0.0
endif
if (dilation ne 0) and (n_elements(dfactor) eq 0) then $
  dfactor = 0.5
curwin = !d.window

;
; determine number of images to read for multi-file mode
;
nfiles = n_elements(lst)
nimages = 0
for i=0,n_elements(lst)-1 do begin
    stis_open,lst(i),fcb,nreads=nreads
    nimages = nimages + nreads
    stis_close,fcb
end
image_number = 0
; If only selected pixel find the smallest square array large enough to include
; all the pixels 
Nselect = N_elements(select)
if Nselect GT 0 then  nx = ceil(sqrt(Nselect))  
;
; loop on files
;
dark_scales = replicate(1.0,nimages)	;optional scale factors for darks
avtemp=0.				; average OCCDHTAV 09aug6 - RCB
for ifile=0,nfiles-1 do begin
    stis_open,lst(ifile),fcb,nreads=nreads
    for iread = 1,nreads do begin
;
; read image
;
        stis_read,fcb,h,data,udl,heng,readout=iread	; rm overscan

	avtemp=avtemp+sxpar(h,'occdhtav')		; 09aug6
        calstis_ref,h       ;get name of ccd parameter table
        calstis_ccdpar,h    ;get readnoise from ccd parameter table.Add ATODGAIN
        eps = byte(0*data)
        if nimages gt 1 then sxaddpar, h, 'errors', '0' $
                        else sxaddpar, h, 'errors', '1' 
        calstis_stat, h, data, err, eps   ; initialize err, eps
;
; Scale darks to 18 degree CCD housing temperature (NOT a calstis default)
;
	if keyword_set(scale_darks) then begin		; NOT usually set - rcb
		get_occdhtav,h,temperature
		dark_scales(image_number) = (18.0-temperature)*0.07 + 1.0
		print,'dark scaled by ',dark_scales(image_number)
		data = data*dark_scales(image_number)
	end
;
; get relevant header info and determine noise from first image
;
        if image_number eq 0 then begin		; BEGIN initialization loop
            detector = strtrim(sxpar(h,'detector'))
            if detector ne 'CCD' then begin
                print,'STIS_CR: ERROR - detector must be CCD'
                retall
            endif
            readnoise = sxpar(h,'readnse')
            gain = sxpar(h,'CCDGAIN')>1
            readnoise = readnoise/gain ;change readnoise to DN
            s = size(data) & ns = s(1) & nl = s(2)
            if Nselect  GT 0 then begin
                     ns = nx  & nl = nx
            endif
;
; just one image?  
;
            if nimages eq 1 then begin
                print,'STIS_CR: Single readout observation: no cosmic ray '+ $
                  'removal done'
                stis_close,fcb
		nused=bytarr(ns,nl)+1b		; fill nused arr w/ unity
                goto,done
            end
;
; set up data cube and exptime vector
;
            cube = fltarr(ns,nl,nimages)
            epscube = bytarr(ns, nl, nimages)
            exptime = fltarr(nimages)
        endif 					; END array initialization loop
        if keyword_set(select) then begin
                 data = data[select]
                 eps = eps[select]
                 Nfill = nx*nx - Nselect
                 if Nfill GT 0 then begin
                        data = [data,fltarr(nfill)]
                        eps = [eps,bytarr(nfill)]
                 endif
                 data = reform( data, nx,nx)
                 eps = reform(eps ,nx,nx)
        endif
        cube(0,0,image_number) = temporary(data)	; sets data to UNDEFINED
        epscube(0,0,image_number) = temporary(eps)
        exptime(image_number) = sxpar(h,'exptime')
        image_number = image_number+1
    end  ; for iread
    input_mask = epscube ne 190 		;flag non-saturated pixels
    stis_close,fcb
end ; for ifile
;
; If all images are sat at a px, then use them all to maintain continuity
;
allbad=total(input_mask,3) eq 0
for i=0,nimages-1 do input_mask(0,0,i) = allbad or input_mask(*,*,i)
;
; combine readouts
;
cr_reject, cube, readnoise, 0.0, gain, mult_noise, data, err, nused, $
  NSIG=nsig, MEDIAN_LOOP=median_loop, MEAN_LOOP=mean_loop, $
  MINIMUM_LOOP=minimum_loop, INIT_MED=init_med, INIT_MIN=init_min, $
  INIT_MEAN=init_mean, EXPTIME=exptime, BIAS=bias, $
  DILATION=dilation, DFACTOR=dfactor, VERBOSE=verbose, $
  NOSKYADJUST=noskyadjust, NOCLEARMASK=noclearmask, XMEDSKY=xmedsky, $
  MASK_CUBE=mask_cube,NULL_VALUE=null_value,INPUT_MASK=input_mask,	$
  WEIGHTING=weighting,SKYBOX=skybox

eps = byte(data*0)
wnu = where(nused lt nimages, cwnu)
if cwnu gt 0 then eps(wnu) = eps(wnu) > (60 + ((nimages-nused(wnu))<5))
wnu0 = where(nused le 0, cwnu0)
if cwnu0 gt 0 then eps(wnu0) = eps(wnu0) > 239
for i=0,nimages-1 do eps = eps > (epscube(*,*,i)*mask_cube(*,*,i)*	$
							input_mask(*,*,i))
if not arg_present(mask_cube) then mask_cube = 0
epscube = 0
sxaddpar,h,'exptime',total(exptime)
sxaddpar,h,'integ',total(exptime) ;add for old timers
hist = strarr(5)
hist(0) = 'STIS_CR version '+string(version,'(F4.1)')+ $
  ': CCD Cosmic Ray Removal'
hist(1) = '   '+strtrim(nimages,2)+' readouts combined with CR_REJECT'
hist(2) = '   '+strtrim(cwnu,2)+' pixels affected by CRs'
hist(3) = '   '+strtrim(cwnu0,2)+' pixels for which all inputs rejected'
if keyword_set(weighting) then wmeth='1' else wmeth='0'
hist(4) = '   '+'Weighting method '+wmeth+' used'
sxaddhist,hist,h
sxaddpar,h,'NCOMBINE',nimages,'Number of CR-SPLITS combined'
if avtemp gt 0 then 							$
	sxaddpar,h,'occdhtav',avtemp/nimages,'Average computed by stis_cr.pro'
if keyword_set(verbose) or (!dump gt 0) then print,hist
;
; display results
;
if keyword_set(display) then begin
    strlist = ''
    for ititle=0,n_elements(lst)-1 do $
      strlist = strlist+','+strtrim(lst(ititle),2)
    null = gettok(strlist,',')
; NG: tvscl, w/ channel # misses RHS w/o window     if !d.name eq 'X' then
    window,curwin+1,xs=ns,ys=nl/2,title='STIS_CR '+strlist
    nlout = nl/2
    nsout = ns/2
    tvscl,alog10(rebin(data(0:nsout*2-1,0:nlout*2-1), $
                       nsout,nlout)>0.1),0
    tvscl,rebin(float(nused(0:nsout*2-1,0:nlout*2-1)), $
                nsout,nlout) lt nimages,1
    plots,[nsout,nsout],[0,nlout],/dev
end
;
; write file
;
done:
if n_elements(outfile) gt 0 then begin
    fits_open,outfile,fcb,/write
    fits_write,fcb,data,h,extname='SCI',extver=1,extlevel=1
    if n_elements(err) gt 1 then $
      fits_write,fcb,err,h,extname='ERR',extver=1,extlev=1
    if n_elements(eps) gt 1 then $
      fits_write,fcb,eps,h,extname='EPS',extver=1,extlev=1
    fits_close,fcb

endif
;
; return user to original window
;
if (curwin ne -1) then wset, curwin

return
end

