pro hot_pixel_table,darklist,table,threshold=threshold
;+
; NAME:
;    HOT_PIXEL_TABLE
; PURPOSE:
;    Routine to create a hot pixel table
;
; CALLING SEQUENCE:
;    hot_pixel_table,darklist,table, [ threshold= ]
;
; INPUTS:
;    darklist - list of dark files to find hot pixels in
;    table - name of the output hot pixel binary table
;		with columns AXIS1, AXIS2, RATE
;
; OPTIONAL KEYWORD INPUTS:
;    threshold - dark rate threshold for bad pixels (electrons/sec)
;                default = 0.01
;
; HISTORY:
;    Version 1  D. Lindler   June 6, 1997
;    Add start time to FITS header   W. Landsman    April 1999
; 2022aug24-fix mixed [i] & (i) - rcb
; 2022aug25-fix reading of BAD files in cdbs.  FOR EXAMPLE:
;  idl-rcblap> fits_read,'/grp/crds/hst/references/hst/68a1901eo_drk.fits',im,hd
;  idl-rcblap> help,im
;	   IM              LONG      =            0

;-
;--------------------------------------------------------------------------
	if n_params() lt 2 then begin
		print,'hot_pixel_table,darklist,table,threshold='
		print,"e.g. darklist = ['ccddrk_oct11_01','ccddrk_oct19_01']"
		print,"hot_pixel_table,darklist,'hpx_oct19_01.fits'"
		return
	endif
	if n_elements(threshold) eq 0 then threshold = 0.01
;
; read files and find max rate from the images.
;
        expstart = strarr(N_elements(darklist) )
	for i=0,n_elements(darklist)-1 do begin
                fdecomp,darklist[i],disk,dir,name,ext
		fname = find_with_def(darklist[i],'SCAL','fits')
;From: Sean Lockwood <lockwood@stsci.edu>
;Subject: Re: Format change of STIS dark ref files
;Date: August 25, 2022 at 11:15:49 AM EDT

;It looks like the zeroth extension of the FITS files got populated with a sparse
;array, whereas before the data in ext=0 was blank.  This probably caused
;fits_read to pick data from ext=0 instead of ext=1.  As a workaround, you can
;specify the fits_read parameter EXTEN_NO=1 in your code.

;I believe the sparse data arrays got populated when I ran the CRDS command to
;generate a DATASUM (checksum of the ext=0 data array) for the extension, to
;avoid new errors on delivery.  I'll reach out to the ReDCaT team to see if we
;can relax the requirements a bit, or if I can deliver with a dummy DATASUM
;keyword in the future (and redeliver those files).
; 2022aug25-rcb		fits_read,fname,d,h0
		fits_read,fname,d,h0,exten=1		; 2022aug25-rcb
		if n_elements(d) lt 1024 then begin	; 2022aug24-rcb
			print,fname,' BAD file. Returning'
			threshold=999
			return
			endif
; can read as im=mrdfits('/grp/crds/hst/references/hst/68a1901eo_drk.fits',1,hd)
;                expstart[i] = sxpar(h0,'EXPSTART') ; No such in hdr-2022aug-rcb
                expstart[i] = sxpar(h0,'USEAFTER') 
		if i eq 0 then begin
			h = h0
			dark = d 
;rcb-keep HOT px from earlier dark images in the current list; BUT the use in
;	makehotpx.pro always calls this subpro for ea image, so i is always =0.:
		endif else dark = dark>d
	endfor
	s = size(dark) & ns = s[1]
;
; find hot pixels greater than the threshold
;
	bad = where(dark ge threshold,nbad)
	print,'ANALYZING:'+name					; rcb 2022aug24
	print,'HOT_PIXEL_TABLE: '+strtrim(nbad,2)+' Hot pixels found above ' + $
		strtrim(threshold,2)+' Electrons/Second'
	x = bad mod ns
	y = bad/ns
	rate = dark(bad)
;
; write table
;
        tstart = min(expstart)
;        sxaddpar,h,'EXPSTART',tstart always zero - rcb
        message,'Use after ' + strtrim(tstart,2),/INF
	sxaddhist,'HOT_PIXEL_TABLE on input files:',h
        ddarklist = '     ' + strtrim(darklist,2)
        sxaddhist,ddarklist, h	
	MAKE_BINARY_TABLE,table,'AXIS1,AXIS2,RATE',x,y,rate,hdr=h
return
end
