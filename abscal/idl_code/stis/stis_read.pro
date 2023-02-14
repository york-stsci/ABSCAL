pro stis_read,file,h,d,udl,heng1,heng2,err,eps,no_abort=no_abort, $
	message=message, time_tag=time_tag, $
	hires=hires, no_overscan=no_overscan,rotate=rotate,nofix=nofix, $
	prelaunch=prelaunch,postlaunch=postlaunch, $
	readout=readout,unrotate=unrotate,nreads=nreads, $
        first=first, last=last, subarray=subarray
;
;+
;	Routine to read stis data
;
; INPUTS:
;	file - id number (entry number in the data catalog), the
;		file name, or the stis_open file control structure
;
; OUTPUTS:
;	h - FITS header
;	d - image (for time-tag files,d will be the undecoded table, users
;		should use routine TIME_TAG to read raw time-tag data)
;	udl - unique data log
;	heng1, heng2 - eng. snapshot headers. (heng2 is not used for
;		postlaunch data)
;	err - statistical error image (if present in input data). If not present
;		err is set to scalar 0.  For post-launch time-tag data, this
;		output will contain undecoded event collection start/stop
;		time binary table data.
;	eps - data quality image (if present in input data). If not present
;		eps is set to all zeros
;
; EXAMPLES:
;	all of the following can be used to read the file for entry number
;	25 in the FSWTEST database
;	
;		stis_read,-25,h,d
;		stis_read,'fsw25.fits',h,d,/prelaunch
;
;
;		stis_open,800,fcb,nreads=nreads
;		stis_read,fcb,h,d,readout=3
;		stis_close,fcb
;
; Optional keyword input:
;	/no_abort - specifies to return instead of retall if error encountered
;	message - output error message if no_abort used.
;	/time_tag - read post-launch raw time tag data
;	/hires - if set 2048x2048 mama data is not automatically
;		converted to lores.
;	/no_overscan - if set CCD is not automatically overscan corrected.
;	/rotate - specify to rotate image so that wavelength increases from
;		the left to the right.
;	/nofix - specify if you do not want the header fixed using info from
;		the FSWTEST data base.  If input ID is not a number
;		the header is not fixed.
;	/prelaunch - search the prelaunch database despite the value
;		of !prelaunch
;	/postlaunch - search the postlaunch database despite the value of
;		!prelaunch
;	readout - number of the readout number in post-launch image files
;		(default = 1)
;	/unrotate - specifies that data should be in the raw readout
;		coordinates.  Default is the ST ScI user coordinates.
;	/rotate - specifies that cross-dispersed modes should be rotated
;		so that wavelengths increase left to right
;	first - first word to read in a timetag dataset
;	last - last word to read in a timetag dataset	
;       /subarray - if data is subarrayed, place it in 1024x1024 array
;
; Optional Keyword Output:
;	nreads - number of image readouts in the file
;
; Notes:
;	Routine will search the users current working directory and then
;	search the data path as defined by logical SDATA or SARC.
;	
;	STIS_READ creates a plot title in !p.title
;
; HISTORY:
;	VERSION 1: D. Lindler  March 5, 1996
;	6-mar-1996	jkf/acc	- added find_with_def to search SDATA path.
;	15-mar-1996	added no_abort and message keywords
;			added hires keyword
;			Added plot title generation to !p.title
;	29-mar-1996	TLB/ACC - add automatic overscan subtraction, plus
;			dubbuging mode.
;	30-mar-1996	DJL - modified to recompute fit when data was
;			rejected from the median overscan.
;	3-apr-96	TLB - Added automatic gunzipping. 	
;	13-apr-96	DJL - made so that other users besides stisdata can
;			gunzip.
;	20 Apr 96	TLB - removed check of image y-size for overscan 
;			correction.
;	23 Apr 96	TLB - Fixed problem with median overscan rejection.
;	10 Jun 96	TLB - Stopped interactive prompting to remove
;			uncompressed file.
;	25 Jul 96	TLB - Changed overscan correction to use columns
;			1046-1062.
;	29 July 96	TLB - Changed overscan subtraction to use
;			floating point math.
;	20 Aug 96	DJL - added call to STIS_FSWFIX and added ROTATE/NOFIX
;			 keyword options. added FSWENTRY keyword to header,
;			added /ttagimage switch to convert time-tagged data
;			to an image.
;	11 Sept 96	TLB - fixed Null filename problem when compressed 
;			filename is given instead of entry number. 
;	25 Sept 96      DJL - modified not to read UDL and eng. snapshots
;				when not in calling sequence.
;	9 Oct 96	TLB - Moved data section to cols 19-1042 & overscan
;				to 1045-1060.
;	24 Jan 97	TLB - CCD Overscan now done by calstis_overscan.pro.
;	27 Jan 97	TLB - routine check to see if Overscan already removed.
;	31 Jan 97	TLB/JLS - added /prelaunch keyword
;	06 Feb 97	JLS - changed code to look at !prelaunch system variable
;	7 Mar 97 	DJL - Additional keywords added to prelaunch data
;			       to match post-launch keywords, now traps 
;				invalid readout number error. 
;	13 Mar 97	DJL - added NREADS keyword output, added appropriate
;				error messages when files not found, added 
;				optional ERR and EPS outputs.
;	22 May 97       DJL - modified to use stis_open, stis_close
;	23 May 97 	DJL - added post-launch time-tag data support
;	27 May 97	DJL - modified to correctly pass /prelaunch
;				parameter to calstis_overscan
;	5 June 97	DJL - added FIRST and LAST keyword inputs.
;				changed to use HREBIN in convertsion
;				from hires to lores.
;	16 June 97	DJL - modified to not rotate timetag vectors
;	24 Nov 97	DJL - modified not to copy DOPPZERO from spt header
;       1  May 98       PP  - added subarray capability
;	13 Oct 98	DJL - added multiple group capability for time
;				tag data
;	20 Oct 98       DJL - modified to correct missing mode_id in raw
;				time-tag data for UDLS after the first one.
;       30 Oct 98       DJL - modified to not require _spt file.
;	24 Apr 99	DJL - modified to change fill data flag from 16 to 255
;				and to remove overscan region from eps array
;	07 Dec 99	DJL - added om1cat and pstrtime keywords to header
;	01 Jul 10	DJL - modified to remove bscale and bzero from header
;-
;-----------------------------------------------------------------------
	if n_params(0) lt 1 then begin
		print,'stis_read,file,h,d,udl,heng1,heng2,err,eps'
		return
	end
;
; defaults
;
	rotation = 1					;stsci coord
	if keyword_set(unrotate) then rotation = 0	;raw readout coord
	if keyword_set(rotate) then rotation = 2	;x-disp rotation
	if n_elements(no_abort) eq 0 then no_abort=0
	if n_elements(hires) eq 0 then hires=0
	if n_elements(no_overscan) eq 0 then no_overscan = 0
	if n_elements(readout) eq 0 then readout = 1
	if keyword_set(nofix) then header_fix=0 else header_fix=1
;
; determine if file is already opened
;
	s = size(file)
	dtype = s(s(0)+1)			;data type
        if dtype eq 8 then begin
                fcb = file
                openfile = 0
           end else begin
		openfile = 1
		stis_open,file,fcb,no_abort=no_abort,message=message, $
			prelaunch=prelaunch,postlaunch=postlaunch, $
			nreads = nreads, time_tag=time_tag
		if !err lt 0 then return
	end
;
; PRELAUNCH DATA
;
	if fcb.prelaunch then begin
		fcb1 = fcb.raw
		FITS_READ,fcb1,d,h
		if strtrim(sxpar(h,'obsmode')) eq 'ACCUM' then $
			if min(d) lt 0 then d = d + (d lt 0)*65536L
		if n_params(0) gt 3 then FITS_READ,fcb1,udl
		if n_params(0) gt 4 then FITS_READ,fcb1,0,heng1
		if n_params(0) gt 5 then FITS_READ,fcb1,0,heng2
		err = 0		;not present in prelaunch data
		eps = 0		;not present in prelaunch data
;
; fix header with FSWTEST info
;
		if (header_fix ne 0) and (fcb.dbentry ne 0) then $
				stis_fswfix,fcb.dbentry,h
		sxaddpar,h,'FSWENTRY',fcb.dbentry
		sxaddpar,heng1,'FSWENTRY',fcb.dbentry
		sxaddpar,heng2,'FSWENTRY',fcb.dbentry
		nreads = 1		;nreads is always 1 for prelaunch data
;
; add postlaunch keywords
;
		stis_pre_to_post,h
	   end else begin
;
; POSTLAUNCH DATA
;
;
; check for valid readout number
;
		extnames = strtrim(fcb.raw.extname)
		if fcb.time_tag then begin
			extname1 = 'EVENTS'
			sci = where(extnames eq extname1,nreads)
			extname2 = 'GTI'
			err_reads = where(extnames eq extname2,nerr_reads)
			ndq_reads = 0
		    end else begin
			extname1 = 'SCI'
			sci = where(extnames eq extname1,nreads)
			extname2 = 'ERR'
			err_reads = where(extnames eq extname2,nerr_reads)
			dq_reads = where(extnames eq 'DQ',ndq_reads)
		end
		if readout gt nreads then begin
			message = 'Invalid readout number requested'
			if keyword_set(no_abort) then begin
				!err = -1
				return
			end
			print,'STIS_READ: '+message
			retall
		endif
;
; read both files
;
		err = 0
		eps = 0
		fits_read,fcb.raw,d,h,extname=extname1, $
				exten_no = sci(readout-1),first=first,last=last
		if (readout le nerr_reads) and (n_params(0) gt 6) and $
			(extname2 eq 'ERR') then begin
			fits_read,fcb.raw,err,hh,extname=extname2, $
				exten_no = err_reads(readout-1) 
			if max(err) eq 0 then err = 0
		end
		if (n_params(0) gt 6) and (extname2 eq 'GTI') then $
			fits_read,fcb.raw,err,hh,extname=extname2
		if (readout le ndq_reads) and (n_params(0) gt 7) then begin
			fits_read,fcb.raw,eps,hh,extname='DQ', $
				exten_no = dq_reads(readout-1) 
			if max(eps) eq 0 then eps = 0
		end
		if fcb.nfiles gt 1 then begin
			fits_read,fcb.spt,udl,heng1,extname='UDL',exten_no=readout
;
; copy some spt header keywords to primary header
;
			hhlist = ['mode_id','slitnum','doppmag','oswabsp', $
				 'omscyl1p','omscyl3p','omscyl4p','pstrtime', $
				 'om1cat']
			for i=0,n_elements(hhlist)-1 do begin
				val = sxpar(heng1,hhlist(i),comment=comment)
				sxaddpar,h,hhlist(i),val,comment
			end
		  end else begin
		  	stis_aper,strtrim(sxpar(h,'aperture')),aper_num=slitnum
			sxaddpar,h,'slitnum',slitnum
		end
		    	
		sxaddpar,h,'integ',sxpar(h,'exptime')
;
; correct mode_id for time-tag groups that were missing it in the raw
; data files.
;
		mode_id = strtrim(sxpar(h,'mode_id'),2)
		if (mode_id eq '') or (mode_id eq '0') then begin
		    opt_elem = strtrim(sxpar(h,'opt_elem'))
		    case opt_elem of
                        'G140L'  : mode_id = '1.1   '
                        'G140M'  : mode_id = '1.2   '
                        'E140M'  : mode_id = '1.3   '
                        'E140H'  : mode_id = '1.4   '
                        'MIRCUV' : mode_id = '1.6   '
                        'MIRFUV' : mode_id = '1.6F  '
			'X140M'  : mode_id = '1.7X3 '
			'X140H'  : mode_id = '1.7X4 '
                        'G230L'  : mode_id = '2.1   '
                        'G230LB' : mode_id = '2.1B  '
                        'G230M'  : mode_id = '2.2   '
                        'G230MB' : mode_id = '2.2B  '
                        'E230M'  : mode_id = '2.3   '
                        'E230H'  : mode_id = '2.4   '
                        'PRISM'  : mode_id = '2.5   '
			'X230M'  : mode_id = '2.7X3 '
			'X230H'  : mode_id = '2.7X4 '
                        'MIRNUV' : mode_id = '2.6   '
                        'G430L'  : mode_id = '3.1   '
                        'G430M'  : mode_id = '3.2   '
                        'MIRVIS' : mode_id = '3.6   '
                        'G750L'  : mode_id = '4.1   '
                        'G750M'  : mode_id = '4.2   '
			else : mode_id = '    '
		    endcase
		    sxaddpar,h,'mode_id',mode_id
		end
	end		
;
; convert mama data to lores
;
	s = size(d) & ns = s(1) & nl = s(2)
	if (hires eq 0) and (ns eq 2048) and (nl eq 2048) then begin
		d = temporary(d)*4L
		hrebin,d,h,out=[1024,1024]
		sxaddhist,'STIS_READ: MAMA data converted to Lo-Res',h
	end

	!err = 0
;
; overscan subtract CCD data if requested. Check header to see if it is already
; done.
;
	det = strtrim(sxpar(h,'detector'),2)
	s = size(d)
	
	hist = sxpar(h,'HISTORY')
	junk = strpos(hist,'CALSTIS_OVERSCAN')
	found = where(junk ne -1,count)
	if (det eq 'CCD' and no_overscan eq 0 and count eq 0) then $
	    calstis_overscan, d, h, prelaunch=fcb.prelaunch, eps = eps

;
; if data is a subarray, insert data into 1024x1024 array
;

        if keyword_set(subarray) then subarray, h,d,err,eps
;
; change flag for fill data from 16 to 255
;
	if n_elements(eps) gt 1 then begin
		fill = where(eps eq 16,n)
		if n gt 0 then eps(fill) = 255b
	end
;
; rotate image if requested
;
	if fcb.prelaunch then sxaddpar,h,'prelaunc',1
	if fcb.time_tag eq 0 then stis_rotate,h,d,rotation=rotation
;
; create a title
;
	if fcb.dbentry eq 0 then name=fcb.fname else name=strtrim(fcb.dbentry,2)
	!p.title = 'ID = '+ name +' '+det+ $
                 ' '+ strtrim(sxpar(h,'mode_id'),2)+' slit#'+ $
                 strtrim(sxpar(h,'slitnum'),2)+'  '+ $
                 strtrim(sxpar(h,'slitsize'),2)+ $
                 ' INTEG = '+strtrim(string(sxpar(h,'integ'),'(F10.1)'),2) + $
		 '  '+strmid(sxpar(h,'expstart'),0,17)
;
; close files if opened by stis_read
;
	if openfile then stis_close,fcb
	sxdelpar,h,['bzero','bscale']
return
end
