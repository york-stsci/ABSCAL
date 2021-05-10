pro wfcREDUCE,FILE,GRAT,aper,STAR,OUT=FILENAME,factor=fudge,merge=merge
;+
;
; PURPOSE:
; 18jun12-COADD wfc IR spectra AND WRITE ASCII OUTPUT FILES
; CALLING SEQUENCE:
;	wfcREDUCE,FILE,GRAT,aper,STAR,OUT=filename,factor=fudge,merge=merge
; INPUT:
;	FILE-NAME OF log file OF wrc3 OBSERVATIONS AS GENERATED BY wfcDIR. MUST
;		start w/ 'dir' or 'DIR', else a list of obs to coadd
;		WITH GRATING NAMES IN THE FILE NAMES.
;		02may16-if file='' then merge existing files
;	GRAT-**VECTOR** list OF GRATING modes TO PROCESS (in merge order)
;		If file is a list of obs, then GRAT must NOT be vector, unless
;		grating name is in the filenames (as in MRGTRANSM.pro).
; OPTIONAL INPUT:
;	aper - restrict aperture selection. Use '' to get all apertures.
;	STAR-RESTRICT PROCESSING TO A CERTAIN STAR NAME *************
;	     OTHERWISE, ALL STARS GET ADDED TOGETHER!!!!!!!!!!!!!!!!!
;		Also, specifying star, makes a nice title, for normal use.
;	out=filename-input/output file name
;	factor=fudge - fudge vector corresponding to grat vector. to fix (k648)
;	merge - keyword to merge spectra from diff gratings. default is NO merge
;		01apr3-if only one grating is spec with /merge, then all name.*
;		get merged.
; OUTPUT:
;	ASCII FILES OF THE COADDS & THE MERGING VIA wfcMRG
; RESTRICTION in wfcobs: all input files must have name dat/spec_rootname.fits
; EXAMPLE:
;	wfcreduce,'dirirstare.log',['G102','G141'],'gd153',		$
;			out='spec/gd153',/merge
; HISTORY
; 	2018jun12 - R. BOHLIN
;-

GRAT=STRUPCASE(GRAT)
if n_params(0) gt 3 then STAR=STRUPCASE(STAR)
st=''
!p.font=-1
!y.style=1  &  !x.style=1
!P.NOCLIP=1
hdr1=["SYS-ERROR is a broadband 1% INTERNAL repeatability of WFC3 fluxes.", $
 "IN ADDITION, THERE IS SOME SYSTEMATIC UNCERTAINTY IN THE ABS CALIB, see ",$
" Bohlin et al.2014,PASP,126,711. BOTH THE STAT-ERR & SYS-ERR ARE 1-SIGMA"]
hdr2=[ '',									     $
 '',                                                                         $
"WL(ang)   COUNT-RATE     FLUX     STAT-ERROR   SYS-ERROR  NPTS "+      $
 "  TIME  QUAL",							     $
 " ###    1"]
nodata=0
file=strlowcase(file)				; dir*.log file
ngrat=N_ELEMENTS(GRAT)
gudgrat=grat
if n_params(0) le 2 then star=''
print,'WFCREDUCE for Star/Gratings: ',star,grat
if file(0) eq '' then goto,mergegrat

FOR IGRAT=0,ngrat-1 DO BEGIN		;MAIN LOOP
	lst=strupcase(file)		;case of an input list instead of 'dir*'
	if strpos(file(0),'dir') lt 0 and ngrat eq 1 then goto,nofidl
; start fiddling for usual case of dir*.log:
	if strpos(file(0),'dir') ge 0 then 				$
		wfcobs,file,lst,GRTS,aps,STARS,GRAT(IGRAT),aper,STAR	$
	      else begin		; input list w/ grating names
		indx=where(strpos(lst,grat(igrat)) ge 0,ngood)
		if ngood gt 0 then lst=lst(indx)
		if ngood le 0 then goto,skipgrat    	; gratings w/ no data  
		endelse
	print,'IN/OUT dir=',filename
nofidl:
	print,'      adding: ',lst
	nodata=1				; data found to process
	lst=strlowcase(lst)
	dum=rem_dup(lst)
	if n_elements(dum) ne n_elements(lst) then stop		; idiot check
	if keyword_set(filename) then begin
		name=strlowcase(filename)
		fdecomp,name,disk,dir,nam
		dum=name
		dum2=gettok(dum,'/')
		input=''
		if strpos(dum,'/') ge 0 then input=gettok(dum,'/')+'/'
	     end else begin
		NAME=strlowcase(strtrim(stars(0),2))
		nam=name
		endelse
; case of special in/out dir, eg noflat:

	wfcadd,lst,HEAD,MT,W,C,F,E,NPTS,TIME,QUAL,input=input

	mt=nam+mt				; for ngc -a, etc star names
	if strpos(name,'.g') gt 1 then outnam=name else 	$	;kludge
				outnam=NAME+'.'+GRAT(IGRAT)
	CLOSE,11 & OPENW,11,strlowcase(outnam)
	PRINTF,11,'FILE WRITTEN BY WFCREDUCE.PRO ON ',!STIME
	last=strpos(lst(0),'.fits')
	lst=strmid(lst,last-9,9)
	printf,11,'coadd lst for '+GRAT(IGRAT)+': ',lst,form='(a/(8(a9,1x)))'
	iepoch=where(strpos(head,'EPOCH:') ge 0,nepoch)
	if nepoch gt 0 then printf,11,replace_char(head(iepoch),'HISTORY ','')
	MX=MAX(f)
	!mtitle=MT
	!ytitle='FLUX(10 !e-15!n erg s!e-1!n cm!e-2!n A!e-1!n)'
	!xtitle='Wavelength (A)'
;	plot,w,f,YR=[0,MX]
;	plotdate,'wfcreduce.pro'
;	if !d.name eq 'X' then read,st
	AP=STRTRIM(SXPAR(HEAD,'aperture'),2)
	det=strtrim(sxpar(head,'detector'),2)
	gwid=sxpar(head,'gwidth')
; check for any fudge factors
	if keyword_set(FUDGE) then begin
		f=f*fudge(igrat)
		if fudge(igrat) ne 1. then printf,11,grat(igrat),              $
		  ' flux multiplied by fudge=',fudge(igrat),form='(2a,f5.3)'
		endif
	ind=where(c ne 0.)
	ferr=f*0  &  FERR(ind)=0.01*abs(F(ind))
	printf,11,'gwidth for '+grat(igrat)+' flux cal=',gwid
	PRINTF,11,hdr1,format='(a)'
	printf,11,hdr2
	printf,11,mt
	FMT='(f8.1,4E12.4,I4,F10.1,I4)'
	good=where(npts gt 0)
	FOR I=good(0),good(-1) DO PRINTF,11,FORMAT=FMT,W(I),C(I),F(I),	$
		E(I),FERR(I),NPTS(I),TIME(I),QUAL(I)
	PRINTF,11,' 0 0 0 0 0 0 0 0'
	close,11
skipgrat:
	ENDFOR		;END MAIN LOOP for 2 gratings
; return unless merging spectra
if (not keyword_set(merge)) or nodata eq 0 then return

; *****************************************************************************
;
; MERGE DATA AND WRITE IN .MRG FILE 
mergegrat:
title=strupcase(nam)						; try 03aug3
IF STAR EQ '' THEN TITLE=strtrim(stars(0),2)
grat=gudgrat
; 01apr3 - merge all NAME.*, if only one grat is specified (for indiv lists)
if ngrat eq 1 then filgrat=findfile(strlowcase(name+'.'+gudgrat(0)))

IF N_ELEMENTS(GRAT) EQ 1 THEN TITLE=TITLE+' '+GRAT(0)		       $
	ELSE							       $
	TITLE=TITLE+' '+GRAT(0)+'+'+GRAT(1)

FIL1=NAME+'.'+GRAT(0) & FIL2=''		;SINGLE GRATING
IF N_ELEMENTS(GRAT) GT 1 THEN FIL2=NAME+'.'+GRAT(1)
ext='.mrg'

wfcmrg,strlowcase(fil1),strlowcase(fil2),strlowcase(name+ext),TITLE=TITLE  

;READ BACK THE MERGED DATA WITH THE ASCII FILE READER
RDF,NAME+ext,1,DAT
W=DAT(*,0)
F=DAT(*,2)
STAT=DAT(*,3)
SYST=DAT(*,4)
;
;full PLOT:
!ytitle='FLUX & ERROR ARRAY (10 !e-15!n erg s!e-1!n cm!e-2!n A!e-1!n)'
flx=f*1.e15
MX=MAX(flx(WHERE((W GT 8000))))*1.2
err=1.e15*Stat
plot,w,flx,YR=[0,MX]
oplot,w,err,linestyle=2,th=2
plotdate,'wfcreduce.pro'
!x.style=0
if !d.name eq 'X' then READ,ST
END
