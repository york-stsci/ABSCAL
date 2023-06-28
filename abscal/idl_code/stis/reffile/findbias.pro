; 09jul13 - find STIS CCD bias files and format for calstis_ccd.txt
; 16feb11 - oops... z* names finished, while new files start w/ zero (0).
; 2018dec24 - see stis/doc/bias. & hot.pixel for some details.
;	In particular, the 'new' files reported by this prog in the 55082-56277
;	mjd range are really OLD obsolete/BAD files that were copied w/o
;	retaining the ancient date. Ignore these ~17 bias files. The proper bias
;	files are already in ~/stisidl/scal/calstis_ccd.txt
;-

;fils=findfile('$oref/*bia.fits')	; 2016feb11-fails as new 0* not in order
; so:
; gives setenv error??? --> no harm tho.
;spawn,'ls -rt /grp/hst/cdbs/oref/*o_bia.fits >! tmp.'	; 2019jun18 - moved
spawn,'ls -rt /grp/crds/hst/references/hst/*o_bia.fits >! tmp.'	; 2019dec5
; Trailing o_ on rootname is STIS (j is ACS, etc..) CDBS convention... 

readcol,'tmp.',fils,form='(a)'

;for i=0,n_elements(fils)-1 do begin
for i=2040,n_elements(fils)-1 do begin			; 2021feb21
;	fits_read,fils(i),im,hd,/header_only -- NG k4b0828lo_bia.fits is BAD
	z=mrdfits(fils(i),0,hd,/silent)
	date=strtrim(sxpar(hd,'useafter'),2)
;	print,date
	mmm=strmid(gettok(date,' '),0,3)
	dd=gettok(date,' ')
	yyyy=gettok(date,' ')
	tim=gettok(date,' ')
	if yyyy lt 2009 then goto,skip
;	print,mmm,dd,yyyy,tim
	amp=strtrim(sxpar(hd,'ccdamp'),2)
	gain=strtrim(sxpar(hd,'ccdgain'),2)
	ccdoff=strtrim(sxpar(hd,'ccdoffst'),2)
	binax1=strtrim(sxpar(hd,'binaxis1'),2)
	binax2=strtrim(sxpar(hd,'binaxis2'),2)
	juldat=jul_date(dd+'-'+mmm+'-'+yyyy+' '+tim)-0.5
	fdecomp,fils(i),disk,dir,name,ext
	print,amp,gain,ccdoff,binax1,binax2,name+'.fits',juldat,form=	$
						'(8x,a1,4a3,4x,a,f10.1)'
	skip:
	endfor
end
	
	
