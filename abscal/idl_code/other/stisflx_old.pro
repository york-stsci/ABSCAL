PRO stisflx_old,file,hdr,wave,flux,net,gross,blower,bupper,epsf,errf,ord,bad=bad, $
			ttcor=ttcor,notime=notime,notemp=notemp,nofile=nofile
;+
;
; PURPOSE:
;  read stis data & calibrate net countrate to flux. Corr both for time &
;	temp changes. Corr flux for slit thruput via APERTAB (ST) APTTAB (IDL) &
;	corr for gwidth extr hgt via PCTTAB (IDL-abscor), both via calstis_abs.
; ***********************************************************************
; *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** 
;     ****GROSS & BKGs are NOT corrected!
;  *** SO: BKG NOT= G(uncor)-Net(corr) !!!!!!!!!!!!!
;	AND do NOT cf uncorr flux from preproc w/ corr flux output here, because
;	preproc has approx flux w/ ctecorr done before defringing, while here
;	flux is from net already defringed & then CTE corr.
;	NO CTEcorr done in preproc to NET; but GAC is applied to net by preproc.
;	{Should put GAC corr here & not in preproc. Otherwise need special 
;		preproc run to make all AGK w/o GAC corr.}
; ***********************************************************************
; INPUT
;	- file-extracted spectral data file
;	- /ttcor, if present make the time, temp, cte corrections & sens* from
;		/scal. Otherwise just return w/ same result as mrdfits,
;		regardless of other keyword setting.
;	- /notime, if present along with /ttcor make only temp correction
;	- /notemp, Do NOT do temp corr. for use in tcorrel.
;	- /nofile. No input file to read for K648 work. Net from geomed image.
; OUTPUT - hdr=header
;	 - wave=wavelengths
;	 - FLUX If /ttcor: corr for cte,time,temp here & also
;			   aper,gwidth here via calstis_abs,_sens
;		gwidth corr uses eg: stisidl/scal/abscor-750.fits
;		halofrac.pro uses stiscal/dat/abscor-g750l, called by ctecorr.
;		BUT 2 files are the same except for added aperture by DJL. See
;		abscor.isr90-01-99update & stis/doc/abscor.pct. 
;		Cked for G750L only.
;	 - EPSF  quality flags
;	 - ERRF  uncertainties in flux units
;	 - GROSS counts/s
;	 - BLOWER lower bkg
;	 - BUPPER upper bkg
;	 - NET counts/s corr for CTE,time,temp, if /ttcor.
;		also E1 corr for GAC, small scale corr, in preproc.pro.
;		NEVER corr for gwidth or aperture, which are done to flux only
;							 via calstis_abs.
;	 - ORD order number
;	- /bad, set bad='B' for late G140L @ +3", where no corr is made
;
; SUBROUTINES CALLED
;	calstis_abs, which also calls calstis_sens to do the aperture and 
;	gwidth corrections to the flux.
; AUTHOR-R.C.BOHLIN
; HISTORY:
; 97aug14 - written
; 97dec10 - remove radial veloc from wavelengths before calibrating;
; 98feb13 - but radial veloc corr in orig wave array is not changed!
; 98apr13 - add echelle capability, as well as order output array
;					of 16-JAN-1998 software update.
; 98jul17 - time,temp corr for MAMA's
; 98sep10 - add sxaddpar,hdr,'pcttab','NONE' and remove my own equiv code.
; 98sep10 - oops, forgot star rad vel is no longer incl in wls, so remove here. 
; 99sep7  - move up time and temp corr to correct net, as well as flux.
; 00jan10 - turn on my good aperture corrections, as fed back from DJL
; 00jan14 - fix bug introduced for spectra processed after 99dec9 w/ new 
;		TIMECORR= 1 header param, which made calstis_abs do corr again.
; 00apr10 - updated errf to get improved error est. after time,temp corr.
; 01apr9  - add CCD CTE correction per ctecorr.pro for net and flux
; 02mar31 - turn off p3m3 corr, as L-flats and new time corr mod are implemented
; 02apr4  - Flag attempts to correct G140L after 1999.2 @ +3"
; 02apr4  - add notemp keyword
; 03jan7 - elim confusing cal keyword and replace w/ bad='B' keyword for 
;					the flagged G140L >1999.2 @+3"
; 03nov7 - add section to do special abscor for wide extract. of sat Vega obs,
; 03dec31	where gross and net get updated to 7px by abscor correction.
; 07aug17 - add /nofile flag to skip reading data and just correct net as input
; 09jul27 - add 1.5 e- for G=4 in ctecorr.pro
; 09jul28 - update G=4/1 ratio from 4.039 to 4.016
; 13feb12 - add G230LB scat lite corr here.
; 14sep29 - add gwidth=11 sens cal for G750L
;-
st=''
if keyword_set(nofile) then begin
	exptm=sxpar(hdr,'exptime')			; IDT geom file only-ok
	goto,skipreading
	endif
z=mrdfits(file,1,hdr,/silent)
indx=where(strpos(hdr,"= 'A2CENTER") ge 0)
if indx(0) gt 0 then begin		; case of reading stsci data
	dum=mrdfits(file,0,hdr,/silent)	; ST header only in extent zero
	exptm=sxpar(hdr,'texptime')
	epsf=z.dq		; However, these have different definitions!
	errf=z.error
	blower=z.background/2	; odd normalization, however. fix for specif app
	bupper=z.background/2
	ord=z.sporder
	sxaddpar,hdr,'pcttab','NONE'		; stsci has no pcttab 02jun17
;2014Jul18-STSCI has pctab, while DJL has pcttab, eg:
;PCTAB   = 'oref$q541740po_pct.fits' / Photometry correction table               
;PCTTAB  = 'abscor-750.fits'    /
    end else begin
	exptm=sxpar(hdr,'exptime')
	epsf=z.epsf
	errf=z.errf
	blower=z.blower
	bupper=z.bupper
	ord=z.order
	endelse
; names common to IDT and STScI:
wave=z.wavelength
flux=z.flux
gross=z.gross
net=z.net
skipreading:
if keyword_set(ttcor) then begin			; else jump to end
	optmode=strtrim(sxpar(hdr,'opt_elem'),2)
	root=strtrim(sxpar(hdr,'rootname'),2)
	msmpos=sxpar(hdr,'OMSCYL1P')	;Mode select cylinder 1 position
	det=strtrim(sxpar(hdr,'detector'),2)
	gwidth=fix(sxpar(hdr,'gwidth'))			; orig gwidth
	targ=strtrim(sxpar(hdr,'targname'),2)
	if targ eq 'HD172167-V6' then targ='HD172167'	;08dec15 - patch
	sens='sens_'
; ff for hz43,vega, etc.
	if gwidth ne 7 and optmode eq 'G750L' and strpos(file,'hgt7')	$
			lt 0 then sens='sens11_'
	sfile=findfile('/internal/1/stisidl/scal/'+sens+		$
						strlowcase(optmode)+'.fits')
	indx=where(strpos(sfile,'PFL') lt 0,ngood)
	if ngood gt 0 then sfile=sfile(indx)
	if sfile(0) eq '' or ngood eq 0 then begin
		if sfile(0) ne '' then goto,OK
		print,'STISFLX: No sensitivity file for ',optmode
		if strpos(optmode,'L') ge 0 then stop
		return					; 98feb13
	    end else begin
OK:
		sxaddpar,hdr,'senstab',sfile(0)		; for calstis_abs
		ind=where(flux ne 0)
; ng-flux is corr	errf(ind)=errf(ind)*net(ind)/flux(ind)	;errf in counts/s
		errf(ind)=errf(ind)/flux(ind)	; 04jan17 - frac error
; IDL intrumental wavelengths for calibrating
		winstr=wave
		helio=strtrim(sxpar(hdr,'helio'),2)
		if helio eq '1' then begin
			evel=sxpar(hdr,'earthvel')	; veloc toward star is +
		        winstr=wave+wave*(-evel)/3e5 ; corr obs to obs. wl frame
;			print,root,evel,'= Heliocentric Rad. Veloc Corr removed'
			endif
; STSCI intrumental wavelengths for calibrating
		helio=strtrim(sxpar(hdr,'helcorr'),2)
		if helio eq 'COMPLETE' then begin
			indx=where(strpos(hdr,'Helio') ge 0)  &  indx=indx(0)
			pos=strpos(hdr(indx),'=')+1
			evel=-float(strmid(hdr(indx),pos,8))
		        winstr=wave+wave*(-evel)/3e5 ; corr obs to obs. wl frame
;			print,root,evel,'= Heliocentric Rad. Veloc Corr removed'
			endif

		oldnet=net
; special processing for heavily saturated Vega & Sirius obs:
		if gwidth ge 35 and (targ eq 'HD172167' or		$
				targ eq 'SIRIUS') then begin
; special abscor correction from Vega for G230LB & AGK for G430L & G750L obs.
; 2013feb-for Sirius use all AGK, as 0.3s G230LB is sat.See sirius/doc.procedure
; Fix for ht=11 per file name to get ttcorr ht=11 corr.
; 2014oct1-new abscor-g750lvega.fits file for corr for using ht=11 sensitivity.
;	increase ihgt count by 1 for added gwidth=11
			newidth='7'
			if sens eq 'sens11_' then newidth='11'
			sxaddpar,hdr,'gwidth',newidth	;now scaled to 7 or 11px
			vegfil='dat/abscor-'+strlowcase(optmode)+'vega.fits'
			ihgt=2		; pick proper extr hgt for the sat. data
			if optmode eq 'G230LB' and targ eq 'HD172167' then begin
				ihgt=1			; no gwidth=11 in file
				vegfil='dat/abscor-g230lbvega3-g230lb.fits'
				endif
; del??	2014oct1	sxaddpar,hdr,'absfile',vegfil		; 03dec28
; del??	2014oct1	sxaddpar,hdr,'pcttab','Special for Saturated Data'
			sxaddpar,hdr,'pcttab',vegfil
			zabs=mrdfits(vegfil,1,hdum)	; vega col1, sirius col2
			wnode=zabs.wavelength
			tnode=zabs.throughput
;corr to def hgt=7 or 11px,even for G750L, where mrgall has made the 7px bit 
;									into 35:
			if targ eq 'SIRIUS' then ihgt=ihgt+1
			corr=cspline(wnode(*,ihgt),tnode(*,ihgt),winstr);cor to7
; CORR G750L from large/7 to large/11px, ie ihgt=1:
			if newidth eq '11' then 			$
			       corr=corr/cspline(wnode(*,1),tnode(*,1),winstr)
			net=net/corr
			errf=errf*sqrt(corr)			; 04jan17
; Now i have 7px net, so must update net & gross in z.structure to make
;	ctecorr work:
			bkg=(gross-oldnet)*float(newidth)/gwidth
			z.net=net
			gross=net+bkg
			z.gross=gross		; bkg is now also per 7px
			print,'Special STISFLX '+targ+		$
					' ABSCOR for gwidth=',gwidth
			endif
; corr net for CTE loss and add any epsf flags for data out of range
; 				reinitializes net & epsf!
;
; ###change
;Print,' *** CTECORR TURNED OFF IN STISFLX'
		if det eq 'CCD' then ctecorr,hdr,z,epsf,net
; ###change end
; 02apr4 - set B flag if G140L obs at +3" after 1999.2
		time=absdate(sxpar(hdr,'pstrtime'))  &  time=time(0)
; 99sep16 - use MSM cyl 1 position to distinguish bwtn +3 and -3 arcsec 
		if msmpos gt 800 and optmode eq 'G140L' and		$
			time gt 1999.2 then begin
				print,time,' BAD time corr for late +3" '+file
				strput,file,'B',0	; flag for stisadd
				bad='B'			; flag for lowsens
				endif
		if keyword_set(ttcor) then ttcorr_old,hdr,wave,net,		$
					notime=notime,notemp=notemp
; 07mar28 - looks like here is the place to make gross=net+orig BKG, sometime??
; reset nocal (251), so calstis_abs can properly update epsf.
		nosens=where(epsf eq 251,nbad)
		if nbad gt 0 then epsf(nosens)=155		; edge flag
; adjust NET at CCD gain=4 for 4.034+-.01 factor of walborn & bohlin
; the 4.034 is used for my 00feb23 flux cal bohlin (2000, AJ, 120, 437)
; adjust NET at CCD gain=4 for 4.039+-.006 (smith, etal. isr00-03) - 00dec26
; .............................4.016    09jul28
; 2017feb16 the orig 4.034 value must be put in by some stisidl program; but 
;	i do not see it. Pipeline headers have 4.015 instead of the 4.016???
; foolproof! see calstis.doc. STSCI-nosuch keyword!!??
		gain=sxpar(hdr,'ATODGAIN')
		if gain gt 3 and gain lt 5 then begin
;			print,'Adj NET (& FLUX) for 4.016 gain ratio from:',gain
			if gain ne 4.034 then stop	; and consider 09jul28
			net=net*4.016/gain  &  endif 	; 00dec26
; 00jan14 - djl says default is to apply time/temp corr to flux (only).
; 	Since i already corrected net, turn off timecorr for calstis_abs
		origtcor=sxpar(hdr,'TIMECORR')
		sxaddpar,hdr,'TIMECORR',0
		dum=errf
; 2013feb11 - put G230LB scat lite corr here & rm from PREPROC.pro:
		if optmode eq 'G230LB' then begin
			scat=sxpar(hdr,'scatlite')
			print,'G230LB scat lite corr=',scat
			net=net-scat
			endif
		calstis_abs,hdr,1,winstr,net,dum,epsf,flux	;1 for first ord
		sxaddpar,hdr,'TIMECORR',origtcor		; replace orig
		errf=errf*flux				; 04jan17
		fracerr=1./sqrt(net*exptm>0.1)		; 0.1 min stis exp 03dec
; Estimate only the new errf values, where orig errf = 0
		if nbad gt 0 then errf(nosens)=flux(nosens)*fracerr(nosens)
		endelse
	endif		; END ttcor
RETURN			; just returns input file values if ttcor is not set
END
