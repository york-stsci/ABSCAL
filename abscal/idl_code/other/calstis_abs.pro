pro calstis_abs,h,m,wave,net,err,eps,flux,pos
;+
;			calstis_abs
;
; Subroutine of CALSTIS to perform flux calibration (sensitivity correction
; and slit throughput correction.
;
; CALLING SEQUENCE:
;	calstis_abs,h,m,wave,net,err,eps,flux,pos
;
; INPUTS:
;	wave - wavelength array
;	net - count rate array
;
; INPUT/OUTPUTS:calstis_
;	h - FITS header
;	m - vector of spectral order numbers (1 for first order gratings)
;	err - error array
;	eps - eps array
;
; OUTPUTS:
;	flux - calibrated flux (ergs/cm2/sec/Ang)
;
; HISTORY:
;	version 1  D. Lindler  6/20/97
;	7/11/97 Lindler, modified to work when error propagation is
;		turned off
;	12/19/97, Lindler, added PCTTAB correction for slit height
;	4/10/98, Lindler, added order dependent sensitivity curves
;	6/11/98, Lindler, now sets flux to zero outside calibrated
;			range of aperture throughputs
;	June 17, 1998, moved throughput correction to calstis_sens.
;	Dec, 2001, Lindler, added echelle blaze shift option
;	June 2004, Lindler, added value of model blaze shift to header
;-
;--------------------------------------------------------------------------
;
	eps_nosens = 251		;data quality for outside of
					; calibrated sens. range
;
; determine data size
;
	s = size(wave) & ns = s(1)
	nspec = n_elements(wave)/ns
;
; get reference file names
;
	senstab = strtrim(sxpar(h,'senstab'))	;sensitivity table
	apttab = strtrim(sxpar(h,'apttab'))	;aperture throughput table
;
; sensitivity correction
;
	flux = net
	if strupcase(senstab) eq 'NONE' then return
	calstis_sens,h,msens,wsens,sens,blaze_model
;
; shift blaze for echelle observations
;
	wsave = wsens
	if blaze_model.mref ne 0 then begin
;
; linearly extrapolate ends of the blaze function
;
		n = n_elements(msens)
		nx = n_elements(wsens(*,0))
		delta = (wsens(1)-wsens(0))*50
		slope1 = (sens(1,*)-sens(0,*))/(wsens(1,*)-wsens(0,*))
		slope2 = (sens(nx-1,*)-sens(nx-2,*))/ $
			 (wsens(nx-1,*)-wsens(nx-2,*))
		sens = [sens(0,*)-delta*slope1,sens, $
				sens(nx-1,*)+delta*slope2]
		wsens = [wsens(0,*)-delta,wsens,wsens(nx-1,*)+delta]
;
; compute blaze shift in pixels
;
		deltat = sxpar(h,'expstart')-blaze_model.mjd
		good = where(m eq blaze_model.mref,ngood)
		if ngood lt 1 then begin
			print,'CALSTIS_ABS: ERRROR - invalid MREF in senstab'
			retall
		end
		j = good(0)
		if ns eq 2048 then hires = 1 else hires = 0
		deltay = pos(ns/2-1,j)/(hires+1)-blaze_model.yref
		deltaw = blaze_model.wref-wave(ns/2-1,j)
		dispersion = (wave(ns/2,j)-wave(ns/2-1,j))*(hires+1)
		deltax = deltaw/dispersion
		bshift = blaze_model.mx*deltax + blaze_model.my*deltay + $
				blaze_model.mt*deltat
		deltamw = dispersion*blaze_model.mref
;
; apply blaze shift to the wavelengths
;
		for i=0,n_elements(msens)-1 do $
			wsens(0,i) = wsens(*,i) + $
				bshift*deltamw/msens(i)
		hist = '  Model blaze shift = '+string(bshift,'$(F5.1)')
		sxaddhist,hist,h
		if !dump gt 0 then print,hist				
	end
;
; apply auto blaze shift
;
	if (max(m) gt 1) and (sxpar(h,'auto_bs') ne 0) then $
			calstis_bs,h,m,wave,net,eps,msens,wsens,sens		
;
; interpolate to observation's wavelengths for each spectral order
;
	for i=0,nspec-1 do begin
	    good = where((msens eq -1) or (msens eq m(i)), ngood)
	    if ngood gt 0 then begin		; ngood=1 for 1st order
		wsens1 = wsens(*,good(0))	;sensitivity for this order
		sens1 = sens(*,good(0))
	    	wmin = min(wsens1)
	    	wmax = max(wsens1)
		wave1 = wave(*,i)		;observed for this order
	    	sens1 = interpol(sens1,wsens1,wave1)
	    	bad = where((wave1 lt wmin) or (wave1 gt wmax) or $ 
			    (sens1 le 0),nbad)			; outside range 
	    	if nbad gt 0 then begin
			sens1(bad) = 1.0
			eps(bad,i) = eps_nosens>eps(bad,i)	; RCB 07may10
			if n_elements(err) gt 1 then err(bad,i) = 0.0
		 endif
;
; divide by sensitivity
;
		 flux(*,i) = net(*,i)/sens1			; the flux cal
		 if nbad gt 0 then flux(bad,i) = 0.0		; no sens cal
		 if n_elements(err) gt 1 then err(*,i) = err(*,i)/sens1
	       end else begin
		 print,'CALSTIS_ABS - Warning: no sensitivities found for '+ $
			'order number '+strtrim(m(i),2)
		 flux(*,i) = 0.0
		 if n_elements(err) gt 1 then err(*,i) = 0.0
		 eps(*,i) = eps_nosens
	    end						; sens for ea order
	end						; each order
return
end
