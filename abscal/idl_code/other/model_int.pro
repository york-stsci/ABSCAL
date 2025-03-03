;+
;      model_int
;
; Routine to interpolation inthe leg_orig data base
;
; CALLING SEQUENCE:
;   model_int, teff, log_g, log_z, wave, flux, scalar_item, scalar_value
;
; INPUTS:
;   teff - effective temperature (K)
;   log_g - Log Surface Gravity
;   Log_z - Log Abundance (x Solar)
;
; OPTIONAL INPUTS:
;   scalar_item - name of a scalar item in the database to interpolate.
;     eg, 'BOL_CORR'
; OUTPUTS:
;   wave - wavelength vector
;   flux - flux vector
;   scalar_value - value of the optional scalar item
;
; OPTIONAL KEYWORD INPUTS:
;   wrange = range of wavelengths [wmin,wmax]
;   /quadratic - perform quadratic interpolation in log_z  
;   /spline - perform cubic spline interpolation in log_z
;   /linear - perform linear interpolation in log_z (default)
;   /castelli - use castelli dbase instead of lej_orig
;   /new_castelli - use castelli 2003 dbase. rcb: ie CK04... see below
;   /d01stis - use d01stis dbase
;   /phoenix - use phoenix dbase
;   vturb - Turbulent velocity for MARCS dbase (default=2)
;   alpha = 'n' or 'a'  if not supplied then alpha='n' log_z>-1 else alpha='a'
;   type = 'p_st','s_st','p_ap', 's_ap', 'hc','mc','gs' for marcs db
; NOTES:
;   if inputs are outside of the grid, !err is set to -1 and the closest
;   gridpoint(s) are used.
;
; HISTORY:
;   Nov 30, 2006 - modified to use closest point when outside of grid
;   09sep22 - this new version from Don does not get the same answer as my
;	old version, so try making the same mods that I made before, viz:
;   09feb5-replace the obsolete get_kur_wave with get_castkur04_wave. rcb
;   09feb5-my castelli DB seems equivalent to Don's new_castelli. Make
;	castelli the default and do NOT use the new_castelli keyword.
;	ie use /castelli to get the castelli & kurucz 2004 database.
;	Also castelli='marcs' gets the marcs database. See ~/models/doc.marcs.
; 09dec8 - Put 657. in zaux/castkur04.wavelengths, because the CK04 grid uses
;	that value on  both the Kurucz and Castelli web sites. (Was 656.125)
;	ie see zaux/castkur04.wavelengths-doc & ~/models/misc-doc.models
; 12May24 - In IDL 8 there is a new intrinsic object function called "list". So
;		change all list --> lst.   rcb
;-
;============================================================================
pro model_int1,teff,log_g_obj,log_z_obj,flux,scalar_item,scalar_value, $
    irange=irange,quadratic=quadratic,spline=spline,vturb=vturb,alpha=alpha, $
    outside=outside, type=type
;
; routine to find flux at a tabulated temperature by
; linear interpolation in Log_g and  linear, quadratic, or spline
; interpolation in Log_z
;

    spar = 'teff='+strtrim(teff,2)
    if (strupcase(db_info('name',0)) ne 'LEJ_ORIG') and $
       (strupcase(db_info('name',0)) ne 'PHOENIX') and $
       (strupcase(db_info('name',0)) ne 'MARCS') then $
    	spar = spar + ',vturb='+strtrim(vturb,2)+',alpha='+strtrim(alpha,2)
    if (strupcase(db_info('name',0)) eq 'PHOENIX') then $
    	spar = spar + ',alpha='+strtrim(alpha,2)
    if (strupcase(db_info('name',0)) eq 'MARCS') then $
    	spar = spar + ',vturb='+strtrim(vturb,2)+',type='+type
    lst = dbfind(spar,/silent)
    dbext,lst,'log_g',log_g
    log_gs = log_g(rem_dup(log_g))

    log_g_int = log_g_obj
    if (log_g_int lt min(log_gs)) or (log_g_int gt max(log_gs)) then begin
       print,'Object is outside of range of tabulated Log_g at '+ $
         'Teff = '+strtrim(teff,2)
       log_g_int = log_g_int < max(log_gs) > min(log_gs)
       outside = 1
    end
    if n_elements(log_gs) gt 1 then tabinv,log_gs,log_g_int,rindex $
    			       else rindex = 0.0

    log_g1 = max(log_gs(fix(rindex[0])))
    log_g2 = min(log_gs(ceil(rindex[0])))
;
; find flux at each log_g and linearly interpolate between them
;

    model_int2,teff,log_g1,log_z_obj,flux1,irange=irange, $
       quadratic=quadratic,spline=spline,vturb=vturb, $
       scalar_item,scalar_value1, alpha=alpha, outside=outside, type=type

    if log_g1 ne log_g2 then begin
       model_int2,teff,log_g2,log_z_obj,flux2,irange=irange, $
         quadratic=quadratic,spline=spline,vturb=vturb, $
         scalar_item,scalar_value2, alpha=alpha, outside=outside, type=type
       frac1 = (log_g2 - log_g_int) / (log_g2-log_g1)
       frac2 = (log_g_int - log_g1) / (log_g2-log_g1)
       flux = frac1*flux1 + frac2*flux2
       if n_elements(scalar_item) eq 1 then $
         scalar_value = frac1*scalar_value1 + $
                 frac2*scalar_value2

        end else begin
       flux = flux1
       if n_elements(scalar_item) eq 1 then $
         scalar_value = scalar_value1
    end
end

pro model_int2,teff,log_g,log_z_obj,flux,scalar_item,scalar_value, $
       irange=irange,quadratic=quadratic,spline=spline,vturb=vturb, $
	alpha=alpha,outside=outside, type=type

;
; routine to find flux at a tabulated temperature and log_g by
; interpolation in Log_z (linear, spline or quadratic)
;
    spar = 'teff='+strtrim(teff,2)+',log_g='+strtrim(log_g,2)
    if (strupcase(db_info('name',0)) ne 'LEJ_ORIG') and $
       (strupcase(db_info('name',0)) ne 'MARCS') and $
    	(strupcase(db_info('name',0)) ne 'PHOENIX') then $
    	spar = spar + ',vturb='+strtrim(vturb,2)+',alpha='+strtrim(alpha,2)
    if (strupcase(db_info('name',0)) eq 'PHOENIX') then $
    	spar = spar + ',alpha='+strtrim(alpha,2)
    if (strupcase(db_info('name',0)) eq 'MARCS') then $
    	spar = spar + ',vturb='+strtrim(vturb,2)+',type='+type
    lst = dbfind(spar,/silent)
;
; I don't remember why I dont use log_z=-2.5
;
    dbext,lst,'log_z',log_z
    if strupcase(db_info('name',0)) eq 'LEJ_ORIG' then begin
       good = where(log_z ne -2.5)
       log_z = log_z(good)
       lst = lst(good)
    end

    sub = sort(log_z)
    log_z = log_z(sub)
    lst = lst(sub)

    if n_elements(scalar_item) ne 1 then dbext,lst,'flux',flux_array $
       else dbext,lst,'flux,'+scalar_item,flux_array,svalues
    log_z_int = log_z_obj
    if (log_z_int lt min(log_z)) or (log_z_int gt max(log_z)) then begin
       print,'Object is outside of range of tabulated Log_z at '+ $
         'Teff = '+strtrim(teff,2)+'    Log_g='+ $
         strtrim(log_g,2)
       log_z_int = log_z_int>min(log_z)<max(log_z)
       outside = 1
    end

    if n_elements(log_z) gt 1 then tabinv,log_z,log_z_int,rindex $
    			      else rindex = 0.0
    i1 = fix(rindex[0]) & i2 = ceil(rindex[0])
    log_z1 = log_z(i1)
    log_z2 = log_z(i2)
;
; interpolate between log_z1 and log_z2
;
    if i1 ne i2 then begin
        if (not keyword_set(spline)) and $
               (not keyword_set(quadratic)) then begin
         flux1 = flux_array(irange[0]:irange[1],i1)
         flux2 = flux_array(irange[0]:irange[1],i2)
         frac1 = (log_z2 - log_z_int) / (log_z2-log_z1)
         frac2 = (log_z_int - log_z1) / (log_z2-log_z1)
         flux = frac1*flux1 + frac2*flux2
         if n_elements(svalues) gt 0 then $
          scalar_value = frac1*svalues(i1) +  $
                     frac2*svalues(i2)
          end else begin
             flux = fltarr(irange[1]-irange[0]+1)
             for i=irange[0],irange[1] do begin
             f = interpol(flux_array(i,*),log_z,log_z_int, $
                   spline=spline,quadratic=quadratic)
             flux[i-irange[0]] = f
         end
         if n_elements(svalues) gt 0 then scalar_value = $
                 interpol(svalues,log_z,log_z_int, $
                   spline=spline,quadratic=quadratic)
        end
      end else begin
       flux = flux_array(irange[0]:irange[1],i1)
       if n_elements(svalues) gt 0 then scalar_value = svalues(i1)
    end
end

pro model_int,teff_obj,log_g_obj,log_z_obj,wave,flux, $
    scalar_item,scalar_value,quadratic=quadratic, $
    spline=spline,linear=linear,wrange=wrange,castelli=castelli, $
    vturb=vturb,new_castelli=new_castelli,d01stis=d01stis, $
    alpha=alpha,phoenix=phoenix, type=type,marcs=marcs

    if n_params(0) lt 1 then begin
       print,'model_int, teff_obj, log_g_obj, log_z_obj, wave, flux
       print,'KEWORDS: /linear, /quadratic (default), /spline
       print,'         wrange=[wmin,wmax],  /castelli, /new_castelli
       print,'         /d01stis, /phoenix, /marcs
       print,'         vturb=vturb, alpha=alpha, type=type
       return
    end
;
; set defaults
;
    if n_elements(quadratic) eq 0 then quadratic = 0
    if keyword_set(linear) then quadratic=0
    if keyword_set(spline) then quadratic=0
    dbase = 'lej_orig'
    dbase = 'castelli'				; 09feb5 - make default rcb
    if keyword_set(castelli) then dbase=castelli  ; 09sep22-for rcb cas='marcs'
;    if keyword_set(castelli) then dbase = 'castelli'
    if n_elements(vturb) eq 0 then vturb=2
    if keyword_set(new_castelli) then $
        dbase = 'castelli_new'
    if keyword_set(d01stis) then dbase = 'd01stis'
    if keyword_set(phoenix) then dbase = 'phoenix'
    if keyword_set(marcs) then dbase = 'marcs'
    if n_elements(alpha) eq 0 then $
    	if log_z_obj le -1 then alpha='a' else alpha='n'
    if n_elements(type) eq 0 then begin
    	if alpha eq 'a' then if log_g_obj gt 3 then type = 'p_st' else type = 's_st'
    	if alpha eq 'n' then if log_g_obj gt 3 then type = 'p_ap' else type = 's_ap'
	if log_z_obj ge 0 then if log_g_obj gt 3 then type = 'p_st' else type = 's_st'
    end

;
; open data base
;
    dbopen,dbase
    if (dbase eq 'marcs') then begin
        lst = dbfind('vturb='+strtrim(vturb,2)+ $
			',type='+strtrim(type,2),/silent) 
      end else begin

    	if (dbase ne 'lej_orig') and (dbase ne 'phoenix') then begin
    		lst = dbfind('vturb='+strtrim(vturb,2)+ $
			',alpha='+strtrim(alpha,2),/silent) 
	      end else begin
		  if dbase eq 'phoenix' then $
		  	lst = dbfind('alpha='+strtrim(alpha,2),/silent) $
		  	else lst = -1
	end
    end
    dbext,lst,'teff',teff
;
; get wavelength vector
;    
    if keyword_set(d01stis) or keyword_set(phoenix) $
    			    then get_d01stis_wave,wave $
; rcb-09sep22    			    else get_kur_wave,wave
    			    else get_castkur04_wave,wave
    if keyword_set(wrange) then begin
       tabinv,wave,wrange,irange
       irange=round(irange)
       wave = wave(irange[0]:irange[1])
       end else begin
           irange = [0,n_elements(wave)-1]
    end
;
; find the temperature above and below the input object's temp.
;
    outside = 0
    teffs = teff(rem_dup(teff))
    good = where(teffs ne 5777)    ;remove sun
    teffs = teffs(good)
    teff_int = teff_obj
    if (teff_int lt min(teffs)) or (teff_int gt max(teffs)) then begin
      print,'Object is outside of range of tabulated temperatures'
       teff_int = teff_int > min(teffs) < max(teffs)
       outside = 1
    end
    tabinv,teffs,teff_int,rindex

    t1 = max(teffs(fix(rindex[0])))
    t2 = min(teffs(ceil(rindex[0])))
;
; find flux at each temperature and linearly interpolate between them
;
    model_int1,t1,log_g_obj,log_z_obj,flux1,irange=irange,spline=spline, $
       quadratic=quadratic,vturb=vturb,scalar_item,scalar_value1, $
       alpha=alpha,outside=outside, type=type
    if t1 ne t2 then begin
       model_int1,t2,log_g_obj,log_z_obj,flux2,irange=irange, $
         spline=spline,quadratic=quadratic,vturb=vturb, $
         scalar_item,scalar_value2, alpha=alpha, outside=outside, type=type
       if n_elements(flux2) eq 1 then begin
         flux = wave*0
         goto,done
       end
       frac1 = (t2 - teff_int) / float((t2-t1))
       frac2 = (teff_int - t1) / float((t2-t1))
       flux = frac1*flux1 + frac2*flux2
       if n_elements(scalar_item) eq 1 then $
         scalar_value = frac1*scalar_value1 + $
                 frac2*scalar_value2
        end else begin
       flux = flux1
       if n_elements(scalar_item) eq 1 then $
         scalar_value = scalar_value1
    end
done:
    if outside eq 1 then !err = -1 else !err = 0
end
