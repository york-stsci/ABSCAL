pro PSFSUB,xc,yc,PEAK,IM,PSF,PSFIM,RESID,frac,ABFRAC
;+
;
; CALLING SEQUENCE:
;      PSFSUB,xc,yc,PEAK,IM,PSF,PSFIM,RESID,frac,ABFRAC
; INPUT
;       XC,YC  - center of STARS in pixels (VECTORS)
;	PEAK   - DN VALUE TO SCALE PSF (VECTOR CORRESPONDING TO XC,YC) 
;	im     - 2-D array containing image
;	PSF    - SQUARE PSF IMAGE WITH SAME PX SIZE AND CENTER AT SIZE(PSF)/2
;
; OUTPUT:
;	PSFIM  - IMAGE OF SAME SIZE AS IM AND CONSISTING OF SCALED PSF SOURCES
;	RESID  - IM-PSFIM
;	frac   - fraction of image remaining after subtracting psf (=resid/im)
;
; ***NOTE - STELLAR PHOTOMETRY IS INDEPENDENTLY SPECIFIED BY THE PEAK VECTOR
;
; HISTORY:
;
; 93JUN19 - RCB  
; 93jun28 - add frac computation
;-
;-----------------------------------------------------------------------------
; check # of parameters
if (n_params(0) eq 0) then begin
  print,'CALLING SEQUENCE: PSFSUB,xc,yc,PEAK,IM,PSF,PSFIM,RESID,frac'
  retall
  endif

Psiz=size(PSF)
CENT=Psiz(1)/2	;PSF CENTER
ISIZ=SIZE(IM)
XSIZ=ISIZ(1)	;IMAGE SIZE
YSIZ=ISIZ(2)
PSFIM=FLTARR(XSIZ,YSIZ)
SIZ=SIZE(XC)
NSTAR=SIZ(1)


FOR I=0,NSTAR-1 DO BEGIN

;COMPUTE X,Y RANGE COVERED BY PSF IN IM
	IXBEG=(XC(I)-CENT)>0
	IXEND=(XC(I)+CENT-1)<(XSIZ-1)
	IYBEG=(YC(I)-CENT)>0
	IYEND=(YC(I)+CENT-1)<(YSIZ-1)
;COMPUTE CORRESPONDING X,Y RANGE IN PSF IMAGE
	PXBEG=CENT-(XC(I)-IXBEG)
	PXEND=CENT+(IXEND-XC(I))
	PYBEG=CENT-(YC(I)-IYBEG)
	PYEND=CENT+(IYEND-YC(I))
PRINT,'RANGE IN IMAGE=',IXBEG,IXEND,IYBEG,IYEND
PRINT,'RANGE IN PSF=',PXBEG,PXEND,PYBEG,PYEND
	PSFIM(IXBEG:IXEND,IYBEG:IYEND)=PSFIM(IXBEG:IXEND,IYBEG:IYEND)+  $
		PSF(PXBEG:PXEND,PYBEG:PYEND)*PEAK(I)/PSF(CENT,CENT)
	ENDFOR

RESID=IM-PSFIM
frac=total(resid)/total(im)
ABfrac=total(ABS(resid))/total(im)
print,'remaining fraction of image=',frac,'  ABS FRAC=',ABFRAC

return
end