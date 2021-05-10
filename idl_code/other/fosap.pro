FUNCTION FOSAP,head,gpar,DETin,GRAT,APER,WL,COEF,WMIN,WMAX,NAV
;+
;
;PURPOSE:
;      93OCT5-GET FOS APER CORRECTIONS FROM DISK$DATA2:[BOHLIN.FOSapr]APCOEF.TAB
;	94MAR10-GET FOS APER CORR FROM DISK$DATA2:[BOHLIN.COapr]APCOEF.TAB
;	94jun8-GET FOS APER CORR FROM disk$data12:[colina.recalib]apcoef.tab
;
; CALLING SEQUENCE:
;	FOSAP(head,gpar,DETIN,GRAT,APER,WL,COEF,WMIN,WMAX,NAV)
; INPUT:
;	head-fos data header
;	gpar  -fos group parmeters---not used and not required to be defined
;	DETin-'BLUE' OR 'AMBER' [or 'RED']
;	GRAT-FOS GRATING, EG. 'H19'
;	APER-FOS SMALL APERTURE, EG. 'B3'.. A1(4.3) {FOSAP=1} ALLOWED-94jan28
;	WL-WAVELENGTH RANGE TO CALC FOSAP, IF THIS RANGE IS GREATER THAN 
;		RANGE COVERED IN APCOEF.TAB SET FOS APER CORR=1000.
;OUTPUT
;	FOSAP-FUNCTION VALUE IS THE APERTURE CORRECTION VECTOR CORRES TO WL
; OPTIONAL OUTPUT:
;	COEF-3X1 OR 3X2 ARRAY OF CORRECTION COEFFICIENTS FROM APCOEF.TAB
;	WMIN-ARRAY OF DIM 1 OR 2 OF MIN WL RANGE FROM APCOEF.TAB
;	WMAX-ARRAY OF DIM 1 OR 2 OF MAX WL RANGE FROM APCOEF.TAB
;	NAV-NUMBER OF INDEPENDENT MEAS THAT ARE AVERAGED TO GET COEF
; HISTORY
;	94MAR10-ADD HEADER (AND GROUP) PARAM TO CALL, TO DISTINGUISH POST-COSTAR
;-

det=strtrim(detin,2)
IF det eq 'RED' then det='AMBER'		;HARMS PATCH
IF det eq 'R' then det='AMBER'	
IF det eq 'B' then det='BLUE'	
IF SXPAR(HEAD,'KYDEPLOY') THEN begin
	SUBDIR='COAPR'  &  print,'post-costar aperture corrections'
      endif ELSE begin
	SUBDIR='FOSAPR'  &  print,'pre-costar aperture corrections'
	endelse
TABLE_EXT,'disk$data12:[colina.recalib]apcoef.tab',			$
	'DETECTOR,FGWA_ID,APER_ID,C0,C1,C2,WMIN,WMAX,NUM_AVG',DETS,     $
	GRATS,APERS,C0,C1,C2,WMINS,WMAXS,NAVS
dets=strtrim(dets,2) & grats=strtrim(grats,2) & apers=strtrim(apers,2)
; TRIM DASH FOR POST COSTAR FILE NAMES
; 94dec7-why??? comment out:
;IF SUBDIR EQ 'COAPR' THEN APERS=STRMID(APERS,0,1)+STRMID(APERS,2,1)

; select desired matrix elements:
; 94FEB1-USE JUST FIRST 3 CHAR OF APER NAME
MODE=WHERE((DET EQ DETS) AND (GRAT EQ GRATS) AND (STRMID(APER,0,3) EQ APERS),  $
								COUNT)

IF COUNT LE 0 THEN BEGIN	;IDIOT CHECK
	PRINT,'IDIOT CK: NON-SUPPORTED MODE IN FOSAP=',DET,GRAT,APER
	STOP
	ENDIF

APCORR=WL*0+1000.		;INITIALLIZE ALL CORRections TO 1000
ISEGS=N_ELEMENTS(MODE)		;NUMBER OF WL SEGMENTS (1 OR 2)
COEF=FLTARR(3,ISEGS)
WMIN=FLTARR(ISEGS)
WMAX=WMIN
NAV=INTARR(ISEGS)

FOR I=0,ISEGS-1 DO BEGIN
	IM=MODE(I)		;INDEX OF SELECTED COEF MODE LINE IN TABLE
	COEF(0,I)=C0(IM)
	COEF(1,I)=C1(IM)
	COEF(2,I)=C2(IM)
	WMIN(I)=WMINS(IM)
	WMAX(I)=WMAXS(IM)
	NAV(I)=NAVS(IM)
	GOOD=WHERE((WL GE WMIN(I)) AND (WL LE WMAX(I)))
	APCORR(GOOD)=C0(IM)+C1(IM)*WL(GOOD)+C2(IM)*WL(GOOD)^2
	ENDFOR

RETURN, APCORR
END
