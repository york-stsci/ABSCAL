PRO FOSFLX,COUNTRATE,FLUX,GRAT,APER
;+
;
; FOSFLX,COUNTRATE,FLUX,GRAT,APER
;
; PURPOSE:
;	CALIBRATE COUNTRATE TO FLUX FOR ONLY A1 IN CYCLE 1-3 AND
;		ALSO INCLUDE APER CORR FOR CYCLE 4-N.
; INPUT - COUNTRATE
;       - GRAT=GRATING, EG. 'AH19' FOR RED H19 OR 'BH13' [RH19 FOR POST-COSTAR]!
;	- APER REQUIRED FOR POST-COSTAR AND FORBIDDEN FOR PRE-COSTAR, EG 'B3'
; OUTPUT - FLUX=calibrated FOS FLUX
;        - approx. pre-costar fluxes are avg IVS of 92mar30 *** USE FOS_PROCESS
;
; AUTHOR-R.C.BOHLIN
; HISTORY:
; 90DEC26-TO MULT.COUNTRATE BY INV. SENS. TO GET FLUX
; 91MAR11-QUICK PATCH TO GET A1 MASTER SENS
; 92dec21-change pointers to jdn master cal of cal/fos 77
; 93feb25-actually achieved above goal
; 94MAR12-MOD FOR COSTAR IVS FILES
; 94AUG12-GO TO lindler AVG A-1 IVS AND APER CORRECTIONS IN .COAPR---but
;	NOT tested... use fos_process, instead.
;-
IF N_PARAMS(0) EQ 0 THEN BEGIN
	PRINT,'FOSFLX,COUNTRATE,FLUX,GRAT'
	PRINT,'typical example: FOSFLX,C,FLX,"BH13","B-3"  POST-COSTAR'
	RETALL
	ENDIF
;
;READ MASTER INVERSE SENSITIV FILES
;

IF N_PARAMS(0) LE 3 THEN BEGIN
; CYCLE 1-3 dir w/ all the NBH AVG. calib files USED FOR PODPS AND ARCHIVE:
	sxopen,1,'DISK$DATA2:[bohlin.30mar92]'+GRAT+'A-1.R2H'
	IVS=SXREAD(1)
      ENDIF ELSE BEGIN
	
sxopen,1,'DISK$DATA10:[lindler.ave_ivs]'+GRAT+'A1.R2H',HEAD
	PRINT,'CYCLE 4 CAL FILE=',GRAT+APER+'.R2H'
	SXADDPAR,HEAD,'KYDEPLOY','T'            ;TO GET POST-COSTAR
	GET_TWAVE,HEAD,WAVE
	IVS=SXREAD(1)/FOSAP(HEAD,GPAR,STRMID(GRAT,0,1),STRMID(GRAT,1,3),       $
				apER,Wave)
	ENDELSE

	FLUX=COUNTRATE*IVS
	CLOSE,1
RETURN
END
