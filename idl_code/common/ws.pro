FUNCTION WS,WAVE,WINDEX
;+
;
; CALLING SEQUENCE:  RESULT=WS(WAVE,WINDEX)
;
; 93JUL16 TO COMPUTE INDICES OF VECTOR WAVE AT WINDEX WL POINTS
; windex may be a vector
;-
TABINV,WAVE,WINDEX,INDEX
RETURN,index
END