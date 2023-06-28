;+
; PREPROC.PRO - rcb
; 97mar17 - run new data thru calstis
; 97apr9  - change from crj input files to raw, because spectrum shifted on re-
;	acq of o3tt21040 2nd spectrum (apr4 obs of gd153) on G750L
; 97aug25 - change to def of extr.slit=11 for g140l, and 7 for ccd
; 97aug26 - install new getflat41 to preserve photometry
; 97nov18 - replace getflat41 by findfringe to 'really?' preserve photometry
; 97nov25 - install my contemporaneous G750L corr proceedure.
; 97dec9  - add stellar radial velocity to header
; 98jan16 - permit default aperture (wrong transm), but good abscor extr hgt cor
; 98mar12 - change echelle extr hgt default from 11 to 7
; 98SEP2  - quit making the stellar radial vel. correction
; 98sep10 - keep default heliocentric wls, except for extracting G750 tung flats
; 99jul15 - don helped me change to the fancy autowavecal that uses before & 
;		after wavecals, unless i have a non-zero stiswlfix, where i 
;		should get same answer as before from using just one wavecal.
;		Resulting WLs for 2 test cases are 0.1-.2 PX different.
; 99nov27 - above NG, as i cannot fix new WL problems. Update to always use 2
;				wavecals as available.
; 99aug26 - update bkg subtr from default to match abscor in 2 main calstis
; 00nov21 - add fix for mult mama readouts. See also [.stis.doc]calstis.doc
; 02jan29 - convert to unix
; 02mar18 - correct +3" G140L WLs by -0.15px & -3" by +0.15px. 
;		because MSM move rotates slit & offsets the spectra from the
;		assumed center of slit.  See stisdoc/dispersion.doc
; 02jun21 - elim median filter for MAMA backgrounds
; 02sep11 - start using targoffs param
; 02nov   - testing CCD L-flats
; 03jul15 - section for special processing of saturated Vega obs. w/o Stis_cr
; 03dec31 - make gross=net+bkg for G750L to simplify ctecorr
; 04jun10 - stray light corr for BD+17d4708 G230LB, shortward of 1802A
; 05sep20 - implement G230LB flat field from proffitt
; 06feb8  - Looks like echscat.pro is preferred for echelle spectra, as that
;	writes also a newspec*, corr for scatt lite. Use Nothing from /direchl
; 06feb16 - new lo-disp trace files installed.
; 06apr13 - new G430L and G750L L-flats installed.
; 06may26 - Mods for E1 etc to fully implement L-flats for 4 L-modes. G230L
;		is ok w/ unit (no) L-flat. See Proffit 2005 Workshop.
; Get rid of E1 standard WL offsets, as the positioning error in the slit
;	decreased to ~0 after ~2003nov
; 06jul31-implement correction to net for GD153 G750L early WIDE o3tt4* obs &
;	remove 14sep16 per gwidth=11 work.
; 09aug6 - update stis_cr to make average CCD temp(OCCDHTAV) from all extents 
;	& rerun all CCD 09sep8
; 11Mar28-Implement GAC file per doc.lpflats to fix 52x2E1 obs for G230LB & ;											G430L
;	Use orig Proffitt corr for other apers. & aug4-fix the Jun20 E1 error
; 11jun20-add scat lite corr for G230LB, 52X2. See /stis/doc/scat.ccdmodes
; *** PROCEDURE: -- See make-tchang.pro ***
; UPDATE-bias,dark, & hot px stuff.see stisdoc/hot.pixel 
; For changes see also: calspec/next.deliv
; *** STATUS ***
; 2011mar29-30 - all low disp
; 2011Aug4 - Run all of G*-E1, ie G230L, & G230LB
; 2013Jan31-gear up for Sirius
; 2013feb5-fix missing comma in G750L fringe flat that made hrepair=-.04 fail,
;			ie in ...ctecorr=0,'+....
; 2013feb7-new abscor=pcttab files w/ more nodes for G750L for Vega & Sirius
; 2013feb11-New SM4 FF. See sirius/doc.procedure.
;           Also update G230LB scat lite corr to use CTE corr Net.
; 2013jul29-New CCD FF w/ Useafter=2011Oct1. See sirius/doc.procedure.
; 2014sep18 - Make simpler G750L corr for early GD153 o3tt4* obs.
; 2014sep29 - Make G750L sens11 cal and implement via  update here.
; 2014Oct8 - Make G750L sens11 cal the default by doing gwidth=11 for dat/spec*
; 2019mar19 - Choose GAC by Gwidth for G230LB & G430L.
; 2019mar26 - Try gwidth=11 GAC for gwidth=7: G230LB-slight improvement - OK.
;	& G430L-NG Go back to gwidth=7 GAC for gwidth=7, which is fine.
; 2019apr25-gwidth=11 the default extr for all low disp, ie change G230LB & 430L
; 2020feb - Mods for 2 stars in each CCD obs of Jesus ... Maiz
; 2020nov10-bkg drops off near E1. Change to bdist=+/-90 for E1 & G430L
;	Gaia star sep=127,442,503,877,138. More detail at calspec/doc/
;	hot * GD153 G430L-E1 change for odud02050 is: +0.1-0.3%, except below
;	~2950A, where the increase shoots up. Rest not yet reprocessed till dec7
; 2020nov23-Try bdist=90 for G230LB, as well. & for anything w/i 300px of edge
;  RESULTS   star    exp(s) type     min     max  change (%) 
;odta14010  FEIGE34   1170   sdO       0       0
;odta60010  HD167060  1100   G3V       0      >+10% at <2000A .001 of peak flux
;odud01030  GD71       580   DA1.5     0        2 at 1650A
;obnl08020  HD37962    850   G2V       0      >-10% at <1800A
;oc3i01010  HD009051  1200   G7III     0      >-10% at <1800A
;oc3i07020  HD185975  1650   G3V       0      >-10% at <1920A flux more Neg!
;  ---> These G230L changes are NOT significant. All small compared to noise, eg
;		the amount of negative flux. KEEP G230LB 300px bkg offsets.
;		See also calspec/doc/changes.deliv 2020nov-dec.
;	Avila G430L SDSS132811 odqg01020 changes by -5% @ ~3100A to +3% @ 3800A,
;		for Bdist 50-->90--OK, mrgpt=3150A.
; 2020dec2-Could apply GAC corr above row ~800, ie e1>(800-511)=289. But there 
;	would be NO changes to Jesus obs.
; 2020dec4-Do G430M-E1 w/ bdist=90, as well.
; 2021feb8 - mods to process CCDflats for CTE revisions 
;-

@../stisidl/stispro/extast.pro	; new one in nidl fails for MAMAs
alfbet=['1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y',',z']
nowavcal=0

tbeg=systime(1)			; start time in sec
!x.style=1
!y.style=0
!p.thick=1
wnodes=[5240,5265,5300,5318,5340,5380,5420,5460,5500,5580,5700,5800,5840,6000, $
        6200,6300,6450,6550,6600,6620,6700,6800,6900,7000,7150,7400,7600,      $
        7800,8000,8300,8500,9100,9500,9800,10200] ; for G750L fringe removal
; NGC6681 - See pre-ngc.pro
; ff are for full reprocessing of all STD star data:
fil63  =findfile('../data/spec/7063/o*raw.fits')	;CCD G..L
fil64=[		'../data/spec/7064/o4dd06020_raw.fits',   $	;GD71,140L
		'../data/spec/7064/o4dd07020_raw.fits'] 	;230L 25mama
stisobs,'dir7064.log',fil64,grat,aper,star,'','',''	     ;all of GRW monitor
fil64='../data/spec/7064/'+(strmid(fil64,9,9))+'_raw.fits'	;-spec+raw
filpg=	findfile('../data/spec/7064/o4pg*_raw.fits')    ;gd71 offsets
;fil64=[fil64,filpg(0:11)]	; elim 52x.05 G230L screwup files

fil96=[	'../data/spec/7096/o3zx08hfm_raw.fits',	$;G140-wavecal
	'../data/spec/7096/o3zx08hhm_raw.fits',	$;G140-gd153
	'../data/spec/7096/o3zx08hjm_raw.fits',	$;g230L-wavecl
	'../data/spec/7096/o3zx08hlm_raw.fits']	;g230L-gd153
;fils=[	'../data/spec/7095/o40301060_raw.fits',     $;G230LB +2"gd153
; 7095 vignetting proposal: Done at +-2" not 0, except for P041C on G750L,M
;	the G750L is -8975, instead of 7751...sigh! But see fil95 below 00dec26
;postarg=+2.0" ff. 750L are NG photometrically, as no short slit fringe flat.
; no flt'../data/spec/7095/O403020J0_raw.fits',  $;G750L pos=2" p041c
;no flt	'../data/spec/7095/O403020K0_raw.fits'] ;G750L  pos=-2" p041c
; 00dec26 incl G430L pos -22 to -6" 
stisobs,'dir7095.log',fils,grat,aper,star,'G430L','52X2','P041C'
;fils='../data/spec/7095/'+(strmid(fils,9,9))+'_raw.fits';-'spec'+'raw
fil95=findfile('../data/spec/7095/o*raw.fits')
; 7097 vignetting proposal: CENTRAL POSITIONS
fil97=[  '../data/spec/7097/O43J01Q2M_raw.fits',	$ ;G140-wavecal
        '../data/spec/7097/O43J01QAM_raw.fits']    ;G140 gd153,pos=3
;	'../data/spec/7097/O43J01QOM_raw.fits',	$ ;g230L-wavecal
;BAD-low'../data/spec/7097/O43J01QYM_raw.fits',   $; 230L g153,pos=.8

;       '../data/spec/7097/O43J02GBM_raw.fits',$ ;g140m-1272 wavecal
;	'../data/spec/7097/O43J02GLM_raw.fits',$ ;g140m gd153,pos=3
;	'../data/spec/7097/O43J03YXM_raw.fits',$ ;g140m-1518 wavecal
;	'../data/spec/7097/O43J03Z7M_raw.fits',$ ;g140m gd153,pos=0
;	'../data/spec/7097/O43J04E2M_raw.fits',$ ;g230m-1851 wavecal
;	'../data/spec/7097/O43J04EEM_raw.fits',$ ;g230m gd153,pos=0
;	'../data/spec/7097/O43J05EGM_raw.fits',$ ;g230m-2659 wavecal
;	'../data/spec/7097/O43J05EUM_raw.fits']  ;g230m gd153,pos=0

stisobs,'dir7100.log',fil00,grat,aper,star,'','',''
fil00='../data/spec/7100/'+(strmid(fil00,9,9))+'_raw.fits';-spec+raw
fil21=findfile('../data/spec/7721/o*_raw.fits')
fil42=findfile('../data/spec/7642/o*_raw.fits')
fil72=findfile('../data/spec/7672/o*_raw.fits')     ;AGK M&L sens mon
stisobs,'direchl.log',filech,grat,aper,star,'','',''		;echelle data
filech='../data/spec/echl/'+(strmid(filech,9,9))+'_raw.fits';-spec+ra
fil32=[findfile('../data/spec/7932/o*_raw.fits')]	;purity, transm
fil32=[	'../data/spec/7932/O4SP03060_raw.fits',  $ ;grw,g430l,52X2
	'../data/spec/7932/O4SP03070_raw.fits',  $ ; grw,g430l,50ccd
	'../data/spec/7932/O4SP01060_raw.fits',  $ ; GD71,g140L,52X2
	'../data/spec/7932/O4SP01070_raw.fits']    ;GD71,g140L,25Mama
fil74=[findfile('../data/spec/7674/o*_raw.fits')]	; IR std's

stisobs,'dir2009.log',fil2009,grat,aper,star,'','',''	; gets G750L flats!!!
fil2009='../data/spec/2009/'+(strmid(fil2009,9,9))+'_raw.fits
filsci=findfile('../data/spec/sci/o*_raw.fits')  	; bianchi,uit340,etc.
filsci=filsci(where(strpos(filsci,'o6dj06') lt 0 and strpos(filsci,'o6dj07') lt 0))		; see pre-m33

stisobs,'dirnewl.log',filnl,grat,aper,star,'','',''
filnl='../data/spec/newl/'+(strmid(filnl,9,9))+'_raw.fits';-spec +raw

; Vega
; 2019apr12 - looks like 7pix files are properly done w/ gwith=11. I could fix
;	the misnomer, however. These .7[ix files are not USED in stisreduce.
fveg=filnl(where(strpos(filnl,'o8i105') ge 0 or strpos(filnl,'o8i106') ge 0))
; call the sat. 430L E1 o8i105050 bad, as noisy, wrong shape, odd overscan bump.
; add std 11px extractions for the 2 g750L + short slit flats:
;  using std. PCTTAB ='abscor-750.fits of 2002 in scal
; The 2 G750L w/ small extr hgt=11. #6 & 12 are the fringe flats:
fveg7=[fveg(5)+'7pix',fveg(6),fveg(11)+'7pix',fveg(12)]	; G750L Vega
filnl=[filnl,fveg7]			; 2014sep12 - fix to match Sirius
; Sirius
fsir=fil2009(where(strpos(fil2009,'obto11') ge 0 or strpos(fil2009,'obto12') ge 0))
; additional std 7px extractions for g750L + short slit flats:
fsir7=[fsir(5)+'7pix',fsir(6),fsir(13)+'7pix',fsir(14)]	; 7pix G750L Sirius
fil2009=[fil2009,fsir7]					; add new

; special narrow slit & Med disp. obs for Finley/Tremblay analysis:
;	See stisdoc/psf-lsf.wd-ccdobs
fils='../data/spec/newl/o6ig010'+['8','9','a','b']+'0_raw.fits' ;GD71 G430L
fils='../data/spec/newl/o6ig100'+['f','g','h','i']+'0_raw.fits' ;G191 g750L
fils='../data/spec/newl/o8v203010_raw.fits' 	;G191 g430L 52x0.1
;stisobs,'dirnewl.log',fils,grat,aper,star,'G430M','','G191B2B' ;few extra files
;stisobs,'dirnewl.log',fils,grat,aper,star,'G750M','','G191B2B'
;fils=['../data/spec/newl/'+(strmid(fils,9,9))+'_raw.fits'];-spec+raw

;fils=[	'../data/spec/7063/O3TT46030_raw.fits',		$
;	'../data/spec/7063/O3TT47030_raw.fits',		$
;	'../data/spec/7063/O3TT48030_raw.fits']  ;Tung flat for SMOV


; special new data & test:
;!path='../data/spec/test:'+!path			; SAVE example

stisobs,'dir2009.log',fils,grat,aper,star,'','',''		; Gaia - jesus
ind=where(strpos(fils,'oe3f') ge 0)			; 15816 ..Maiz
;	ind=where(strpos(fils,'oe3f06') ge 0)			; 587_0224,_8560
fils=fils(ind)
grat=grat(ind)  &  aper=aper(ind)  &  star=star(ind)
fils=(strmid(fils,9,9))
fils='../data/spec/2009/'+fils+'_raw.fits'		; end Gaia

;fils=fil2009						; last done ..
fils=[filnl,fil2009]					; ..........2021mar22

; everything:
;fils=[filnl,fil63,fil64,fil96,fil95,fil97,fil00,fil21,fil42,fil72,fil74,$
;	fil32,fil2009,filsci]	;last done 19aug6-7. 21jan

good=where(strpos(fils,'ocmv') lt 0 and strpos(fils,'ocy5') lt 0)		;omit Massa & Worthy 0.5" slit
fils=fils(good)
; ###change - fix star & aper for Gaia cases. Comment ff to do just Gaia
apertmp=strarr(n_elements(fils))  &  startmp=apertmp
ind=where(strpos(fils,'oe3f') gt 0,ngaia)

;*****
if ngaia ne n_elements(aper) then stop		; idiot ck
;*****

apertmp(ind)=aper  &  startmp(ind)=star		; insert gaia
aper=apertmp  &  star=startmp			; end everything

;fils=findfile('../data/spec/newl/o5i0010*0_raw.fits'); KEEP BASELINE
;fils=findfile('../data/spec/newl/o5i008020_raw.fits') ; G191 M Ha-em
;fils='../data/spec/2009/obvp080'+['2','3']+'0_raw.fits'
; HZ43
fil57=findfile('../data/spec/*/o57t*_raw.fits')
fil07=findfile('../data/spec/*/o69u07*_raw.fits')
fil08=findfile('../data/spec/*/o69u08*_raw.fits')
;	 fils=[fil57(0:3),fil57(10:13),fil07(0:3),fil08(0:3)]
; HZ43B
fil17=findfile('../data/spec/7674/o49x17*_raw.fits')
fil18=findfile('../data/spec/7674/o49x18*_raw.fits')
;	fils=[fil17(0:1),fil18(0:1)]

; small aper data for finley... hmmm... how does the fringe flat work here??:
;  stisobs,'dirnewl.log',fil19sm,grat,aper19sm,star,'G750L','','G191B2B'
;  stisobs,'dirnewl.log',fil71sm,grat,aper71sm,star,'G750L','','GD71'
;  fils=[fil19sm,fil71sm]  &  aper=[aper19sm,aper71sm]
;  good=where(strmid(aper,0,4) ne '52X2')
;  fils=fils(good)
;  fils=['../data/spec/newl/'+(strmid(fils,9,9))+'_raw.fits']
; gets G750L flats for star name =''   ONLY!!!:
; start new agk or grw   ##############################################
;stisobs,'dir2009.log',fils,grat,aper,star,'','','' ;new agk&grw gets fring w/ ''
;	fils=(strmid(fils,9,9))
;	fils=fils(where(strpos(fils,'oe36') ge 0 or strpos(fils,'oe5c') ge 0  $
;		or strpos(fils,'oefsl') ge 0))	;AGK&GRW
;	fils='../data/spec/2009/'+fils+'_raw.fits
; end new agk or grw   ##############################################
;fils='../data/spec/2009/o'+['ceim10b','ceil406']+'0_raw.fits'	; TEST for _flc
;fils='../data/spec/2009/o'+['ceim10b','ceil406']+'0_flt.fits'	; TEST for _flc
;fils='../data/spec/2009/o'+['ceim10b','ceil406']+'0_flc.fits'
;;fils='../data/spec/2009/o'+['ceim10b','ceil406']+'0_cte.fits'
;	fils=[fils,'../data/spec/2009/oceil4070_raw.fits']	; add fringe flt

;stisobs,'dirccdflat.log',fils,grat,aper,star,'','0.3X0.09','CCDFLAT'  ;CTE study
;	fils='../data/spec/ccdflat/'+(strmid(fils,9,9))+'_raw.fits'	
;fils=[findfile('../data/spec/ccdflat/ocmv18030_raw.fits')]	; gen purpose

;fils=[findfile('../data/spec/2009/oefrl2*_raw.fits')]		 ; gen purpose
;fils=[findfile('../data/spec/2009/oefrl1*_raw.fits')]		 ; gen
;fils=[findfile('../data/spec/2009/oehj02*_raw.fits'),		$ ; gen
;      findfile('../data/spec/2009/oehj03*_raw.fits')] ;,          $;purpose
;	findfile('../data/spec/2009/oe7q06020_raw.fits')]
fils= [findfile('../data/spec/2009/oddg01030_raw.fits'), '../data/spec/2009/oddg01040_raw.fits']
;fils=['../data/spec/newl/o8kh030'+['2','3']+'0_raw.fits',$
;	findfile('../data/spec/2009/ob6h*_raw.fits'),	$	; G750M
;	findfile('../data/spec/2009/oc33010*_raw.fits')]	; G430L
fils=findfile('../data/spec/2009/oe7q02020_raw.fits')
fils=[fils,'../data/spec/2009/oe7q02030_raw.fits']	; Fringe flat
;fils=[findfile('../data/spec/newl/o8kh03010_raw.fits'),		$
;	findfile('../data/spec/newl/o8kh03030_raw.fits')]	; fringe flat
;fils=findfile('../data/spec/2009/odck01050_raw.fits')
;fils=findfile('../data/spec/7932/o4sp010b0_raw.fits')		; g140M GD71
;fils='../data/spec/2009/'+['odta14010','odta60010','odud01030','obnl08020', $
;	'oc3i01010','oc3i07020']+'_raw.fits'
;fils=findfile('../data/spec/newl/o8v204*_raw.fits')		; gen purpose
;fils=findfile('../data/spec/2009/ocy5*_raw.fits')		; gaia hot stars

;stisobs,'dir2009.log',fils,grat,aper,star,'','',''  ;misses 0.3X0.09 w/ starname
;  fils='../data/spec/2009/'+(strmid(fils,9,9))+'_raw.fits
;  good=where(strpos(fils,'obnk01') gt 0 or strpos(fils,'oe7q01') gt 0	$
;	or strpos(fils,'oe7q02') gt 0) 				;wd0308
;  good=where(strpos(fils,'o6il01') gt 0 or strpos(fils,'o8h107') gt 0	$
;	or strpos(fils,'o8h108') gt 0 or strpos(fils,'o8h109') gt 0	$
;	or strpos(fils,'obbm01') gt 0)   			; LDS749B
;fils=fils(good)
;fils=[fveg,fsir] ; sat data

fils=strlowcase(fils)
print,'processing files:',fils

; 2021jan15 - set aper array for all runs (aper always defined, but not proper)
ind=where(strpos(fils,'oe3f') ge 0,nlow)		; 15816 Gaia
if nlow le 0 then aper=strarr(n_elements(fils))		; if NOT Gaia

st=''
; ###change:
subdir='cte'  &  gwidir='11'				; rare use	
;;;if strpos(fils(ifil),'_raw') lt 0 then subdir='cte'
subdir='dat'  &  gwidir='11'				; gwidth for this run
;subdir='test'  &  gwidir='11'
;subdir='hgt22'  &  gwidir='22'
;subdir='hgtp4'
;subdir='carlos/g430l-11'  &  gwidir='11'
;subdir='g750lhgt11'  &  gwidir='11'
;subdir='g750lhgt9'  &  gwidir='9'
;subdir='g750lhgt15'  &  gwidir='15'

; 2019mar26-make gac11 the default for G230LB. gac7 best for G430L,gwidth=7:
gac7=mrdfits('../stisidl/scal/ccd_gac7.fits',0,hgac70)	;2011Mar28 CCD aper corr
gac7=mrdfits('../stisidl/scal/ccd_gac7.fits',1,hgac7)	;2011Mar28 CCD aper corr
; 2019apr1-fix old err in make-gac.pro:
sxaddpar,hgac70,'filename','scal/ccd_gac7.fits'
gacwl7=gac7.wavelength
gacth7=gac7.throughput
gac=mrdfits('../stisidl/scal/ccd_gac11.fits',0,hgac0)	;2019Mar26 CCD aper corr
gac=mrdfits('../stisidl/scal/ccd_gac11.fits',1,hgac)	;2019Mar26 CCD aper corr
gacaper=strtrim(gac.aperture,2)
gacmode=strtrim(gac.opt_elem,2)
gacwl=gac.wavelength
gacth=gac.throughput
; GAC summary:
; ccd_gac7 - original Proffitt results, as updated by RCB for 52x2E1 28Mar2011
; ccd_gac11.fits-updated for 52x2E1 21Mar2019-Gen. use,except above G430L,gwid=7
; ccd_gac11.fits-updated for 52x2E1 2020dec1 for G230LB & G430L w/ <0.1% 
;	change, even for change to bdist=90px for G430L-E1


lastfil=''  &  crskip=''
; ###change:
for ifil=0,n_elements(fils)-1 do begin
;for ifil=1660,n_elements(fils)-1 do begin
	flg7px=0		; special for G750L vega to get 7px hgt.
	if strpos(fils(ifil),'7pix') ge 0 then begin
		flg7px=1			; 7px extract for vega & sirius
		fils(ifil)=replace_char(fils(ifil),'7pix','')	; .fits still OK
		endif
	sptpos=strpos(fils(ifil),'.fits')
	sptfil=fils(ifil)
	strput,sptfil,'spt',sptpos-3
	fits_read,sptfil,im,hspt,/header_only
	pstrtim=strtrim(sxpar(hspt,'pstrtime'),2)	;yyyy.ddd... format time
	date=strmid(pstrtim,0,8)
;fits_read does not pop all keywds!! ff triggers subarr msg for T/ACQ images.
	stis_read,fils(ifil),hd,infil,udl,heng1,heng2,errin,/hires
	
	mode=strtrim(sxpar(hd,'opt_elem'),2)
	gwidth=gwidir				; reset from the odd cases below
	if strpos(mode,'M') gt 0 then gwidth='7'; 202Aug31 for Lennon 16079 keep
	otemp='OM1CAT' & if strpos(mode,'230') gt 0 then otemp='OM2CAT'
	OMCAT=sxpar(hspt,otemp,comment=comment)	; 20011 raw data has om1cat
; 02sep11 use defaults - bkgfit='3'  &  if mode eq 'G140L' then bkgfit='-1'
; poly_fit bkg
	cenwave=strtrim(sxpar(hd,'cenwave'),2)
	targ=strtrim(sxpar(hd,'targname'),2)
	if targ eq 'HD172167-V6' then targ='HD172167'	;19mar28 - patch
	root=strlowcase(strtrim(sxpar(hd,'rootname'),2))
	if root eq 'o49x17010' or root eq 'o49x18010' then targ='HZ43B'  ;patch
	obsmode=strtrim(sxpar(hd,'obsmode'),2)
	postarg=sxpar(hd,'postarg2')
	aperture=strtrim(sxpar(hd,'propaper'),2)	; 02mar7 - was aperture
	gain = sxpar(hd,'ccdgain')
	if aperture eq '0' or aperture eq '' then aperture=strtrim(sxpar(hd,'aperture'),2);02mar26
	if targ eq 'NONE' then lastfil=fils(ifil)	; MAMA wavecals
	svel=starvel(targ)				; add to headers only.
; Lennon LS-V-+22-25 binary w/ diff phase obs (default=26.2km/s):
	if root eq 'oe9l01010' or root eq 'oe9l01020' then svel=+65.3
; ###change:
	if strpos(mode,'MIR') ge 0  or cenwave eq '8975'    		$
;		or (targ ne 'GD153' and targ ne 'GD71' and targ ne 'HZ43'    $
;		and targ ne 'G191B2B' and strpos(targ,'WD') lt 0)	$
;		and targ ne 'AGK+81D266')				$
		or targ eq 'NONE'					$
; ###change temp star:
;				or targ ne 'BD+75D325'			$
;				or targ ne 'GRW+70D5824' 		$
;				or targ ne 'AGK+81D266' 		$
;				or strpos(targ,'LDS749B') lt 0		$
;				or (targ ne 'VB-8' and strpos(targ,'2M') lt 0) $
;				or targ ne 'WD1657+343' 		$
;				or strpos(targ,'P') ne 0		$
; 03jul9 .31 changed to .5 for some Vega obs
; 16feb8 - Carlos HD146233 has postarg=0.525
; 16apr28 - must extract all postargs once to make nowavcal work. chg 1 to 199:
; 18jul   - SDSS132811, G750L has postarg=1
; 21jan21 - must extract all postargs once to make nowavcal work. chg 1 to 199:
; ###change temp 100, usual=1:
		or (abs(postarg) gt 1 and (postarg ne 3 or mode ne 'G140L') $
; 2014sep16-Allow P041C postargs ne 0, as best avail:
				and targ ne 'P041C' and			$
; Gaia 405_1056,6912 - uniq special test to cf center to off-center G750L
				root ne 'oe3f01030')	$	; permanent
			or obsmode eq 'ACQ/PEAK'			$
			or (targ eq 'CCDFLAT' and strpos(fils(ifil),'ccdflat') lt 0)	$
;			or mode ne 'PRISM'				$
;			or gain lt 1					$ ; CCD
;			or mode ne 'G140L'				$
; 2019mar - ff could be used to go back to incl. G750L E1;
			or (mode eq 'G750L' and targ ne '2M0036+18'	and strpos(targ,'GAIA') lt 0 and strpos(aperture,'E1') gt 0)		$
;			or mode ne 'G750M'				$
; Special for testing new E1 bdist=90
; del			or (mode ne 'G430L' or strpos(aperture,'E1') lt 0) $
			or (strpos(mode,'M') gt 0 and strmid(root,0,4) ne 'oe9l')	$ ; Gets Prism 
;			or strpos(mode,'G230') lt 0 			$
; allow few 0.5arcsec slit, otherwise consider only photometric 2" aper:
			or (strpos(aperture,'.') gt 0 and $
			    root ne 'obto03010' and strmid(root,0,4) ne 'oc8c' $ ;LCB
			    and strmid(root,0,4) ne 'od6j'		$ ;Clayt
			    and strmid(root,0,4) ne 'odi9'		$ ;Clayt
			    and strmid(root,0,4) ne 'ocy5'		$ ;Gaia
			    and strmid(root,0,4) ne 'oe9l'		$;Lennon
			    and strpos(fils(ifil),'ccdflat') lt 0)	$
; 02mar26 - 6x6 early aper not used for std *, just transm. Findfringe fails??
			or aperture eq '6X6'			$
			or date lt '1997.138'  $ ;pstrtime MAMA turn-on yr.ddd
			or fils(ifil) eq crskip		$ ;odd 14141 G.Worthy
; ##change - special run
;			or date lt '2018.180'  $ 		
	then begin
		print,'skipping:',root,' ',fils(ifil),' ',mode,' targ=',targ,' postarg=',postarg,' aper=',aperture,' at ',date
		crskip=''
		goto,skipit
		endif

	print,'processing file,mode=',fils(ifil),' ',mode,' targ=',targ,' postarg=',postarg,' aper=',aperture
; do wavecals when there was a targ/acq
	fdecomp,fils(ifil),disk,dir,name,qual
	if name eq 'o49x01010r5_raw' then begin
		name='o49x01010_raw'	;5th sat readout
		root='o49x01010r5'
		endif
; 98apr29 - No wavecal for .2x.09 E140H-1416 BD+75 - o4rq01gfq
;;;;;;;;if aperture ne '6X6' and strpos(fils(ifil),'o4rq01gfq') lt 0 then begin

	wfil=fils(ifil)
; 97nov10 - wavecal for tungsten for 7063 & 7805
	if targ eq 'NONE' then wfil=replace_char(fils(ifil),'30_','40_')
	if targ eq 'NONE' and strpos(wfil,'o4d1') ge 0 then	wfil=replace_char(fils(ifil),'40_','30_')
	if targ eq 'NONE' and strpos(wfil,'o49q') ge 0 then	wfil='../7642/o49q02060_raw.fits'
;2017apr	strput,wfil,'wav',strpos(fils(ifil),'raw')
	strput,wfil,'wav',strpos(fils(ifil),'.fits')-3
; patches for wavecal not being assoc w/ early MAMA data !!##$$%%:
	if strpos(fils(ifil),'o3zx09') ge 0 or 				$
		strpos(fils(ifil),'o4vt') ge 0 or 			$
		strpos(fils(ifil),'o5i002') ge 0 or 			$ ;01feb
		(strpos(fils(ifil),'o6ig') ge 0 and gain eq 0) or 	$ ;02jul
		strpos(fils(ifil),'o4pg')   ge 0 
	then wfil=fils(ifil+1) $
	else if gain eq 0 and name lt 'o45901010_raw' 	   $
	then wfil=lastfil ;early MAMA wavecals
	
	if strpos(fils(ifil),'o3zx08i9m') ge 0 then wfil='o3zx08ifm_raw.fits' ;98sep1 for Ed
; patches for "orphan" wavecal=old TRANS error in 2 cases
	if strpos(fils(ifil),'o4dd05lgq') ge 0 then	wfil='../data/spec/echl/o4dd05080_wav.fits'
	if strpos(fils(ifil),'o4pg02qcq') ge 0 then	wfil='../data/spec/echl/o4pg02qiq_raw.fits'
	if strpos(fils(ifil),'o4pg02qkq') ge 0 then	wfil='../data/spec/echl/o4pg02qmq_raw.fits'
	if strpos(fils(ifil),'o6ig01040') ge 0 then	wfil='../data/spec/newl/o6ig01bkq_raw.fits'
	if strpos(fils(ifil),'o6ig01050') ge 0 then	wfil='../data/spec/newl/o6ig01bnq_raw.fits'
	if strpos(fils(ifil),'o6ig01060') ge 0 then	wfil='../data/spec/newl/o6ig01bsq_raw.fits'
	if strpos(fils(ifil),'o6ig01b7q') ge 0 then	wfil='../data/spec/newl/o6ig01baq_raw.fits' ;2019may6
; propid=10039. SOME parts are hand tooled to have NON-assoc wavecals:
; GD71
	if strpos(fils(ifil),'o8v201040') ge 0 or strpos(fils(ifil),'o8v201050') ge 0 then wfil='../data/spec/newl/o8v201gmq_raw.fits' ;G430L
	if strpos(fils(ifil),'o8v201090') ge 0 then wfil='../data/spec/newl/o8v201h4q_raw.fits' ;G750L
	if strpos(fils(ifil),'o8v201g6q') ge 0 then wfil='../data/spec/newl/o8v201geq_raw.fits' ;G140L
	if strpos(fils(ifil),'o8v201guq') ge 0 then	wfil='../data/spec/newl/o8v201gzq_raw.fits' ;G230L
; GD153
	if strpos(fils(ifil),'o8v202030') ge 0 or strpos(fils(ifil),'o8v202040') ge 0 then wfil='../data/spec/newl/o8v202fzq_raw.fits' ;G230LB
	if strpos(fils(ifil),'o8v202050') ge 0 or strpos(fils(ifil),'o8v202060') ge 0 then wfil='../data/spec/newl/o8v202g8q_raw.fits' ; G430L
	if strpos(fils(ifil),'o8v202070') ge 0 then	wfil='../data/spec/newl/o8v202gcq_raw.fits' ;G750L
	if strpos(fils(ifil),'o8v202f7q') ge 0 then	wfil='../data/spec/newl/o8v202fgq_raw.fits' ;G140L
	if strpos(fils(ifil),'o8v202fhq') ge 0 then	wfil='../data/spec/newl/o8v202flq_raw.fits' ;G230L
	if strpos(wfil,strmid(name,0,6)) lt 0 then wfil='none'
;AGK, propid=9265, dir7672, no wavecals. Meas and put 16 offsets in stiswlfix!!
	if strpos(wfil,'7672/o6il02') ge 0 or strpos(wfil,'o6ig010c') ge 0 or strpos(wfil,'o6b6010') ge 0 then wfil='none'	; 02dec10 8891 
	if strpos(wfil,'6ig0108') ge 0 or strpos(wfil,'6ig0109') ge 0 or strpos(wfil,'6ig010a') ge 0 or strpos(wfil,'6ig010b') ge 0 then wfil='none'	;
	if targ eq 'CCDFLAT' then begin			;21feb-non-assoc wavcal
		pos=strpos(wfil,'_')-2
		pick=strmid(wfil,pos,1)
		indx=where(pick eq alfbet,n)
		if n le 0 then stop
		if n eq 1 then new=alfbet(indx(0)-1)
;?		if pick eq '1' then new='2'
		strput,wfil,new,pos
		wfilck=findfile(wfil)  &  wfilck=wfilck(0)
		if wfilck eq '' then begin
			new=alfbet(indx(0)-2)		; 2 eaarlier
			strput,wfil,new,pos
			wfilck=findfile(wfil)  &  wfilck=wfilck(0)
			endif
		if wfilck eq '' then begin
			new=alfbet(indx(0)+1)		; 1 later, eg. o6h3020g0
			strput,wfil,new,pos
			if root eq 'o6ig01070' then strput,wfil,'g',pos
			wfilck=findfile(wfil)  &  wfilck=wfilck(0)
			endif
		if wfilck eq '' then begin		; skip and/or odd cases
;			read,st
			nowavcal=nowavcal+1
			goto,skipit
			endif
		endif					; End CCDflats
	print,'#'
	print,'wavecal file=',wfil
	print,'#'

	if subdir eq 'hgtp4' then begin
		gwidth='15'
		if mode eq 'G230LB' or mode eq 'G430L' or mode eq 'G230MB' or mode eq 'G430M' then gwidth='11'
		endif

;eg o57t02020 w/ M star@+45px. 2020mar16-change from all=7px:
	if targ eq 'HZ43' and mode eq 'G750L' then gwidth='7'
	if root eq 'o69u07030' then gwidth='11'		; v.faint M* HZ43B F750L

;Massa dbl stars
	if targ eq 'HD281159' then gwidth='15'
	if targ eq 'HD73882' then gwidth='23'
;	gwidth='40'			; for the wavecal
;	fils(ifil)=wfil			;to process a wavecal
; Saturated Vega G230LB:
	if root eq 'o8i105020' or root eq 'o8i106030' then gwidth='84'
; Saturated Vega G430L:
	if root eq 'o8i105050' or root eq 'o8i105060' or root eq 'o8i106010' or root eq 'o8i106040' then gwidth='54'
;*****BY*****
; Saturated Vega G750L:
	if flg7px eq 0 and (root eq 'o8i105070' or root eq 'o8i106050')	then gwidth='40'
; Saturated Sirius G230LB:
	if root eq 'obto11010' or root eq 'obto11020' or root eq 'obto11030' or root eq 'obto12010' or root eq 'obto12020' or root eq 'obto12030' then gwidth='206'
; Saturated Sirius G430L:
	if root eq 'obto11040' or root eq 'obto11050' or root eq 'obto12040' or root eq 'obto12050' then gwidth='182'
; Saturated Sirius G750L:
	if flg7px eq 0 and (root eq 'obto11060' or root eq 'obto12060')	then gwidth='148'
	fixpar=''
; ###change:   97nov6 Echelle quick added kludge
	if strpos(wfil,'o40p01d4m') ge 0 or strmid(mode,0,1) eq 'E' then fixpar=',mincounts=50,minline=10'    ;weak mama1 wavecal
	toffix=stiswlfix(root,cenwave)		;special WL fixes.

; do fancy wavecal that uses before and after wavecals:
;	toffset=0				; for cases of NO WAVECAL below
	toffset={mjd:0,offsets:0,history:strarr(4)}	; 2010may5
	wfilck=findfile(wfil)  &  wfilck=wfilck(0)
	if wfilck eq '' and targ ne 'CCDFLAT' then wfil='none'	; 2020feb5
	if wfil ne 'none' then begin
; ff. must produce a set of 'calstis...' printout lines per each readout
		stis_woffset,wfil,trace=0,/noplot,toffset,hspc	; calls calstis
	end else begin
	     print,'*******WARNING: NO WAVECAL',format='(///25x,a////)'
; no wavecal ok for HD209458 and for 'none'
; 09jul27-a few other cases, like agk above, where none is OK:
;		if targ ne 'HD209458' and targ ne 'BPM16274' 	$
;			and wfilck eq '' then stop
	endelse
	print,'stis_woffset offset in px from template=',toffset
	mulnoise=0.03			; stis_cr default
	if float(date) ge 2009. then mulnoise=0.     ; new data post-SM4 default
; CCD IMAGES W/ SHIFTS BETWEEN THE CR-SPLITS:
;  Watch the image from stis_cr to see if the rejected px track the spectrum.
;###change all CR-splits w/ shifts must be explicitly called out as below.
;  for best noise rejection, use min mult_noise needed to fix photom.
; HD209458 zillion NRPTobs data seems v. sensitive to rejection. set all to 0.1
	if targ eq 'HD209458' or					$
		strpos(fils(ifil),'o3wy020a0') ge 0 or  		$
		strpos(fils(ifil),'o3tt21040') ge 0 or  		$
		strpos(fils(ifil),'o3tt40040') ge 0 or  		$
		strpos(fils(ifil),'o3tt42040') ge 0 or  		$
		strpos(fils(ifil),'o3tt43040') ge 0 or  		$
		strpos(fils(ifil),'o3tt44040') ge 0 or  		$
		strpos(fils(ifil),'o3tt45040') ge 0 or  		$
		strpos(fils(ifil),'o3tt46040') ge 0 or  		$
		strpos(fils(ifil),'o3tt47040') ge 0 or  		$
		strpos(fils(ifil),'o49q02050') ge 0 or  		$
		strpos(fils(ifil),'o49q02070') ge 0 or  		$
		strpos(fils(ifil),'o4sp04010') ge 0 or  		$
		strpos(fils(ifil),'o4sp040d0') ge 0 or  		$
		strpos(fils(ifil),'o49x06010') ge 0 or  		$
		strpos(fils(ifil),'o49x07010') ge 0 or  		$
		strpos(fils(ifil),'o49x16010') ge 0 or  		$
		strpos(fils(ifil),'o49x17010') ge 0 or  		$
		strpos(fils(ifil),'o49x18010') ge 0 or  		$
		strpos(fils(ifil),'o49x19010') ge 0 or  		$
		strpos(fils(ifil),'o49x20010') ge 0 or  		$
		strpos(fils(ifil),'o49x28010') ge 0 or  		$
		strpos(fils(ifil),'o57t02020') ge 0 or  		$
		strpos(fils(ifil),'o61002030') ge 0 or  		$
		strpos(fils(ifil),'o61003030') ge 0 or  		$
		strpos(fils(ifil),'o61004030') ge 0 or  		$
		strpos(fils(ifil),'o3tt41040') ge 0 then begin
		 mulnoise=0.1
		 if strpos(fils(ifil),'o57t02020') ge 0 then mulnoise=0.4
		 if strpos(fils(ifil),'o49x16010') ge 0 then mulnoise=0.06
		 print,'Shifted ccd cr-split. mult_noise set to ',mulnoise
; 2016mar7	if targ eq 'HD209458' then nosky=1	; no sky in subarrays
		endif
	nosky=0				;default=DO sky adjustment
; for HD189733 (g750M) ck 18sco-No sig. diff, so not always required:
	if sxpar(hd,'subarray') then nosky=1	; no sky in subarrays
; after 09jul16 - for cases of change from the default of 0:
	if strpos(fils(ifil),'oa9j020l0') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'oa9j020m0') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'obau01020') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'obbc10020') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'obbc07040') ge 0 then mulnoise=.03	;10may4
	if strpos(fils(ifil),'obau02020') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'obau02040') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'obau04020') ge 0 then mulnoise=.03
	if strpos(fils(ifil),'obc401060') ge 0 then mulnoise=.03 ;G191,G430L
	if strpos(fils(ifil),'obc404010') ge 0 then mulnoise=.03 ;p330e,G430L
	if strpos(fils(ifil),'obc404020') ge 0 then mulnoise=.03 ;p330e,G750L
	if strpos(fils(ifil),'obc405030') ge 0 then mulnoise=.03 ;hd165459 230lb
	if strpos(fils(ifil),'obc405050') ge 0 then mulnoise=.03 ;hd165459 g430l
	if strpos(fils(ifil),'obc405060') ge 0 then mulnoise=.03 ;hd165459 g50l
	if strpos(fils(ifil),'obmzl1020') ge 0 then mulnoise=.03 ;agk g230lb
	if strpos(fils(ifil),'obc402050') ge 0 then mulnoise=.03 ;gd153 g430l
	if strpos(fils(ifil),'obnf06010') ge 0 then mulnoise=.03 ;P330E g430l-E1
	if strpos(fils(ifil),'obnl09010') ge 0 then mulnoise=.10 ;HD38949 g750L
	if strpos(fils(ifil),'obnl09030') ge 0 then mulnoise=.10 ;HD38949 g430L
	if strpos(fils(ifil),'obnl08010') ge 0 then mulnoise=.03 ;HD37962 g430L
	if strpos(fils(ifil),'obnl08020') ge 0 then mulnoise=.01 ;HD37962 g230LB
	if strpos(fils(ifil),'obnl08030') ge 0 then mulnoise=.03 ;HD37962 g750L
	if strpos(fils(ifil),'obnl10010') ge 0 then mulnoise=.03 ;HD106252 g430L
	if strpos(fils(ifil),'obnl10020') ge 0 then mulnoise=.03 ;HD106252g230LB
	if strpos(fils(ifil),'obnl10030') ge 0 then mulnoise=.03 ;HD106252 g750L
	if strpos(fils(ifil),'obnl11010') ge 0 then mulnoise=.03 ;HD205905 g430L
	if strpos(fils(ifil),'obnl11020') ge 0 then mulnoise=.02 ;HD205905g230LB
	if strpos(fils(ifil),'obnl11030') ge 0 then mulnoise=.03 ;HD205905 g750L
	if strpos(fils(ifil),'obnl03020') ge 0 then mulnoise=.03 ;HD116405G230LB
	if strpos(fils(ifil),'obnl03040') ge 0 then mulnoise=.01 ;HD116405 G750L
	if strpos(fils(ifil),'obnl06020') ge 0 then mulnoise=.03 ;HD180609G230LB
	if strpos(fils(ifil),'obnl07030') ge 0 then mulnoise=.03 ;BD+60D1753 G430L
	if strpos(fils(ifil),'obnl12010') ge 0 then mulnoise=.05 ;HD27836 G430L
	if strpos(fils(ifil),'obnl12030') ge 0 then mulnoise=.06 ;HD27836 G430L
	if strpos(fils(ifil),'obvp06060') ge 0 then mulnoise=.03 ;gd71 G750L
	if strpos(fils(ifil),'obto03020') ge 0 then mulnoise=.01;HD158485 g230LB
	if strpos(fils(ifil),'obto09010') ge 0 then mulnoise=.03 ;hd60753 g430L
	if strpos(fils(ifil),'obto09020') ge 0 then mulnoise=.01 ;hd60753 g430LE
	if strpos(fils(ifil),'obto09030') ge 0 then mulnoise=.03;hd60753 g230LBE
	if strpos(fils(ifil),'obto09040') ge 0 then mulnoise=.01 ;hd60753 g230LB
	if strpos(fils(ifil),'obto09050') ge 0 then mulnoise=.03 ;hd60753 g750L
;sat	if strpos(fils(ifil),'obto05010') ge 0 then mulnoise=.03 ;lamlep g430L
;sat	if strpos(fils(ifil),'obto05020') ge 0 then mulnoise=.1 ;lamlep g430LE1
;sat	if strpos(fils(ifil),'obto05030') ge 0 then mulnoise=.03 ;lamlep g430LE1
	if strpos(fils(ifil),'obto05050') ge 0 then mulnoise=.03 ;lamlep g750L
	if strpos(fils(ifil),'obto07050') ge 0 then mulnoise=.03 ;mu col g750L
; few extra pts rejected & no mulnoise will fix - OK:
;SAT!	if strpos(fils(ifil),'obto03040') ge 0 then mulnoise=.03 ;HD158485 g430L
;SAT!	if strpos(fils(ifil),'obto03050') ge 0 then mulnoise=.01 ;HD158485 g750L
;SAT!	if strpos(fils(ifil),'obto01010') ge 0 then mulnoise=.2 ;159222 g430L
	if strpos(fils(ifil),'obto01020') ge 0 then mulnoise=.03 ;159222 g230LB
	if strpos(fils(ifil),'obto01030') ge 0 then mulnoise=.1 ;159222 g750L
	if strpos(fils(ifil),'obto08010') ge 0 then mulnoise=.1 ;ksi2Cet g430L
;SAT!!	if strpos(fils(ifil),'obto08020') ge 0 then mulnoise=.2 ;ksi2Cet g430LE1
	if strpos(fils(ifil),'obto08030') ge 0 then mulnoise=.03 ;ksi2Cet g430L
	if strpos(fils(ifil),'obto08040') ge 0 then mulnoise=.03 ;ksi2Cet g430L
	if strpos(fils(ifil),'obto08050') ge 0 then mulnoise=.03 ;ksi2Cet g430L
	if strpos(fils(ifil),'obto04020') ge 0 then mulnoise=.03 ;Hd163466 230lb
	if strpos(fils(ifil),'obto04030') ge 0 then mulnoise=.03 ;Hd163466 430LE
	if strpos(fils(ifil),'obto04050') ge 0 then mulnoise=.03 ;Hd163466 G750L
	if strpos(fils(ifil),'obto02040') ge 0 then mulnoise=.03 ;Hd14943 G430L
	if strpos(fils(ifil),'obto02050') ge 0 then mulnoise=.03 ;Hd14943 G750L
	if strpos(fils(ifil),'obto06010') ge 0 then mulnoise=.03 ;10Lac G430L
	if strpos(fils(ifil),'obto06020') ge 0 then mulnoise=.03 ;10Lac G430LE1
	if strpos(fils(ifil),'obto06030') ge 0 then mulnoise=.06 ;10Lac G230LBE1
	if strpos(fils(ifil),'obto06040') ge 0 then mulnoise=.06 ;10Lac G230LB
	if strpos(fils(ifil),'obto06050') ge 0 then mulnoise=.06 ;10Lac G750L
	if strpos(fils(ifil),'oc3i01010') ge 0 then mulnoise=.04 ;hd9051 G230lbE
	if strpos(fils(ifil),'oc3i01020') ge 0 then mulnoise=.03 ;hd9051 G430lE1
	if strpos(fils(ifil),'oc3i01030') ge 0 then mulnoise=.03 ;hd9051 G430l
	if strpos(fils(ifil),'oc3i01040') ge 0 then mulnoise=.03 ;hd9051 G750l
	if strpos(fils(ifil),'oc3i02010') ge 0 then mulnoise=.06 ;hd31128G230lbE
	if strpos(fils(ifil),'oc3i02020') ge 0 then mulnoise=.03 ;...    G430LE1
	if strpos(fils(ifil),'oc3i02040') ge 0 then mulnoise=.03 ;...    G750L
	if strpos(fils(ifil),'oc3i03010') ge 0 then mulnoise=.02 ;hd74000G230lbE
	if strpos(fils(ifil),'oc3i04040') ge 0 then mulnoise=.03 ;hd111980 G750L
	if strpos(fils(ifil),'oc3i05010') ge 0 then mulnoise=.03 ;h160617G230lbE1
	if strpos(fils(ifil),'oc3i05020') ge 0 then mulnoise=.03 ;h160617G430lE1
	if strpos(fils(ifil),'oc3i05040') ge 0 then mulnoise=.03 ;h160617G750L
	if strpos(fils(ifil),'oc3i06020') ge 0 then mulnoise=.01 ;h200654G430lE1
	if strpos(fils(ifil),'oc3i06030') ge 0 then mulnoise=.03 ;h200654G430l
	if strpos(fils(ifil),'oc3i07020') ge 0 then mulnoise=.02 ;h185975G230lbE1
	if strpos(fils(ifil),'oc3i07040') ge 0 then mulnoise=.01 ;h185975G750lE1
	if strpos(fils(ifil),'oc3i08010') ge 0 then mulnoise=.02 ;b21d0607G230lbE1
	if strpos(fils(ifil),'oc3i09010') ge 0 then mulnoise=.03 ;bd54 G230lb E1
	if strpos(fils(ifil),'oc3i09030') ge 0 then mulnoise=.03 ;bd54 G750l
	if strpos(fils(ifil),'oc3i10020') ge 0 then mulnoise=.01 ;bd29 G230lb E1
	if strpos(fils(ifil),'oc3i12020') ge 0 then mulnoise=.03 ;bd02 G430l E1
	if strpos(fils(ifil),'oc3i12030') ge 0 then mulnoise=.02 ;bd02 G750L
	if strpos(fils(ifil),'oc3i13010') ge 0 then mulnoise=.01 ;gj75 G230lb E1
	if strpos(fils(ifil),'oc3i14040') ge 0 then mulnoise=.03 ;G191 G750L
	if strpos(fils(ifil),'oc3i11010') ge 0 then mulnoise=.03 ;BD26 G430L E1
	if strpos(fils(ifil),'oc3i11030') ge 0 then mulnoise=.03 ;BD26 G750L
	if strpos(fils(ifil),'ocga04060') ge 0 then mulnoise=.05 ;GD71 G750L
	if strpos(fils(ifil),'oceil4060') ge 0 then mulnoise=.15 ;AGK G750L
	if strpos(fils(ifil),'ocmv01010') ge 0 then mulnoise=.03 ;HD29647 G430L
	if strpos(fils(ifil),'ocmv26010') ge 0 then mulnoise=.03 ;HD281159 G430L
	if strpos(fils(ifil),'ocmv32010') ge 0 then mulnoise=.03 ;HD46660 G430L
	if strpos(fils(ifil),'ocmv0q010') ge 0 then mulnoise=.35 ;HD91983 G430L
	if strpos(fils(ifil),'ocmv05020') ge 0 then mulnoise=.2;HD73882 G750Ldbl
	if strpos(fils(ifil),'ocmv0w010') ge 0 then mulnoise=.03 ;HD142096 G430L
	if strpos(fils(ifil),'ocmv0w020') ge 0 then mulnoise=.02 ;HD142096 G750L
	if strpos(fils(ifil),'ocmv59010') ge 0 then mulnoise=.03 ;HD198781 G430L
	if strpos(fils(ifil),'ocmv66010') ge 0 then mulnoise=.03 ;HD199216 G430L
	if strpos(fils(ifil),'ocmv92020') ge 0 then mulnoise=.02 ;HD239745 G750L
	if strpos(fils(ifil),'ocmv13010') ge 0 then mulnoise=.02 ;HD210121 G430L
	if strpos(fils(ifil),'ocmv0x010') ge 0 then mulnoise=.01 ;HD142165 G430L
	if strpos(fils(ifil),'ocmv0x020') ge 0 then mulnoise=.01 ;HD142165 G750L
	if strpos(fils(ifil),'ocmv53020') ge 0 then mulnoise=.01 ;HD142165 G750L
	if strpos(fils(ifil),'ocmv85020') ge 0 then mulnoise=.03 ;HD217086 G750L
	if strpos(fils(ifil),'ocmv69010') ge 0 then mulnoise=.03 ;HD18352 G430L
	if strpos(fils(ifil),'ocmv22010') ge 0 then mulnoise=.01 ;HD197702 G430L
	if strpos(fils(ifil),'ocmv0j010') ge 0 then mulnoise=.04 ;HD146285 G430L
	if strpos(fils(ifil),'ocmv0j020') ge 0 then mulnoise=.01 ;HD146285 G750L
	if strpos(fils(ifil),'ob6h01010') ge 0 then mulnoise=.1 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h01020') ge 0 then mulnoise=.1 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h01030') ge 0 then mulnoise=.1 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h01040') ge 0 then mulnoise=.1 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h03010') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h03020') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h03030') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h03040') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h52010') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h52020') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h52030') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h52040') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h54010') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h54020') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h54030') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h54040') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h55010') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h55020') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h55030') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'ob6h55040') ge 0 then mulnoise=.2 ;HD189733 G750M
	if strpos(fils(ifil),'odd707010') ge 0 then mulnoise=.05 ;KF08T3 G430LE1
	if strpos(fils(ifil),'odck01060') ge 0 then mulnoise=.1  ;GD71 G750L
	if strpos(fils(ifil),'odd708020') ge 0 then mulnoise=.05 ;KF08T3 G750L
	if strpos(fils(ifil),'odck02060') ge 0 then mulnoise=.05 ;GD153 G4750L
	if strpos(fils(ifil),'odd709020') ge 0 then mulnoise=.03 ;KF08T3 G750L
	if strpos(fils(ifil),'odbvl3010') ge 0 then mulnoise=.03 ;AGK G230LBE1
	if strpos(fils(ifil),'odohl2020') ge 0 then mulnoise=.03 ;AGK G230L
	if strpos(fils(ifil),'odohm10b0') ge 0 then mulnoise=.03 ;AGK G430L
	if strpos(fils(ifil),'odohm10a0') ge 0 then mulnoise=.1 ;AGK G430L -10.4
	if strpos(fils(ifil),'odohm10d0') ge 0 then mulnoise=.03 ;AGK G430L 16.6
	if strpos(fils(ifil),'odohl3010') ge 0 then mulnoise=.03 ;AGK G230LBE1
	if strpos(fils(ifil),'odta05020') ge 0 then mulnoise=.01 ;HD2811g230lbE1
	if strpos(fils(ifil),'odta05030') ge 0 then mulnoise=.02 ;HD2811 g430lE1
	if strpos(fils(ifil),'odta05040') ge 0 then mulnoise=.07  ;HD2811 g750l
	if strpos(fils(ifil),'odta08010') ge 0 then mulnoise=.02 ;16cygBgyroBad?
	if strpos(fils(ifil),'odta08020') ge 0 then mulnoise=.04 ;16cygBgyroBad?
	if strpos(fils(ifil),'odta08030') ge 0 then mulnoise=.1 ;16cygB gyroBad?
	if strpos(fils(ifil),'odta08040') ge 0 then mulnoise=.3 ;16cygB gyroBad?
	if strpos(fils(ifil),'odta09020') ge 0 then mulnoise=.04 ;HD142331 g430l
	if strpos(fils(ifil),'odta09030') ge 0 then mulnoise=.04 ;HD142331 g750l
	if strpos(fils(ifil),'odta13010') ge 0 then mulnoise=.03 ;feige110g230lb
	if strpos(fils(ifil),'odta13020') ge 0 then mulnoise=.05 ;feige110 g430l
	if strpos(fils(ifil),'odta03010') ge 0 then mulnoise=.06 ;HD128998g230lb
	if strpos(fils(ifil),'odta03020') ge 0 then mulnoise=.04 ;HD128998g2lbE1
	if strpos(fils(ifil),'odta03030') ge 0 then mulnoise=.03 ;HD128998g4lbE1
	if strpos(fils(ifil),'odta03040') ge 0 then mulnoise=.04 ;HD128998 g430l
	if strpos(fils(ifil),'odta03050') ge 0 then mulnoise=.27 ;HD128998 g750l
	if strpos(fils(ifil),'odta03060') ge 0 then mulnoise=.03 ;HD128998 g750l
	if strpos(fils(ifil),'odta51010') ge 0 then mulnoise=.09 ;delumi g230lb
	if strpos(fils(ifil),'odta51020') ge 0 then mulnoise=.09;delumi g230lbE1
	if strpos(fils(ifil),'odta51030') ge 0 then mulnoise=.10 ;delumi g430lE1
	if strpos(fils(ifil),'odta51040') ge 0 then mulnoise=.09 ;delumi g430l
	if strpos(fils(ifil),'odta51050') ge 0 then mulnoise=.25 ;delumi g750l
	if strpos(fils(ifil),'odta51060') ge 0 then mulnoise=.30 ;delumi g750l
	if strpos(fils(ifil),'odud03010') ge 0 then mulnoise=.01 ; g191 g230lb
	if strpos(fils(ifil),'odud03030') ge 0 then mulnoise=.10 ; g191 g230lb
	if strpos(fils(ifil),'odud03040') ge 0 then mulnoise=.03 ; g191 g230lb
	if strpos(fils(ifil),'odta12010') ge 0 then mulnoise=.12 ; etauma g230lb
; sat	if strpos(fils(ifil),'odta12030') ge 0 then mulnoise=.05 ; etauma g430l
	if strpos(fils(ifil),'odta12050') ge 0 then mulnoise=.06 ; etauma g750l
	if strpos(fils(ifil),'odvkl1030') ge 0 then mulnoise=.02 ; agk g430l
	if strpos(fils(ifil),'odvkl1040') ge 0 then mulnoise=.05 ; agk g430l E1
	if strpos(fils(ifil),'odta14010') ge 0 then mulnoise=.08;feige34 230lbE1
	if strpos(fils(ifil),'odta14020') ge 0 then mulnoise=.02 ;feige34 430lE1
	if strpos(fils(ifil),'odta18010') ge 0 then mulnoise=.01 ;HZ44 430lE1
	if strpos(fils(ifil),'odta18020') ge 0 then mulnoise=.02 ;HZ44 G230LB E1
	if strpos(fils(ifil),'odta19010') ge 0 then mulnoise=.02 ;109vir G230LB
	if strpos(fils(ifil),'odta19020') ge 0 then mulnoise=.06;109vir G230LBE1
	if strpos(fils(ifil),'odta19030') ge 0 then mulnoise=.04;109vir G430LBE1
	if strpos(fils(ifil),'odta19040') ge 0 then mulnoise=.05 ;109vir G430LB
	if strpos(fils(ifil),'odta19050') ge 0 then mulnoise=.11; 109vir G750L
	if strpos(fils(ifil),'odta19060') ge 0 then mulnoise=.09 ;109vir G750L
	if strpos(fils(ifil),'odta11010') ge 0 then mulnoise=.04;HD115169-230lbE
	if strpos(fils(ifil),'odta11020') ge 0 then mulnoise=.01;HD115169-430l E
	if strpos(fils(ifil),'odta11030') ge 0 then mulnoise=.02 ;HD115169-750l
	if strpos(fils(ifil),'odta07010') ge 0 then mulnoise=.04  ;18sco-230lbE1
	if strpos(fils(ifil),'odta07020') ge 0 then mulnoise=.02  ;18sco-G430E1
	if strpos(fils(ifil),'odta07030') ge 0 then mulnoise=.12  ;18sco-G750L
	if strpos(fils(ifil),'odta07040') ge 0 then mulnoise=.07  ;18sco-G750L
	if strpos(fils(ifil),'odta17020') ge 0 then mulnoise=.04  ;HZ4-G430E1
	if strpos(fils(ifil),'odta15010') ge 0 then mulnoise=.04  ;HD93521-230lb
	if strpos(fils(ifil),'odta15020') ge 0 then mulnoise=.04  ;93521-230lbE1
	if strpos(fils(ifil),'odta15030') ge 0 then mulnoise=.06  ;93521-G430lE1
	if strpos(fils(ifil),'odta15040') ge 0 then mulnoise=.15  ;93521-G430L
	if strpos(fils(ifil),'odta15050') ge 0 then mulnoise=.01  ;93521-G750L
	if strpos(fils(ifil),'odta06030') ge 0 then mulnoise=.01  ;HD55677-G430l
	if strpos(fils(ifil),'odta06040') ge 0 then mulnoise=.02  ;HD55677-G750L
	if strpos(fils(ifil),'odud01060') ge 0 then mulnoise=.05  ;GD71-G750L
	if strpos(fils(ifil),'odta72010') ge 0 then mulnoise=.04  ;eta1dor-230lb
	if strpos(fils(ifil),'odta72020') ge 0 then mulnoise=.06  ;eta1dor230lbE
	if strpos(fils(ifil),'odta72030') ge 0 then mulnoise=.26  ;eta1dor430lE1
	if strpos(fils(ifil),'odta72040') ge 0 then mulnoise=.04  ;eta1dor g430l
	if strpos(fils(ifil),'odta72050') ge 0 then mulnoise=.07  ;eta1dor g750l
	if strpos(fils(ifil),'odta72060') ge 0 then mulnoise=.08  ;eta1dor g750l
	if strpos(fils(ifil),'odvkl2040') ge 0 then mulnoise=.15  ;AGK g430l
	if strpos(fils(ifil),'odta04020') ge 0 then mulnoise=.03  ;hd101452-lbE1
; SAT.  if strpos(fils(ifil),'odta04030') ge 0 then mulnoise=.5  ;hd101452-43E1
	if strpos(fils(ifil),'odta04040') ge 0 then mulnoise=.02  ;hd101452 G750
	if strpos(fils(ifil),'odta60010') ge 0 then mulnoise=.04  ;hd167060-lbE1
	if strpos(fils(ifil),'odta60020') ge 0 then mulnoise=.01  ;hd167060-43E1
	if strpos(fils(ifil),'odta60030') ge 0 then mulnoise=.12 ;hd167060-G750L
	if strpos(fils(ifil),'odud02060') ge 0 then mulnoise=.05 ;G153-G750L
	if strpos(fils(ifil),'odta16020') ge 0 then mulnoise=.02  ;hz21-g430-E1
	if strpos(fils(ifil),'oe3653010') ge 0 then mulnoise=.08  ;agk-g230lb-E1
	if strpos(fils(ifil),'oe3653020') ge 0 then mulnoise=.09  ;agk-g230lb
	if strpos(fils(ifil),'oe3653030') ge 0 then mulnoise=.12  ;agk-g430l
	if strpos(fils(ifil),'oe3653040') ge 0 then mulnoise=.10  ;agk-g430l-E1
	if strpos(fils(ifil),'oe3653060') ge 0 then mulnoise=.33  ;agk-g750l
	if strpos(fils(ifil),'ocy525f8q') ge 0 then mulnoise=.10  ;HD93250 g750l
	if strpos(fils(ifil),'ocy527i8q') ge 0 then mulnoise=.25;HD164794 g230lb
	if strpos(fils(ifil),'ocy527ieq') ge 0 then mulnoise=.05 ;HD164794 g430l
	if strpos(fils(ifil),'ocy527iiq') ge 0 then mulnoise=.03 ;HD164794 g750l
	if strpos(fils(ifil),'ocy531thq') ge 0 then mulnoise=.3 ;HD217086 g750l
	if strpos(fils(ifil),'ocy533hvq') ge 0 then mulnoise=.05 ;NGC6611- g750l
	if strpos(fils(ifil),'ocy533hxq') ge 0 then mulnoise=.15 ;NGC6611- g750l
	if strpos(fils(ifil),'ocy533hzq') ge 0 then mulnoise=.05 ;NGC6611- g750l
	if strpos(fils(ifil),'ocy538alq') ge 0 then mulnoise=.15 ;HD207198 g430l
	if strpos(fils(ifil),'ocy539qcq') ge 0 then mulnoise=.1  ;HD210839 g430l
	if strpos(fils(ifil),'odvkm1090') ge 0 then mulnoise=.1  ;AGK g430l
	if strpos(fils(ifil),'oehj01060') ge 0 then mulnoise=.05  ;AGK g750l

	if mulnoise ne 0 and float(date) gt 2009 then print,fils(ifil)+' CR-split. mult_noise set to ',mulnoise
	dist1=string(100+abs(postarg)/.05,'(i3)')	;vignetting data-00dec26
; 2012Feb7 - OOPS, I saturated the HD158485 G430L & 750L too much for cr rej.
; 2012apr25 - and also lam lep, mu Col G230Lb & 430L. 750L!!
; 2012May21 - and also HD159222 430L.
; 2012Aug14 - and still ksi2 Ceti, but only G430L-E1
; 2019apr26 - once again. G430L for HD101452
; ****BUT**** gwidth=7 is enough for these slight saturations.

	if gain gt 0 then begin
	    ; CCD reduction path
	
        ; no DQ flagging. even 7pxhgt has no problems. Data w/ saturation:
		hotpix='0,ccdpar=ccd_nosat.fits,hotthresh=10,warmthresh=10,hpxtab=none,bpxtab=none,'
		nobs=max([sxpar(hd,'crsplit'),sxpar(hd,'nrptexp')])

        ; SATURATED OBS:
		if  (strpos(targ,'HD172167') ge 0 and flg7px eq 0) or	$
		    (targ eq 'SIRIUS' and flg7px eq 0) or 		$
		    (targ eq 'HD158485' and (mode eq 'G430L' or mode eq 'G750L')) or	$
		    (targ eq 'LAMLEP' and (mode eq 'G430L' or mode eq 'G230LB')) or	$
		    (targ eq 'HD159222' and mode eq 'G430L') or 	$
		    (targ eq 'KSI2CETI' and root eq 'obto08020') or $
		    (strpos(targ,'189733') gt 0 and mode eq 'G430L') or $
            ; MUCOL 3 modes satur:
		    (targ eq 'MUCOL') or $
		    root eq 'odta12020' or	$
            ; odta12030 is eta UMa & flgged as BAD. odta04030 is HD101452
		    root eq 'odta12030' or $
		    root eq 'odta04030' or	$
		    root eq 'odta12040' or $
		    root eq 'odta12060' then begin		
		        ; odta1212040 is etaUMa,which is too sat for gwid=7, but going to gwid=25
                ;	gets only 0.26% more avg flux. See /plots/modcf-etauma-orig.ps
                ;	But .1s/.3s shows only a bump peaking @ ~1% high at 4000A, ie the 
                ;	dominant long odta1212040 is ~1% low at 4000A. I see no obvious fix.
	    		if root eq 'odta12030' or root eq 'odta12040' then gwidth='11' ; Sat.
                ; Saturated data - crej is summed data-overscan:
                ; Vega all images satur... even 0.9s g230lb at wl>~2900A. 
                ; flg7px=1 only for 7px high extraction of the 2 G750L obs.
                ;   Comment mrgall to avoid making *.7pix files.  See sirius/doc.procedure
                ; pcttab=calstis gets same (uncorr) NET w/ or w/o pcttab default or none; but
                ; calstis Flux is closer w/ default pcttab, so run calstis w/ default.
			    crerr=0			; see vega-test.pro
			    eps=0			; 2016feb9 for subarray.pro
			    crej=0.
			    print,'***'+root+' '+targ+' sat. NO CR-rej gwidth='+gwidth
                ; 09Aug5 & see corresponding temp average changes in stis_cr
			    avtemp=0.
			    for i=1,nobs do begin	; readout 0 is just a header
                    ; see overscan.doc: (Yes, this reads the counts. DQ & err arrays are a mystery.)
				    stis_read,fils(ifil),hdcr,infil,readout=i
				    avtemp=avtemp+sxpar(hdcr,'occdhtav')
			 	    crej=crej+infil		; sum of all crsplit img
			 	endfor
			    sxaddpar,hdcr,'occdhtav',avtemp/nobs,'Average computed by preproc.pro'
			    sxaddpar,hdcr,'integ',sxpar(hdcr,'integ')*nobs
			    sxaddpar,hdcr,'exptime',sxpar(hdcr,'exptime')*nobs
			    sxaddpar,hdcr,'ncombine',nobs 	; subtract Bias N times
                ; no stis_cr to  pop. senstab for gwidth=7.( =11 fixed below)
			    sxaddpar,hdcr,'senstab','sens_'+strlowcase(mode)+'.fits'	; 2014Oct7
			    sxaddhist,'All '+targ+' have saturation. Use Standard overscan & NO cr-rej.',hdcr
                ; APPROX: no dark subtr=negligible error.
                ;	BIAS is subtracted by calstis.

		end else begin			; END Satur. processing
		    
                ; case of no Sat. data. Hot pix files made w/ .05 threshold, so .04 is equiv to
                ;	.05, Minus sign is for interpol along spectrum. See calstis.doc &
                ;	calstis_bpx.pro
			    hotpix='-0.04,'
                ; stis_cr must populate senstab keyword & subtract darks, but meandark not pop.
			    print,'** Do stis_cr **'
			    crlist=fils(ifil)
                ; Fix screwy 14141 G. Worthy obs:
			    if nobs eq 1 and strpos(fils(ifil),'ocy5') ge 0 then begin ;crsplit=1
				    crlist=[fils(ifil),fils(ifil+1)]
				    stis_read,fils(ifil),hdcr
				    fdecomp,fils(ifil),dum,dum,f1
				    fdecomp,fils(ifil+1),dum,dum,f2
				    f1=gettok(f1,'_')
				    f2=gettok(f2,'_')
				    crskip=fils(ifil+1)
				    print,'stis_cr for ',crlist
			    endif
                ; stis_cr Does NOT rm any HOT px! Cannot cf. to _flt, as  pipeline must rm hot
                ;	px by another method in _flt,_crj files!
			    stis_cr,crlist,hdcr,crej,crerr,eps,nused,noskyadjust=nosky,mult_noise=mulnoise

                ; hdcr is OUTPUT of stis_cr
			    if nobs eq 1 and strpos(fils(ifil),'ocy5') ge 0 then begin
				    sxaddhist,'Do STIS_CR for '+f1+' '+f2,hdcr
				    sxaddpar,hdcr,'crsplit',2
				endif
			    sxaddhist,'MULT_NOISE='+string(mulnoise,'(f5.2)')+' IDL STIS_CR Noise mult. factor',hdcr
		endelse				; End non-sat processing

        ; special cases of big CR hit on spectrum, requiring "hand patch"
		if strpos(fils(ifil),'o3tt21040') ge 0 then crej(889,494:496) = avg([crej(888,494:496),crej(890,494:496)])
		if strpos(fils(ifil),'o4a502020') ge 0 then crej(731,510)=1142
		if strpos(fils(ifil),'o4a502030') ge 0 then crej(731,510)=4660
    	if strpos(fils(ifil),'o4a502040') ge 0 then crej(731,510)=9130
	    if strpos(fils(ifil),'o49x07010') ge 0 then crej(363,510)=28000
		if strpos(fils(ifil),'obbc09010') ge 0 then begin ; KF06T2-G430L
		    crej(196,892)=160.
    		crej(196,893)=440.
	    	crej(196,894)=260.
			crej( 91,895)= 50.
		endif
    	if strpos(fils(ifil),'o4a520040') ge 0 then begin
	    	crej(329,510)=7600
	    	crej(330,511)=7000
			crej(968,509)=9500
			crej(968,510)=11300
		endif
    	if strpos(fils(ifil),'o61002030') ge 0 then begin
	    	crej(352:353,510)=11200
	    	crej(350:353,511)=9700
			crej(603,511)=2800
		endif
        ; 00nov21 - P330E, G750L:
		if strpos(fils(ifil),'o49x28010') ge 0 then begin
		    crej(658,510)=6700
		    crej(658,511)=3600
		    crej(703,510)=6200
		    crej(703,511)=3600
		endif
    	if strpos(fils(ifil),'o4d101020') ge 0 then begin
	    	crej(287:289,509)=4400
	    	crej(628,509)=5250
	    endif
		if strpos(fils(ifil),'o53002030') ge 0 then begin  ;No CR split!
		    crej(822,509)=1650
		    crej(881,511)=238-150
		endif
		if strpos(fils(ifil),'o8h105030') ge 0 then crej(362,895)=0
		if strpos(fils(ifil),'o8h105030') ge 0 then	crej(362:364,894)=860
    	if strpos(fils(ifil),'o8vj11020') ge 0 then crej(568,512)=100
	    if strpos(fils(ifil),'obc402030') ge 0 then crej(39:40,514)=30
		if strpos(fils(ifil),'obc406040') ge 0 then crej(603,511)=4700
		if strpos(fils(ifil),'obc406040') ge 0 then crej(603,516:518)=110
		if strpos(fils(ifil),'obc406040') ge 0 then crej(348,502:507)=380
    	if strpos(fils(ifil),'obc407040') ge 0 then crej(603,511)=2000
	    if strpos(fils(ifil),'obc408040') ge 0 then crej(603,511)=2300
		if strpos(fils(ifil),'obc410040') ge 0 then crej(603,511)=6200
    	if strpos(fils(ifil),'obc411040') ge 0 then crej(603,511)=5200
	    if strpos(fils(ifil),'obmzl3010') ge 0 then crej(168,895)=1000
		if strpos(fils(ifil),'obbc09020') ge 0 then crej(39:40,514)=420
    	if strpos(fils(ifil),'obto03050') ge 0 then crej(939,512:513)=3550
    	if strpos(fils(ifil),'oc3i13010') ge 0 then crej(97,892)=50
	    if strpos(fils(ifil),'oc3i13010') ge 0 then crej(363,890)=118
		if strpos(fils(ifil),'obau0306') ge 0 then crej(763,508)=900
	    if strpos(fils(ifil),'obto10030') ge 0 then crej(325:327,508)=417
    	if strpos(fils(ifil),'obto10030') ge 0 then crej(326,509)=1300
	    if strpos(fils(ifil),'obmzl3040') ge 0 then crej(168,895)=2100
		if strpos(fils(ifil),'oceil4060') ge 0 then begin
		    crej(206:207,510)=23840
		    crej(207,511)=21550
		    crej(429,510)=13000
		    crej(428:429,511)=10840
		    crej(433,510:511)=11200
		endif
		if strpos(fils(ifil),'obc407040') ge 0 then crej(348,504:505)=400
    	if strpos(fils(ifil),'obc408040') ge 0 then crej(348,504:505)=170
		if strpos(fils(ifil),'obc410040') ge 0 then crej(348,504:505)=400
		if strpos(fils(ifil),'obc411040') ge 0 then crej(348,504:505)=600
		if strpos(fils(ifil),'obnl09030') ge 0 then crej(348,505)=2750
		if strpos(fils(ifil),'o4pz01030') ge 0 then crej(394,504)=15
		if strpos(fils(ifil),'o4pz01030') ge 0 then crej(394:395,505:506)=25
		if strpos(fils(ifil),'obto07050') ge 0 then crej(512,505)=1100
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(260:261,506)=65
		if strpos(fils(ifil),'o8vj13010') ge 0 then crej(693,508:509)=50
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(157:159,510)=2200
        ; ff was 599, but ~5 rows that are ~100 too high:
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(261,512)=0
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(636:637,508)=180
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(636:637,509)=300
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(660:661,514)=0
		if strpos(fils(ifil),'o8vj12020') ge 0 then crej(745:746,512)=150
		if strpos(fils(ifil),'oddg02010') ge 0 then begin ; grw-G430L E1
    		crej(168,895)=1500.
	    	crej(296,890)=290.
			crej(887,893)=4170.
		    crej(944,895)= 470.
		endif
		if strpos(fils(ifil),'oddg02020') ge 0 then begin
		    crej(594,506)=465.			; grw-G430L
    		crej(603,511)=1986.
		endif
	    if strpos(fils(ifil),'oddg02030') ge 0 then begin ; grw-G750L
			crej(200,505)=200.
		    crej(348,505)=140.
    		crej(603,511)=1275.
		endif
	    if strpos(fils(ifil),'odd707020') ge 0 	then crej(594,506)=2080.
		if strpos(fils(ifil),'odd708010') ge 0 	then crej(152,892)=350.
    	if strpos(fils(ifil),'odck02050') ge 0 then begin ; GD153 G430L
	    	crej(538,510)=450.
	    	crej(538:539,511)=150.
	    endif
    	if strpos(fils(ifil),'o8h104020') ge 0 	then crej(200,505)=60.
        ; 2017jul17-intermit. HOT px shows mostly in long GRW exp:
    	if strpos(fils(ifil),'obc402060') ge 0 then crej(40,514)=275
	    if strpos(fils(ifil),'odck02060') ge 0 then crej(40,514)=355
		if strpos(fils(ifil),'oddg03030') ge 0 then crej(40,514)=232
    	if strpos(fils(ifil),'obc403020') ge 0 then crej(40,514)=163
	    if strpos(fils(ifil),'obnk01040') ge 0 then crej(40,514)=55
		if strpos(fils(ifil),'obto05040') ge 0 then crej(20,509)=4060
    	if strpos(fils(ifil),'oddg04030') ge 0 then crej(40,514)=200
	    if strpos(fils(ifil),'oddg04030') ge 0 then crej(40,514)=200
		if strpos(fils(ifil),'odqg01030') ge 0 then crej(40,514)=25
    	if strpos(fils(ifil),'odqg01030') ge 0 then crej(40,513)=40
	    if strpos(fils(ifil),'odqg01030') ge 0 then crej(49,506:520)=crej(49,506:520)-20
    	if strpos(fils(ifil),'odqg02030') ge 0 then crej(40,514)=30
	    if strpos(fils(ifil),'odqg03020') ge 0 then crej(168,895)=100
		if strpos(fils(ifil),'odud03010') ge 0 then crej(508,506)=700
	    if strpos(fils(ifil),'obc406030') ge 0 then crej(881:882,898)=240
		if strpos(fils(ifil),'obc408030') ge 0 then crej(881:882,898)=90
	    if strpos(fils(ifil),'obc410030') ge 0 then crej(881:882,898)=280
	    if strpos(fils(ifil),'o49x16010') ge 0 then crej(482:483,510)=870
	    if strpos(fils(ifil),'o49x16010') ge 0 then crej(581,511)=600
	    if strpos(fils(ifil),'odud01030') ge 0 then crej(168,895)=630
	    if strpos(fils(ifil),'odud01040') ge 0 then crej(168,895)=335
	    if strpos(fils(ifil),'odta60010') ge 0 then crej(305:306,888:889)=16
	    if strpos(fils(ifil),'odud03040') ge 0 then crej(508,506)=315
	    if strpos(fils(ifil),'odqg03030') ge 0 then crej(129,530)=73
	    if strpos(fils(ifil),'odqg03030') ge 0 then crej(130,531)=110
	    if strpos(fils(ifil),'odqg03030') ge 0 then crej(129,533)=26
	    if strpos(fils(ifil),'odqg03030') ge 0 then crej(130,536)=16
	    if strpos(fils(ifil),'oe7q02020') ge 0 then crej(841,506)=195
	    if strpos(fils(ifil),'oe7q02020') ge 0 then crej(841,507)=330
	    if strpos(fils(ifil),'oe7q02020') ge 0 then crej(841,508)=450
	    if strpos(fils(ifil),'oe7q02020') ge 0 then crej(841,509)=1020
	    if strpos(fils(ifil),'oe7q06020') ge 0 then crej(841,506)=120
	    if strpos(fils(ifil),'oe7q06020') ge 0 then crej(841,507)=220
	    if strpos(fils(ifil),'oe7q06020') ge 0 then crej(841,508)=260
	    if strpos(fils(ifil),'oe7q06020') ge 0 then crej(841,509)=660
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,504)=0
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,505)=25
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,506)=25
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,507)=50
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,508)=95
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,509)=130
	    if strpos(fils(ifil),'oe7q08010') ge 0 then crej(841,510)=170
	    if strpos(fils(ifil),'oe7q01010') ge 0 then crej(145,893)=2710
	    if strpos(fils(ifil),'oe7q01010') ge 0 then	crej(145:146,894)=1710
	    if strpos(fils(ifil),'oe3f01010') ge 0 then	crej(307:308,767)=20 ; Gaia 405_1056
	    if strpos(fils(ifil),'oe3f03010') ge 0 then crej(99:100,891)=0
	    if strpos(fils(ifil),'oe3f03010') ge 0 then crej(518,887)=0
	    if strpos(fils(ifil),'oe3f03020') ge 0 then crej(99:100,891)=15
	    if strpos(fils(ifil),'oe3f05010') ge 0 then crej(975,887)=0
	    if strpos(fils(ifil),'oe3f06010') ge 0 then crej(99:100,891)=0
	    if strpos(fils(ifil),'oe3f06010') ge 0 then crej(518,887)=0
	    if strpos(fils(ifil),'odud02030') ge 0 then crej(39:40,514)=46
	    if strpos(fils(ifil),'ocy514030') ge 0 then crej(882,56)=2140
	    if strpos(fils(ifil),'ocy514040') ge 0 then crej(882,56)=2140
        ; 2020dec9 - subarr starting at y=834 (895-834=61
	    if strpos(fils(ifil),'ocy541010') ge 0 then crej(168,61)=1010
	    if strpos(fils(ifil),'ocy541020') ge 0 then crej(168,61)=1010
	    if strpos(fils(ifil),'oefrl1060') ge 0 then crej(762,510)=1561
	    if strpos(fils(ifil),'oehj02030') ge 0 then crej(40,513:514)=42
	    if strpos(fils(ifil),'oehj02030') ge 0 then crej(39,504:505)=20
	    if strpos(fils(ifil),'oehj03010') ge 0 then crej(40,513:514)=100
	    if strpos(fils(ifil),'oehj03010') ge 0 then crej(39,504:505)=70

        ; ff is NG, as the dark has a feature here that gets subtr.
        ;;	       if strpos(fils(ifil),'oe3f06010') ge 0 then crej(491,748:763)=50
        ;stop

    	nextend=sxpar(hdcr,'nextend')	; 2017jan30-Clarify Nextend:
	    sxaddpar,hdcr,'nextend',nextend,'Orig value. All IDL spec_* files have nextend=1.'
		sxaddpar,hdcr,'pstrtime',pstrtim,'predicted obs. start time from SPT file'
		sxaddpar,hdcr,otemp,OMCAT,comment ; 20011 raw data has om1cat !
        sxaddpar,hdcr,'rad_vel',svel,'Star Velocity (km/s) NOT applied'
		sxaddpar,hdcr,'date',!stime
		
		if targ eq 'HD209458' and wfilck eq '' then sxaddhist,'No Wavecal for HD209458',hdcr

        ; For CCD subimages, insert into full size image
    	subkg='1'
	    if sxpar(hdcr,'subarray') eq 1 then begin
            ; do not subtract any background... (bkgr still computed and written)
		    subkg='0'				; for calstis	
		    subarray,hdcr,crej,crerr,eps		;04sep-djl fancy
		    sxaddhist,'Sub-array inserted into 1024x1024 by preproc.pro',hdcr
		    sxaddhist,'No background subtraction.',hdcr
		endif

    	e1='0'
        ; E1 normally at row 511+383=894.
		if strpos(aperture,'E1') ge 0 then e1='383'	; rows of offset, ie 511+383=894
        ; G750L GAIA405_1056: oe3f01020 @ center, oe3f01030 Postarg=+6.4":
		if root eq 'oe3f01030' then e1='126'		; G750L
        ; G750L GAIA405_6912: oe3f01020 @ low, oe3f01030 @ center:
		if strpos(aper(ifil),'low') gt 0 then begin	; 2020Feb2
		    dist1='20'
    	    if root eq 'oe3f01010' then e1='256'	;G430L 766-511
	    	if root eq 'oe3f01020' then e1='-126'	; G750L @ -6.4"
			if root eq 'oe3f01030' then e1='0'	; G750L @ center
    		if strpos(root,'oe3f030') eq 0 then e1='-60'  ;G587_3008
    		if strpos(root,'oe3f040') eq 0 then e1='-222' ;G593_9680
	    	if strpos(root,'oe3f050') eq 0 then e1='-495' ;G588_7712
			if root eq 'oe3f06010' then e1='243'  ;G587_8560 G430@E1
		    if root eq 'oe3f06020' then e1='-138' ;G587_8560G750cent
		endif
    	;Gaia Jesus
	    if e1 ne '0' then print,targ,' ',root,' ',mode,' ',aperture,' @ E1 offset=',e1

		sxaddhist,'Wavecal file: '+wfil,hdcr
    	if toffix ne 0 then sxaddhist,'STISWLFIX OFFSET(px)='+string(toffix,'(f5.2)'),hdcr
    	sxaddhist,'Written by preproc.pro '+!stime,hdcr
	endif			; end gain gt 0 CCD loop

	sxaddpar,hdcr,'targname',targ	; 2020mar17 eg HZ43B

    ; ###change:
	ext='.fits'
    ;	ext='.trace'			; for special trace files 06jan27
    ;	ext='.noflat'			; to test non-default flat fields
    ;	pfltfil='none'
	pfltfil='def'

	det=strlowcase(strmid(sxpar(hd,'detector'),0,3))

    ; G750L fringe FLATS:
	shrtflg=0		; NOT a G750L short slit CONTEMPORANEOUS flat
	if targ eq 'CCDFLAT' then goto,skipfring

	if strpos(mode,'G750L') ge 0 and cenwave eq '7751' and ext ne '.noflat' then begin
        ;1-D short slit G750L processing - bohlin standard process:
        ; exceptions: {o49q=hr7615 7642, o3wy=bd75 o403 7095, o408 7100}
		if name ge 'o3tt42040_raw' and strpos(name,'o49q') lt 0	$
		   and strpos(name,'o3wy') lt 0 and strpos(name,'o403') lt 0 $	
		   and strpos(name,'o408') lt 0 and strpos(name,'o4sp') lt 0 $
		   and strpos(name,'oe9l02040') lt 0  $		; Lennon 16079
		   and strpos(name,'ocy5') lt 0  $		; Gaia hot stars
           ; ###change? NO-the small slit does not align w/ spectrum @ E1:
		   and e1 eq '0'		 	$	; 20jun26
		then begin   ;do normal short slit flat
			shrtflg=1		; flag for short slit flat
			if name ge 'o3tt42040_raw' and name le 'o3tt48040_raw' then begin
				print,'SPECIAL AVG CONTEMPORANEOUS FLAT USED'
				sxaddhist,'SPECIAL AVG CONTEMPORANEOUS FLAT USED',hdcr
                ; Read un-norm tung flat from preproc and mrgall
				name='tung-noflat.g750l'	;2014sep14 updat
				if gwidth ge '11' then rdf,'dat/tung-noflat11.g750l',1,d else rdf,'dat/tung-noflat.g750l',1,d	; 7px
				wd=d(*,0)  &  nd=d(*,1)
				goto,skipcal
			endif
			orignam=name
            ; find right short slit flat file:
			name=replace_char(name,'40_','30_')
			if name ge 'o45a01030_raw' then	name=replace_char(name,'30_','40_')
			if strpos(name,'o4a5') ge 0 then name=replace_char(name,'20_','10_')
            ; Later data - GENERAL Purpose:
            ;newagk,gd71
			if strpos(name,'o6') ge 0 and strpos(name,'o69u') lt 0 then name=replace_char(name,'40_','50_')
			
            ; DEFAULT name DEFAULT DEFAULT. GENERAL purpose: ########################
            ; FF is for switch of file order. Note : the GE and files+1 in ff. line:
			if orignam ge 'o6ig01060_raw' and targ ne 'HD209458' then fdecomp,fils(ifil+1),dum,dum,name	;DEFAULT

            ; special cases G750L fringe flat files:
			if orignam eq 'o61001040_raw' then name='o61001060_raw'
 			if orignam eq 'o6ig01060_raw' then name='o6ig01070_raw'
			if orignam eq 'o4a505010_raw' then name=replace_char(name,'10_','20_') ; odd!
			if strpos(orignam,'o49x') ge 0 then name=replace_char(name,'10_','20_') ;IRSTD
			if strpos(orignam,'o57t') ge 0 then name=replace_char(name,'20_','30_') ;IRSTD
			if orignam eq 'o49x03020_raw' then name=replace_char(name,'20_','10_')
			if orignam eq 'o6i903060_raw' then name=replace_char(name,'60_','70_')
			if orignam eq 'o6i904060_raw' then name=replace_char(name,'60_','70_');50s
			if strpos(orignam,'o49x130') ge 0 then name='o49x13050_raw' ;IRSTD
			if strpos(orignam,'o49x140') ge 0 then name='o49x14050_raw' ;IRSTD
			if orignam eq 'o5ja040i0_raw' then name='o5ja040j0_raw'	; fring Flat
            ; HD209458
			if strpos(orignam,'o6d2') ge 0 then name='o6d2a4040_raw'
			if strpos(orignam,'o6d22') ge 0 or strpos(orignam,'o6d2b') ge 0 then name='o6d2b2050_raw'
			if strpos(orignam,'o6n3') ge 0 then name='o6n3a2050_raw'
			if strpos(orignam,'o6n3a4') ge 0 or	strpos(orignam,'o6n304') ge 0 then name='o6n3a4050_raw'
			if strpos(orignam,'o6ig100d0') ge 0 then name='o6ig100k0_raw'
			if strpos(orignam,'o5i') ge 0 then name=replace_char(name,'50_','60_') ;IRSTD
			if orignam eq 'o6i904050_raw' then name='o6i9040a0_raw'	;  0.1X0.09 E1
			if strpos(orignam,'o8h2010') ge 0 then name='o8h2010i0_raw' ;02dec10
			if strpos(orignam,'oa9j01090') ge 0 then name='oa9j010b0_raw' ;09jul
            ; 18Nov - two G750L obs:
			if strpos(orignam,'odta08030') ge 0 then name='odta08050_raw' ; 18nov 16CYGB
			if strpos(orignam,'odta51050') ge 0 then name='odta51070_raw' ; 18nov DELUMI
			if strpos(orignam,'odta12') ge 0 then name='odta12070_raw' ; 18DEC27 eta uma
			if strpos(orignam,'odta19') ge 0 then name='odta19070_raw' ; 19jan11 109vir
			if strpos(orignam,'odta07') ge 0 then name='odta07050_raw' ; 19jun19 18sco
			if strpos(orignam,'odta72') ge 0 then name='odta72070_raw' ; 19jun19 ETA1DOR
			if strpos(orignam,'odta03') ge 0 then name='odta03070_raw' ; 19jun19 HD128998
			if name eq orignam then begin		; idiot check
				print,'IDIOT ck: g750l flat in preproc is data'
				stop
			endif
			flat=disk+dir+name+'.'+qual
			sxaddhist,'G750L fringe flat & gwidth='+name+' '+string(gwidth),hdcr
			print,'G750L flat='+flat

            ; ###change
            ; stis_cr,flat,hp,pflat,/displ
			stis_cr,flat,hp,pflat

			if sxpar(hp,'subarray') eq 1 then begin	; for HD209458
                ; do not subtract any background... (bkgr still computed and written)
				subkg='0'		
				offst=abs(sxpar(hp,'ltv2'))
				siz=size(pflat)  &  siz2=siz(2)
				tmp=fltarr(1024,1024)
				tmp(*,offst:offst+siz2-1)=pflat	;insert data
				pflat=tmp
				sxaddpar,hp,'naxis2',1024
				sxaddpar,hp,'subarray','F'
				sxaddhist,'Sub-array inserted into 1024x1024 by preproc.pro',hp
				sxaddhist,'No background subtraction.',hp
			endif

            ; obc407050=fringe flat for obc407040,etc have hot pixel:
			if strpos(fils(ifil),'obc407040') ge 0 then pflat(348,505)=400
			if strpos(fils(ifil),'obc408040') ge 0 then pflat(348,505)=240
			if strpos(fils(ifil),'obc411040') ge 0 then pflat(348,505)=360
			if strpos(fils(ifil),'obnl09030') ge 0 then pflat(348,505)=675
			if strpos(fils(ifil),'obc406040') ge 0 then pflat(348,505)=240
			if strpos(fils(ifil),'obc406040') ge 0 then pflat(603,511)=6000
			if strpos(fils(ifil),'obc406040') ge 0 then pflat(603,517)=200

            ;short slit G750L fringe flat (nd=1D flat). NO CTE corr to * net, so no flat CTE
            ;	corr. pflat=cr-rejected RAW Tungsten 1024x1024 image , ie. w/ NO FF. 
            ; pflfil='DEF' means tch16087o_pfl.fits for post-sm4? nope x6417094o_pfl is used
			calstis,9999,data=pflat,header=hp,'apttab=NONE,pfltfile='+pfltfil+',trace=0,soffset='+e1+ $
			    ',targoffs='+string(-toffix*.0507)+',gwidth='+gwidth+',ctecorr=0,hrepair=-.04,' $
			    'outspec=none',autowave=toffset,hspc,orders,wd,fd,eps,err,gross,back,nd

skipcal:	; skip normal extraction of fringe flat.
        	rnodes=wnodes*0+1
			nd(1020:1023)=avg(nd(1017:1020))	; last 4 pts bad
            yfit=SPLINEFIT(wd,nd,wd*0+1,wnodes,rnodes)
	        flat=nd/yfit            ;norm to unity
            ;			plot,wd,nd  &  oplot,wd,yfit  &  stop ; for bad px
	        ind=where(wd gt 6620)	; 97dec3
            ;reduce tung fringes by 11%. 2014jun19-No! only by ~2% for .8 or 1.2 high fringe
	        flat(ind)=flat(ind)*.89+.11
			print,'G750L SHORT SLIT pflat derived from ',name,' gwidth=',gwidth
; END skipcal
		
		end else begin		; no short slit fringe flat:
            ; CASE OF LONG SLIT FLATS --> result=pflat.fits:
            ; extract wave,flux from stellar data to make findfringe work.
            ; 2019apr4 - mod from gwidth=7 to current gwidth:
			print,'No short Slit Fringe Flat. Do Calstis, findfinge.pro at row offset=',e1

	   	    calstis,9999,data=crej,header=hdcr,errin=crerr,'apttab=NONE,darkfile=def,gwidth='+ $
	   	        gwidth+',bdist=300,trace=0,soffset='+e1+',dist1='+dist1+',pfltfile=none,hrepair=-.04,'+	$
			    'targoffs='+string(-toffix*.0507)+',outspec=tmp.fits',autowave=toffset,		$
                ; Fails.		'senstab='+senstab+','+ 			$ ;20jul
                ;			',soffset=0,subback=0',	$	; for wavecal
			    hspc,orders,wave,flux,eps,err,gross,back,net,pos,image

			print,'Run findfringe.pro'
		    sxaddhist,'PFLAT from FINDFRINGE w/ CROSS_CORRELATE',hdcr

            ; median ,wid=5 and filt both data and flat. Prob not important.
            ; width default = 15, zero if no good xcorel...
			findfringe,flux,wave,fringe,fringeh,exttab='tmp.fits',/median,file='pflat.fits', $
			    /verbose,/shift,/CROSS_CORRELATE
			pfltfil='pflat.fits'		; FINDFRINGE flat
			print,'G750L LIBRARY pflat-fringeflat written'
		endelse			; end library fringe flat
	endif				;END G750L fringeflat processing
skipfring:

	bdist='300'				; 97jun6
	if gain gt 0 then begin		; CCD processing:
        ; HZ43 in G750L, where red star is brighter & contrib sig to bkg. Use
        ;	bdist='88', as bupper is 44px from hz43b, which is 44 px from HZ43!
		if targ eq 'HZ43' and mode eq 'G750L' and root ne 'o69u07030' then begin
			dist1='10'
			bdist='88'
		endif ; 3 spectra

		if targ eq '2M0559-14' then dist1='25'
        ; 02dec13-L-flat for E1 position
        ; 05sep21 - omit for G230*B modes, as there is now a real L-flat
        ; 06may26 - omit for all MODES ................................. 
        ;;;		if strpos(aperture,'E1') gt 0 and strpos(mode,'B') 	$
        ;;;								lt 0 then begin
        ;;;			if mode ne 'G750L' then	begin
        ;;;				rdf,'lflat-e1.'+strlowcase(mode),1,lflt
        ;;;				lflat=lflt(*,1)
        ;;;			     end else lflat=0.974	; temp constant
        ; make the L-flat correction to lines 850:950, in region of E1 psuedo aperture:
        ;;;			for i=850,950 do crej(*,i)=crej(*,i)/lflat
        ;;;			sxaddhist,'PREPROC applied special L-flat at E1',hdcr
        ;;;			endif
		if flg7px eq 1 then ext='.7pix'
		senstab=sxpar(hdcr,'senstab')
        ; 2019mar28-new g230lb & g430l defaults. Prefer gwidth=11 for sat wide gwidths
		if long(gwidth) ge 11 and mode eq 'G230LB' then	senstab='sens11_g230lb.fits'
		if long(gwidth) ge 11 and mode eq 'G430L' then senstab='sens11_g430l.fits'
		if long(gwidth) ge 11 and mode eq 'G750L' then senstab='sens11_g750l.fits'
		sxaddpar,hdcr,'senstab',senstab			; 2019mar26

        ; ###change
        ;bdist='19'
        ; crej NOT modified by calstis. Calstis corrects Flux, NOT net for CTE, etc.
		fltparam=''				; for raw files
		fitspos=strpos(fils(ifil),'.fits')
		filtyp=strmid(fils(ifil),fitspos-3,3)

        ; pixel based cte corr test:
		if strmid(filtyp,0,2) eq 'fl' or filtyp eq 'cte' then begin
			fltparam=',biasfile=none,darkfile=none,pfltfile=none,lfltfile=none'
			if filtyp eq 'flc' or filtyp eq 'cte' then sxaddpar,hdcr,'pixelcte','YES',after='RAD_VEL','From px-based CTE corr flcfile'
		endif

        ; CCD processing: MAIN CALL---- 
		if strpos(root,'ocy5') eq 0 then begin		;2020dec9
			ext='-'+mode+'-'+targ+'.fits' 		;14141 Gaia hot
			if targ eq 'GJ-894.3' then ext='-'+mode+'-FEIGE110.fits'
		endif
		if strpos(root,'oe3f') eq 0 and targ ne 'CCDFLAT'then ext='-'+mode+'-'+star(ifil)+'.fits' ; Gaia 15816

        ; 2020nov10 - Reduce bdist to 90 for G430L-E1 & most Gaia G430L G430M:
        ; 2020nov23 - trace too close to edge for std bdist=300:
		if strpos(mode,'G430') ge 0 and	abs(e1) gt 203 then bdist='90'

		print,'Main CCD CALSTIS: START.'
        ; Calstis bug w/ trace=1 --> first CCD stiscr /display window is erased for a
        ;	fresh IDL session. Just rerun to see the /displ.

        ; defaults:   'bwidth=11,b_order='+bkgfit+',b_mean1=1,b_mean2=1,b_median=9,'+$
        ; ###change - special test: NOTE: after 06jul31, corr for o3tt below,
        ;	invalidates this special test:
        ;		 'exttab=o3tt43040_ext.fits,'+	$	; trace table
        ;		 ',soffset=0,subback=0',	$	;for wavecal
		calstis,9999,data=crej,header=hdcr,errin=crerr,'subback='+subkg+',dist1='+dist1+','+ $
	        'gwidth='+gwidth+',bdist='+bdist+',senstab='+senstab+',pfltfile='+pfltfil+ $
	        ',trace=0,soffset='+e1+',targoffs='+string(-toffix*.0507)+',hrepair='+hotpix+  $
	        'outspec='+subdir+'/spec_'+root+ext+fltparam,autowave=toffset, $
		    hspc,orders,wave,flux,eps,err,gross,back,net,pos,image

        ;above image is calibrated to ct/s & flat fielded. Use crej to edit input above.
		if targ eq 'CCDFLAT' then goto,skipfixes

        ; Special Gaia bkg adjustments:
		if (strpos(fils(ifil),'oe3f06010') gt 0 or 		$
			strpos(fils(ifil),'oe3f05010') gt 0 or 	$
			strpos(fils(ifil),'oe3f03010') gt 0)	$
			and aper(ifil) eq '52X2low' 
		then begin
			z=mrdfits(subdir+'/spec_'+root+ext,1,hdr)
			net=z.net
			bkg=z.gross-net
			if strpos(fils(ifil),'oe3f06010') gt 0 then net(491)=0.05	; fix bad dark subtr
			if strpos(fils(ifil),'oe3f03010') gt 0 then	net=net-0.02	;bad bkg drop @ spectrum
			if strpos(fils(ifil),'oe3f05010') gt 0 then begin
				net=net-0.02	; at Row=15
				net(180:181)=0	; bad dark subtr.
			endif
			z.gross=net+bkg
			z.net=net
			mwrfits,z,subdir+'/spec_'+root+ext,hdr,/create  ;fix net
			stisflx,subdir+'/spec_'+root+ext,hdr,wdum,flux,/ttcor
			z.flux=flux
			mwrfits,z,subdir+'/spec_'+root+ext,hdr,/create ;fix flux
		endif

        ; 2011mar28-add GAC=grating, aperture corr. see doc.lpflats and 
        ; 2019mar-E1 GAC revised. See stis/doc/abscor.pct
        ; skip ONLY 52x2 cent. Orig GAC file does all narrow apertures:
        ; 2020feb3 - also skip J..Maiz 2nd targets below E1 & flagged as 52X2low:
		if aperture ne '52X2' and aper(ifil) ne '52X2low' then begin
			z=mrdfits(subdir+'/spec_'+root+ext,1,hdr)
			net=z.net
			bkg=z.gross-net
            ; Gen. purpose gac11:
			ind=where(gacaper eq aperture and gacmode eq mode)
			ind=ind(0)
			if ind eq -1 then ind=27

			print,'GAC corr for ',gacaper(ind),' ',gacmode(ind),' gwidth=',gwidth
            ; 7px GAC needed for G430L, gwidth=7; but 11px GAC good for G230LB 7&11px gwidth
            ;	Wide sat. obs. at E1 are NO GOOD. 2019mar
            ; Special case of using orig Proffit as updated for G230LB&G430L 52X2E1 2011mar:
			if gwidth eq '7' and mode eq 'G430L' then linterp,gacwl7(*,ind),gacth7(*,ind),	$
					z.wavelength,thcorr	$
			    else					$
                ; General case of using my gwidth=11 E1, G230LB & G430L (only) update of 2019mar
                ;	ie G750L E1 gets orig proffitt result.
				    linterp,gacwl(*,ind),gacth(*,ind),z.wavelength,thcorr
        	
        	z.net=net/thcorr
		    z.gross=z.net+bkg
			if gwidth eq '7' and mode eq 'G430L' then sxaddhist,'GAC Corr of NET & FLUX per '+$
					sxpar(hgac70,'filename'),hdr	$
		    	else					$
				    sxaddhist,'GAC Corr  of NET & FLUX per '+sxpar(hgac0,'filename'),hdr

			mwrfits,z,subdir+'/spec_'+root+ext,hdr,/create  ;fix net
            ; Main flux cal is here, 2019aug28-hdum-->hdr, ie write stisflx hdr updates.
			stisflx,subdir+'/spec_'+root+ext,hdr,wdum,flux,/ttcor
			z.flux=flux
			mwrfits,z,subdir+'/spec_'+root+ext,hdr,/create ;fix flux
		endif				; END E1 GAC corr.

        ; 2011jun20-add red scat lite corr for G230LB, 52x2. See /stis/doc/scat.ccdmodes
        ; 2019apr1 - Mv all scat lite calc to stisflx, but do not write corr. NET.
        ;	Add stisflx corr for all CCD modes, so above CALSTIS flux is irrelevant.

		z=mrdfits(subdir+'/spec_'+root+ext,1,hdr)

        ; scat lite is subtracted in stisflx for G230LB
		stisflx,subdir+'/spec_'+root+ext,hdr,wdum,flux,/ttcor
		z.flux=flux		; scat lite corr to flux only
		sxaddhist,'NET is NOT corrected in this file, except for GAC at E1.',hdr
		mwrfits,z,subdir+'/spec_'+root+ext,hdr,/create ;fix flux

skipfixes:
	end else begin				;end CCD, start MAMA

        ; MAMA details for CALSTIS:
		print,'Main MAMA CALSTIS: START.'

        ; 02mar26 - G140L on repeller wire at center before ~jul1, 1997;
		if mode eq 'G140L' and float(date) gt 1997.181 then begin		; See dispersion.doc

            ; 02mar18 corr. for -3, bigger wls. 02sep11 - should switch to targoffs??:
			if sxpar(hd,'OMSCYL1P') lt 800 then	toffset.offsets=toffset.offsets-0.15	$
                ; 02mar18 corr. for +3, smaller wls:
			    else toffset.offsets=toffset.offsets+0.15
			
			sxaddhist,'G140L wavelengths corr by 0.15px for MSM offset of 3 arcsec',hd
		endif

        ; defaults		gwidth='11'				; MAMA default.
        ; apparently		if strmid(mode,0,1) eq 'E' then gwidth='7'
        ; ##change  for hires extractions:
        ;		hires='pfltfile=none,hires=1,exttab=spec_'+root+'.lores,'
		hires='pfltfile='+pfltfil+',hires=0,'
		sxaddpar,hd,'pstrtime',pstrtim,'predicted obs. start time from SPT file'
		sxaddpar,hd,otemp,OMCAT,comment  ; 20011 raw data has om1cat
        sxaddpar,hd,'rad_vel',svel,'Star Velocity (km/s) NOT applied'
		sxaddpar,hd,'date',!stime
		sxaddhist,'Wavecal file: '+wfil,hd
		if toffix ne 0 then sxaddhist,'STISWLFIX OFFSET(px)='+string(toffix,'(f5.2)'),hd
		sxaddhist,'Written by preproc.pro '+!stime,hd
		
        ; trace=1 fails sometimes for hires w/ err about bpos. uses lores pos for extr
        ; 00nov21 - fix for mult MAMA readout, eg. o47z01060,40 P330E:
        fits_open,fils(ifil),fcb
            extnames = strtrim(fcb.extname)
            crsplit = long(total(extnames eq 'SCI'))
        fits_close,fcb
		if crsplit ne 1 then begin
			print,' *** WARNING *** MAMA CRSPLIT=',CRSPLIT,' FOR ',FILS(ifil),' ********************'
			infil=0
			exptime=0
			integ=0
			for ii=0,crsplit-1 do begin
				stis_read,fils(ifil),hdum,d,readout=ii+1
				infil=infil+d
				exptime=exptime+sxpar(hdum,'exptime')
				integ=integ+sxpar(hdum,'integ')
			endfor
			print,' *** TOTAL EXPTIME = ',EXPTIME
			sxaddpar,hd,'exptime',exptime
			sxaddpar,hd,'integ',integ
		endif			; end special fix
		soffst='0'			; normal default
        ; soffsets for M33 - UIT339: G140L=-137
		if root eq 'o5co02020' or root eq 'o5co02030' then soffst='-137,dist1=20'
        if root eq 'obrr02030' then infil(1784:1785,1142:1144)=0	;a hit or hot px
        if root eq 'oc8c20020' then infil(344,753)=19		;a hit or hot px
		bdist='300'
        ;2018aug31-Special Bkg test for Clayton (May be good for any faint exp.)
		if strpos(root,'od6j') ge 0 then begin
			bdist='50'
			print,'SPECIAL MAMA BKG TEST FOR bdist='+bdist
		endif
		mincnt='500'  &  dist1='200'
        ; 2021mar22 - weak GD71 G140L vignetting obs:
		if strpos(root,'o5300') eq 0 then  begin
			mincnt='50'
			dist1='500'
		endif
        ; MAIN MAMA EXTRACTion
		calstis,9999,data=infil,header=hd,errin=errin,udl=udl,'targoffs='+string(-toffix*.0246)+ $
		    ',pfltfile='+pfltfil+',bdist='+bdist+',soffset='+soffst+',gwidth='+gwidth+','+ $
            ; ###change - normally use dist1=200 default
            ; hires+     					$
			'dist1='+dist1+',mincounts='+mincnt+',outspec='+subdir+'/spec_'+root+ext+','+ $
            ; 'gwidth=26,'+ 		$
            ; 02aug5-new defaults:'bwidth=11,b_order='+bkgfit+',b_mean1=1,b_mean2=1,'+$
	   		'trace=0',autowave=toffset,hspc,orders,wave,flux,eps,err,gross,back,net,pos,image
	endelse			; data are now extracted for CCD & MAMAs

	print,'Main CALSTIS: FINISH: ',ifil,' of ',n_elements(fils)-1,form='(a,i4,a,i4,//)'

    ; 07may2 - fix the spike at 3100+A in VB8, g430l:		
    ;Pixels in column 91 are very hot in the dark and are over-subtracted.  because
    ;of hot pixels both in columns 90 and 92, the column 91 pixels are un-repairable
    ;(The current algorithm requires at least one good left or right neighbor).
	if root eq 'o8ui03010' then begin		; VB8 - G430L
		z=mrdfits(subdir+'/spec_o8ui03010.fits',1,hdr)
		bkg=z.gross-z.net
        ; 2011aug-not		bad=where(z.wavelength lt 2950)		; big noise
        ;   so bad now?		z.epsf(bad)=254
		z.net(91:92)=0.				; the big neg spike
		z.flux(91:92)=0.
		z.gross=z.net+bkg	; maintain smooth bkg
		sxaddhist,'Big negative spike at 3135-3138A set to zero.',hdr
		mwrfits,z,subdir+'/spec_o8ui03010.fits',hdr,/create;fix flux
	endif

    ; 2018jul23-spotted the bad neg. glitch in o8kh03020 G750L 2m0036 that is NOT
    ;	negl. in crej image & must get screwed by attempted auto-fix: (eps=170)
	if root eq 'o8kh03020' then begin		; G750L 2m0036...
		z=mrdfits(subdir+'/spec_o8kh03020.fits',1,hdr)
		bkg=z.gross-z.net
		z.net(180)=.07
		z.gross=z.net+bkg
		z.flux(180)=2e-17
		sxaddhist,'Big negative spike at 6140A fixed.',hdr
		mwrfits,z,subdir+'/spec_o8kh03020.fits',hdr,/create
	endif

    ;
    ;bohlin special G750L 1-D Tungsten short slit flat application:
    ; 
	if shrtflg eq 1 and targ ne 'CCDFLAT' then begin	; Defringe G750L
		z=mrdfits(subdir+'/spec_'+root+ext,1,hdr)
		wave=z.wavelength
        ; remove heliocentric radial vel corr w/ backwrd djl conven. on earthvel
        ;	CALSTIS converts to heliocentric WLs by default.
		evel=sxpar(hdr,'earthvel')
		wavno=wave/(1+evel/3e5)			; Rm just for flux cal
		net=z.net
		bkg=z.gross-net				;Calstis OK net not corr
        ; HZ43 in G750L, where red star is brighter & contrib sig to bkg:
		if strpos(fils(ifil),'o57t01020') ge 0 or		$
		   strpos(fils(ifil),'o57t02020') ge 0 or		$
           ;;;		   strpos(fils(ifil),'o69u07030') ge 0 or	$  ; HZ43B faint
		   strpos(fils(ifil),'o69u08030') ge 0 $
		then begin
            ; hz43B is above hz43 by 44 px, use BU @ +88 for bkg to estim effect of B on A
			bu=z.bupper 
			bl=z.blower
			bkgup=median(bu,9)
			bfit=smooth(smooth(bkgup,39),39)
			bkg=bfit*gwidth/11.		; bwidth=11 default
			net=z.gross-bkg
			sxaddhist,'Special HZ43 Bkg. subtraction',hdr
			plot_io,wave,net,yr=[.01,100]
			oplot,wave,bkgup*gwidth/11 & oplot,wave,bkg,thic=2,lin=2
			oplot,wave,bl*gwidth/11.,lin=1
			XYouts,.15,.7,'HZ43',/norm
            ;			read,st
		endif
        ; norm data to unity to make cross-correl work;
		tmp=net/SPLINEFIT(wavno,net,wavno*0+1,wnodes,rnodes)
        ; - need larger search range for E1. automate 02nov12. 06jul31-backwards. Fix:
		if strpos(aperture,'E1') lt 0 then hrs_offset,tmp(650:935),flat(650:935),offset,0,5	$ ;normal
		   else	hrs_offset,tmp(650:935),flat(650:935),offset,0,16	 ; E1

		print,'Tung. flat offset in px=',offset
        ; flat is the 1-D fringe corr from above
	    linterp,indgen(1024)+offset,flat,indgen(1024),linflat
	    ind=where(wd le 6620)
        linflat(ind)=1            	; do not add noise at short WL
        ;		save,WAVE,FLUX,LINFLAT  &  stop	; for running fringtweak
	    net=net/linflat			; no CTE corr to numer or denom
        ; 06jul28 - fix early gd153 g750L for wide FWHM w/ results of o3tt4-corr.pro:
        ; 2014sep18 -Do More minor fix for gwidth=11, instead. see plots/14sep750gd153*,
        ;			which do NOT include ff fix;
		if strpos(root,'o3tt4') ge 0 then begin
			wfix=[5000,9600,10230]
			ffix=[1,1.,0.95]
			linterp,wfix,ffix,wave,linfix
			net=net/linfix
			z.flux=z.flux/linfix
			z.errf=z.errf/linfix
			print,root +' ***FIXED for wide spectrum problem***' 
		endif

        ; then replace defringed net in output:
		z.net=net		; corr for flat only
		z.gross=net+bkg		; 03dec31
        ; flatten linearised flx, while stisflx fluxes flattened GROSS. But systematics
        ;	should come out in calib wash, as this written flux is not used.
		z.flux=z.flux/linflat
		z.errf=z.errf/linflat	; 04jan2
		sxaddhist,'Flat at WL>6620A is .89+.11 of Tungsten fringes',hdr
		mwrfits,z,subdir+'/spec_'+root+ext,hdr,/create
		flux=z.flux  &  gross=z.gross	; for plot below.
	endif				; end G750L fringe flats

    ; ###change
    ;	goto,skipit				; no Plots

    ; DISPLAY THE SPECTRUM
	!ytitle='FLUX'
	!mtitle=mode+' '+cenwave+' '+targ+' '+aperture+' '+root
	ngood=n_elements(wave)-3
	window,1			; because,0 kills the trace window.
	plot,wave,flux,yr=[0,max(flux(25:-20))]
	back=gross-net
    ;	oplot,wave,back
	if flg7px eq 1 then xyouts,.15,.5,'Special 7 pix high extraction',/norm
	plotdate,'preproc'
	read,st			; stop and wait for <cr> before contin.
skipit:
endfor

;print,nowavcal,' CCDflats w/ NO wavecal.'
print,'Total Clock time in Minutes',(systime(1)-tbeg)/60.
end
