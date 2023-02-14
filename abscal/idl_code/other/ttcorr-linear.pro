pro ttcorr,hd,wav,flux,notime=notime,version=version,tcorr=tcorr,silent=silent,$
	wcorr=wcorr,time=time,notemp=notemp
;+
;
; 98jul17 - time and temperature correction for STIS low dispersion
; INPUT:
;	hd - stis header
;	wav - Instrumental. obs wavelength vector w/ heliocentric vel removed.
;	flux - vector to correct, eg NET or FLUX
;	notime - keyword: do not make time corr. (just temp corr)
;	notemp - do not make temp correction for tcorrel.pro
; OUTPUT:
;	flux - as corrected for time and/or temp changes
;	tcorr - keyword w/ value of temp correction
;	wcorr - the time correction, added by DJL
;	time - keyword: for time of obs.
; REFERENCE - Bohlin, R. 1999,Instrument Science Report, STIS 99-07.
; HISTORY
;	99dec7 - DJL, added version, corrections, time, and silent keyword
;		parameters, modified to work for Medium and Echelle
;	Jan 21, 2000 - DJL, modfied to use EXPSTART keyword when PSTRTIME is
;		missing
;	01mar - rcb, tidy up comments,etc.
;	02apr4  - use two discontinuous segments for G140L
;	02apr4  - add notemp keyword for tcorrel.
;	02jul23 - L-flat for G140L & 0.36 PX WL fix for G230L
;	02nov17 - add 2 WL bins for G230L & one for G430L
;	02nov18 - go to 3 line segments for G230L
;	03mar 1 - add temp corrections for G230LB and G430L
;	03may1  - new flats, new WLs, new bkg for CTE corr
;	03May15 - added prism support
;	03Dec   - new time nodes for G140L, G230Ls
;	05may   - new time, temp corr coef.
;	05oct3	- G750L is solid 3sig result for temp change... add this code.
;	05dec   - new time, temp corr coef. per new CTE
;	06mar2  - .............................new Trace files (slight changes)
;	06apr13 - install G430L and G750L L-flats (~negligible diffs)
;	06jul31 - 6th digit update for G750L per early GD153 fix
;	08nov24 - few 4th digit changes to G140L
;	09jul20 - new 2009 pts
; 09sep11-'big' fix of G=4 ratio & avg occdhav ~09aug7 & ctecorr tweak
;	10May   - new break points
;	11Aug   - new segments & and endpts for G230L & LB
;	14Feb - new G230L 2002.6, instead of 2002.5 
;		& new 2011.5 endpt for G230*)
;	14jul - add new 1999.3 pt for G430 & 750L, plus 2009.5 for G750L
;	14sep18 - add gwidth=11 time corrections.
;	14sep22 - add gwidth=11 temperture corrections.
;	14dec23 - New monitor obs in 2014Nov.
;	2017jan27 - fix ttcorr for g750L tchange for ht=7 HZ43. Vega?
; 2018jun8-11 - add nodes G140L-2001.8,2004.7,2009.5 & G230L,B 2013.3,2015.6
; 2019mar20 - add gwidth=11 corr for g230lb and g430L & use for all gwid>=11
; 2019mar25 - change a fit endpoint for G430L to 1999.0
;-
	VERSION = 'Bohlin/March, 2019 ttcorr '

;G140L time change array
; 1175. 1225. 1275. 1325. 1375. 1425. 1475. 1525. 1575. 1625. 1675.
tc140l=[								$
;   1997.38   1999.20   2001.80   2009.50   20120.00
;  1.000000,  0.966988,  0.964673,  0.901645,  0.772346,  0.850827,  0.779294, $
;  1.000000,  0.971823,  0.972107,  0.950810,  0.899910,  0.938397,  0.911760, $
;  1.000000,  0.985747,  0.983394,  0.956951,  0.914237,  0.949700,  0.930063, $
;  1.000000,  0.990678,  0.984822,  0.951582,  0.906370,  0.941707,  0.923816, $
;  1.000000,  0.985150,  0.985647,  0.945475,  0.880665,  0.927204,  0.907394, $
;  1.000000,  0.982486,  0.973249,  0.925921,  0.848003,  0.899682,  0.871722, $
;  1.000000,  0.968782,  0.961380,  0.907873,  0.814958,  0.871342,  0.846710, $
;  1.000000,  0.953481,  0.950865,  0.889327,  0.767418,  0.845561,  0.822117, $
;  1.000000,  0.952365,  0.949290,  0.879353,  0.766653,  0.835694,  0.827303, $
;  1.000000,  0.955784,  0.948308,  0.869363,  0.765240,  0.829082,  0.819714, $
;  1.000000,  0.982713,  0.980205,  0.892585,  0.795794,  0.853005,  0.837508]
; 18jun-Nodes revised
;  1.000000,  0.966988,  0.964673,  0.901645,  0.772346,  0.850897,  0.779132, $
;  1.000000,  0.971823,  0.972107,  0.950810,  0.899910,  0.938993,  0.910380, $
;  1.000000,  0.985747,  0.983394,  0.956951,  0.914237,  0.949747,  0.929956, $
;  1.000000,  0.990678,  0.984822,  0.951582,  0.906370,  0.941673,  0.923895, $
;  1.000000,  0.985150,  0.985647,  0.945475,  0.880665,  0.927550,  0.906594, $
;  1.000000,  0.982486,  0.973249,  0.925921,  0.848003,  0.899956,  0.871088, $
;  1.000000,  0.968782,  0.961380,  0.907873,  0.814958,  0.871970,  0.845256, $
;  1.000000,  0.953481,  0.950865,  0.889327,  0.767418,  0.846278,  0.820459, $
;  1.000000,  0.952365,  0.949290,  0.879353,  0.766652,  0.836916,  0.824475, $
;  1.000000,  0.955784,  0.948308,  0.869363,  0.765240,  0.830052,  0.817470, $
;  1.000000,  0.982713,  0.980205,  0.892585,  0.795795,  0.853771,  0.835736]
; 18aug
;  1.000000,  0.966959,  0.964659,  0.901629,  0.772281,  0.850277,  0.773086, $
;  1.000000,  0.971874,  0.972131,  0.950832,  0.900032,  0.938878,  0.907763, $
;  1.000000,  0.985748,  0.983387,  0.956950,  0.914204,  0.948803,  0.930234, $
;  1.000000,  0.990674,  0.984827,  0.951584,  0.906374,  0.940724,  0.924401, $
;  1.000000,  0.985154,  0.985651,  0.945478,  0.880672,  0.926798,  0.906273, $
;  1.000000,  0.982483,  0.973251,  0.925915,  0.848038,  0.899321,  0.869625, $
;  1.000000,  0.968782,  0.961393,  0.907874,  0.814934,  0.871448,  0.843741, $
;  1.000000,  0.953488,  0.950870,  0.889336,  0.767486,  0.845972,  0.818541, $
;  1.000000,  0.952350,  0.949287,  0.879342,  0.766629,  0.836921,  0.823123, $
;  1.000000,  0.955786,  0.948294,  0.869355,  0.765218,  0.830289,  0.815560, $
;  1.000000,  0.982697,  0.980182,  0.892572,  0.795820,  0.854060,  0.833051]
; 19jan
  1.000000,  0.966959,  0.964659,  0.901629,  0.772281,  0.850832,  0.771768, $
  1.000000,  0.971874,  0.972131,  0.950832,  0.900033,  0.939034,  0.907391, $
  1.000000,  0.985748,  0.983387,  0.956950,  0.914204,  0.949430,  0.928745, $
  1.000000,  0.990674,  0.984827,  0.951584,  0.906374,  0.941103,  0.923501, $
  1.000000,  0.985154,  0.985651,  0.945478,  0.880672,  0.927184,  0.905357, $
  1.000000,  0.982483,  0.973251,  0.925915,  0.848038,  0.899760,  0.868581, $
  1.000000,  0.968782,  0.961393,  0.907874,  0.814934,  0.872073,  0.842254, $
  1.000000,  0.953487,  0.950870,  0.889336,  0.767486,  0.846979,  0.816146, $
  1.000000,  0.952350,  0.949286,  0.879342,  0.766629,  0.838160,  0.820178, $
  1.000000,  0.955786,  0.948294,  0.869355,  0.765218,  0.831793,  0.811987, $
  1.000000,  0.982697,  0.980182,  0.892572,  0.795820,  0.855817,  0.828876]
; 19may

; G230L time change array 
;   1997.38   1998.84   1999.50   2002.70   2012.00  -- old
;   1997.38   1998.12   1999.05   2002.70   2012.00  -- old
;new 14feb17 1997.38   1998.40   1999.30   2002.60   2009.50   2011.20   2018.00
;   1997.38   1998.40   1999.30   2002.60   2009.50   2011.20   2013.30
;   2015.60   2020.00	2019march (& earlier)

; 1650. 1750. 1850. 1950. 2050. 2150. 2250. 2350. 2450. 2550. 2650. 2750. 2850.
; 2950. 3050.
tc230l=[								$
;  1.000000,  1.011024,  1.020722,  0.960947,  0.966704,  0.930445,  0.939458, $
;  0.959502,  1.000000,  1.012274,  1.021680,  0.947009,  0.915203,  0.919745, $
;  0.933370,  0.928962,  1.000000,  1.015462,  1.032093,  0.948955,  0.929294, $
;  0.922227,  0.944088,  0.911288,  1.000000,  1.015680,  1.032157,  0.940112, $
;  0.927317,  0.917809,  0.935681,  0.884035,  1.000000,  1.019259,  1.029316, $
;  0.936455,  0.923811,  0.916174,  0.923955,  0.871838,  1.000000,  1.022032, $
;  1.025069,  0.945976,  0.922642,  0.924716,  0.925151,  0.883298,  1.000000, $
;  1.031599,  1.024362,  0.957351,  0.926651,  0.933735,  0.929186,  0.896951, $
;  1.000000,  1.024863,  1.014345,  0.957684,  0.921375,  0.926651,  0.929405, $
;  0.895156,  1.000000,  1.011182,  1.006599,  0.959320,  0.913929,  0.919319, $
;  0.929316,  0.895696,  1.000000,  1.007585,  1.003458,  0.963206,  0.923110, $
;  0.924155,  0.936771,  0.902525,  1.000000,  1.014336,  1.007400,  0.969297, $
;  0.948667,  0.935792,  0.944876,  0.916090,  1.000000,  1.013165,  1.010371, $
;  0.974767,  0.953537,  0.947354,  0.948456,  0.922605,  1.000000,  1.013200, $
;  1.008319,  0.975553,  0.961113,  0.952952,  0.952729,  0.930675,  1.000000, $
;  1.018072,  1.009831,  0.981213,  0.946860,  0.954944,  0.960392,  0.932131, $
;  1.000000,  1.017931,  1.012488,  0.984003,  0.952612,  0.959825,  0.963944, $
;  0.937494]				; 18jun
;  1997.38   1998.40   1999.30   2002.60   2009.50   2011.20   2013.30   2015.60
;   2019.00
;  1.000000,  1.011024,  1.020722,  0.960947,  0.968983,  0.930445,  0.937014, $
;  0.948501,  0.952065,  0.957343,  1.000000,  1.012274,  1.021680,  0.947009, $
;  0.918011,  0.919744,  0.931475,  0.937401,  0.928564,  0.931130,  1.000000, $
;  1.015462,  1.032093,  0.948955,  0.933658,  0.922227,  0.945503,  0.935358, $
;  0.912941,  0.927097,  1.000000,  1.015680,  1.032157,  0.940112,  0.929955, $
;  0.917809,  0.939009,  0.919412,  0.889832,  0.906894,  1.000000,  1.019259, $
;  1.029315,  0.936455,  0.922984,  0.916174,  0.926952,  0.910145,  0.881302, $
;  0.891185,  1.000000,  1.022032,  1.025069,  0.945976,  0.922560,  0.924716, $
;  0.926056,  0.917769,  0.894965,  0.892858,  1.000000,  1.031599,  1.024362, $
;  0.957351,  0.926221,  0.933735,  0.927758,  0.929323,  0.906491,  0.904222, $
;  1.000000,  1.024863,  1.014345,  0.957684,  0.923117,  0.926651,  0.929276, $
;  0.922458,  0.906813,  0.899380,  1.000000,  1.011182,  1.006599,  0.959320, $
;  0.918563,  0.919319,  0.929327,  0.919158,  0.907327,  0.897443,  1.000000, $
;  1.007585,  1.003458,  0.963206,  0.929377,  0.924155,  0.937302,  0.923276, $
;  0.913848,  0.903722,  1.000000,  1.014336,  1.007400,  0.969297,  0.954672, $
;  0.935792,  0.944853,  0.936321,  0.926546,  0.916030,  1.000000,  1.013165, $
;  1.010371,  0.974767,  0.959363,  0.947354,  0.949170,  0.941560,  0.933249, $
;  0.923901,  1.000000,  1.013200,  1.008319,  0.975552,  0.964000,  0.952952, $
;  0.950179,  0.953586,  0.937852,  0.932491,  1.000000,  1.018072,  1.009831, $
;  0.981213,  0.950554,  0.954944,  0.960436,  0.951793,  0.942116,  0.934141, $
;  1.000000,  1.017931,  1.012488,  0.984003,  0.951647,  0.959825,  0.962809, $
;  0.960401,  0.947423,  0.938434]		; 18jun new nodes
;  1.000000,  1.011024,  1.020722,  0.960947,  0.968983,  0.930445,  0.937014, $
;  0.948501,  0.953704,  0.950072,  1.000000,  1.012274,  1.021680,  0.947009, $
;  0.918011,  0.919744,  0.931475,  0.937401,  0.929242,  0.928121,  1.000000, $
;  1.015462,  1.032093,  0.948955,  0.933658,  0.922227,  0.945503,  0.935358, $
;  0.913693,  0.923760,  1.000000,  1.015680,  1.032157,  0.940112,  0.929955, $
;  0.917809,  0.939008,  0.919412,  0.890507,  0.903897,  1.000000,  1.019259, $
;  1.029316,  0.936455,  0.922984,  0.916174,  0.926952,  0.910145,  0.881856, $
;  0.888731,  1.000000,  1.022032,  1.025069,  0.945976,  0.922560,  0.924716, $
;  0.926056,  0.917769,  0.895145,  0.892059,  1.000000,  1.031599,  1.024362, $
;  0.957351,  0.926220,  0.933735,  0.927758,  0.929323,  0.906773,  0.902974, $
;  1.000000,  1.024863,  1.014345,  0.957684,  0.923117,  0.926651,  0.929276, $
;  0.922458,  0.907051,  0.898325,  1.000000,  1.011182,  1.006599,  0.959320, $
;  0.918563,  0.919319,  0.929327,  0.919158,  0.907699,  0.895793,  1.000000, $
;  1.007585,  1.003458,  0.963206,  0.929377,  0.924155,  0.937303,  0.923276, $
;  0.914207,  0.902131,  1.000000,  1.014336,  1.007400,  0.969297,  0.954672, $
;  0.935792,  0.944853,  0.936321,  0.926787,  0.914959,  1.000000,  1.013165, $
;  1.010371,  0.974767,  0.959363,  0.947354,  0.949170,  0.941560,  0.933673, $
;  0.922018,  1.000000,  1.013200,  1.008319,  0.975552,  0.963999,  0.952952, $
;  0.950179,  0.953586,  0.938003,  0.931821,  1.000000,  1.018072,  1.009831, $
;  0.981213,  0.950554,  0.954944,  0.960436,  0.951793,  0.942308,  0.933291, $
;  1.000000,  1.017931,  1.012488,  0.984003,  0.951647,  0.959825,  0.962809, $
;  0.960401,  0.947892,  0.936354]		 ; 18aug
  1.000000,  1.011024,  1.020726,  0.960933,  0.969058,  0.930421,  0.937038, $
  0.948507,  0.954899,  0.941445,  1.000000,  1.012270,  1.021676,  0.947018, $
  0.917989,  0.919774,  0.931470,  0.937408,  0.929559,  0.925514,  1.000000, $
  1.015467,  1.032099,  0.948956,  0.933649,  0.922229,  0.945519,  0.935370, $
  0.914087,  0.923669,  1.000000,  1.015686,  1.032156,  0.940115,  0.929983, $
  0.917822,  0.939019,  0.919418,  0.890776,  0.905570,  1.000000,  1.019245, $
  1.029309,  0.936446,  0.922974,  0.916173,  0.926948,  0.910146,  0.882087, $
  0.888468,  1.000000,  1.022040,  1.025072,  0.945982,  0.922503,  0.924737, $
  0.926056,  0.917781,  0.895498,  0.888353,  1.000000,  1.031584,  1.024356, $
  0.957339,  0.926266,  0.933723,  0.927759,  0.929321,  0.906940,  0.900253, $
  1.000000,  1.024878,  1.014359,  0.957698,  0.923127,  0.926667,  0.929301, $
  0.922475,  0.906200,  0.900179,  1.000000,  1.011183,  1.006604,  0.959324, $
  0.918574,  0.919332,  0.929333,  0.919161,  0.907833,  0.891079,  1.000000, $
  1.007586,  1.003463,  0.963206,  0.929303,  0.924164,  0.937311,  0.923288, $
  0.914440,  0.897023,  1.000000,  1.014323,  1.007387,  0.969290,  0.954616, $
  0.935788,  0.944843,  0.936316,  0.927038,  0.909500,  1.000000,  1.013172, $
  1.010366,  0.974761,  0.959432,  0.947360,  0.949177,  0.941569,  0.934089, $
  0.916058,  1.000000,  1.013200,  1.008328,  0.975572,  0.963918,  0.952959, $
  0.950195,  0.953598,  0.939110,  0.923943,  1.000000,  1.018073,  1.009836, $
  0.981214,  0.950520,  0.954978,  0.960442,  0.951803,  0.941956,  0.932830, $
  1.000000,  1.017938,  1.012476,  0.984000,  0.951610,  0.959855,  0.962820, $
  0.960419,  0.946969,  0.937520]		;19may

; G230LB time change array
;   1997.38   1999.00   2000.50   2002.50   2009.50   2011.20   2018.00
; 1750. 1850. 1950. 2050. 2150. 2250. 2350. 2450. 2550. 2650. 2750. 2850. 2950.
tc230lb=[								$
;  1997.38   1999.00   2000.50   2002.50   2009.50   2011.20   2013.30   2015.60
;   2020.00
;  1.000000,  1.011789,  0.986396,  0.931878,  0.882166,  0.895944,  0.901659, $
;  0.899747,  0.896614,  0.882799,  1.000000,  1.028070,  1.002305,  0.939842, $
;  0.901619,  0.907495,  0.917215,  0.903099,  0.884204,  0.891070,  1.000000, $
;  1.028542,  0.997674,  0.933212,  0.894926,  0.904819,  0.913949,  0.890092, $
;  0.861423,  0.876450,  1.000000,  1.020855,  0.984892,  0.924484,  0.890529, $
;  0.897345,  0.898400,  0.878698,  0.852213,  0.856514,  1.000000,  1.017688, $
;  0.986462,  0.935272,  0.899001,  0.900801,  0.903275,  0.887802,  0.868248, $
;  0.860721,  1.000000,  1.016530,  0.988986,  0.945879,  0.906692,  0.906312, $
;  0.908931,  0.896586,  0.884092,  0.868331,  1.000000,  1.013696,  0.990371, $
;  0.953302,  0.907764,  0.910852,  0.913218,  0.902778,  0.891331,  0.872748, $
;  1.000000,  1.012690,  0.990044,  0.956787,  0.917377,  0.918191,  0.917497, $
;  0.907433,  0.898116,  0.879589,  1.000000,  1.011390,  0.993423,  0.964764, $
;  0.930912,  0.929004,  0.926571,  0.918262,  0.908627,  0.893646,  1.000000, $
;  1.008544,  0.992342,  0.968741,  0.935763,  0.935167,  0.933062,  0.923549, $
;  0.915229,  0.901549,  1.000000,  1.007472,  0.992683,  0.971942,  0.943764, $
;  0.939158,  0.935419,  0.928423,  0.920638,  0.906081,  1.000000,  1.002845, $
;  0.989209,  0.970376,  0.934759,  0.937781,  0.934360,  0.926872,  0.919243, $
;  0.905364,  1.000000,  0.998710,  0.987053,  0.968946,  0.938349,  0.938371, $
;  0.933682,  0.927399,  0.921172,  0.906905]		 ;19jan
;  1997.38   1999.00   2000.50   2002.50   2009.50   2011.20   2013.30   2015.60 ;				 2020.00
;  1.000000,  1.011789,  0.986396,  0.931878,  0.882166,  0.895944,  0.901659, $
;  0.899747,  0.896614,  0.882799,  1.000000,  1.028069,  1.002304,  0.939842, $
;  0.901619,  0.907495,  0.917215,  0.903099,  0.884204,  0.891070,  1.000000, $
;  1.028542,  0.997674,  0.933212,  0.894926,  0.904819,  0.913949,  0.890092, $
;  0.861423,  0.876450,  1.000000,  1.020855,  0.984892,  0.924484,  0.890529, $
;  0.897345,  0.898400,  0.878698,  0.852213,  0.856514,  1.000000,  1.017688, $
;  0.986462,  0.935272,  0.899000,  0.900801,  0.903275,  0.887802,  0.868248, $
;  0.860721,  1.000000,  1.016530,  0.988986,  0.945879,  0.906692,  0.906312, $
;  0.908931,  0.896586,  0.884092,  0.868331,  1.000000,  1.013696,  0.990371, $
;  0.953302,  0.907764,  0.910852,  0.913218,  0.902778,  0.891331,  0.872748, $
;  1.000000,  1.012690,  0.990044,  0.956787,  0.917377,  0.918191,  0.917497, $
;  0.907433,  0.898116,  0.879589,  1.000000,  1.011390,  0.993423,  0.964764, $
;  0.930912,  0.929004,  0.926571,  0.918262,  0.908627,  0.893646,  1.000000, $
;  1.008544,  0.992342,  0.968741,  0.935763,  0.935167,  0.933062,  0.923549, $
;  0.915229,  0.901549,  1.000000,  1.007472,  0.992683,  0.971941,  0.943764, $
;  0.939158,  0.935419,  0.928423,  0.920638,  0.906081,  1.000000,  1.002845, $
;  0.989209,  0.970376,  0.934759,  0.937781,  0.934360,  0.926872,  0.919243, $
;  0.905364,  1.000000,  0.998709,  0.987053,  0.968946,  0.938349,  0.938371, $
;  0.933682,  0.927399,  0.921172,  0.906905]		 ;18feb
  1.000000,  1.011796,  0.986602,  0.931094,  0.881353,  0.895115,  0.900832, $
  0.898924,  0.895803,  0.881994,  1.000000,  1.028074,  1.002862,  0.937646, $
  0.899290,  0.905152,  0.914842,  0.900770,  0.881924,  0.888776,  1.000000, $
  1.028546,  0.997928,  0.932219,  0.893876,  0.903756,  0.912875,  0.889048, $
  0.860412,  0.875425,  1.000000,  1.020857,  0.985365,  0.922612,  0.888538, $
  0.895336,  0.896385,  0.876728,  0.850303,  0.854602,  1.000000,  1.017691, $
  0.986775,  0.934035,  0.897684,  0.899479,  0.901949,  0.886502,  0.866980, $
  0.859462,  1.000000,  1.016532,  0.989294,  0.944659,  0.905396,  0.905016, $
  0.907630,  0.895304,  0.882828,  0.867090,  1.000000,  1.013697,  0.990474, $
  0.952898,  0.907344,  0.910427,  0.912793,  0.902358,  0.890915,  0.872344, $
  1.000000,  1.012692,  0.990347,  0.955588,  0.916101,  0.916913,  0.916222, $
  0.906172,  0.896869,  0.878369,  1.000000,  1.011392,  0.993549,  0.964275, $
  0.930388,  0.928483,  0.926052,  0.917748,  0.908119,  0.893147,  1.000000, $
  1.008545,  0.992362,  0.968672,  0.935690,  0.935095,  0.932990,  0.923478, $
  0.915160,  0.901483,  1.000000,  1.007474,  0.992680,  0.971966,  0.943793, $
  0.939188,  0.935450,  0.928454,  0.920669,  0.906112,  1.000000,  1.002846, $
  0.989174,  0.970526,  0.934921,  0.937944,  0.934524,  0.927035,  0.919403, $
  0.905525,  1.000000,  0.998711,  0.987123,  0.968674,  0.938057,  0.938080, $
  0.933395,  0.927114,  0.920889,  0.906629]	; 19apr-gwidth=11 below, new GAC

; 2019mar20 - begin the gwidth=11 for G230LB. see g750L example below.
; HGT=7 above. Do NOT change
;   1997.38   1999.00   2000.50   2002.50   2009.50   2011.20   2013.30   
; 2015.60   2020.00
; 1750. 1850. 1950. 2050. 2150. 2250. 2350. 2450. 2550. 2650. 2750. 2850. 2950.
tc230lb11=[								$ 
  1.000000,  1.010201,  0.984554,  0.931272,  0.880254,  0.894074,  0.898542, $
  0.898042,  0.893783,  0.879007,  1.000000,  1.028316,  1.002354,  0.938983, $
  0.899993,  0.905734,  0.913849,  0.901401,  0.882312,  0.883609,  1.000000, $
  1.027759,  0.996999,  0.932307,  0.895282,  0.903358,  0.911910,  0.888822, $
  0.860238,  0.870886,  1.000000,  1.020315,  0.984423,  0.923365,  0.888010, $
  0.895721,  0.896065,  0.876952,  0.850882,  0.850743,  1.000000,  1.016141, $
  0.985197,  0.933894,  0.897743,  0.899374,  0.901555,  0.886504,  0.867182, $
  0.856542,  1.000000,  1.014624,  0.987563,  0.944295,  0.905276,  0.904785, $
  0.907329,  0.895351,  0.883498,  0.864102,  1.000000,  1.011045,  0.988401, $
  0.951453,  0.907410,  0.909560,  0.912119,  0.901623,  0.890975,  0.870389, $
  1.000000,  1.009510,  0.987695,  0.954366,  0.915284,  0.916473,  0.916269, $
  0.905708,  0.897806,  0.877065,  1.000000,  1.009049,  0.991714,  0.963056, $
  0.929496,  0.928051,  0.926590,  0.917386,  0.909808,  0.890775,  1.000000, $
  1.006876,  0.990815,  0.967817,  0.935492,  0.934718,  0.933986,  0.923508, $
  0.917496,  0.899655,  1.000000,  1.005888,  0.991689,  0.970551,  0.943808, $
  0.938382,  0.936414,  0.928081,  0.922756,  0.903520,  1.000000,  1.002471, $
  0.988987,  0.970006,  0.937086,  0.938287,  0.936684,  0.927829,  0.922832, $
  0.903757,  1.000000,  0.999061,  0.987580,  0.968366,  0.940161,  0.939014, $
  0.936191,  0.928444,  0.924871,  0.906465]	; 2019may

; G430L time change array:
;   1997.38   1999.0   2009.50   2020.00
; 3000. 3200. 3400. 3600. 3800. 4000. 4200. 4400. 4600. 4800. 5000. 5200. 5400.
; 5600.
tc430l=[								$ 
; 13jul - new flats
; 2014aug18 - fix end WL bins.
;  1.000000,  1.007145,  0.966155,  0.977309,  0.955550,  1.000000,  1.002927, $
;  0.962389,  0.971147,  0.945976,  1.000000,  0.999258,  0.963078,  0.969435, $
;  0.945944,  1.000000,  0.996105,  0.966316,  0.969758,  0.947118,  1.000000, $
;  0.995857,  0.963318,  0.967430,  0.946867,  1.000000,  0.994707,  0.966432, $
;  0.968568,  0.948429,  1.000000,  0.993260,  0.964523,  0.966549,  0.946923, $
;  1.000000,  0.991668,  0.963592,  0.965935,  0.945611,  1.000000,  0.992788, $
;  0.963918,  0.965843,  0.944886,  1.000000,  0.991561,  0.964219,  0.965411, $
;  0.942834,  1.000000,  0.989544,  0.962687,  0.963915,  0.941730,  1.000000, $
;  0.990577,  0.963764,  0.963965,  0.942495,  1.000000,  0.988830,  0.962733, $
;  0.962064,  0.939790,  1.000000,  0.990247,  0.965698,  0.963686,  0.942636]
; 19jan
  1.000000,  1.011540,  0.968359,  0.979760,  0.957998,  1.000000,  1.007673, $
  0.964487,  0.973601,  0.948417,  1.000000,  1.003951,  0.965426,  0.971754, $
  0.948258,  1.000000,  1.000947,  0.968168,  0.972026,  0.949384,  1.000000, $
  1.000290,  0.965459,  0.969471,  0.948915,  1.000000,  0.998853,  0.968623, $
  0.970438,  0.950311,  1.000000,  0.997139,  0.966361,  0.968168,  0.948559, $
  1.000000,  0.995831,  0.965282,  0.967634,  0.947324,  1.000000,  0.996127, $
  0.964859,  0.967067,  0.946133,  1.000000,  0.995810,  0.965967,  0.967161, $
  0.944592,  1.000000,  0.992828,  0.964069,  0.965013,  0.942853,  1.000000, $
  0.993334,  0.964263,  0.964749,  0.943310,  1.000000,  0.991155,  0.962859, $
  0.962506,  0.940271,  1.000000,  0.992538,  0.965825,  0.964199,  0.943186]
; 19mar  - use gwidth=7 GAC

; 2019mar20 - begin the gwidth=11 for G430L. see g750L example below.
; HGT=7 above. Do NOT change
;   1997.38   1999.0   2009.50   2020.00
; 3000. 3200. 3400. 3600. 3800. 4000. 4200. 4400. 4600. 4800. 5000. 5200. 5400.
; 5600.

tc430l11=[								$ 
  1.000000,  1.005380,  0.960585,  0.975039,  0.951959,  1.000000,  1.002295, $
  0.952954,  0.967114,  0.941036,  1.000000,  0.998944,  0.956802,  0.966464, $
  0.941728,  1.000000,  0.996370,  0.959147,  0.966591,  0.942805,  1.000000, $
  0.996689,  0.957856,  0.965016,  0.943369,  1.000000,  0.995632,  0.961821, $
  0.966399,  0.944927,  1.000000,  0.994804,  0.960167,  0.964856,  0.943879, $
  1.000000,  0.993677,  0.962460,  0.965909,  0.944833,  1.000000,  0.993980, $
  0.961542,  0.965113,  0.944133,  1.000000,  0.993989,  0.963258,  0.965775, $
  0.944103,  1.000000,  0.991288,  0.962988,  0.964778,  0.943736,  1.000000, $
  0.992163,  0.962444,  0.964620,  0.945012,  1.000000,  0.990410,  0.962599, $
  0.963643,  0.944448,  1.000000,  0.992387,  0.964898,  0.965654,  0.948266]
; 2019may

; G750L time change array:
;   1997.38   2016.0
;   1997.38   1999.30   2009.50   2020.00
; 5700. 6100. 6500. 6900. 7300. 7700. 8100. 8500. 8900. 9300. 9700.
tc750l=[								$ 
;   1997.38   1999.30   2009.50   2015.00
  1.000000,  0.997970,  0.984762,  0.985602,  0.978012,  1.000000,  0.998179, $
  0.982551,  0.982932,  0.972493,  1.000000,  0.996686,  0.980973,  0.981388, $
  0.970433,  1.000000,  0.995702,  0.979386,  0.977641,  0.967323,  1.000000, $
  0.993745,  0.978995,  0.976174,  0.962808,  1.000000,  0.993273,  0.981307, $
  0.974183,  0.961352,  1.000000,  0.992655,  0.987198,  0.973595,  0.963048, $
  1.000000,  0.985828,  0.966214,  0.959103,  0.946583,  1.000000,  0.987191, $
  0.971510,  0.963247,  0.949702,  1.000000,  0.981083,  0.964509,  0.954902, $
  0.941289,  1.000000,  0.974496,  0.930921,  0.940577,  0.922219]	;14sep23
; #######!!!!!!!!!!!
; HGT=7 above. Do NOT change, but would need to temp. edit back to 2015 end time
;	to use.	So 2017jan27-corr for hz43 use, see farther below. 11 just below
;   1997.38   1999.30   2009.50   2018.00
; gwidth=11:
tc750l11=[								$ 
;  1.000000,  0.997072,  0.982275,  0.983376,  0.976385,  1.000000,  0.996913, $
;  0.979514,  0.980997,  0.971013,  1.000000,  0.995396,  0.977747,  0.979960, $
;  0.970677,  1.000000,  0.994755,  0.976468,  0.977090,  0.968281,  1.000000, $
;  0.993492,  0.975941,  0.976277,  0.963731,  1.000000,  0.994233,  0.978724, $
;  0.976039,  0.962966,  1.000000,  0.994998,  0.985130,  0.977610,  0.963791, $
;  1.000000,  0.989631,  0.964573,  0.964992,  0.948104,  1.000000,  0.992040, $
;  0.970246,  0.969813,  0.953241,  1.000000,  0.987383,  0.964740,  0.962798, $
;  0.946162,  1.000000,  0.980980,  0.933573,  0.948447,  0.927047]	; 17aug
;  1.000000,  0.997051,  0.982478,  0.983343,  0.976314,  1.000000,  0.996892, $
;  0.979717,  0.980823,  0.971022,  1.000000,  0.995375,  0.977950,  0.979935, $
;  0.970316,  1.000000,  0.994734,  0.976670,  0.976961,  0.968273,  1.000000, $
;  0.993471,  0.976143,  0.976184,  0.963175,  1.000000,  0.994212,  0.978926, $
;  0.976036,  0.962297,  1.000000,  0.994978,  0.985333,  0.978025,  0.961759, $
;  1.000000,  0.989611,  0.964774,  0.965282,  0.945953,  1.000000,  0.992020, $
;  0.970448,  0.970036,  0.950430,  1.000000,  0.987362,  0.964941,  0.963297, $
;  0.943355,  1.000000,  0.980959,  0.933771,  0.948687,  0.923947]	; 18jun
;  1.000000,  0.997051,  0.982478,  0.983519,  0.975849,  1.000000,  0.996892, $
;  0.979717,  0.980962,  0.970753,  1.000000,  0.995375,  0.977950,  0.980046, $
;  0.969979,  1.000000,  0.994734,  0.976670,  0.977089,  0.967891,  1.000000, $
;  0.993471,  0.976143,  0.976187,  0.963078,  1.000000,  0.994212,  0.978926, $
;  0.976014,  0.961954,  1.000000,  0.994978,  0.985333,  0.977964,  0.961884, $
;  1.000000,  0.989611,  0.964774,  0.965159,  0.946186,  1.000000,  0.992020, $
;  0.970448,  0.970236,  0.950610,  1.000000,  0.987362,  0.964941,  0.963287, $
;  0.942960,  1.000000,  0.980959,  0.933771,  0.948452,  0.924552]	; 18aug
  1.000000,  0.997068,  0.982270,  0.983507,  0.974492,  1.000000,  0.996916, $
  0.979517,  0.980961,  0.969118,  1.000000,  0.995398,  0.977751,  0.980079, $
  0.968245,  1.000000,  0.994759,  0.976467,  0.977128,  0.966242,  1.000000, $
  0.993483,  0.975934,  0.976093,  0.961360,  1.000000,  0.994237,  0.978736, $
  0.976052,  0.959799,  1.000000,  0.995000,  0.985118,  0.978008,  0.959461, $
  1.000000,  0.989631,  0.964577,  0.965038,  0.943995,  1.000000,  0.992046, $
  0.970261,  0.970119,  0.948342,  1.000000,  0.987384,  0.964737,  0.963552, $
  0.939473,  1.000000,  0.980948,  0.933573,  0.948345,  0.921772]	; 19may

opt_elem=strtrim(sxpar(hd,'opt_elem'),2)
case opt_elem of	
  'G140L' : mode = 'G140L'
  'G140M' : mode = 'G140L'
  'E140M' : mode = 'G140L'
  'E140H' : mode = 'G140L'
  'G230L' : mode = 'G230L'
  'G230M' : mode = 'G230L'
  'E230M' : mode = 'G230L'
  'E230H' : mode = 'G230L'
  'PRISM' : mode = 'G230L'
  'G230LB': Mode = 'G230LB'
  'G230MB': Mode = 'G230LB'
  'G430L' : Mode = 'G430L'
  'G430M' : Mode = 'G430L'
  'G750L' : Mode = 'G750L'
  'G750M' : Mode = 'G750L'
  else: return
  endcase

; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
;###change to new yr, (BUT not before getting the new fits from make-tchang??)
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
endpts=[1997.38,2020]
gwid=fix(sxpar(hd,'gwidth'))
if gwid ge 11 then tc230lb=tc230lb11		; 2019apr8 was eq
if gwid ge 11 then tc430l=tc430l11
if gwid ge 11 then tc750l=tc750l11

; 2019mar20-no-op?? Move down??
; ??? if gwid eq 7 and mode eq 'G750L' then endpts=[1997.38,2015]

;if mode eq 'G140L' then endpts=[endpts(0),1999.2,2002.7,endpts(1)]
if mode eq 'G140L' then endpts=[endpts(0),1999.2,2001.80,2009.50,endpts(1)]
;if mode eq 'G230L' then endpts=[endpts(0),1998.4,1999.3,2002.6,2009.5,	$
;			2011.2,endpts(1)]
if mode eq 'G230L' then endpts=[endpts(0),1998.4,1999.3,2002.6,2009.5,	$
			2011.2,2013.30,2015.60,endpts(1)]
;if mode eq 'G230LB' then endpts=[endpts(0),1999.0,2000.5,2002.5,2009.5,$
;			2011.2,endpts(1)]
if mode eq 'G230LB' then endpts=[endpts(0),1999.0,2000.5,2002.5,2009.5,	$
			2011.2,2013.30,2015.60,endpts(1)]
;if mode eq 'G430L' then endpts=[endpts(0),1999.30,2009.5,endpts(1)] ; 2014jul22
if mode eq 'G430L' then endpts=[endpts(0),1999.0,2009.5,endpts(1)] ; 2019mar25
if mode eq 'G750L' then endpts=[endpts(0),1999.30,2009.5,endpts(1)] ; 2014jul22

; ******** END SECTION FOR UPDATES (except temp coeff below) ***************

time=absdate(sxpar(hd,'pstrtime'))  &  time=time(0)	; change yyddd to yr.xxx
if time eq 0 then begin					;use MJD from EXPSTART
        mjd = sxpar(hd,'texpstrt')			; 02sep23
 	caldat,mjd+2400000.5d0,month,day,year,hour	; Julian to calendar
	time = year + (ymd2dn(year,month,day)+hour/24.0)/365.25d0
	endif
tcorr = 1.0				; default, applies to G230L
wcorr = 1.0

case mode of
; TEMPERATURE corrections from temp only fits. Time changes SHOULD be removed!
  'G140L' : begin
 	temp=sxpar(hd,'om1cat')
	if temp le 20 or temp ge 50 then begin			; idiot check
		print,'WARNING: NO om1cat for temp corr in TTCOR',	$
				format='(///8x,a///)'
	    end else begin
;		tcorr=1-(temp-36)*.0032		; 13jan  +- .00036 3 segments
;		tcorr=1-(temp-36)*.0033		; 13jun  +- .00037 3 segments
;		tcorr=1-(temp-36)*.0033		; 14feb  +- .00036 3 segments
;		tcorr=1-(temp-36)*.0033		; 14aug  +- .00037 3 segments
;		tcorr=1-(temp-36)*.0033		; 16dec  +- .00037 3 segments
;		tcorr=1-(temp-36)*.0035		; 17jul  +- .00038 3 segments
		tcorr=1-(temp-36)*.0036		; 19may  +- .00037 4 segments

		if n_elements(flux) gt 1 and not keyword_set(notemp) then begin
			flux=flux/tcorr
	    		if not keyword_set(silent) then print,mode,	$
							' temp. corr=',tcorr
			endif
		endelse
	endcase
  'G230LB': begin 
 	temp=sxpar(hd,'OCCDHTAV')
 	if temp gt 0 and time gt 2001.5 then begin
; 03may1 - ff are temp only fits and have ~0 change w/ the new flats and WLs
;		tcorr=1+(temp-19)*.0033	    ;11jun20 +- .00024 w/ scat lite corr
;		tcorr=1+(temp-19)*.0031			;11aug10 +- .00015
;		tcorr=1+(temp-19)*.0032			;13sep +- .00018
;		tcorr=1+(temp-19)*.0031			;14feb17 +- .00015
;		tcorr=1+(temp-19)*.0031			;14aug08 +- .00014
;		tcorr=1+(temp-19)*.0032			;17aug +- .00017
;		tcorr=1+(temp-19)*.0031			;18jun +- .00019
		tcorr=1+(temp-19)*.0032			;18jun +- .00017more seg
		tcorr=1+(temp-19)*.0033			;19jan +- .00016
		if gwid eq 11 then tcorr=1+(temp-19)*.0034 	;19may+/-.00014

		if n_elements(flux) gt 1 and not keyword_set(notemp) then begin
			flux=flux/tcorr
			if not keyword_set(silent) then print,mode,	$
							' temp. corr=',tcorr
			endif
		endif
	endcase
  'G430L' : begin 
 	temp=sxpar(hd,'OCCDHTAV')
 	if temp gt 0 and time gt 2001.5 then begin
;		tcorr=1+(temp-19)*.0020			;04jun5 +- .0005
;		tcorr=1+(temp-19)*.0021			;05may13 +- .00045
;		tcorr=1+(temp-19)*.0019			;05dec21  +- .00061
;		tcorr=1+(temp-19)*.0021			;09jul20  +- .00051
;		tcorr=1+(temp-19)*.0023		;09aug13 +/-.00034per fix of g=4
;		tcorr=1+(temp-19)*.0026			;09oct2 +/-.00032
;		tcorr=1+(temp-19)*.0027			;10may7 +/-.00023
;		tcorr=1+(temp-19)*.0023			;11aug11 +/-.00015 w/ E1
;		tcorr=1+(temp-19)*.0021			;12jan12+/-.14 w/3nodes
;		tcorr=1+(temp-19)*.0021			;13jan+/-.14 w/3nodes
;		tcorr=1+(temp-19)*.0021			;14feb+/-.00016 w/3nodes
;		tcorr=1+(temp-19)*.0022			;16dec+/-.00019 w/3nodes
;		tcorr=1+(temp-19)*.0022			;18jun+/-.00019 w/3nodes
		tcorr=1+(temp-19)*.0023			;19mar+/-.00019 w/3nodes
		if gwid eq 11 then tcorr=1+(temp-19)*.0022 	;19may+/-.00013

		if n_elements(flux) gt 1 and not keyword_set(notemp) then begin
			flux=flux/tcorr
			if not keyword_set(silent) then print,mode,	$
							' temp. corr=',tcorr
			endif
		endif
	endcase
  'G750L' : begin 
 	temp=sxpar(hd,'OCCDHTAV')
 	if temp gt 0 and time gt 2001.5 then begin
;		tcorr=1+(temp-19)*.0011			;09aug13  +- .00020
;		tcorr=1+(temp-19)*.0008			;09oct2  +- .00022
;		tcorr=1+(temp-19)*.0006			;10may7  +- .00017
;		tcorr=1+(temp-19)*.0007			;11aug11  +- .00017
;		tcorr=1+(temp-19)*.0006			;12sep  +- .00016
;		tcorr=1+(temp-19)*.0005			;14feb  +- .00023
		tcorr=1+(temp-19)*.0006			;14jul  +- .00022
;		if gwid eq 11 then tcorr=1+(temp-19)*.0012 	;14sep+/-.00015
;		if gwid eq 11 then tcorr=1+(temp-19)*.0013 	;16sep+/-.00014
;		if gwid eq 11 then tcorr=1+(temp-19)*.0012 	;16dec+/-.00015
;		if gwid eq 11 then tcorr=1+(temp-19)*.0013 	;17aug+/-.00015
;		if gwid eq 11 then tcorr=1+(temp-19)*.0012 	;18jun+/-.00015
		if gwid eq 11 then tcorr=1+(temp-19)*.0013 	;19may+/-.00015
		if n_elements(flux) gt 1 and not keyword_set(notemp) then begin
			flux=flux/tcorr
			if not keyword_set(silent) then print,mode,	$
							' temp. corr=',tcorr
			endif
		endif
	endcase
  else:
  endcase
if keyword_set(notime) then return

; TIME correction

case mode of	; reform: eg reformat 5 endpts*11wl-bins to 5x11 array
  'G140L' : begin & timcor=reform(tc140l,7,11) & wch=indgen(11)*50.+1175 & end
  'G230L' : begin & timcor=reform(tc230l,10,15) & wch=indgen(15)*100.+1650 & end
  'G230LB': begin & timcor=reform(tc230lb,10,13) & wch=indgen(13)*100.+1750 & end
  'G430L' : begin & timcor=reform(tc430l,5,14) & wch=indgen(14)*200.+3000 & end
  'G750L' : begin & timcor=reform(tc750l,5,11) & wch=indgen(11)*400.+5700 & end
  endcase

if time le 1997 or time ge max(endpts) then begin	; idiot check
	print,'ttcorr: Invalid Observation Date in Header = ',time
	stop
	endif
corr=timefit(endpts,timcor,time,mode=mode)	; get corr in ea WL bin @ t=time
;orig:linterp,wch,corr,wav,wcorr;05may truncates to endpt val. Best 4 AGK @ 1mic
wcorr=interpol(corr,wch,wav)		; 14Jul24 - extrap beyond wch ends
; 2014Jul24-As the corr can be a bit lumpy, 5 node spline fit makes worse lumps!
; 2018jan - ie the lumps look real. Leave as is.
;	extrap w/ interpol makes GD71 resid shoot up after 9900A, so extrap only
;	to 9900= last node @ 9700A+dlam/2:

del=wch(1)-wch(0)
bad=where(wav lt wch(0)-del/2)  &  wcorr(bad)=wcorr(bad(-1))
bad=where(wav gt wch(-1)+del/2)  &  wcorr(bad)=wcorr(bad(0))
if n_elements(flux) gt 1 then flux=flux/wcorr
if not keyword_set(silent) then $
    if min(wcorr) lt 1 or max(wcorr) gt 1 then 				$
		print,version,mode,' Time corr. made to NET or FLUX @ ',time
return
END