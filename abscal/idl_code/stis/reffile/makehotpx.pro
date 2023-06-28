; 09sep28 - RCB routine to make hot pixel tables.
; First edit this program to get the list of new drk files from $oref and
;	edit results into stisidl/scal/calstis_dark.txt
; Doc is at stisdoc/hot.pixel. eg: run findbias.pro
;-

; ###change - enter the input new dark file names:
darks=[	'taj2055ko','taj2055lo','taj20569o','taj2056no','tau1346eo',	$
	'tau1346go','tau1346io']	; 09nov12.
; MAMA NUV dark:  'tcl1742bo'
darks=[ 'tbn1517io','tce1631go','tce1631io','tc420049o','tc42004bo',	$
	'tcg1733ko','tcg1733mo','u211310eo','u211310go',		$
	'u221550fo','u221550ho','u221550jo','u2g1858co','u2g1858eo',	$
	'u321936bo','u321936do','u3u1152oo','u3u1152qo','u3u1152so',	$
	'u3u11530o','u4j1725so','u4j1725to','u4j17263o','u4j17266o',	$
	'u4r2025go','u4r2025io','u4r2025ko','u4r2025mo']	; 2010May3
; 2010jun9
darks=['u611612ho', 'u611612io', 'u611612jo', 'u611612ko', 'u681759co',	      $
'u681759do', 'u681759eo', 'u681759fo', 'u681759go', 'u681759ho', 'u681759io', $
'u681759jo', 'u681759ko', 'u681759lo', 'u681759mo', 'u681759no', 'u681759oo', $
'u681759po', 'u681759qo', 'u681759ro', 'u681759so', 'u681759to', 'u6818000o', $
'u6818001o', 'u6818002o', 'u6818003o', 'u6818004o', 'u6818005o', 'u6818006o', $
'u6818007o', 'u6818008o', 'u6818009o', 'u681800ao', 'u681800bo', 'u681800co', $
'u681800do', 'u681800eo', 'u681800fo', 'u681800go', 'u681800ho', 'u681800io', $
'u681800jo', 'u681800ko', 'u681800lo', 'u681800mo', 'u681800no', 'u681800oo', $
'u681800po', 'u681800qo', 'u681800ro', 'u681800so', 'u681800to', 'u6818010o', $
'u6818011o', 'u6818012o', 'u6818013o', 'u6818014o', 'u6818015o', 'u6818016o', $
'u6818017o', 'u6818018o', 'u6818019o', 'u681801ao', 'u681801bo', 'u681801co', $
'u681801do', 'u681801eo', 'u681801fo', 'u681801go', 'u681801ho', 'u681801io', $
'u681801jo', 'u681801ko', 'u681801lo', 'u681801mo', 'u691334fo', 'u691334go', $
'u691334ho', 'u691334io', 'u691334jo', 'u691334ko', 'u691334lo', 'u691334mo', $
'u691334no', 'u691334oo', 'u691334po', 'u691334qo', 'u691334ro', 'u691334so', $
'u691334to', 'u6913350o', 'u6913351o', 'u6913352o', 'u6913353o', 'u6913354o', $
'u6913355o', 'u6913356o', 'u6913357o', 'u6913358o', 'u6913359o', 'u691335ao', $
'u691335bo', 'u691335co', 'u691335do', 'u691335eo', 'u691335fo', 'u691335go', $
'u691335ho', 'u691335io', 'u691335jo', 'u691335ko', 'u691335lo', 'u691335mo', $
'u691335no', 'u691335oo', 'u691335po', 'u691335qo', 'u691335ro', 'u691335so', $
'u691335to', 'u6913360o', 'u6913361o', 'u6913362o', 'u6913363o', 'u6913364o', $
'u6913365o', 'u6913366o', 'u6913367o', 'u6913368o', 'u6913369o', 'u691336ao', $
'u691336bo', 'u691336co', 'u691336do', 'u691336eo', 'u691336fo', 'u691336go']
darks=['u691805po','u691805qo','u691805ro','u691805so','u691805to',	      $
'u6918060o','u6918061o','u6918062o','u6918063o','u6918064o','u6918065o',      $
'u6918066o','u6918067o','u6918068o','u6918069o','u691806ao','u691806bo',      $
'u691806co','u691806do','u691806eo','u691806fo','u6o1317po','u6o1317qo',      $
'u6o1317ro','u6o1317so']				; 2010jul6
; 2011jan25
darks=[									$
'u7u1516no','u7u1516oo','u7u1516po','u7u1516qo','u7u1516ro','u8r1931oo',  $
'u8r1931po','u8r1931qo','u8r1931ro','u9m1955to','u9m19560o','u9m19561o',  $
'uas1836ro','uas1836so','uas1836to','uas18370o','ubj15315o','ubj15316o',  $
'ubj15317o','ubj15318o','ucm1604qo','ucm1604ro','ucm1604so','ucm1604to',  $
'v1514376o','v1514377o','v1514378o','v1514379o']
darks=[									$
'v2320397o','v2320398o','v2320399o','v232039ao','v372106so','v372106to',	$
'v3721070o','v3721071o']
darks=[					$			; 2011Aug4
'v471332ro','v471332so','v471332to','v4713330o','v4c0057no',              $
'v4c0057oo','v4c0057po','v4c0057qo','v5b1216ao','v5b1216bo','v5b1216co',  $
'v5b1216do','v661427bo','v661427co','v661427do','v661427eo','v7i1156po',  $
'v7i1156qo','v7i1156ro']
darks=[			$					; 2011Aug10
'v891810ao','v891810bo','v891810co','v891810do']
darks=[			$					; 2011Oct21
'v911541ro','v911541qo','v911541po','v911541so','v9j1616ro','v9j1616so',  $
'v9j1616to','v9j16170o']					; 2011Oct21
darks=[			$					; 2011Oct21
'vap1253ho','vap1253io','vap1253jo','vap1253ko','vbf18120o','vbf18121o',  $
'vbf18122o','vcd1735co','vcd1735do','vcd1735eo','vcd1735fo']	; 2012Jan2
darks=[			$					; 2012sep6
'vcd1735fo','w1j1436io','w1j1436jo','w1j1436ko','w1j1436lo','w2g1626ho',  $
'w2g1626io','w2g1626jo','w2g1626ko','w3k1237po','w3k1237qo','w3k1237ro',  $
'w3k1237so','w3n1413po','w3n1413qo','w3n1413ro','w3n1413so','w521315io',  $
'w521315jo','w521315ko','w521315lo','w6k1722lo','w6k1722mo','w6k1722no',  $
'w6k1722oo','w7r1710io','w7r1710jo','w7r1710ko','w7r1710ro','w7r1710so',  $
'w7r1710to','w7r17110o','w8n1954eo','w8n1954fo','w8n1954go','w8n1954ho']
darks=[			$					; 2013Jan7
'wb91745eo','wb91745fo','wb91745go','wbe1830ho','wbe1830io','wbe1830jo',  $
'wbe1830ko','wbk1745fo','wbk1745go','wbk1745ho','x1218198o','x1218199o',  $
'x121819ao','x121819bo','x1314591o','x1314592o','x1314593o']
darks=[			$					; 2013Jun6
'x4b14103o','x4b14104o','x4b14105o','x4b14106o','x4b14107o','x4b14108o',  $
'x4b14109o','x4t15459o','x4t1545ao','x4t1545bo','x4t1545co']
darks=[			$					; 2013Jul26
'x6e1332to','x6e13330o','x6e13331o','x6i1153bo','x6i1153co','x6i1153do',  $
'x6i1153eo']
darks=['x8c1626co','x8c1626do','x8c1626eo','x8c1626fo']		; 2013Sep9
darks=['x9c18021o','x9c18022o','x9c18023o','x9j1854ko','x9j1854lo',	$
			'x9j1854mo','x9j1854no']		; 2013Sep20
darks=[								$ 2014feb10
'xcg2017ao','xcg2017co','xcg2017eo','xcg2017go','xcg2017ko','xcg2017mo',  $
'xcg2017oo','xcg2017so','xcg20180o','xcg20182o','xcg20184o','y2519101o',  $
'y2519102o','y2519103o','y2519104o','y261955go','y261955ho','y261955io',  $
'y2716123o','y2716124o','y2716125o','y2716126o']
darks=[								$ 2014jul28
'y5k1338oo','y5k1338po','y5k1338qo','y5k1338ro','y6316570o','y6316571o',  $
'y6316572o','y6316573o','y6913480o','y6913481o','y6913482o','y6g2048ho',  $
'y6g2048io','y6g2048jo','y6g2048ko']
darks=[								$ 2014aug6
'y851541ho','y851541io','y851541jo','y851541ko']
darks=[								$ 2014sep19
'y9i1617lo]','y9i1617mo','y9i1617no','y9i1617oo','y9i1617po','y9i1617qo']
darks=[								$ 2014nov25
'yal1700oo','yal1700po','yal1700qo','yal1700ro','yal1700so',	$
'yal1700to','yal17010o','yb41544po','yb41544qo','yb41544ro','yb41544so']
darks=[								$ 2015may3
'z3b1918oo','z3b1918po','z3b1918qo','z3b1918ro','z3b1918so','z3b1918to', $
'z3b19190o']
darks=[								$ 2015jun1
'z5c20075o','z5c20076o','z5c20077o','z5c20078o','z5c20079o','z5c2007ao', $
'z5c2007bo','z5j1735so','z5j1735to','z5j17360o','z5j17361o','z5j17362o', $
'z5j17363o','z5j17364o','z5j17365o','z5j17366o']
darks=[								$ 2015oct12
'z7t1652no','z7t1652oo','z7t1652po','z7t1652qo','z7t1652ro','z7t1652so', $
'z7t1652to','z7t16530o','z7t16531o','z7t16532o','z8o1932do','z8o1932eo', $
'z8o1932fo','z8o1932go']
darks=[								$ 2016feb5
'zbi1324bo','zbi1324co','zbi1324do','zbi1324eo','zbi1324fo','zbi1324go']
darks=[								$ 2016feb11
'02a19206o','02a19207o','02a19208o','02a19209o','02a1920ao','02a1920bo', $
'02a1920co','02a1920do','02a1920eo']
darks=[								$ 2016jun7
'0521335no','0521335po','0521335to','05213362o','05213363o','05213364o', $
'05213365o','05314034o','05314039o','0531403ao','0531403bo','0531403co', $
'0531403go']
darks=[								$ ;2016sep20
'08o14511o','08o14512o','08o14513o','08o14514o','08o14515o','08o14516o',  $
'08o14517o','08o14518o','08o14519o','08o1451ao','08o1451bo','08o1451co',  $
'08o1451do','09k14504o','09k14505o','09k14506o','09k14507o']
darks=[								$ ;2016nov10
'0ac1554io','0ac1554jo','0ac1554ko','0ac1554lo','0ba1437ko','0ba1437lo',  $
'0ba1437mo']
darks=[								$ ;2016dec1
'0bs1534ho','0bs1534io','0bs1534jo']
darks=[								$ ;2016dec29
'0cl1918fo','0cl1918ho','0cl1918jo']
darks=[								$ ;2017mar16
'11r2137eo','11r2137go','11r2137io','11r2137jo','1371350ao','1371350co',  $
'1371350eo']
darks=[								$ ;2017apr22
'13o1804no','13o1804po','13o1804ro','14j2024qo','14j2024ro','14j2024to', $
'14j20251o'] 
darks=[								$ ;2017jun22
'15h1342so','15h13430o','15h13431o','16j1349ko','16j1349lo','16j1349no']
darks=[								$ ;2017jul14
'17e1456oo','17e1456qo','17e1456ro','17e1456to']
darks=[								$ ;2018jun1
'17e14590o','17e1459to','17e15003o','17e15008o','1871743no','1871743oo', $
'1871743qo','1a314267o','1a314269o','1a31426go','1a314273o','1a314275o', $
'1a31427eo','1a31427po','1a31427ro','1a31427so','1a314280o','1a31428mo', $
'1a31428oo','1a31428to','1a31428go','1a31428ho','1a31428jo','1a31428lo', $
'1a31428po','1a31428ro','1a314290o','1a314292o','1a314294o','1a314295o', $
'1a314297o','1a314266o','1aj2103jo','1aj2103no','1aj2103lo','1aj2102so', $
'1aj2103ao','1aj2103co','1aj2102po','1b717088o','1b71708co','1b717087o', $
'1bt1510io','1bt1510go','1bt1510jo','2191719go','2191719oo','2191719mo', $
'2191719po','2271950ro','2271950no','2271950so','22j2042lo','22j2042mo', $
'22j2042jo']
darks=[								$ ;2018jul20
'26i12005o','26i12003o','26i1159co','26i1200fo','26i11595o','26i1200io', $
'26i1200ko','26i12000o','26i1159qo','26i1159so','26i12001o','26i11599o', $
'26i11592o','27j14002o','27j1359no','27j1359oo','27j1359ro']
darks=[								$ ;2018aug11
'2831116to','28311172o','28311174o']
darks=[								$ ;2018dec27
'29b1320jo','29b1320fo','29b1320ko','29b1320no','2ab1300fo','2ab13007o', $
'2ab1300ho']
darks=[								$ ;2019jan7
'31317115o','3131711bo','31317113o','31316051o','31316053o','31316056o']
darks=[								$ ;2019may1
'31p2135go','31p2135io','31p2135po','31p2135no','32r2042co',	$
'32r20428o','32r2042go','33l1706oo','33l1706ro','33l1706ho','33l1706lo']
darks=[								$ ;2020mar19
'43i20166o','43i20168o','43i2016ko','43i2016qo','43i20175o','43i20178o', $
'43i2017ao','43i2017fo','43i2017ho','43i2017po','43i20180o','43i20189o', $
'43i2018no','43i2018ro','43i2018to','43i20192o','43i20199o','43i2019bo', $
'43i2019mo','43i2020eo','43i2020lo','43i2020so','43i20210o','43i20213o', $
'43i20218o','43i2021co','43i2021ko','43i2021mo','43i20222o','43i20225o', $
'43i20229o','43i2022co','43i2022so','43i20231o','43i20233o','43i2023fo', $
'43i2023ko','43i2023oo','43i2023ro','43i20244o','43i20246o','43i20249o', $
'43i2024bo','43i2024eo','43i2024jo','43i2024lo','43i2024oo','43i2024qo', $
'43i2024to','43i2025bo','43i2025mo','43i20263o','43i20268o']
darks=[								$ ;2020mar19
'4491915mo','44919161o','44919163o','44919167o','4491916co','4491916io', $
'4491916no','4491916po']
darks=[								$ ;2020nov10
'48d2029mo','48d2030ao','48d2031bo','48d2031go','48d2031lo','48d2032fo', $
'48d2033oo','48d20344o','48d2034fo','48d2034mo','48d2034to','48d20357o', $
'48d2035ao','48d20372o','48d2037co','48d2037mo','48d20380o','48d20389o']
darks=[								$ ;2020nov13
'4bc1439ro','4bc1439so','4bc14401o','4bc1440do','4bc1440io','4bc1440lo', $
'4bc1440mo','4bc1440ro','4bc14412o','4bc14415o','4bc14417o','4bc14419o']
darks=[								$ ;2021feb3
'4c71826co','4c71826io','4c71826ko','4c71826no','4c71826ro','5162002do', $
'5162002fo','5162002io','5162002ko','51s2043po','51s2043so','51s20443o', $
'51s20448o']
darks=[								$ ;2021mar19
'52o1652bo','52o1652do','52o1652eo','52o1652jo','53i1606ao','53i1606bo',  $
'53i1606do','53i1606jo']
darks=[								$ ;2021aug9
'5491743ko','54917441o','5491744co','5491744go','55a20445o','55a20447o',  $
'55a2044fo','55a2044jo']
darks=[								$ ;2021aug20
'58c1857eo','58c1857po]','58c1857so','58c18582o','58c18587o','58c1858io',  $
'58c1858ko','58c1858ro']
darks=['58p16347o','58p1634bo','58p1634do','58p1634go']		  ;2021sep24
darks=[									$
'59o1636oo','59o1636qo','59o16370o','59o16372o','59o16374o','5at1826oo',   $
'5at1826ro','5at18271o','5at18276o','5at18279o','5c717329o']	;2022jan25
darks=[									$
'62i21226o','62i2122fo','62i2122ho','62i2122jo','62i2122mo','62i2122oo',$
'62i2122ro','62i21236o','62i2123bo']				; 2022feb22
darks=['63f1625jo','63f1625lo','63f1625qo','63f16263o','63f16266o']  ; 2022mar29
darks=['64c1930io','64c1930no','64c1930so','64c1930to']		; 2022apr18
darks=['65o18246o','65o1824co','65o1824do','65o1824fo']		; 2022may25
darks=['6691532co','6691532eo','6691532go','6691532no','6691532ro']  ; 2022jun11
darks=[									$
'6781450ho','6781450lo','6781450oo','6781450qo','68a1901eo','68a1901jo',$
'68a1901no','68a1901po']
darks=[									$
'6ae1520po','6ae1520so','6ae15216o','6ae15218o','6ae15219o','6ae1521eo',$
'6ae1521go','6ae1521ho']
darks=[									$
'6bs2140do','6bs2140go','6bs2140io','6bs2140jo','6bs21417o','6bs2141ao',$
'6bs2141co','6bs2141do','6bs2141ho','6bs2141io','6bs2141ko','6bs2141no',$
'6bs2141to']
darks=['6cm19370o','6cm19376o','6cm19379o','6cm1937co']

;  ck mjd order of results
	  
outtabl=darks
; fails on laptop, works on desktop???? darkfils='$oref/'+darks+'_drk.fits'
darkfils='/grp/crds/hst/references/hst/'+darks+'_drk.fits'

; ###change - special case fix from drkfix.pro
;darks=['drkfix']					  ;17may1 
;darkfils=['../data/spec/dark/2017Jan12_drk.fits']

for i=0,n_elements(darks)-1 do begin
	fits_read,darkfils(i),dum,hd
	date=sxpar(hd,'useafter')
	mo=strlowcase(strmid(gettok(date,' '),0,3))
	day=gettok(date,' ')
	yr=gettok(date,' ')
	hr=gettok(date,':')
	outfil='hpx'+darks(i)+'_'+day+mo+yr+'.fits'
;	if hr ne 0 then stop	see hot.pixel email of 2010Jun10
	juldat=jul_date(day+'-'+mo+'-'+yr+' '+hr)
	print,'juldate=',juldat
	thresh=0.05
	hot_pixel_table,darkfils(i),'../stisidl/scal/'+outfil,threshold=thresh
	if thresh eq 999 then outtabl(i)=				$
	     darks(i)+' BAD cdbs file '+string(juldat,'(f10.2)') else	$
	     outtabl(i)='     CCD             '+outfil+string(juldat,'(f10.2)')
	endfor
print,outtabl,form='(a)'
end
