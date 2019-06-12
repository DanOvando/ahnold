
data fish;
set final.ucsbfish;
where level ne 'CAN'  and campus='UCSB' and method in ('SBTL_FISH' 'SBTL_FISH_NPS' 'SBTL_FISH_CRANE' 'SBTL_FISH_VRG') ;  
*if site = 'SCI_YELLOWBANKS' and side='E' then side='CEN';
if site = 'SCI_PELICAN' and side='FAR WEST' then delete;
/*if fish_tl lt 10 then do; count=0; fish_tl=.; end;*/

*size at 1st maturity - females;
mature=0;
if classcode='SATR' and fish_tl ge 15.8	then mature=1; *Romero 1988;
if classcode='SCHR' and fish_tl ge 24.3 then mature=1; *Zaitlin 1986;
if classcode='SCAU' and fish_tl ge 29.5 then mature=1; *Lea et al 1999;
if classcode='SMIN' and fish_tl ge 36.5 then mature=1; *Lea et al 1999;
if classcode='OYT' and fish_tl ge 28.5 then mature=1; *Lea et al 1999;
if classcode='SMAR' and fish_tl ge 44.5 then mature=1;* Fitch and Lavenberg 1971;
if classcode='HDEC' and fish_tl ge 31.6 then mature=1; *Rothrock 1973;
if classcode='SCAR' and fish_tl ge 20.7 then mature=1; *Lea et al 1999;
if classcode='SMEL' and fish_tl ge 30 then mature=1; *Wyllie Echeverria 1987;
if classcode='SMYS' and fish_tl ge 19.6 then mature=1; *Lea et al 1999;
if classcode='SPUL' and fish_tl ge  20/.831 then mature=1; *JENN;
if classcode='PCLA' and fish_tl ge  19 then mature=1; *Lea et al 1999;
if classcode='SDAL' and fish_tl ge  9 then mature=1; *VRG;
if classcode='SSAX' and fish_tl ge  10 then mature=1; *VRG;
if classcode='SSEM' and fish_tl ge  11 then mature=1; *VRG;
if classcode='OPIC' and fish_tl ge  14 then mature=1; *VRG;
if classcode='SHOP' and fish_tl ge  14 then mature=1; *VRG;
if classcode='SROS' and fish_tl ge  15 then mature=1; *VRG;
if classcode='SGUT' and fish_tl ge  18 then mature=1; *VRG;
if classcode='EJAC' and fish_tl ge  19 then mature=1; *VRG;
if classcode='HRUB' and fish_tl ge  21.7 then mature=1; *VRG;
if classcode='EJAC' and fish_tl ge  19 then mature=1; *VRG;
if classcode='SAUR' and fish_tl ge  26 then mature=1; *VRG;
if classcode='RVAC' and fish_tl ge  29 then mature=1; *VRG;
if classcode='SRUB' and fish_tl ge  34 then mature=1; *VRG;
if classcode='SENT' and fish_tl ge  35 then mature=1; *VRG;
if classcode='SSER' and fish_tl ge  35 then mature=1; *VRG;
if classcode='SPAU' and fish_tl ge  36 then mature=1; *VRG;
if classcode='SPIN' and fish_tl ge  44 then mature=1; *VRG;
if classcode='OELO' and fish_tl ge  55.7 then mature=1; *VRG;

*Length specific fecundity - length (when ge size at first maturity) sometimes converted to SL first;
if mature then do;
if classcode='SATR' and mature then LSF=0.000021*(fish_tl*10/1.233)**4.123; *Romero 1988,Lea et al 1999;
if classcode='SCHR' and mature then LSF=0.0000000136*(0.822*(fish_tl*10)-0.123)**5.59; *Zaitlin 1986, Lea et al 1999 
					- coeff. is 1.36e-5 in Gregors DB;
if classcode='SCAU' and mature then LSF=0.0000000026*(fish_tl*10)**5.347; *DeLacy et al 1964;
if classcode='SMIN' and mature then LSF=0.000002*(fish_tl)**5.0226; *Love et al 1990;
if classcode='OYT' and mature then LSF=0.006*(fish_tl)**4.619; *Love & Westphal 1981;
if classcode='SMAR' and mature then LSF= 0.29*fish_tl*10 - 87.54; *Lauth 1988,(range: 500 - 775 mm TL);
if classcode='SPUL' and mature then LSF=(0.00131*(0.831*fish_tl)**2.95)*5377; *Caselle & Hamilton;
if classcode='PCLA' and mature then  LSF=10**(3.02*(log10(fish_tl*10)) - 3.13); *Caselle & Hamilton;
end;
*if classcode='SMEL' and fish_tl ge 30
	then LSF=64098*(fish_tl)-2083531; *Wallace & Tagart 1994, Wyllie Echeverria 1987;

*if classcode='SPAU' and fish_tl ge ??? 
	then LSF=0.00116*(fish_tl)**3.2696; *Love et al 1990;

*length-weight calculation - updated 4_7_14;
if classcode='AFLA' then weight= 0.0004*(fish_tl)**3.4302; *Fishbase;
if classcode='AAFF' then weight= 0.00000886*(fish_tl*10)**3.03574; *Lea et al 1999;
if classcode='ACAL' then weight= 0.0079*(fish_tl)**3.03748; *Fishbase;
if classcode='ADAV' then weight= 0.0168*(fish_tl)**2.9856; *Fishbase;
if classcode='AHOL' then weight= 0.0180*(fish_tl)**2.7566; *Fishbase;
if classcode='AOCE' then weight= 0.0038*(fish_tl)**3.1581; *Fishbase;
if classcode='ASAN' then weight= 0.0062*(fish_tl)**3.249; *Fishbase;
if classcode='ATHE' then weight= 0.0004*(fish_tl)**2.519; *Lea et al 1999;
if classcode='AVUL' then weight= 0.0188*(fish_tl*0.526)**3.4302; *Fishbase, converted TL to FL;
if classcode='BAITBALL' then weight= 0.0076*(fish_tl)**3.150; *Fishbase for Sardines;
if classcode='BFRE' then weight= 0.0250*(fish_tl)**3.0; *Fishbase for family;
if classcode='BOTH' then weight= 0.000038442*(fish_tl*10)**3.25755; *Lea et al 1999;
if classcode='BPOL' then weight= 0.0268*(fish_tl)**2.8983; *Fishbase;
if classcode='BRAY' then weight= 0.0089*(fish_tl)**3.1133; *Fishbase;
if classcode='CAGG' then weight= 0.00407*(fish_tl)**3.45; *Lea et al 1999, TL in cm;
if classcode='CFAL' then weight= 0.0256*(fish_tl)**3.006; *Fishbase;
if classcode='CITH' then weight= 0.000003934*(fish_tl/1.12*10)**3.25755; *Chamberlin 1979, W=g, SL=mm;
if classcode='CLUP' then weight= 0.0076*(fish_tl)**3.150; *Fishbase for Sardines;
if classcode='CLIN' then weight= 0.00310456*(fish_tl)**3.243; *Steipen 1986, for HROS, assume TL in cm;
if classcode='CNIC' or classcode='RNIC' then weight= 0.01*(fish_tl)**3.0247; *Fishbase;
if classcode='COTT' then weight= 0.00514*(fish_tl)**3.31; *Bayer 1985, for Staghorn Sculpin;
if classcode='CPRI' then weight= 0.00000304*(fish_tl*10)**3.22; *Unknown ref;
if classcode='CPUG' then weight= 0.00514*(fish_tl)**3.31; *Bayer 1985, for Staghorn Sculpin;
if classcode='CPUN' then weight= 0.0000185*(fish_tl*10/1.213)**3.09; *Quast 1968b, SL = mm, L-L from fishbase for Garibaldi;
if classcode='CSAT' then weight= 0.0109*(fish_tl)**3.024; *Lea et al 1999, TL in cm;
if classcode='CSOR' then weight= 0.000003934*(fish_tl/1.12*10)**3.25755; *Chamberlin 1979, W=g, SL=mm;
if classcode='CSTI' then weight= 0.0000019188*(fish_tl/1.12*10)**3.44488; *Lea et al 1999, W=g, SL=mm;
if classcode='CVEN' then weight= 0.0027*(fish_tl)**3.0477; *Fishbase, TL in cm;
if classcode='DVAC' or classcode='RVAC' then weight= 0.208462*(fish_tl/2.54)**3.039; *DeMartini et al 1994, TL inches, converted to cm, corrected value from error in LH DB;
if classcode='EJAC' then weight= 0.230156*(fish_tl/2.54)**3.137; *DeMartini et al 1994, TL inches, converted to cm, corrected value from error in LH DB;
if classcode='ELAT' then weight= 0.000021*(fish_tl*0.925*10)**3.007; *Gnose 1967, FL in mm;
if classcode='EMBI' then weight= 0.025*(fish_tl)**3.0; *fishbase;
if classcode='EMOR' then weight= 0.017*(fish_tl*0.853)**2.95; *fishbase, SL = cm;
if classcode='GBY' then  weight= 0.00001299*(fish_tl*10)**3.077; *used SCAR, Lea et al 1999;
if classcode='GGAL' then weight= 0.0008*(fish_tl)**3.350; *fishbase;
if classcode='GIBB' then weight= 0.0095*(fish_tl)**3.0829; *fishbase;
if classcode='GMOR' then weight= 0.00000011*(fish_tl*10*0.992)**3.43; *Quast 1968, SL = mm;
if classcode='GNIG' then weight= 0.0000528*(fish_tl*10*0.851)**2.93; *Quast 1968, SL = mm;
if classcode='HANA' then weight= 0.000032*(fish_tl*10/1.179)**2.99467; *Anderson 1969, SL = mm, eq. for HARG;
if classcode='HARG' then weight= 0.000032*(fish_tl*10/1.179)**2.99467; *Anderson 1969, SL = mm;
if classcode='HCAR' then weight= 0.000002*(fish_tl*10*0.818)**3.56; *Quast 1968, SL = mm;
if classcode='HDEC' then weight= 0.00000418*(fish_tl*10)**3.198; *Rothrock 1973;
if classcode='HELL' then weight= 0.00000783*(fish_tl*10)**3.131; *Wydoski and Bennett 1973, for males;
if classcode='HFRA' then weight= 0.0027*(fish_tl)**3.0477; *Fishbase for swell shark, no info on this spp or genus;
if classcode='HLAG' then weight= 0.00000418*(fish_tl*10)**3.198; *Rothrock 1973, for kelp greenling HDEC;
if classcode='HROS' then weight= 0.00310456*(fish_tl)**3.243; *Steipen 1986, assume TL = cm;
if classcode='HRUB' then weight= 0.383276*(fish_tl/2.54)**3.077; *DeMartini et al 1994, TL inches, converted to cm, corrected value from error in LH DB;
if classcode='HSEM' then weight= 0.0848392*(fish_tl/2.54)**3.412; * DeMartini et al 1994, TL inches, converted to cm, corrected value from error in LH DB;
if classcode='HSUP' then weight= 0.00000418*(fish_tl*10)**3.198; *Rothrock 1973, for kelp greenling HDEC;
if classcode='KGB' then weight= 0.00000629*(fish_tl*10)**3.172; *Lea et al 1999 for SATR;
if classcode='LCON' then weight= 0.01*(fish_tl)**3.0247; *fishbase for Family values;
if classcode='LDAL' then weight= 0.0120*(fish_tl)**2.744; *fishbase for Genus values;
if classcode='LHIR' then weight= 0.0122*(fish_tl)**3.0436; *fishbase for Family values;
if classcode='LMUC' then weight= 0.0302*(fish_tl)**3.0; *fishbase for Genus values;
if classcode='MCAL' then weight= 0.0000083*(fish_tl*10*0.92)**3.26; *Quast 1968, SL = mm, L-L from fishbase;
if classcode='MCEP' then weight= 0.0109*(fish_tl)**3.089; *fishbase;
if classcode='NSTE' then weight= 0.0077*(fish_tl)**2.9620; *fishbase for Family;
if classcode='NUNI' then weight= 0.0077*(fish_tl)**2.9620; *fishbase for Family;
if classcode='OCAL' then weight= 0.00000128*(fish_tl*10*0.9)**3.5; *Demartini et al 1994, SL = mm, L-L guess no info in fishbase;
if classcode='OELO' then weight= 0.0000007125*(fish_tl*10)**3.405; *Lea et al 1999;
if classcode='OPIC' then weight= 0.0127*(fish_tl)**3.0264; *fishbase for Family;
if classcode='OYB' then weight= 0.00978*(fish_tl)**3.09; *Love 1978 for Olive RF;
if classcode='OYT' or classcode='SSER' or classcode='SFLA' then weight= 0.00978*(fish_tl)**3.09; *Love 1978;
if classcode='PATR' then weight= 0.000019*(fish_tl*10*0.7405)**3.065; *Antrim 1981 for White surfperch PFUR, SL = mm;
if classcode='PCAL' then weight= 0.0000939*(fish_tl*10*0.843)**3.0088; *Haacker 1975, SL = mm;
if classcode='PCLA' then weight= 0.0000066*(fish_tl*10/1.18)**3.26; *Quast 1968, SL = mm;
if classcode='PCOE' then weight= 0.00000956*(fish_tl*10*0.84)**3.19554; *Quast 1968, SL = mm;
if classcode='PFUR' then weight= 0.000019*(fish_tl*10*0.7405)**3.065; *Antrim 1981, SL = mm;

if classcode='PHOL' then weight= 0.0062*(fish_tl)**3.249; *Fishbase for Family;
if classcode='PLUE' or classcode='PLEU' then weight= 0.01917*(fish_tl)**2.9541; *Forrester et al 1969 for starry flounder;
if classcode='PMAC' then weight= 0.000026*(fish_tl*10/1.169)**3.0187; *Allen et al 1995, SL = mm;
if classcode='PNEB' then weight= 0.000026*(fish_tl*10/1.169)**3.0187; *Allen et al 1995, SL = mm, used eq for PMAC because DeMartini et al 1994 was weird;
if classcode='PTRI' then weight= 0.0056*(fish_tl)**2.9088; *Fishbase for Family;
if classcode='RHYP' then weight= 0.0091*(fish_tl)**3; *Fishbase for body shape;
if classcode='RJOR' then weight= 0.0091*(fish_tl)**3; *Fishbase for body shape;
if classcode='RTOX' then weight= 0.00000457*(fish_tl*10*0.853)**3.36; *Quast 1968, SL = mm;
if classcode='RYOY' then weight= 0.00000629*(fish_tl*10)**3.172; *Lea et al 1999 for SATR;
if classcode='SARG' then weight= 0.003962*(fish_tl)**2.983; *Walford 1932;
if classcode='SATR' then weight= 0.000006291*(fish_tl*10)**3.172; *Lea et al 1999;
if classcode='SAUR' then weight= 0.00000172*(fish_tl*10)**3.02; *Lea et al 1999;
if classcode='SCAL' then weight= 0.0117*(fish_tl)**2.934; *Fishbase for Genus;
if classcode='SCAR' then weight= 0.00001299*(fish_tl*10)**3.077; *Lea et al 1999;
if classcode='SCAU' then weight= 0.000008976*(fish_tl*10)**3.132; *Lea et al 1999;
if classcode='SCHI' then weight= 0.0117*(fish_tl)**3.03; *Fishbase;
if classcode='SCHR' then weight= 0.00001117*(fish_tl*10)**3.114; *Lea et al 1999;
if classcode='SDIP' then weight= 0.0115*(fish_tl)**3.0; *Fishbase;
if classcode='SEBSPP' then weight= 0.000006291*(fish_tl*10)**3.172; *Lea et al 1999, using values for SATR;
if classcode='SENT' then weight= 0.0159*(fish_tl)**3.0056; *Fishbase for Genus;
if classcode='SGIG' then weight= 0.0221*(fish_tl)**3.0; *Fishbase;
if classcode='SGUT' then weight= 0.0196*(fish_tl)**3.0; *Love et al 1987;
if classcode='SJAP' then weight= 0.0014*(fish_tl)**3.4; *Fishbase;
if classcode='SLAL' then weight= 0.0432 *(fish_tl/1.116)**2.85; *Fishbase;
if classcode='SMAR' then weight= 0.000005498*(fish_tl*10)**3.185; *Lea et al 1999;
if classcode='SMEL' then weight= 0.00000581*(fish_tl*10)**3.187; *Lea et al 1999;
if classcode='SMIN' then weight= 0.00001458*(fish_tl*10)**3.041; *Lea et al 1999;
if classcode='SMYS' then weight= 0.000009774*(fish_tl*10)**3.09; *Lea et al 1999;
if classcode='SNEB' then weight= 0.000007789*(fish_tl*10)**3.177; *Lea et al 1999;
if classcode='SPAU' then weight= 0.0132*(fish_tl)**3.0; *Fishbase;
if classcode='SPIN' then weight= 0.000006883*(fish_tl*10)**3.147; *Lea et al 1999;
if classcode='SPUL' then weight= 0.0238*(fish_tl)**2.8847; *Hamilton and Caselle unpublished;
if classcode='SRAS' then weight= 0.000007310*(fish_tl*10)**3.187; *Lea et al 1999;
if classcode='SSAG' then weight= 0.0076*(fish_tl*0.84)**3.150; *Fishbase;
if classcode='SSEM' then weight= 0.0127*(fish_tl)**3.016; *Love et al 1990 for females;
if classcode='STRE' then weight= 0.0000182*(fish_tl*10*0.88)**3.07; *Quast 1968, SL = mm, L-L in fishbase for blue rockfish;
if classcode='TREE' then weight= 0.0000182*(fish_tl*10*0.88)**3.07; *Quast 1968, SL = mm, L-L in fishbase for blue rockfish;
if classcode='SYNG' then weight= 0.00000144*(fish_tl*10*0.9)**2.65; *Mahan 1985, SL = mm;
if classcode='TCAL' then weight= 0.0487*(fish_tl)**2.7479; *Fishbase;
if classcode='TSEM' then weight= 453.59237*0.00000575*(fish_tl)**3.1044; *Ackerman 1971, TL = cm (?), W = lbs;
if classcode='TSYM' then weight= 0.0044*(fish_tl)**3.212; *Fishbase;


*CHECK ALL THESE UNITS;

if classcode='CVIO' then weight= 0.01289*(fish_tl)**2.9	; *NOTE THIS IS FOR SL, THERE IS NO CONVERSION FROM SL TO TL ;
if classcode='HAZU' then weight= 0.0000107*(fish_tl*10)**3.207435	;
if classcode='MMIN' then weight= 0.010398*(fish_tl)**3.207435	; *VRG table ;
if classcode='MMOL' then weight= 0.0454*(fish_tl)**3.0496	; *VRG table ;
if classcode='PNOT' then weight= 0.004163*(fish_tl)**3.282277	; *VRG table ;
if classcode='SACA' then weight= 0.0031*(fish_tl)**3.106; *VRG table ;
if classcode='SDAL' then weight= 0.00945*(fish_tl)**3.21542; *VRG table ;
if classcode='SDIP' then weight= 0.0195*(fish_tl)**2.927; *VRG table ;
if classcode='SHOP' then weight= 0.01464*(fish_tl)**2.96355;	*VRG table ;
if classcode='SLUC' then weight= 0.005551*(fish_tl)**3.028658;	;
if classcode='SROS' then weight= 0.00000989*(fish_tl*10)**3.068	;
if classcode='SRUB' then weight= 0.0000234675*(fish_tl*10)**2.9431;
if classcode='SSAX' then weight= 0.02479*(fish_tl)**2.8048	;
if classcode='XCAL' then weight= 0.0000107*(fish_tl*10)**2.91	;

*Check these equations;
if classcode='PGAL' or classcode = 'PGLA' then weight= 0.0032*(fish_tl*0.822)**3.131; *Fishbase, FL = cm SHOULD THIS BE PGLA= BLUE SHARK ??? ;
if classcode='ACOR' then weight= .0004*(.8414*fish_tl-.10798)**2.7971; *VRG table ;
if classcode='AINE' then weight=0.003145*(fish_tl)**3.030433	;  *from VRG table, Love unpublished data;
if classcode='BATH' or classcode='RATH' then weight= 0.01*(fish_tl)**3.04;  *fishbase for rathbunnella hypoplecta;
if classcode='PSTE' or classcode='FLAT' then weight=0.00000163*(fish_tl*10)**3.34  ; *starry founder Platichthys stellatus, from VRG table;
if classcode='AGON' then weight=0.003145*(fish_tl)**3.030433 ; *VRG for Agonidae MMS - Xeneretmus latifrons as L-W substitute; 
if classcode='HGUT' then weight=0.009057*(fish_tl)**3.172961 	; *diamond turbot Hypsopsetta guttulata, USED eqn for  hornyhead turbot Pleuronichthys verticalis from VRG	table ;
if classcode='XLIO' then weight=0.012048*(fish_tl)**3.050394;  *fantail sole Xystreurys liolepis, from VRG table;
if classcode='HGRI' then weight=0.01350*(fish_tl)**3 	; *BIOMASS BUSTER !!  bluntnose sixgill shark, Hexanchus griseus;
if classcode='NBLA' then weight=0.00931*(fish_tl)**2.99	; *USED A TROPICAL JAWFISH Opistognathus whitehursti, DUSKY JAWFISH, FROM FISHBASE, THIS NEED TO BE IMPROVED;
if classcode='ANOB' then weight=  0.00004*(fish_tl*10/1.2444)**2.7991 ;
if classcode='UHAL' then weight= 0.00734*(fish_tl)**3 	; *from FISHBASE CHECK THIS	;
if classcode='ZROS' then weight= 0.01532*(fish_tl)**2.92796	; *from VRG table ;
if classcode='AOCE' then weight= 0.00000503*(fish_tl)**3.091;
if classcode='OTRI' then weight=  0.00007*(fish_tl)**2.6479	;
if classcode='SUMB' then weight= 0.0067*(fish_tl)**3.319; *Chen 1971 from Love rockfish book ;
if classcode='RALL' then weight= 0.00000114*(fish_tl)**3.442696	; *VRG ;
if classcode='RPRO' then weight= 	0.00000343*(fish_tl)**3.0125;  *VRG	;
if classcode='ZEXA' then weight= 0.000001*(fish_tl*10)**3.4226	;  *Blanco-Parra et al. 2009 from VRG;
if classcode='LZEB' then weight= 0.00000343*(fish_tl)**3.101332;*(fish_tl)**3.0125; *VRG	;
if classcode='AGUA' then weight= 0.01034*(fish_tl*10)**3.07492	; *VRG;
if classcode='AARG' then weight= 0.03*(fish_tl*0.801*10)**3; *VRG and fishbase	;
if classcode='NCEP' then weight= 0.000000861*(fish_tl)**3.2152	; *VRG Sevengill

*separate out YOY categories,  can we get better l/w equations for YOY sizes?;
if fish_tl le 5 then do;
	if classcode = 'CPUN' then classcode = 'CPUNYOY';
    if classcode = 'BFRE' then classcode = 'BFREYOY';
    if classcode = 'OCAL' then classcode = 'OCALYOY';
    if classcode = 'HSEM' then classcode = 'HSEMYOY';
end;
if fish_tl le 10 then do;
    if classcode = 'PCLA' then classcode = 'PCLAYOY';
    if classcode = 'SMYS' then classcode = 'SMYSYOY';
    if classcode = 'SPAU' then classcode = 'SPAUYOY';
    if classcode = 'SPIN' then classcode = 'SPINYOY';
    if classcode = 'SPUL' then classcode = 'SPULYOY';
    if classcode  in ('SATR' 'SCAR' 'SCAU' 'SCHR' 'GBY') then classcode ='KGB';
    if classcode  in ('SFLA' 'SSER' 'SMEL' 'OYT') then classcode ='OYB';
    if classcode in ('SMIN' 'SNEB' 'SRAS' 'STRE' 'SSAX' 'SDAL' 'SDIP' 'SEBSPP')
        then classcode = 'RYOY';
end;
run;

data missingweight;
set fish;
where weight = . and classcode ne 'NO_ORG';
run;

proc sort data=fish; by classcode; run;
%excelin(&fileserver\UCSB\CI_reserves\PISCO_Fish_size_trophic_table.xls,trophic_group);
proc sort data=trophic_group; by classcode; run;
data fishforsize fishfordiversity;
merge fish (in=data) trophic_group;
by classcode;
if data;
*remove transect placeholder and cryptic spp that were not looked for in all years;
if classcode in ('NO_ORG' 'LDAL' 'CNIC') then delete;
*determine if fish are legal size;
legal=0;
if legal_size and fish_tl ge legal_size then legal=1;
output fishforsize;
*remove non-specific species obs for diversity calculations;
if for_diversity ne '' then classcode=for_diversity;
if classcode in ('DELETE') then delete;
output fishfordiversity;
run;
proc sort data=fishforsize;
by year site side classcode mature;
run;
proc univariate data=fishforsize noprint;
by year site side classcode;
var fish_tl;
output out=lengths mean=meanlen median=medianlen;
run;
*calculate weight and eggproduction ;
proc means noprint data=fishforsize;
by year site side classcode;
weight count;
var lsf weight;
output out=size sum=/autoname;
run;
*caluclate density and biomass of legal sized fish;
proc means noprint data=fishforsize ;
where legal;
by year site side classcode;
weight count;
var weight;
output out=legal n=legalcount sum=legalweight;
run;
*calculate density and biomass of mature sized fish;
proc means noprint data=fishforsize ;
where mature;
by year site side classcode;
weight count;
var weight;
output out=mature n=maturecount sum=matureweight;
run;
proc sql;
create table sizewt as select classcode,min(fish_tl) as minlength, min(weight) as minweight, max(fish_tl) as maxlength, max(weight) as maxweight
	from fishforsize group by classcode;

create table trans as select distinct year,site,side,zone,transect from fish;
create table sumtrans as select distinct year,site,side,count(transect) as numtrans from  trans 
	group by year,site,side;

create table propmature as select year,site,side,classcode,sum(mature*count)/sum(count) as propmature from fishforsize 
	group by year,site,side,classcode;
quit;
data size2;
merge size mature legal;
by year site side classcode;
run;
data size3;
merge size2 (in=indata) sumtrans;
by year site side;
eggprod=lsf_sum/(2*numtrans);
biomass=(weight_sum)/1000000/(numtrans*60/10000);  *NOTE IS THIS IN GRAMS/M^2??? ;
legalbiomass=(legalweight)/1000000/(numtrans*60/10000);
legaldens=legalcount/numtrans;
maturebiomass=(matureweight)/1000000/(numtrans*60/10000);
maturedens=maturecount/numtrans;

run;

proc sql;

select classcode from trophic_group where Legal_Size;

create table richness as select year,site,side,count(distinct classcode) as richness from fishfordiversity 
	group by year,site,side;
create table sitetotal as select year,site,side,sum(count) as total from fishfordiversity group by year,site,side;
create table sitespptotal as select year,site,side,classcode,sum(count) as count from fishfordiversity group by
	year,site,side,classcode;
quit;
data sppcounts;
merge sitespptotal sitetotal;
by year site side;
sw=-count/total*log(count/total);
simp=count*(count-1);
run;
proc means noprint data=sppcounts;
by year site side;
var sw;
output out=shannon sum=shannon;
proc means noprint data=sppcounts;
by year site side total;
var simp;
output out=simpson sum=;
run;
data simpson2;
set simpson;
simpsons=1-simp/(total*(total-1));
drop total simp;
run;
data diversity;
merge richness shannon simpson2;
by year site side;
evenness=shannon/log(richness);
drop _type_ _freq_;
run;

data combined;
merge propmature size3 lengths;
by year site side classcode;
run;
proc transpose data=combined out=combined2;
by year site side classcode;
var  propmature eggprod biomass maturedens maturebiomass legaldens legalbiomass meanlen medianlen;
run;
data combined3;
set combined2;
newname=compress(_name_)||'_'||compress(classcode);
*clean up meaningless lines;
if classcode not in ('SATR' 'SCHR' 'SCAU' 'SMIN' 'OYT' 'SMAR' 'HDEC' 'SCAR' 'SMEL' 'SMYS' 'SPUL' 'PCLA''CAGG'
'SDAL' 'SSAX' 'SSEM' 'OPIC' 'SHOP' 'SROS' 'SGUT' 'EJAC' 'HRUB' 'SAUR' 'RVAC' 'SRUB'
'SENT' 'SSER' 'SPAU' 'SPIN' 'OELO') 
and _name_ in ('propmature' 'maturedens' 'maturebiomass') then delete;

if classcode not in ('HDEC' 'OELO' 'PCAL' 'PCLA' 'PMAC' 'PNEB' 'SARG' 'SCHI' 'SGUT' 'SLAL' 'SMAR' 'SPIN' 'SPUL' 'TSEM' 'SARG' 'SCAR' 'SRAS'  'ANOB' 'SCAL') 
and _name_ in ('legaldens' 'legalbiomass') then delete;

if classcode not in ('SATR' 'SCHR''SCAU' 'SMIN' 'OYT' 'SMAR' 'SPUL' 'PCLA') 
and _name_ in ('eggprod') then delete;
run;

proc transpose data=combined3 out=combined4;
by year site side;
id newname;
var col1;
run;
data combined5;
merge combined4 diversity;
by year site side;
array all(*) biomass_: eggprod_: legaldens_: legalbiomass_: maturedens_: maturebiomass_:;
do i=1 to dim(all);
if all(i)=. then all(i)=0;
end;
run;
data mpasites;
set ref_tabl.Final_Site_Table_UCSB;
run;
proc sort data=mpasites; by site side; run;	*use Final Site Table UCSB.xls or new site table.xls in MLPA 2013	;
proc sort data=combined5; by site side year; run;
data ci_reserve_data_final ;
set CIMPA.ci_reserve_data_final;
run;
proc sort data=ci_reserve_data_final; by site side year; run;  *MUST MAKE THIS FILE FROM program ci-reserve_data which puts it in CIMPA13 ;
data fishsizedata mismatch;
merge combined5 (in=data) mpasites (in=sites);
by site side;
if (not data or not sites) and (campus='UCSB' or year) then output mismatch;
*if site='ANACAPA_ADMIRALS' and year=1999 then delete;
if year lt 2004  and site ne 'ANACAPA_EAST_ISLE' then period='BEFORE'; *this assumes all reserves started in 2003 and 2004 was first year we could have seen effects, ANA_EAST_ISLE;
	else period='AFTER';		      *must be kicked out of before/after analysis because all data is after;
if  
	/*and (reserve='OUT' or year ge year_mpa)  this takes only years after MPA started for reserve sites and all years for reference sites*/ 
	region in ('ANA' 'SCI' 'SRI' 'SBI' 'SMI')  and data
		then output fishsizedata;
run;
data biomassdata (keep=year site side eggprod: biomass: mpagroup--reserve)
	legalsize (keep=year site side legaldens: legalbiomass: mpagroup--reserve)
	maturesize (keep=year site side  maturedens: maturebiomass: mpagroup--reserve)
	lengthdata (keep=year site side meanlen: medianlen: mpagroup--reserve)
	diversity (keep=year site side richness evenness shannon simpsons mpagroup--reserve);
set fishsizedata;
run;
data ci_reserve_data_final final_mismatch;
merge ci_reserve_data_final (in=data) biomassdata (in=biomass);
by site side year;
if not data or not biomass then output final_mismatch;
output ci_reserve_data_final;
run;

data CIMPA.ci_reserve_data_final2;
set  ci_reserve_data_final;
run;


%excelout(biomassdata,&fileserver\UCSB\2013_MLPA_processing\CI reserve fishsizedata.xls );
%excelout(legalsize,&fileserver\UCSB\2013_MLPA_processing\CI reserve fishsizedata.xls );
%excelout(maturesize,&fileserver\UCSB\2013_MLPA_processing\CI reserve fishsizedata.xls );
%excelout(lengthdata,&fileserver\UCSB\2013_MLPA_processing\CI reserve fishsizedata.xls );	*THIS WONN"T EXPORT	;
%excelout(diversity,&fileserver\UCSB\2013_MLPA_processing\CI reserve fishsizedata.xls );

%excelout(ci_reserve_data_final,&fileserver\UCSB\2013_MLPA_processing\CI reserve data all.xls );


* cannot do this right now because I took out the clusters, we would want to rerun the clusters;
title 'univariate tests using fish cluster groups - all fish size vars';
proc glm data=fishsizedata (where=(region in ('A' 'B' 'C' 'D' 'E' 'F')));
class   year region reserve site;
model richness shannon evenness propmature_: eggprod_: biomass_:  maturedens_: maturebiomass_: legaldens_: legalbiomass_:= 
		year region|reserve site(region*reserve) ; 
*lsmeans reserve*region / pdiff;
test h=reserve e=region*reserve;
test h=year region region*reserve e=site(region*reserve);
quit;

title 'univariate tests using islands as groups - all fish size vars';
proc glm data=fishsizedata ;
class   year region reserve site;
model richness shannon evenness propmature_: eggprod_: biomass_:  maturedens_: maturebiomass_: legaldens_: legalbiomass_:= 
	year region|reserve site(region*reserve) ; 
*lsmeans reserve / pdiff;
test h=reserve e=region*reserve;
test h=region region*reserve e=site(region*reserve);
quit;

title 'univariate tests using fish cluster groups (before and after) - all fish size vars';
proc glm data=fishsizedata (where=(region in ('A' 'B' 'C' 'D' 'E' 'F') and site ne 'ANACAPA_EAST_ISLE'));
class region reserve site period;
model richness shannon evenness propmature_: eggprod_: biomass_:  maturedens_: maturebiomass_: legaldens_: legalbiomass_:= 
	region reserve period region*reserve reserve*period region*reserve*period site(region*reserve) ; *left year out as it was non-sig;
*lsmeans reserve / pdiff;
test h=reserve e=region*reserve;
test h=region region*reserve e=site(region*reserve);
test h=period e=reserve*period;
test h= reserve*period e=region*reserve*period;
quit;
*log response ratios;
proc sort data=fishsizedata ; by year region reserve site side; run;
proc transpose data=fishsizedata out=fishsizedata2 (rename=(_name_=classcode col1=biomass));
by year region reserve site side;
var biomass_:;
run;
proc sort data=fishsizedata2 
(where=(classcode in ('biomass_AFLA' 'biomass_BFRE' 'biomass_CPRI' 'biomass_CPUN' 'biomass_DVAC' 'biomass_EJAC' 'biomass_ELAT'
'biomass_GNIG' 'biomass_HRUB' 'biomass_HDEC' 'biomass_MCAL' 'biomass_OCAL' 'biomass_OPIC' 'biomass_PFUR' 'biomass_RTOX' 'biomass_PCLA'
'biomass_OELO' 'biomass_OYT' 'biomass_SATR' 'biomass_SAUR' 'biomass_SCAR' 'biomass_SCAU' 'biomass_SCHR' 'biomass_SMAR' 'biomass_SGUT'
'biomass_SMEL' 'biomass_SMIN' 'biomass_SMYS' 'biomass_SNEB' 'biomass_SPAU' 'biomass_SPIN' 'biomass_SPUL' 'biomass_SRAS' 'biomass_STRE')))
out=select; 
by year region reserve classcode; run;
proc means data=select noprint;
by year region reserve classcode; 
var biomass;
output out=means mean=;
run;
proc sort data=means; by year region classcode reserve; run;
proc transpose data=means out=ratios;
by year region classcode;
var biomass;
id reserve;
run;
data ratios2;
set ratios;
if in=0 then in=0.00001;
if out=0 then out=0.00001;
ratio=log(in/out);
rename _name_=classcode;
run;
proc sort data=ratios2; by classcode;run;
title 'univariate tests using year/island log response ratios - individual spp biomass';
proc glm data=ratios2; 
by classcode;
class   year region;
model ratio  = year region; 
lsmeans year / pdiff;
quit;
proc sort data=select;by classcode; run;

ods pdf file = "&fileserver\analysis\ci_reserves\indiv. spp. biomass by year.pdf";
ods graphics on;
 goptions reset=all;

title 'fish spp biomass response ratios - island replicates';
proc gchart data=ratios2;
by classcode;
vbar year/ midpoints=(1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012)
				patternid=by
					sumvar=ratio
					errorbar=top
					type=mean; 
run;
quit;
title 'fish spp biomass - sides as replicates';
proc gchart data=select;
by classcode;
vbar reserve/ group=year
			midpoints=(1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012)
			errorbar=top 	
				patternid=midpoint
             
			sumvar=biomass
					type=mean; 

run;
quit; 

proc gchart data=select (where=(year ge 2004));
by classcode;
title 'fish spp biomass - After MPAs';
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=biomass
					errorbar=top
					type=mean; 
run;
quit;

ods pdf close;
   ods graphics off;
title;


/**%excelout(fish6,&fileserver\analysis\ci_reserves_sitemeans.xls);
%macro graph;
%let spp = HDEC OELO OYT SATR SCAR SCAU SCHR SMAR SMEL SMIN SMYS SNEB SPAU SPIN SPUL SRAS STRE;
%do i=1 %to 17;
%let sp=%scan(&spp,&i);
proc gchart data=fishsizedata;
*axis1  origin=(,20);
title ;
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=propmature_&sp
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=biomass_&sp
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=eggprod_&sp
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=maturedens_&sp
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=maturebiomass_&sp
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=legaldens_&sp
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=legalbiomass_&sp
					errorbar=top
					type=mean; 

run;
quit;
%end;
%mend graph;
%graph
proc gchart data=fishsizedata (where=(period='BEFORE'));
*axis1  origin=(,20);
title 'before';
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=richness
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=shannon
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=evenness
					errorbar=top
					type=mean; 
run;
quit;
proc sort data=fishsizedata; by year; run;
proc gchart data=fishsizedata;* (where=(period='AFTER'));
*axis1  origin=(,20);
title 'after';
*by year;
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=richness
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=shannon
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=evenness
					errorbar=top
					type=mean; 
vbar reserve/ group=region 
				patternid=midpoint
					sumvar=simpsons
					errorbar=top
					type=mean; 
run;
quit;

title1 'Richness by year & region';

symbol1 interpol=mean
        value=dot
        height=3;
proc gplot data=fishsizedata;
plot richness*(year)= region;
run;
quit;
proc sql;
create table sites as select site,side,year,reserve,period from fishsizedata;
quit;
proc transpose data=sites out=sites2;
by site side;
id year;
var reserve;
run;
%excelout(sites2,ci reserve sites.xls);


