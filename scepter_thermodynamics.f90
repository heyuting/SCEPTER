module scepter_thermodynamics
    use scepter_constants
    use scepter_variables
    use scepter_concentration
    use scepter_findloc
    implicit none
    private
    public :: sld_therm, k_arrhenius, K_q10, calc_gamma_davies
    real(kind=8), parameter :: cal2j = 4.184d0

contains

    subroutine sld_therm( &
        & rg,tc,tempk_0,ss_x,ss_y,ss_z &! input
        & ,mineral &! input
        & ,therm &! output
        & ) 
        implicit none

        real(kind=8),intent(in)::rg,tc,tempk_0,ss_x,ss_y,ss_z
        real(kind=8),intent(out):: therm
        character(5),intent(in):: mineral
        real(kind=8) tc_ref,ha,therm_ref,delG
        real(kind=8) tc_ref_1,ha_1,therm_ref_1,therm_1,delG_1
        real(kind=8) tc_ref_2,ha_2,therm_ref_2,therm_2,delG_2
        real(kind=8) tc_ref_3,ha_3,therm_ref_3,therm_3,delG_3
        real(kind=8) tc_ref_4,ha_4,therm_ref_4,therm_4,delG_4
        real(kind=8) tc_ref_5,ha_5,therm_ref_5,therm_5,delG_5
        real(kind=8) tc_ref_6,ha_6,therm_ref_6,therm_6,delG_6
        real(kind=8) tc_ref_7,ha_7,therm_ref_7,therm_7,delG_7
        real(kind=8) tc_ref_8,ha_8,therm_ref_8,therm_8,delG_8

        ! real(kind=8) k_arrhenius

        therm = 0d0

        select case(trim(adjustl(mineral))) 
            case('ka') 
                ! Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
                therm_ref = 10d0**(7.435d0)
                ha = -35.3d0*cal2j
                tc_ref = 25d0
                ! from PHREEQC.DAT 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            ! case('ab')
                ! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
                ! therm_ref = 10d0**3.412182823d0
                ! ha = -54.15042876d0
                ! tc_ref = 15d0
                ! from Kanzaki and Murakami 2018
                ! therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('kfs')
                ! K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
                therm_ref = 10d0**0.227294204d0
                ha = -26.30862098d0
                tc_ref = 15d0
                ! from Kanzaki and Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('sdn')
                ! Sanidine_high: KAlSi3O8 +4.0000 H+  =  + 1.0000 Al+++ + 1.0000 K+ + 2.0000 H2O + 3.0000 SiO2
                therm_ref = 10d0**(0.9239d0)
                ha = -35.0284d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('leu')
                ! Leucite: KAlSi2O6 + 2H2O + 4H+ = 2H4SiO4 + Al+3 + K+
                therm_ref = 10d0**(6.423d0)
                ha = -22.085d0*cal2j
                tc_ref = 25d0
                ! from minteq.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('anl')
                ! NaAlSi2O6*H2O  + 5 H2O  = Na+  + Al(OH)4-  + 2 Si(OH)4(aq)
                therm_ref = 10d0**(-16.06d0)
                ha = 101d0
                tc_ref = 25d0
                ! from Wilkin and Barnes 1998
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('nph')
                ! Nepheline  + 4 H+  = 2 H2O  + SiO2(aq)  + Al+++  + Na+
                therm_ref = 10d0**(14.93646757d0)
                ha = -130.8197467d0
                tc_ref = 15d0
                ! from Kanzaki and Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('fo')
                ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O
                therm_ref = 10d0**29.41364324d0
                ha = -208.5932252d0
                tc_ref = 15d0
                ! from Kanzaki and Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('fa')
                ! Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O
                therm_ref = 10d0**19.98781342d0
                ha = -153.7676621d0
                tc_ref = 15d0
                ! from Kanzaki and Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            ! case('an')
                ! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
                ! therm_ref = 10d0**28.8615308d0
                ! ha = -292.8769275d0
                ! tc_ref = 15d0
                ! from Kanzaki and Murakami 2018
                ! therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('cc')
                ! CaCO3 = Ca2+ + CO32-
                therm_ref = 10d0**(-8.43d0)
                ! therm_ref = therm_ref*10d0 ! testing higher solubility
                ha = -8.028943471d0
                tc_ref = 15d0
                ! from Kanzaki and Murakami 2015
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('arg')
                ! CaCO3 = Ca2+ + CO32-
                therm_ref = 10d0**(-8.3d0)
                ha = -12d0
                tc_ref = 25d0
                ! from minteq.v4
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('dlm') ! disordered
                ! CaMg(CO3)2 = Ca+2 + Mg+2 + 2CO3-2
                therm_ref = 10d0**(-16.54d0)
                ha = -46.4d0
                tc_ref = 25d0
                ! from minteq.v4
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('amal')
                ! Al(OH)3 + 3 H+ = Al+3 + 3 H2O
                therm_ref = 10d0**(10.8d0)
                ha = -26.500d0*cal2j
                tc_ref = 25d0
                ! from PHREEQC.DAT 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('gb')
                ! Al(OH)3 + 3 H+ = Al+3 + 3 H2O
                therm_ref = 10d0**(8.11d0)
                ha = -22.80d0*cal2j
                tc_ref = 25d0
                ! from PHREEQC.DAT 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('amsi','sio2')
                ! SiO2 + 2 H2O = H4SiO4
                therm_ref = 10d0**(-2.71d0)
                ha = 3.340d0*cal2j
                tc_ref = 25d0
                ! from PHREEQC.DAT 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('phsi')
                ! SiO2 + 2 H2O = H4SiO4
                therm_ref = 10d0**(-2.74d0)
                ha = 10.85d0
                tc_ref = 25d0
                ! from Fraysse et al. 2006 for banboo phytolith 
                ! (see Fraysse et al. 2009 for different values for different plant phytoliths) 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('qtz')
                ! SiO2 + 2H2O = H4SiO4
                therm_ref = 10d0**(-4d0)
                ha = 22.36d0
                tc_ref = 25d0
                ! from minteq.v4 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('amfe3')
                ! Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
                therm_ref = 10d0**(4.891d0)
                ha = 0d0*cal2j
                tc_ref = 25d0
                ! from PHREEQC.DAT 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('gt')
                ! Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
                therm_ref = 10d0**(0.5345d0)
                ha = -61.53703d0
                tc_ref = 25d0
                ! from Sugimori et al. 2012 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('hm')
                ! Fe2O3 + 6H+ = 2Fe+3 + 3H2O
                therm_ref = 10d0**(-1.418d0)
                ha = -128.987d0
                tc_ref = 25d0
                ! from minteq.v4
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('ct')
                ! Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
                therm_ref = 10d0**(32.2d0)
                ha = -46.800d0*cal2j
                tc_ref = 25d0
                ! from PHREEQC.DAT 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('mscv')
                ! KAl2(AlSi3O10)(OH)2 + 10 H+  = 6 H2O  + 3 SiO2(aq)  + K+  + 3 Al+++
                therm_ref = 10d0**(15.97690572d0)
                ha = -230.7845245d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('plgp')
                ! KMg3(AlSi3O10)(OH)2 + 10 H+  = 6 H2O  + 3 SiO2(aq)  + Al+++  + K+  + 3 Mg++
                therm_ref = 10d0**(40.12256823d0)
                ha = -312.7817497d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018 
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('cabd')
                ! Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
                therm_ref = 10d0**(7.269946518d0)
                ha = -157.0186168d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('mgbd')
                ! Beidellit-Mg  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Mg++
                therm_ref = 10d0**(7.270517113d0)
                ha = -160.1864268d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('nabd')
                ! Beidellit-Na  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .33 Na+
                therm_ref = 10d0**(7.288837383d0)
                ha = -150.7328834d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('kbd')
                ! Beidellit-K  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .33 K+
                therm_ref = 10d0**(6.928086412d0)
                ha = -145.6776905d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('casp')
                ! Ca.165Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.1650 Ca++ + 0.3300 Al+++ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O
                therm_ref = 10d0**(26.2900d0)
                ha = -207.971d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('ksp')
                ! K.33Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 K+ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O
                therm_ref = 10d0**(26.0075d0)
                ha = -196.402d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('nasp')
                ! Na.33Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 Na+ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O
                therm_ref = 10d0**(26.3459d0)
                ha = -201.401d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('mgsp')
                ! Mg3.165Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 3.1650 Mg++ + 3.6700 SiO2 + 4.6600 H2O
                therm_ref = 10d0**(26.2523d0)
                ha = -210.822d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('ill')
                ! Illite  + 8 H+  = 5 H2O  + .6 K+  + .25 Mg++  + 2.3 Al+++  + 3.5 SiO2(aq)
                therm_ref = 10d0**(10.8063184d0)
                ha = -166.39733d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            ! case('dp')
                ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
                ! therm_ref = 10d0**(21.79853309d0)
                ! ha = -138.6020832d0
                ! tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                ! therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            ! case('hb')
                ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
                ! therm_ref = 10d0**(20.20981116d0)
                ! ha = -128.5d0
                ! tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                ! therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('jd')
                ! Jadeite  NaAl(SiO3)2 +4.0000 H+  =  + 1.0000 Al+++ + 1.0000 Na+ + 2.0000 H2O + 2.0000 SiO2
                therm_ref = 10d0**(8.3888d0)
                ha = -84.4415d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                ! $Id: llnl.dat 12776 2017-08-02 20:02:16Z dlpark $
                ! Data are from 'thermo.com.V8.R6.230' prepared by Jim Johnson at
                ! Lawrence Livermore National Laboratory, in Geochemist's Workbench
                ! format. Converted to Phreeqc format by Greg Anderson with help from
                ! David Parkhurst. A few organic species have been omitted.  
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('wls')
                ! Wollastonite  CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2
                therm_ref = 10d0**(13.7605d0)
                ha = -76.5756d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('tm')
                ! Tremolite  + 14 H+  = 8 H2O  + 8 SiO2(aq)  + 2 Ca++  + 5 Mg++
                therm_ref = 10d0**(61.6715d0)
                ha = -429.0d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('antp')
                ! Anthophyllite (Mg2Mg5(Si8O22)(OH)2) + 14 H+  = 8 H2O  + 7 Mg++  + 8 SiO2(aq)
                therm_ref = 10d0**(70.83527792d0)
                ha = -508.6621624d0
                tc_ref = 15d0
                ! from Kanzaki & Murakami 2018
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('splt')
                ! Sepiolite: Mg4Si6O15(OH)2:6H2O +8.0000 H+  =  + 4.0000 Mg++ + 6.0000 SiO2 + 11.0000 H2O
                therm_ref = 10d0**(30.4439d0)
                ha = -157.339d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('ep')
                ! Epidote: Ca2FeAl2Si3O12OH +13.0000 H+  =  + 1.0000 Fe+++ + 2.0000 Al+++ + 2.0000 Ca++ + 3.0000 SiO2 + 7.0000 H2O
                therm_ref = 10d0**(32.9296d0)
                ha = -386.451d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('clch')
                ! Clinochlore-14A: Mg5Al2Si3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Mg++ + 12.0000 H2O
                therm_ref = 10d0**(67.2391d0)
                ha = -612.379d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('cdr')
                ! Cordierite_anhyd: Mg2Al4Si5O18 +16.0000 H+  =  + 2.0000 Mg++ + 4.0000 Al+++ + 5.0000 SiO2 + 8.0000 H2O
                therm_ref = 10d0**(52.3035d0)
                ha = -626.219d0
                tc_ref = 25d0
                ! from llnl.dat in Phreeqc
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('gps')
                ! CaSO4*2H2O = Ca+2 + SO4-2 + 2H2O
                therm_ref = 10d0**(-4.61d0)
                ha = 1d0
                tc_ref = 25d0
                ! from minteq.v4
            case('caso4')
                ! CaSO4 = Ca+2 + SO4-2 + 2H2O
                therm_ref = 10d0**(-4.36d0)
                ha = -7.2d0
                tc_ref = 25d0
                ! from minteq.v4
            case('fe2o')
                ! FeO +2.0000 H+  =  + 1.0000 Fe++ + 1.0000 H2O
                therm_ref = 10d0**(13.5318d0)
                ha = -106.052d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('mgo')
                ! MgO +2.0000 H+  =  + 1.0000 H2O + 1.0000 Mg++
                therm_ref = 10d0**(21.3354d0)
                ha = -150.139d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('k2o')
                ! K2O +2.0000 H+  =  + 1.0000 H2O + 2.0000 K+
                therm_ref = 10d0**(84.0405d0)
                ha = -427.006d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('cao')
                ! Lime
                ! CaO +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O
                therm_ref = 10d0**(32.5761d0)
                ha = -193.832d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('na2o')
                ! Na2O +2.0000 H+  =  + 1.0000 H2O + 2.0000 Na+
                therm_ref = 10d0**(67.4269d0)
                ha = -351.636d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('al2o3')
                ! Corundum
                ! Al2O3 +6.0000 H+  =  + 2.0000 Al+++ + 3.0000 H2O
                therm_ref = 10d0**(18.3121d0)
                ha = -258.626d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('nacl')
                ! Halite
                ! NaCl  =  + 1.0000 Cl- + 1.0000 Na+
                therm_ref = 10d0**(1.5855d0)
                ha = 3.7405d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('kcl')
                ! Sylvite
                ! KCl  =  + 1.0000 Cl- + 1.0000 K+
                therm_ref = 10d0**(0.8459d0)
                ha = 17.4347d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('cacl2')
                ! Hydrophilite
                ! CaCl2  =  + 1.0000 Ca++ + 2.0000 Cl-
                therm_ref = 10d0**(11.7916d0)
                ha = -81.4545d0
                tc_ref = 25d0
                ! from llnl.dat
                therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
            case('la','ab','an','by','olg','and')
                ! CaxNa(1-x)Al(1+x)Si(3-x)O8 + (4x + 4) = xCa+2 + (1-x)Na+ + (1+x)Al+++ + (3-x)SiO2(aq) 
                ! obtaining Anorthite 
                therm_ref_1 = 10d0**28.8615308d0
                ha_1 = -292.8769275d0
                tc_ref_1 = 15d0
                ! from Kanzaki and Murakami 2018
                therm_ref_1 = 10d0**(-19.714d0 - 2d0* (-22.7d0))
                ha_1 = (11.580d0 - 2d0* (42.30d0))*cal2j
                tc_ref_1 = 25d0
                ! from PHREEQC.DAT
                therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg) ! rg in kJ mol^-1 K^-1
                delG_1 = - rg*(tc+tempk_0)*log(therm_1) ! del-G = -RT ln K  now in kJ mol-1
                ! Then albite 
                therm_ref_2 = 10d0**3.412182823d0
                ha_2 = -54.15042876d0
                tc_ref_2 = 15d0
                ! from Kanzaki and Murakami 2018
                therm_ref_2 =  10d0**(-18.002d0 - 1d0* (-22.7d0))
                ha_2 = (25.896d0 - 1d0* (42.30d0))*cal2j
                tc_ref_2 = 25d0
                ! from PHREEQC.DAT
                therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg)
                delG_2 = - rg*(tc+tempk_0)*log(therm_2) ! del-G = -RT ln K  now in kJ mol-1
                
                if (ss_x == 1d0) then 
                    delG = delG_1 ! ideal anorthite
                elseif (ss_x == 0d0) then 
                    delG = delG_2 ! ideal albite
                elseif (ss_x > 0d0 .and. ss_x < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x))
                endif 
                therm = exp(-delG/(rg*(tc+tempk_0)))
            case('cpx','hb','dp')
                ! FexMg(1-x)CaSi2O6 + 4 H+  = Ca++  + 2 H2O  + xFe++ + (1-x)Mg++  + 2 SiO2(aq)
                ! obtaining hedenbergite 
                therm_ref_1 = 10d0**(20.20981116d0)
                ha_1 = -128.5d0
                tc_ref_1 = 15d0
                ! from Kanzaki & Murakami 2018
                therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg) ! rg in kJ mol^-1 K^-1
                delG_1 = - rg*(tc+tempk_0)*log(therm_1) ! del-G = -RT ln K  now in kJ mol-1
                ! Then diopside 
                therm_ref_2 = 10d0**(21.79853309d0)
                ha_2 = -138.6020832d0
                tc_ref_2 = 15d0
                ! from Kanzaki and Murakami 2018
                therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg)
                delG_2 = - rg*(tc+tempk_0)*log(therm_2) ! del-G = -RT ln K  now in kJ mol-1
                
                if (ss_x == 1d0) then 
                    delG = delG_1 ! ideal hedenbergite
                elseif (ss_x == 0d0) then 
                    delG = delG_2 ! ideal diopside
                elseif (ss_x > 0d0 .and. ss_x < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x))
                endif 
                therm = exp(-delG/(rg*(tc+tempk_0)))
            case('opx','en','fer')
                ! FexMg(1-x)SiO3 + 2 H+  = xFe++ + (1-x)Mg++  +  SiO2(aq)
                ! obtaining ferrosilite
                therm_ref_1 = 10d0**(7.777162795d0)
                ha_1 = -60.08612326d0
                tc_ref_1 = 15d0
                ! from Kanzaki & Murakami 2018
                therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg) ! rg in kJ mol^-1 K^-1
                delG_1 = - rg*(tc+tempk_0)*log(therm_1) ! del-G = -RT ln K  now in kJ mol-1
                ! Then enstatite 
                therm_ref_2 = 10d0**(11.99060855d0)
                ha_2 = -85.8218778d0
                tc_ref_2 = 15d0
                ! from Kanzaki and Murakami 2018
                therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg)
                delG_2 = - rg*(tc+tempk_0)*log(therm_2) ! del-G = -RT ln K  now in kJ mol-1
                
                if (ss_x == 1d0) then 
                    delG = delG_1 ! ideal ferrosilite
                elseif (ss_x == 0d0) then 
                    delG = delG_2 ! ideal enstatite
                elseif (ss_x > 0d0 .and. ss_x < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x))
                endif 
                therm = exp(-delG/(rg*(tc+tempk_0)))
            case('agt')
                ! obtaining opx (FexMg(1-x)SiO3 + 2 H+  = xFe++ + (1-x)Mg++  +  SiO2(aq))
                ! obtaining ferrosilite
                therm_ref_1 = 10d0**(7.777162795d0)
                ha_1 = -60.08612326d0
                tc_ref_1 = 15d0
                ! from Kanzaki & Murakami 2018
                therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg) ! rg in kJ mol^-1 K^-1
                ! converting to the formula Fe2Si2O6 + 4 H+  = 2Fe++ + 2SiO2(aq)
                therm_1 = therm_1**2d0
                delG_1 = - rg*(tc+tempk_0)*log(therm_1) ! del-G = -RT ln K  now in kJ mol-1
                
                ! Then enstatite 
                therm_ref_2 = 10d0**(11.99060855d0)
                ha_2 = -85.8218778d0
                tc_ref_2 = 15d0
                ! from Kanzaki and Murakami 2018
                therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg)
                ! converting to the formula Mg2Si2O6 + 4 H+  = 2Mg++ + 2SiO2(aq)
                therm_2 = therm_2**2d0
                delG_2 = - rg*(tc+tempk_0)*log(therm_2) ! del-G = -RT ln K  now in kJ mol-1
                
                if (ss_x == 1d0) then 
                    delG_3 = delG_1 ! ideal ferrosilite
                elseif (ss_x == 0d0) then 
                    delG_3 = delG_2 ! ideal enstatite
                elseif (ss_x > 0d0 .and. ss_x < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG_3 = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x))
                endif 
                therm_3 = exp(-delG_3/(rg*(tc+tempk_0)))
                
                ! obtaining cpx (FexMg(1-x)CaSi2O6 + 4 H+  = Ca++  + 2 H2O  + xFe++ + (1-x)Mg++  + 2 SiO2(aq))
                ! obtaining hedenbergite 
                therm_ref_4 = 10d0**(20.20981116d0)
                ha_4 = -128.5d0
                tc_ref_4 = 15d0
                ! from Kanzaki & Murakami 2018
                therm_4 = k_arrhenius(therm_ref_4,tc_ref_4+tempk_0,tc+tempk_0,ha_4,rg) ! rg in kJ mol^-1 K^-1
                delG_4 = - rg*(tc+tempk_0)*log(therm_4) ! del-G = -RT ln K  now in kJ mol-1
                ! Then diopside 
                therm_ref_5 = 10d0**(21.79853309d0)
                ha_5 = -138.6020832d0
                tc_ref_5 = 15d0
                ! from Kanzaki and Murakami 2018
                therm_5 = k_arrhenius(therm_ref_5,tc_ref_5+tempk_0,tc+tempk_0,ha_5,rg)
                delG_5 = - rg*(tc+tempk_0)*log(therm_5) ! del-G = -RT ln K  now in kJ mol-1
                
                if (ss_x == 1d0) then 
                    delG_6 = delG_4 ! ideal hedenbergite
                elseif (ss_x == 0d0) then 
                    delG_6 = delG_5 ! ideal diopside
                elseif (ss_x > 0d0 .and. ss_x < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG_6 = ss_x*delG_4 + (1d0-ss_x)*delG_5 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x))
                endif 
                therm_6 = exp(-delG_6/(rg*(tc+tempk_0)))
                
                ! finally mixing opx and cpx 
                if (ss_y == 1d0) then 
                    delG_7 = delG_3 ! ideal opx
                elseif (ss_y == 0d0) then 
                    delG_7 = delG_6 ! ideal cpx
                elseif (ss_y > 0d0 .and. ss_y < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG_7 = ss_y*delG_3 + (1d0-ss_y)*delG_6 + rg*(tc+tempk_0)*(ss_y*log(ss_y)+(1d0-ss_y)*log(1d0-ss_y))
                endif 
                therm_7 = exp(-delG_7/(rg*(tc+tempk_0)))
                
                ! obtaining Jadeite  
                therm_ref_8 = 10d0**(8.3888d0)
                ha_8 = -84.4415d0
                tc_ref_8 = 25d0
                ! from llnl.dat in Phreeqc
                therm_8 = k_arrhenius(therm_ref_8,tc_ref_8+tempk_0,tc+tempk_0,ha_7,rg)
                delG_8 = - rg*(tc+tempk_0)*log(therm_5) ! del-G = -RT ln K  now in kJ mol-1
                
                ! finally mixing non-Na pyroxene (delG_7) and Na-bearing pyroxene (delG_8) 
                if (ss_z == 1d0) then 
                    delG = delG_8 ! ideal Na-bearlng pyroxene
                elseif (ss_z == 0d0) then 
                    delG = delG_7 ! ideal no-Na-bearlng pyroxene
                elseif (ss_z > 0d0 .and. ss_z < 1d0) then  ! solid solution 
                    ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
                    delG = ss_z*delG_8 + (1d0-ss_z)*delG_7 + rg*(tc+tempk_0)*(ss_z*log(ss_z)+(1d0-ss_z)*log(1d0-ss_z))
                endif 
                therm = exp(-delG/(rg*(tc+tempk_0)))
                
            case('g1')
                therm = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
            case('g2')
                therm = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
                ! therm = 0.121d-1 ! mo2 Michaelis, Davidson et al. (2012) x 0.1
                ! therm = 0.121d-2 ! mo2 Michaelis, Davidson et al. (2012) x 0.01
                ! therm = 0.121d-3 ! mo2 Michaelis, Davidson et al. (2012) x 0.001
                ! therm = 0.121d-6 ! mo2 Michaelis, Davidson et al. (2012) x 1e-6
            case('g3')
                therm = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
            case('amnt')
                therm = 0.121d0 ! mo2 Michaelis for fertilizer
            case('gac','mesmh','ims','teas','naoh','naglp')
                therm = 1d0 ! reacting in any case
            case('inrt')
                therm = 1d0 ! not reacting in any case
            case default 
                therm = 0d0
        endselect 

    endsubroutine sld_therm



    subroutine calc_gamma_davies( &
        & nz,iosx,tc,charge &
        & ,gamma,dgamma_dis &
        & )
        implicit none

        integer,intent(in)::nz
        real(kind=8),intent(in)::iosx(nz),tc,charge
        real(kind=8),intent(out)::gamma(nz),dgamma_dis(nz)
        real(kind=8) epsiron,a,b

        epsiron = 87.74d0 - 0.40008d0*tc+0.0009398d0*tc**2d0 - 0.00000141d0*tc**3d0 ! dielectric constant  Malmberg & Maryott 1956
        a = 1.824d6*( epsiron*(tc + 273.15d0 ) )**(-3d0/2d0)
        b = 0.3d0
        ! b = 0.2d0
        gamma = -a*charge**2d0*(iosx**0.5d0/(1d0+iosx**0.5d0) -b*iosx)
        dgamma_dis = -a*charge**2d0*( &
            & 0.5d0*iosx**(-0.5d0)/(1d0+iosx**0.5d0)   & 
            & + iosx**0.5d0*(-1d0)/(1d0+iosx**0.5d0)**2d0*0.5d0*iosx**(-0.5d0)  &
            & - b &
            & )
        ! dgamma_dis = d( 10d0**gamma ) /d (gamma) * d( gamma )/d( is ) = (log10)*10**gamma * dgamma_dis
        dgamma_dis = log(10d0) * 10d0**gamma * dgamma_dis
        gamma = 10d0**gamma

    endsubroutine calc_gamma_davies

    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    function k_arrhenius(kref,tempkref,tempk,eapp,rg)
        implicit none
        real(kind=8) k_arrhenius,kref,tempkref,tempk,eapp,rg
        k_arrhenius = kref*exp(-eapp/rg*(1d0/tempk-1d0/tempkref))
    endfunction k_arrhenius
    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    function k_q10(kref,tc,tc_ref,q10)
        implicit none
        real(kind=8) k_q10,kref,tc_ref,tc,q10
        k_q10 = kref*q10**( (tc - tc_ref)/10d0  )
    endfunction k_q10
    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

end module scepter_thermodynamics 