!-----------------------------------------------------------------------
! Module: scepter_eq_coefs
! Purpose: Calculate various equilibrium constants, diffusion coefficients,
!          and reaction rates for aqueous, gaseous, and solid species
!-----------------------------------------------------------------------    
module scepter_eq_coefs
    implicit none
    private
    public :: coefs_v2

    ! Constants
    real(kind=8), parameter :: cal2j = 4.184d0
    
    contains

    !-----------------------------------------------------------------------
    !Subroutine: coefs_v2
    !Purpose: Main coefficient calculation subroutine, called by scepter_equilibrium.f90
    !-----------------------------------------------------------------------
    subroutine coefs_v2( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pro,cec_pH_depend,mcec_all,logkhaq_all &! input
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
        & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
        & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,mwt_all,staq_all &!input
        & ,nsp_aq,nsp_aq_cnst,chraq,chraq_cnst,maq,maqc &!input
        & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 &! output
        & ,keqaq_oxa,keqaq_cl &! output
        & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
        & ,keqcec_all,keqiex_all &! output 
        & ) 
        
        implicit none

        integer,intent(in)::nz
        real(kind=8),intent(in)::rg,rg2,tc,sec2yr,tempk_0
        real(kind=8),dimension(nz),intent(in)::pro
        real(kind=8),dimension(nz)::oh,po2,kin,dkin_dmsp
        real(kind=8) kho,po2th,mv_tmp,therm,ss_x,ss_y,ss_z,ss_tmp,therm_tmp,mwt_tmp,visc
        real(kind=8),intent(out)::ucv,kw

        integer,intent(in)::nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        real(kind=8),dimension(nsp_aq_all),intent(out)::daq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all
        real(kind=8),dimension(nsp_gas_all),intent(out)::dgasa_all,dgasg_all
        real(kind=8),dimension(nsp_gas_all,3),intent(out)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(out)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_c
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc

        integer ieqaq_co3,ieqaq_hco3
        data ieqaq_co3,ieqaq_hco3/1,2/

        integer ieqaq_so4,ieqaq_so42
        data ieqaq_so4,ieqaq_so42/1,2/

        integer ieqaq_no3,ieqaq_no32
        data ieqaq_no3,ieqaq_no32/1,2/

        integer ieqaq_nh3,ieqaq_nh32
        data ieqaq_nh3,ieqaq_nh32/1,2/

        integer ieqaq_oxa,ieqaq_oxa2
        data ieqaq_oxa,ieqaq_oxa2/1,2/

        integer ieqaq_cl,ieqaq_cl2
        data ieqaq_cl,ieqaq_cl2/1,2/

        integer isps,ispss,ispa

        ! real(kind=8)::thon = 1d0
        real(kind=8)::thon = -1d100
        character(5) mineral,ssaq,sssld,aqsp
        character(5),dimension(7):: chrss_gbas_aq,chrss_cbas_aq
        character(5),dimension(7):: chrss_gbas_sld,chrss_cbas_sld
        logical zero_cec

        real(kind=8),parameter::mcec_threshold = 0d0
        real(kind=8) keqiex_prona

        chrss_gbas_aq  = (/'si   ','al   ','na   ','k    ','mg   ','ca   ','fe2  '/)
        chrss_cbas_aq  = (/'si   ','al   ','na   ','k    ','mg   ','ca   ','fe2  '/)
        chrss_gbas_sld = (/'amsi ','al2o3','na2o ','k2o  ','mgo  ','cao  ','fe2o '/)
        chrss_cbas_sld = (/'qtz  ','al2o3','na2o ','k2o  ','mgo  ','cao  ','fe2o '/)

        ucv = 1.0d0/(rg2*(tempk_0+tc))

        ! water viscosity (Pa.s) from Likhachev 2003
        visc = 0.000024152d0*EXP(4.7428d0/(tc+tempk_0-139.86d0)/rg)
        ! converting Pa.s to centipoise
        visc = 1d3*visc

        ! Aq species diffusion from Li and Gregory 1974 except for Si which is based on Rebreanu et al. 2008
        daq_all(findloc(chraq_all,'fe2',dim=1)) = k_arrhenius(1.7016d-2    , 15d0+tempk_0, tc+tempk_0, 19.615251d0, rg)
        daq_all(findloc(chraq_all,'fe3',dim=1)) = k_arrhenius(1.5664d-2    , 15d0+tempk_0, tc+tempk_0, 14.33659d0 , rg)
        daq_all(findloc(chraq_all,'so4',dim=1)) = k_arrhenius(2.54d-2      , 15d0+tempk_0, tc+tempk_0, 20.67364d0 , rg)
        daq_all(findloc(chraq_all,'no3',dim=1)) = k_arrhenius(4.6770059d-2 , 15d0+tempk_0, tc+tempk_0, 18.00685d0 , rg)
        daq_all(findloc(chraq_all,'na' ,dim=1)) = k_arrhenius(3.19d-2      , 15d0+tempk_0, tc+tempk_0, 20.58566d0 , rg)
        daq_all(findloc(chraq_all,'k'  ,dim=1)) = k_arrhenius(4.8022699d-2 , 15d0+tempk_0, tc+tempk_0, 18.71816d0 , rg)
        daq_all(findloc(chraq_all,'mg' ,dim=1)) = k_arrhenius(1.7218079d-2 , 15d0+tempk_0, tc+tempk_0, 18.51979d0 , rg)
        daq_all(findloc(chraq_all,'si' ,dim=1)) = k_arrhenius(2.682396d-2  , 15d0+tempk_0, tc+tempk_0, 22.71378d0 , rg)
        daq_all(findloc(chraq_all,'ca' ,dim=1)) = k_arrhenius(1.9023312d-2 , 15d0+tempk_0, tc+tempk_0, 20.219661d0, rg)
        daq_all(findloc(chraq_all,'al' ,dim=1)) = k_arrhenius(1.1656226d-2 , 15d0+tempk_0, tc+tempk_0, 21.27788d0 , rg)
        daq_all(findloc(chraq_all,'cl' ,dim=1)) = k_arrhenius(4.9363501d-2 , 15d0+tempk_0, tc+tempk_0, 18.948983d0, rg)

        ! organic acid 
        ! oxalic acid (value at 25 oC from Wen et al. 2008; activation just assumed)
        daq_all(findloc(chraq_all,'oxa',dim=1)) = k_arrhenius(3.114735d-02, 25d0+tempk_0, tc+tempk_0, 20.00000d0  , rg)
        ! acetic acid (Schulz and Zabel, 2006)
        daq_all(findloc(chraq_all,'ac',dim=1)) = k_arrhenius(2.51198D-02 , 15d0+tempk_0, tc+tempk_0, 21.569542d0 , rg)
        ! MES ( Othmer and Thakar (1953) equation and molecular volume from solid mes monohydrate )
        daq_all(findloc(chraq_all,'mes',dim=1))  &
            & = 14d-5 /( visc**1.1d0 * (mv_all(findloc(chrsld_all,'mesmh',dim=1)))**0.6d0 ) * sec2yr *1d-4 ! sec2yr*1d-4 converting cm2/s to m2/yr
        ! Imidazole ( Othmer and Thakar (1953) equation and molecular volume from La-Scalea et al. 2005 )
        daq_all(findloc(chraq_all,'im',dim=1))  = 14d-5 /( visc**1.1d0 * (75.3d0)**0.6d0 ) * sec2yr *1d-4 ! sec2yr*1d-4 converting cm2/s to m2/yr
        ! TriEthanolAmine ( Othmer and Thakar (1953) equation and molecular volume from La-Scalea et al. 2005 )
        daq_all(findloc(chraq_all,'tea',dim=1)) = 14d-5 /( visc**1.1d0 * (177.3d0)**0.6d0 ) * sec2yr *1d-4 ! sec2yr*1d-4 converting cm2/s to m2/yr
        ! Glycerophosphate (value for glycerol from Schramke et al. 1999 for now)
        daq_all(findloc(chraq_all,'glp',dim=1)) = 0.93d-5 * sec2yr *1d-4 ! sec2yr*1d-4 converting cm2/s to m2/yr

        #ifdef disp_cnst
        daq_all=disp_cnst
        #endif 

        ! --------------------------------- gas diff

        ! values used in Kanzaki and Murakami 2016 for oxygen 
        dgasa_all(findloc(chrgas_all,'po2',dim=1)) = k_arrhenius(5.49d-2 , 15d0+tempk_0, tc+tempk_0, 20.07d0 , rg)
        dgasg_all(findloc(chrgas_all,'po2',dim=1)) = k_arrhenius(6.09d2  , 15d0+tempk_0, tc+tempk_0, 4.18d0  , rg)

        ! assuming a value of 0.14 cm2/sec (e.g., Pritchard and Currie, 1982) and O2 gas activation energy for CO2 gas 
        ! and CO32- diffusion from Li and Greogy 1974 for aq CO2 
        dgasa_all(findloc(chrgas_all,'pco2',dim=1)) = k_arrhenius(2.2459852d-2, 15d0+tempk_0, tc+tempk_0, 21.00564d0, rg)
        dgasg_all(findloc(chrgas_all,'pco2',dim=1)) = k_arrhenius(441.504d0   , 15d0+tempk_0, tc+tempk_0, 4.18d0    , rg)

        ! NH4+ diffusion for aqueous diffusion from Schulz and Zabel 2005
        ! NH3 diffusion in air from Massman 1998
        dgasa_all(findloc(chrgas_all,'pnh3',dim=1)) = k_arrhenius(4.64d-02    , 15d0+tempk_0, tc+tempk_0, 19.15308d0, rg)
        dgasg_all(findloc(chrgas_all,'pnh3',dim=1)) = 0.1978d0*((tc+tempk_0)/(0d0+tempk_0))**1.81d0 * sec2yr *1d-4 ! sec2yr*1d-4 converting cm2 to m2 and sec-1 to yr-1

        ! assuming the same diffusion as CO2 diffusion (e.g., Pritchard and Currie, 1982) for gaseous N2O 
        ! N2O(aq) diffusion from Schulz and Zabel 2005
        dgasa_all(findloc(chrgas_all,'pn2o',dim=1)) = k_arrhenius(4.89d-02    , 15d0+tempk_0, tc+tempk_0, 20.33417d0, rg)
        dgasg_all(findloc(chrgas_all,'pn2o',dim=1)) = k_arrhenius(441.504d0   , 15d0+tempk_0, tc+tempk_0, 4.18d0    , rg)

        #ifdef disp_cnst
        dgasa_all=disp_cnst
        #endif 

        kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc**2d0+0.000000555d0*tc**3d0-0.0000000007581d0*tc**4d0  ! Murakami et al. 2011
        kw = k_arrhenius(10d0**(-14.35d0), tempk_0+15.0d0, tempk_0+tc, 58.736742d0, rg) ! from Kanzaki and Murakami 2015

        oh = kw/pro


        keqgas_h = 0d0

        ! kho = k_arrhenius(10.0d0**(-2.89d0), tempk_0+25.0d0, tempk_0+tc, -13.2d0, rg)
        keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0) = &
            & k_arrhenius(10d0**(-2.89d0), tempk_0+25.0d0, tempk_0+tc, -13.2d0, rg)
        kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)

        keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0) = &
            & k_arrhenius(10d0**(-1.34d0), tempk_0+15.0d0, tempk_0+tc, -21.33183d0, rg) ! from Kanzaki and Murakami 2015
        keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1) = &
            & k_arrhenius(10d0**(-6.42d0), tempk_0+15.0d0, tempk_0+tc, 11.94453d0, rg) ! from Kanzaki and Murakami 2015
        keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2) = &
            & k_arrhenius(10d0**(-10.43d0), tempk_0+15.0d0, tempk_0+tc, 17.00089d0, rg) ! from Kanzaki and Murakami 2015

        keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h0) = &
            & k_arrhenius(10d0**(1.770d0), tempk_0+25.0d0, tempk_0+tc, -8.170d0*cal2j, rg) ! from WATEQ4F.DAT 
        keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h1) = &
            & k_arrhenius(10d0**(-9.252d0), tempk_0+25.0d0, tempk_0+tc, 12.48d0*cal2j, rg) ! from WATEQ4F.DAT (NH4+ = NH3 + H+)

        keqgas_h(findloc(chrgas_all,'pn2o',dim=1),ieqgas_h0) = &
            & k_arrhenius(0.033928709d0, tempk_0+15.0d0, tempk_0+tc, -22.21661d0, rg) ! ! N2O solubility from Weiss & Price 1980 MC assuming 0 salinity

            
        keqaq_c     = 0d0
        keqaq_h     = 0d0
        keqaq_s     = 0d0
        keqaq_no3   = 0d0
        keqaq_nh3   = 0d0
        keqaq_oxa   = 0d0
        keqaq_cl    = 0d0

        ! SO4-2 + H+ = HSO4- 
        ! keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4) = &
            ! & k_arrhenius(10d0**(1.988d0),25d0+tempk_0,tc+tempk_0,3.85d0*cal2j,rg) ! from PHREEQC.DAT
        keqaq_h(findloc(chraq_all,'so4',dim=1),ieqaq_h1) = &
            & k_arrhenius(10d0**(1.988d0),25d0+tempk_0,tc+tempk_0,3.85d0*cal2j,rg) ! from PHREEQC.DAT
        ! SO4-2 + NH4+ = NH4SO4-
        keqaq_nh3(findloc(chraq_all,'so4',dim=1),ieqaq_nh3) = &
            & k_arrhenius(10d0**(1.03d0),25d0+tempk_0,tc+tempk_0,0d0,rg) ! from MINTEQV4.DAT 

        ! H+ + NO3- = HNO3 
        ! keqaq_no3(findloc(chraq_all,'no3',dim=1),ieqaq_no3) = 1d0/35.5d0 ! from Levanov et al. 2017 
        ! keqaq_no3(findloc(chraq_all,'no3',dim=1),ieqaq_no3) = 1d0/(10d0**1.3d0) ! from Maggi et al. 2007 
        keqaq_h(findloc(chraq_all,'no3',dim=1),ieqaq_h1) = 1d0/35.5d0 ! from Levanov et al. 2017 
        keqaq_h(findloc(chraq_all,'no3',dim=1),ieqaq_h1) = 1d0/(10d0**1.3d0) ! from Maggi et al. 2007 
        ! (temperature dependence is assumed to be 0) 

        ! 1.0000 H+ + 1.0000 Cl-  =  HCl
        keqaq_h(findloc(chraq_all,'cl',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(-0.67d0),25d0+tempk_0,tc+tempk_0,0d0,rg) ! from LLNL.DAT saying 'Not possible to calculate enthalpy of reaction'

        ! Oxa= + H+ = OxaH-
        ! keqaq_h(findloc(chraq_all,'oxa',dim=1),ieqaq_h1) = 1d0/(10d0**-4.266d0) ! from Lawrence et al., GCA, 2014
        ! Oxa= + 2H+ = OxaH2
        ! keqaq_h(findloc(chraq_all,'oxa',dim=1),ieqaq_h2) = 1d0/(10d0**-5.516d0) ! from Lawrence et al., GCA, 2014

        !  OxaH- = Oxa= + H+ 
        keqaq_h(findloc(chraq_all,'oxa',dim=1),ieqaq_h1) = (10d0**-4.266d0) ! from Lawrence et al., GCA, 2014
        !  OxaH- + H+ = OxaH2
        keqaq_h(findloc(chraq_all,'oxa',dim=1),ieqaq_h2) = 1d0/(10d0**-1.25d0) ! from Lawrence et al., GCA, 2014

        ! Sikora buffer (consts from Goldberg et al., 2002)
        #ifdef Goldberg_Sikora
        ! AcO- + H+ = AcOH  
        keqaq_h(findloc(chraq_all,'ac',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(4.756d0),25d0+tempk_0,tc+tempk_0,0.41d0,rg)  
            
        ! MESH- + H+ = MES  
        keqaq_h(findloc(chraq_all,'mes',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(6.270d0),25d0+tempk_0,tc+tempk_0,-14.8d0,rg)  
            
        ! Imidazole + H+ = Imidazole+  
        keqaq_h(findloc(chraq_all,'im',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(6.993d0),25d0+tempk_0,tc+tempk_0,-36.64d0,rg) 
            
        ! TEA + H+ = TEA+  
        keqaq_h(findloc(chraq_all,'tea',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(7.762d0),25d0+tempk_0,tc+tempk_0,-33.6d0,rg)  


        ! Sikora buffer (consts at 25oC from Sikora 2006 with temperature dependence remaining from Goldberg)
        #else
        ! AcO- + H+ = AcOH  
        keqaq_h(findloc(chraq_all,'ac',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(4.48d0),25d0+tempk_0,tc+tempk_0,0.41d0,rg)  
            
        ! MESH- + H+ = MES  
        keqaq_h(findloc(chraq_all,'mes',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(6.18d0),25d0+tempk_0,tc+tempk_0,-14.8d0,rg)  
            
        ! Imidazole + H+ = Imidazole+  
        keqaq_h(findloc(chraq_all,'im',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(7.10d0),25d0+tempk_0,tc+tempk_0,-36.64d0,rg) 
            
        ! TEA + H+ = TEA+  
        keqaq_h(findloc(chraq_all,'tea',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(8.09d0),25d0+tempk_0,tc+tempk_0,-33.6d0,rg)  
        #endif 
            
        ! Mehlich buffer (consts from Goldberg et al., 2002)

        !  GlpH- = Glp= + H+ 
        keqaq_h(findloc(chraq_all,'glp',dim=1),ieqaq_h1) = &! (10d0**-4.266d0)  &
            & k_arrhenius(10d0**(-6.650d0),25d0+tempk_0,tc+tempk_0,-1.85d0,rg) 
        !  GlpH- + H+ = GlpH2
        keqaq_h(findloc(chraq_all,'glp',dim=1),ieqaq_h2) = &! 1d0/(10d0**-1.25d0) 
            & k_arrhenius(10d0**(1.329d0),25d0+tempk_0,tc+tempk_0, 12.2d0,rg) 

        ! Mehlich buffer (consts from Hoskins and Erich 2008)

        ! ----------------

        ! Al3+ + H2O = Al(OH)2+ + H+
        keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h1) = &
            & k_arrhenius(10d0**(-5d0),25d0+tempk_0,tc+tempk_0,11.49d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al3+ + 2H2O = Al(OH)2+ + 2H+
        keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h2) = &
            & k_arrhenius(10d0**(-10.1d0),25d0+tempk_0,tc+tempk_0,26.90d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al3+ + 3H2O = Al(OH)3 + 3H+
        keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h3) = &
            & k_arrhenius(10d0**(-16.9d0),25d0+tempk_0,tc+tempk_0,39.89d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al3+ + 4H2O = Al(OH)4- + 4H+
        keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h4) = &
            & k_arrhenius(10d0**(-22.7d0),25d0+tempk_0,tc+tempk_0,42.30d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al3+ + SO4-2 = AlSO4+
        keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4) = &
            & k_arrhenius(10d0**(3.5d0),25d0+tempk_0,tc+tempk_0,2.29d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al3+ + 2SO4-2 = Al(SO4)2-
        ! ignoring for now
        ! keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42) = &
            ! & k_arrhenius(10d0**(5.0d0),25d0+tempk_0,tc+tempk_0,3.11d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al3+ + OxaH- = AlOxa+ + H+ (Al3+ + Oxa= = AlOxa+  plus OxaH- = Oxa= + H+ )
        ! keqaq_oxa(findloc(chraq_all,'al',dim=1),ieqaq_oxa) = 1d0/(10d0**-7.26d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        ! Al3+ + H2O + Oxa= = Al(OH)Oxa + H+ (dominant reaction according to Lawrence et al., GCA, 2014)
        ! Al3+ + H2O + HOxa- = Al(OH)Oxa + 2H+ <----> Al3+ + H2O + Oxa= = Al(OH)Oxa + H+  plus  OxaH- = Oxa= + H+
        keqaq_oxa(findloc(chraq_all,'al',dim=1),ieqaq_oxa) = 1d0/(10d0**-2.57d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        ! Al3+ + 2H2O + Oxa= = Al(OH)2Oxa- + 2H+ (dominant reaction according to Lawrence et al., GCA, 2014)
        ! Al3+ + 2H2O + HOxa- = Al(OH)2Oxa- + 3H+ <----> Al3+ + 2H2O + Oxa= = Al(OH)2Oxa + 2H+  plus  OxaH- = Oxa= + H+
        keqaq_oxa(findloc(chraq_all,'al',dim=1),ieqaq_oxa2) = 1d0/(10d0**3.12d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999

        ! H4SiO4 = H3SiO4- + H+
        keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h1) = &
            & k_arrhenius(10d0**(-9.83d0),25d0+tempk_0,tc+tempk_0,6.12d0*cal2j,rg) ! from PHREEQC.DAT 
        ! H4SiO4 = H2SiO4-2 + 2 H+
        keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h2) = &
            & k_arrhenius(10d0**(-23d0),25d0+tempk_0,tc+tempk_0,17.6d0*cal2j,rg) ! from PHREEQC.DAT 


        ! Mg2+ + H2O = Mg(OH)+ + H+
        keqaq_h(findloc(chraq_all,'mg',dim=1),ieqaq_h1) = &
            & k_arrhenius(10d0**(-11.44d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Mg2+ + CO32- = MgCO3 
        keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_co3) = &
            & k_arrhenius(10d0**(2.98d0),25d0+tempk_0,tc+tempk_0,2.713d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Mg2+ + H+ + CO32- = MgHCO3
        keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_hco3) = & 
            & k_arrhenius(10d0**(11.399d0),25d0+tempk_0,tc+tempk_0,-2.771d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Mg+2 + SO4-2 = MgSO4
        keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4) = & 
            & k_arrhenius(10d0**(2.37d0),25d0+tempk_0,tc+tempk_0, 4.550d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Mg2+ + OxaH- = MgOxa + H+ (Mg2+ + Oxa= = MgOxa  plus OxaH- = Oxa= + H+ )
        keqaq_oxa(findloc(chraq_all,'mg',dim=1),ieqaq_oxa) = 1d0/(10d0**-3.43d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        ! 1.0000 Mg++ + 1.0000 Cl-  =  MgCl+
        keqaq_cl(findloc(chraq_all,'mg',dim=1),ieqaq_cl) =  &
            & k_arrhenius(10d0**(-0.1349d0),25d0+tempk_0,tc+tempk_0,-0.58576d0,rg) ! from LLNL.DAT 


        ! Ca2+ + H2O = Ca(OH)+ + H+
        keqaq_h(findloc(chraq_all,'ca',dim=1),ieqaq_h1) =  &
            & k_arrhenius(10d0**(-12.78d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ! from PHREEQC.DAT 
        ! (No delta_h is reported so used the same value for Mg)
        ! Ca2+ + CO32- = CaCO3 
        keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_co3) = &
            & k_arrhenius(10d0**(3.224d0),25d0+tempk_0,tc+tempk_0,3.545d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Ca2+ + H+ + CO32- = CaHCO3
        keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_hco3) = &
            & k_arrhenius(10d0**(11.435d0),25d0+tempk_0,tc+tempk_0,-0.871d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Ca+2 + SO4-2 = CaSO4
        keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4) = &
            & k_arrhenius(10d0**(2.25d0),25d0+tempk_0,tc+tempk_0,1.325d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Ca+2 + NO3- = CaNO3+
        keqaq_no3(findloc(chraq_all,'ca',dim=1),ieqaq_no3) = &
            & k_arrhenius(10d0**(0.5d0),25d0+tempk_0,tc+tempk_0,-5.4d0,rg) ! from MINTEQV4.DAT           
        ! Ca+2 + NH4+ = CaNH3+2 + H+
        keqaq_nh3(findloc(chraq_all,'ca',dim=1),ieqaq_nh3) = &
            & k_arrhenius(10d0**(-9.144d0),25d0+tempk_0,tc+tempk_0,0d0,rg) ! from MINTEQV4.DAT 
        ! Ca+2 + 2NH4+ = Ca(NH3)2+2 + 2H+
        ! ignoring for now
        ! keqaq_nh3(findloc(chraq_all,'ca',dim=1),ieqaq_nh32) = &
            ! & k_arrhenius(10d0**(-18.788d0),25d0+tempk_0,tc+tempk_0,0d0,rg) ! from MINTEQV4.DAT 
        ! Ca2+ + Oxa= = CaOxa
        ! keqaq_oxa(findloc(chraq_all,'ca',dim=1),ieqaq_oxa) = 1d0/(10d0**-3.19d0) ! from Prapaipong et al., GCA, 1999
        ! Ca2+ + OxaH- = CaOxa + H+ (Ca2+ + Oxa= = CaOxa  plus OxaH- = Oxa= + H+ )
        keqaq_oxa(findloc(chraq_all,'ca',dim=1),ieqaq_oxa) = 1d0/(10d0**-3.19d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        ! 1.0000 Cl- + 1.0000 Ca++  =  CaCl+
        keqaq_cl(findloc(chraq_all,'ca',dim=1),ieqaq_cl) =  &
            & k_arrhenius(10d0**(-0.6956d0),25d0+tempk_0,tc+tempk_0,2.02087d0,rg) ! from LLNL.DAT 

            
        ! Fe2+ + H2O = Fe(OH)+ + H+
        keqaq_h(findloc(chraq_all,'fe2',dim=1),ieqaq_h1) = &
            & k_arrhenius(10d0**(-9.51d0),25d0+tempk_0,tc+tempk_0, 40.3d0,rg) ! from Kanzaki and Murakami 2016
        ! Fe2+ + CO32- = FeCO3 
        keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_co3) = &
            & k_arrhenius(10d0**(5.69d0),25d0+tempk_0,tc+tempk_0, -45.6d0,rg) ! from Kanzaki and Murakami 2016
        ! Fe2+ + H+ + CO32- = FeHCO3
        keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_hco3) = &
            & k_arrhenius(10d0**(1.47d0),25d0+tempk_0,tc+tempk_0, -18d0,rg) &! from Kanzaki and Murakami 2016 
            & /keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2) 
        ! Fe+2 + SO4-2 = FeSO4
        keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4) = &
            & k_arrhenius(10d0**(2.25d0),25d0+tempk_0,tc+tempk_0,3.230d0*cal2j,rg) ! from PHREEQC.DAT 
        ! 1.0000 Fe++ + 1.0000 Cl-  =  FeCl+
        keqaq_cl(findloc(chraq_all,'fe2',dim=1),ieqaq_cl) =  &
            & k_arrhenius(10d0**(-0.1605d0),25d0+tempk_0,tc+tempk_0,3.02503d0,rg) ! from LLNL.DAT 


        ! Fe3+ + H2O = Fe(OH)2+ + H+
        keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h1) = &
            & k_arrhenius(10d0**(-2.19d0),25d0+tempk_0,tc+tempk_0,10.4d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Fe3+ + 2H2O = Fe(OH)2+ + 2H+
        keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h2) = &
            & k_arrhenius(10d0**(-5.67d0),25d0+tempk_0,tc+tempk_0,17.1d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Fe3+ + 3H2O = Fe(OH)3 + 3H+
        keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h3) = &
            & k_arrhenius(10d0**(-12.56d0),25d0+tempk_0,tc+tempk_0,24.8d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Fe3+ + 4H2O = Fe(OH)4- + 4H+
        keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h4) = &
            & k_arrhenius(10d0**(-21.6d0),25d0+tempk_0,tc+tempk_0,31.9d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Fe+3 + SO4-2 = FeSO4+
        keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4) = &
            & k_arrhenius(10d0**(4.04d0),25d0+tempk_0,tc+tempk_0,3.91d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Fe+3 + 2 SO4-2 = Fe(SO4)2-
        ! ignoring for now
        ! keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42) = &
            ! & k_arrhenius(10d0**(5.38d0),25d0+tempk_0,tc+tempk_0,4.60d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Fe+3 + NO3- = FeNO3+2
        keqaq_no3(findloc(chraq_all,'fe3',dim=1),ieqaq_no3) = &
            & k_arrhenius(10d0**(1d0),25d0+tempk_0,tc+tempk_0,-37d0,rg) ! from MINTEQV4.DAT 
        ! Fe+3 + OxaH- = FeOxa+ + H+ (Fe+3 + Oxa= = FeOxa+  plus OxaH- = Oxa= + H+ )
        ! keqaq_oxa(findloc(chraq_all,'fe3',dim=1),ieqaq_oxa) = 1d0/(10d0**-9.33d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        keqaq_oxa(findloc(chraq_all,'fe3',dim=1),ieqaq_oxa) = 1d0/(10d0**-9.15d0)*(10d0**-4.266d0) ! from Perez-Fodich and Derry GCA, 2019
        ! Fe+3 + 2 OxaH- = Fe(Oxa)2- + 2 H+ (Fe+3 + 2Oxa= = Fe(Oxa)2-  plus 2  {OxaH- = Oxa= + H+} )
        ! keqaq_oxa(findloc(chraq_all,'fe3',dim=1),ieqaq_oxa) = 1d0/(10d0**-15.45d0)*(10d0**(-4.266d0*2d0)) ! from Perez-Fodich and Derry, 2019
        ! 1.0000 Fe+++ + 1.0000 Cl-  =  FeCl++
        keqaq_cl(findloc(chraq_all,'fe3',dim=1),ieqaq_cl) =  &
            & k_arrhenius(10d0**(-0.8108d0),25d0+tempk_0,tc+tempk_0,36.6421d0,rg) ! from LLNL.DAT 



        ! Na+ + CO3-2 = NaCO3-
        keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3) = & 
            & k_arrhenius(10d0**(1.27d0),25d0+tempk_0,tc+tempk_0, 8.91d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Na+ + H + CO3- = NaHCO3
        keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3) = & 
            & k_arrhenius(10d0**(-0.25d0),25d0+tempk_0,tc+tempk_0, -1d0*cal2j,rg) &! from PHREEQC.DAT for Na+ + HCO3- = NaHCO3
            & /keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)  ! HCO3- = CO32- + H+
        ! Na+ + SO4-2 = NaSO4-
        keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4) = & 
            & k_arrhenius(10d0**(0.7d0),25d0+tempk_0,tc+tempk_0, 1.120d0*cal2j,rg) ! from PHREEQC.DAT 
        ! +1.000Na+     +1.000NO3-          = Na(NO3)
        keqaq_no3(findloc(chraq_all,'na',dim=1),ieqaq_no3) = & 
            & k_arrhenius(10d0**(0.4d0),25d0+tempk_0,tc+tempk_0, 0d0,rg) ! from SIT.DAT (no enthalpy data)
        ! Na+ + Oxa= = NaOxa-
        ! keqaq_oxa(findloc(chraq_all,'na',dim=1),ieqaq_oxa) = 1d0/(10d0**-0.86d0) ! from Prapaipong et al., GCA, 1999
        ! Na+ + OxaH- = NaOxa- + H+ (Na+ + Oxa= = NaOxa-  plus OxaH- = Oxa= + H+ )
        keqaq_oxa(findloc(chraq_all,'na',dim=1),ieqaq_oxa) = 1d0/(10d0**-0.86d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        ! 1.0000 Na+ + 1.0000 Cl-  =  NaCl
        keqaq_cl(findloc(chraq_all,'na',dim=1),ieqaq_cl) =  &
            & k_arrhenius(10d0**(-0.777d0),25d0+tempk_0,tc+tempk_0,5.21326d0,rg) ! from LLNL.DAT 



        ! K+ + SO4-2 = KSO4-
        keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4) = & 
            & k_arrhenius(10d0**(0.85d0),25d0+tempk_0,tc+tempk_0, 2.250d0*cal2j,rg) ! from PHREEQC.DAT 
        ! +1.000K+     +1.000NO3-          = K(NO3)
        keqaq_no3(findloc(chraq_all,'k',dim=1),ieqaq_no3) = & 
            & k_arrhenius(10d0**(-0.15d0),25d0+tempk_0,tc+tempk_0, 0d0,rg) ! from SIT.DAT (no enthalpy data)
        ! K+ + OxaH- = KOxa- + H+ (K+ + Oxa= = KOxa-  plus OxaH- = Oxa= + H+ )
        keqaq_oxa(findloc(chraq_all,'k',dim=1),ieqaq_oxa) = 1d0/(10d0**-0.80d0)*(10d0**-4.266d0) ! from Prapaipong et al., GCA, 1999
        ! 1.0000 K+ + 1.0000 Cl-  =  KCl
        keqaq_cl(findloc(chraq_all,'k',dim=1),ieqaq_cl) =  &
            & k_arrhenius(10d0**(-1.4946d0),25d0+tempk_0,tc+tempk_0,14.1963d0,rg) ! from LLNL.DAT 


        ! keqaq_s = 0d0

        !!! ----------- Ion exchange ------------------------!!
        keqcec_all = 0d0
        keqiex_all = 0d0

        call get_base_charge( &
            & nsp_aq_all & 
            & ,chraq_all & 
            & ,base_charge &! output 
            & )

        do isps = 1, nsp_sld_all
            ! exchange capacity in (eq)mol / mineral mol 
            
            mwt_tmp = mwt_all(isps)
            zero_cec = .false.
            
            if ( mcec_all(isps) >= mcec_threshold ) then 
                keqcec_all(isps) = (    &
                    & mcec_all(isps)    &! cmol/kg
                    & *1d-2/1d3         &! mol/g
                    & *mwt_tmp          &! mol/mol with molar weight g/mol
                    & )
            else 
                zero_cec = .true.
            endif 
            
            if (.not.zero_cec) then
                ! half exchange reaction (cf. Turner et al. GCA 1996 Sect 4.1) 
                ! Na+(aq) + X- = NaX 
                ! keqiex_all(findloc(chrsld_all,'ka',dim=1),findloc(chraq_all,'na',dim=1)) = 1d10
                ! 
                keqiex_prona = 5.9d0 ! log KH\Na default from Appelo 1994
                ! keqiex_prona = 6.7d0 ! log KH\Na random test
                keqiex_prona = 8.0d0 ! log KH\Na random test
                keqiex_prona = 7.8d0 ! log KH\Na random test
                keqiex_prona = 7.0d0 ! log KH\Na random test
                ! keqiex_prona = 7.2d0 ! log KH\Na random test
                ! keqiex_prona = 6.5d0 ! log KH\Na random test
                ! keqiex_prona = 6.0d0 ! log KH\Na random test
                keqiex_prona = 5.5d0 ! log KH\Na random test
                ! keqiex_prona = keqiex_prona - keqiex_prona *0.5d0* mcec_all(isps)/100d0 ! assume log KH\Na decreases with CEC
                if (mcec_all(isps) < 100d0) then
                    ! keqiex_prona = keqiex_prona + keqiex_prona *0.5d0* mcec_all(isps)/100d0 ! assume log KH\Na decreases with CEC
                    ! keqiex_prona = keqiex_prona + keqiex_prona *1.0d0* mcec_all(isps)/100d0 ! assume log KH\Na decreases with CEC
                    keqiex_prona = keqiex_prona + keqiex_prona *1.5d0* mcec_all(isps)/100d0 ! assume log KH\Na decreases with CEC
                    ! keqiex_prona = keqiex_prona + keqiex_prona *0.5d0* mcec_all(isps)/200d0 ! assume log KH\Na decreases with CEC
                endif 
                
                ! keqiex_prona = logkhna_all(isps)
                
                ! 
                do ispa=1,nsp_aq_all
                    aqsp = chraq_all(ispa)
                    selectcase(trim(adjustl(aqsp)))
                        
                        case('na')
                            if (cec_pH_depend(isps)) then
                                ! Na+(aq) + X-H = X-Na + H+
                                ! K = [X-Na]/CEC * [H+] / ([Na+] * [X-H]/CEC)
                                ! [X-Na]/CEC = K * [Na+] * [X-H]/CEC / [H+]
                                ! [X-Na] = K * [Na+] * [X-H] / [H+]
                                keqiex_all(isps,ispa) = 10d0**(-5.883d0)
                                ! from log KH\Na = 5.883 in Appeolo (1994)
                                ! keqiex_all(isps,ispa) = 10d0**(-5.9d0)
                                ! keqiex_all(isps,ispa) = 10d0**(-keqiex_prona)
                                keqiex_all(isps,ispa) = 10d0**(-logkhaq_all(isps,ispa))
                                ! from log KH\Na = 5.9 - 3.4 * f[X-H] in Appeolo (1994)
                                ! log KNa\H = 3.4 * f[X-H] - 5.9 
                                ! the 'activity coefficient' term 10**(3.4 * f[X-H]) will be added when calculating f[X-H]
                                ! dependence on f[X-Na] is not accounted for 
                                ! and parameterization adopted for freshwater by Appelo (1994) is used
                                ! To include the dependency on f[X-Na], two equations must be solved instead of one 
                            else 
                                keqiex_all(isps,ispa) = 10d0**(0.00d0)
                            endif 
                        
                        case('k')
                            if (cec_pH_depend(isps)) then
                                ! K+(aq) + X-H = X-K + H+
                                ! K = [X-K]/CEC * [H+] / ([K+] * [X-H]/CEC)
                                ! [X-K]/CEC = K * [K+] * [X-H]/CEC / [H+]
                                ! [X-K] = K * [K+] * [X-H] / [H+]
                                keqiex_all(isps,ispa) = 10d0**(-4.783d0)
                                ! from log KK\Na = 1.10 and log KH\Na = 5.883 in Appeolo (1994)
                                ! K+(aq)  + X-H  = X-K  + H+    (log K = X      ) ----- (A)
                                ! K+(aq)  + X-Na = X-K  + Na+   (log K = 1.10   ) ----- (B)
                                ! Na+(aq) + X-H  = X-Na + H+    (log K = -5.883 ) ----- (C)
                                ! (A) = (B) + (C)
                                ! X = 1.10 - 5.883 = -4.783
                                
                                ! case adopting log KK\Na = 0.902 and log KH\Na = 5.9 - 3.4 * f[X-H]
                                ! X = 0.902 - 5.9 = -4.998
                                keqiex_all(isps,ispa) = 10d0**(-4.998d0)
                                ! the 'activity coefficient' term 10**(3.4 * f[X-H]) will be added when calculating f[X-H]
                                ! according to PHREEQC.DAT log KK\Na = 0.7 and X = 0.7 - 5.9 = -5.2
                                keqiex_all(isps,ispa) = 10d0**(-5.2d0)
                                keqiex_all(isps,ispa) = 10d0**(-5.9d0) ! assuming Na
                                keqiex_all(isps,ispa) = 10d0**(-6.9d0) ! assuming Na in seawater
                                keqiex_all(isps,ispa) = 10d0**(-4.8d0) ! log KK\Na = 1.10 
                                keqiex_all(isps,ispa) = 10d0**(1.1d0 - keqiex_prona)  ! log KK\Na = 1.10 
                                keqiex_all(isps,ispa) = 10d0**(-logkhaq_all(isps,ispa))
                                ! keqiex_all(isps,ispa) = 10d0**(0.9d0 - keqiex_prona)  ! log KK\Na = 0.90 
                                ! keqiex_all(isps,ispa) = 10d0**(0.7d0 - keqiex_prona)  ! log KK\Na = 0.70 From phreeqc.dat
                            else 
                                keqiex_all(isps,ispa) = 10d0**(1.10d0)
                                keqiex_all(isps,ispa) = 10d0**(0.902d0)
                            endif 
                        
                        case('mg')
                            if (cec_pH_depend(isps)) then
                                ! Mg++(aq) + 2X-H = X2-Mg + 2H+
                                ! K = (2[X2-Mg]/CEC) * [H+]2 / ([Mg++] * ([X-H]/CEC)^2 )
                                ! 2[X2-Mg]/CEC = K * [Mg++] * ([X-H]/CEC)^2 / [H+]^2
                                ! [X2-Mg] = CEC/2 * K * [Mg++] * ([X-H]/CEC)^2 / [H+]^2
                                ! [X2-Mg] = 1/(2*CEC) * K * [Mg++] * [X-H]^2 / [H+]^2
                                keqiex_all(isps,ispa) = 10d0**(-10.752d0)
                                ! from log KMg\Na = 0.507 and log KH\Na = 5.883 in Appeolo (1994)
                                ! Mg++(aq)  + 2 X-H  = X2-Mg  + 2 H+    (log K = X       ) ----- (A)
                                ! Mg++(aq)  + 2 X-Na = X2-Mg  + 2 Na+   (log K =  0.507*2) ----- (B)
                                ! 2 Na+(aq) + 2 X-H  = 2 X-Na + 2 H+    (log K = -5.883*2) ----- (C)
                                ! (A) = (B) + (C)
                                ! X = 0.507* 2 - 5.883 * 2 = -10.752
                                
                                ! case adopting log KMg\Na = 0.307 and log KH\Na = 5.9 - 3.4 * f[X-H]
                                ! X = 0.307*2 - 5.9*2 = -11.186
                                keqiex_all(isps,ispa) = 10d0**(-11.186d0)
                                ! the 'activity coefficient' term 10**(2*3.4 * f[X-H]) will be added when calculating f[X-H]
                                keqiex_all(isps,ispa) = 10d0**(-10.786d0)  ! log KMg\Na = 0.507 
                                keqiex_all(isps,ispa) = 10d0**( (0.507d0 - keqiex_prona)*2d0)  ! log KMg\Na = 0.507 
                                keqiex_all(isps,ispa) = 10d0**(-logkhaq_all(isps,ispa))
                            else 
                                keqiex_all(isps,ispa) = 10d0**(1.014d0)
                                keqiex_all(isps,ispa) = 10d0**(0.614d0)
                            endif 
                        
                        case('ca')
                            if (cec_pH_depend(isps)) then
                                ! Ca++(aq) + 2X-H = X2-Ca + 2H+
                                ! K = (2[X2-Ca]/CEC) * [H+]2 / ([Ca++] * ([X-H]/CEC)^2 )
                                ! 2[X2-Ca]/CEC = K * [Ca++] * ([X-H]/CEC)^2 / [H+]^2
                                ! [X2-Ca] = CEC/2 * K * [Ca++] * ([X-H]/CEC)^2 / [H+]^2
                                ! [X2-Ca] = CEC/2 * K * [Ca++] * ([X-H]/CEC)^2 / [H+]^2
                                keqiex_all(isps,ispa) = 10d0**(-10.436d0)
                                ! from log KCa\Na = 0.665 and log KH\Na = 5.883 in Appeolo (1994)
                                ! Ca++(aq)  + 2 X-H  = X2-Ca  + 2 H+    (log K = X       ) ----- (A)
                                ! Ca++(aq)  + 2 X-Na = X2-Ca  + 2 Na+   (log K =  0.665*2) ----- (B)
                                ! 2 Na+(aq) + 2 X-H  = 2 X-Na + 2 H+    (log K = -5.883*2) ----- (C)
                                ! (A) = (B) + (C)
                                ! X = 0.507* 2 - 5.883 * 2 = -10.436
                                
                                ! case adopting log KCa\Na = 0.465 and log KH\Na = 5.9 - 3.4 * f[X-H]
                                ! X = 0.465*2 - 5.9*2 = -10.87
                                keqiex_all(isps,ispa) = 10d0**(-10.87d0)
                                ! the 'activity coefficient' term 10**(2*3.4 * f[X-H]) will be added when calculating f[X-H]
                                keqiex_all(isps,ispa) = 10d0**(-10.47d0) ! log KCa\Na = 0.665  
                                keqiex_all(isps,ispa) = 10d0**( (0.665d0 - keqiex_prona)*2d0)  ! log KCa\Na = 0.665  
                                ! keqiex_all(isps,ispa) = 10d0**( (0.0d0 - keqiex_prona)*2d0)  ! log KCa\Na = 0.0 (just for testing)  
                                ! keqiex_all(isps,ispa) = 10d0**( (0.4d0 - keqiex_prona)*2d0)  ! log KCa\Na = 0.4 From phreeqc.dat 
                                keqiex_all(isps,ispa) = 10d0**(-logkhaq_all(isps,ispa)) 
                            else 
                                keqiex_all(isps,ispa) = 10d0**(1.33d0)
                                keqiex_all(isps,ispa) = 10d0**(0.93d0)
                            endif 
                        
                        case('al')
                        ! case('XXXXX')
                            if (cec_pH_depend(isps)) then
                                ! Al+++(aq) + 3X-H = X3-Al + 3H+
                                ! K = (3[X3-Al]/CEC) * [H+]3 / ([Al++] * ([X-H]/CEC)^3 )
                                ! 3[X3-Al]/CEC = K * [Al+++] * ([X-H]/CEC)^3 / [H+]^3
                                ! [X3-Al] = CEC/3 * K * [Al+++] * ([X-H]/CEC)^3 / [H+]^3
                                ! [X3-Al] = CEC/3 * K * [Al+++] * ([X-H]/CEC)^3 / [H+]^3
                                keqiex_all(isps,ispa) = 10d0**(-16.08d0)
                                ! from KNa\Al = 0.3 (log KAl\Na = 0.523) from Tang et al. 2013 and log KH\Na = 5.883 in Appeolo (1994)
                                ! Al+++(aq)  + 3 X-H  = X3-Al  + 3 H+    (log K = X       ) ----- (A)
                                ! Al+++(aq)  + 3 X-Na = X3-Al  + 3 Na+   (log K =  0.523*3) ----- (B)
                                ! 3 Na+(aq)  + 3 X-H  = 3 X-Na + 3 H+    (log K = -5.883*3) ----- (C)
                                ! (A) = (B) + (C)
                                ! X = 0.523* 3 - 5.883 * 3 = -16.08
                                
                                ! case adopting log KAl\Na = 0.41 (PHREEQC.DAT) and log KH\Na = 5.9 - 3.4 * f[X-H]
                                ! X = 0.41 - 5.9*3 = -17.29
                                keqiex_all(isps,ispa) = 10d0**(-17.29d0)
                                ! the 'activity coefficient' term 10**(3*3.4 * f[X-H]) will be added when calculating f[X-H]
                                keqiex_all(isps,ispa) = 10d0**( (0.41d0 - keqiex_prona)*3d0 )  ! log KAl\Na = 0.41 
                                keqiex_all(isps,ispa) = 10d0**(-logkhaq_all(isps,ispa)) 
                            else 
                                ! keqiex_all(isps,ispa) = 10d0**(1.569d0)
                                keqiex_all(isps,ispa) = 10d0**(0.41d0)
                            endif 
                        
                        case default 
                            ! do nothing 
                            continue 
                            
                    endselect
                enddo
            endif 
        enddo 

        !!! ----------- Solid phases ------------------------!!
        ksld_all = 0d0 
        keqsld_all = 0d0

        ! call get_mgasx_all( &
            ! & nz,nsp_gas_all,nsp_gas,nsp_gas_cnst &
            ! & ,chrgas,chrgas_all,chrgas_cnst &
            ! & ,mgas,mgasc &
            ! & ,mgas_loc  &! output
            ! & )

        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maq,mgas,maqc,mgasc &
            & ,maqf_loc,mgas_loc  &! output
            & )

        do isps = 1, nsp_sld_all
            mv_tmp = mv_all(isps)
            mineral = chrsld_all(isps)
            
            call sld_kin( &
                & nz,rg,tc,sec2yr,tempk_0,pro,kw,kho,mv_tmp &! input
                & ,nsp_gas_all,chrgas_all,mgas_loc &! input
                & ,nsp_aq_all,chraq_all,maqf_loc &! input
                & ,mineral,'xxxxx' &! input 
                & ,kin,dkin_dmsp &! output
                & ) 
            ksld_all(isps,:) = kin
            
            ! check for solid solution 
            select case (trim(adjustl(mineral))) 
                case('la','ab','an','by','olg','and')
                    ss_x = staq_all(isps, findloc(chraq_all,'ca',dim=1))
                    ss_y = 0d0 
                    ss_z = 0d0 
                case('cpx','hb','dp') 
                    ss_x = staq_all(isps, findloc(chraq_all,'fe2',dim=1))
                    ss_y = 0d0 
                    ss_z = 0d0 
                case('opx','en','fer') 
                    ss_x = staq_all(isps, findloc(chraq_all,'fe2',dim=1))
                    ss_y = 0d0 
                    ss_z = 0d0 
                case('agt') 
                    ss_z = staq_all(isps, findloc(chraq_all,'na',dim=1))
                    ss_y = 1d0 - staq_all(isps, findloc(chraq_all,'ca',dim=1))/(1d0 - ss_z)
                    ss_x = staq_all(isps, findloc(chraq_all,'fe2',dim=1))/(1d0+ ss_y )/(1d0 -ss_z)
                case default 
                    ss_x = 0d0 
                    ss_y = 0d0 
                    ss_z = 0d0 
            endselect 
            
            select case(trim(adjustl(mineral))) 
                case('gbas','cbas') 
                    ! doing rather complicated solid solution though simplified
                    ! following Pollyea and Rimstidt 2017; Aradóttir et al. 2012
                    ss_x = 0d0 
                    ss_y = 0d0 
                    ss_z = 0d0 
                    
                    therm = 0d0
                    
                    do ispss=1,7
                        if (trim(adjustl(mineral)) =='gbas') then
                            ssaq = chrss_gbas_aq(ispss)
                            sssld = chrss_gbas_sld(ispss)
                        elseif (trim(adjustl(mineral)) =='cbas') then
                            ssaq = chrss_cbas_aq(ispss)
                            sssld = chrss_cbas_sld(ispss)
                        endif 
                        call sld_therm( &
                            & rg,tc,tempk_0,ss_x,ss_y,ss_z &! input
                            & ,sssld &! input
                            & ,therm_tmp &! output
                            & ) 
                        ! required oxide mole per basalt mole
                        ss_tmp = staq_all(isps, findloc(chraq_all,trim(adjustl(ssaq)),dim=1)) &! contained cation per basalt
                            & /staq_all(findloc(chrsld_all,trim(adjustl(sssld)),dim=1), findloc(chraq_all,trim(adjustl(ssaq)),dim=1)) ! divided by cation per oxide
                        therm = therm + ss_tmp*log(therm_tmp) + ss_tmp*log(ss_tmp) 
                    enddo
                    
                    therm = exp(therm)
                    
                case default
                    call sld_therm( &
                        & rg,tc,tempk_0,ss_x,ss_y,ss_z &! input
                        & ,mineral &! input
                        & ,therm &! output
                        & ) 
            endselect
            
            ! correction of thermodynamic data wrt primary species
            select case (trim(adjustl(mineral))) 
                case('anl')
                    ! replacing Al(OH)42- with Al+++ as primary Al species
                    therm = therm/keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h4)
                case default 
                    ! do nothing
            endselect 
            
            keqsld_all(isps) = therm
        enddo


        !--------- other reactions -------------! 
        krxn1_ext_all = 0d0
        krxn2_ext_all = 0d0

        krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) = &
            & 1d0
            ! & max(8.0d13*60.0d0*24.0d0*365.0d0*(kw/pro)**2.0d0, 1d-7*60.0d0*24.0d0*365.0d0) &   
            ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)excluding the term (c*po2)
            ! & *merge(0d0,1d0,po2<po2th*thon)
            
        krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = 0.71d0 ! vmax mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
        krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = &
            & krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) !*1d1 ! reducing a bit to be fitted with modern soil pco2

        krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
            
            
        krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:) = 0.01d0*24d0*365d0 ! mg C mg-1 MBC yr-1
        ! converted from 0.01 mg C mg-1 MBC hr-1 Georgiou et al. (2017)

        krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:) = 250d0 ! mg C g-1 soil  Georgiou et al. (2017)
            
            
        krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:) = 0.00028d0*24d0*365d0 ! mg C mg-1 MBC yr-1
        ! converted from 0.00028 mg C mg-1 MBC hr-1 Georgiou et al. (2017)

        krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:) = 2d0 ! beta value Georgiou et al. (2017)




        krxn1_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1),:) = & 
            & 10.0d0**(-6.07d0)*60.0d0*60.0d0*24.0d0*365.0d0  !! excluding the term (fe3**0.93/fe2**0.40)  
            
        krxn1_ext_all(findloc(chrrxn_ext_all,'g2k',dim=1),:) = ksld_all(findloc(chrsld_all,'g2',dim=1),:)
        krxn1_ext_all(findloc(chrrxn_ext_all,'g2ca',dim=1),:) = ksld_all(findloc(chrsld_all,'g2',dim=1),:)
        krxn1_ext_all(findloc(chrrxn_ext_all,'g2mg',dim=1),:) = ksld_all(findloc(chrsld_all,'g2',dim=1),:)


    endsubroutine coefs_v2

    
endmodule scepter_eq_coefs
