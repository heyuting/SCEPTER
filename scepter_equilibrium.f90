module scepter_equilibrium
    implicit none
    private
    public :: coefs_v2, calc_pH_v7_4, calc_charge_balance, calc_charge_balance_point

    ! Constants
    real(kind=8), parameter :: cal2j = 4.184d0

contains

    ! Main coefficient calculation subroutine
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

        ! real(kind=8) k_arrhenius
        real(kind=8) :: cal2j = 4.184d0 

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
        real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_s
        real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_no3
        real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_nh3
        real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_oxa
        real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_cl
        real(kind=8),dimension(nsp_sld_all,nz),intent(out)::ksld_all
        real(kind=8),dimension(nsp_sld_all),intent(in)::mv_all,mwt_all,mcec_all
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::logkhaq_all
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::staq_all
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(out)::keqiex_all
        real(kind=8),dimension(nsp_sld_all),intent(out)::keqsld_all,keqcec_all
        real(kind=8),dimension(nrxn_ext_all,nz),intent(out)::krxn1_ext_all
        real(kind=8),dimension(nrxn_ext_all,nz),intent(out)::krxn2_ext_all

        integer,intent(in)::nsp_gas,nsp_gas_cnst,nsp_aq,nsp_aq_cnst
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maq
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all

        logical,dimension(nsp_sld_all),intent(in)::cec_pH_depend

        real(kind=8),dimension(nsp_gas_all,nz)::mgas_loc
        real(kind=8),dimension(nsp_aq_all,nz)::maqf_loc
        real(kind=8),dimension(nsp_aq_all)::base_charge

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

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
        ! Al+3 + SO4-2 = AlSO4+
        keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4) = &
            & k_arrhenius(10d0**(3.5d0),25d0+tempk_0,tc+tempk_0,2.29d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al+3 + 2SO4-2 = Al(SO4)2-
        ! ignoring for now
        ! keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42) = &
            ! & k_arrhenius(10d0**(5.0d0),25d0+tempk_0,tc+tempk_0,3.11d0*cal2j,rg) ! from PHREEQC.DAT 
        ! Al+3 + OxaH- = AlOxa+ + H+ (Al+3 + Oxa= = AlOxa+  plus OxaH- = Oxa= + H+ )
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

    ! pH calculation subroutine
    subroutine calc_pH_v7_4( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input
        & ,poro,sat,tc &! input  
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 &! input
        & ,keqaq_oxa,keqaq_cl &! input
        & ,print_cb,print_loc,z,act_ON &! input 
        & ,dprodmaq_all,dprodmgas_all &! output
        & ,iosx,diosdmaq_all,diosdmgas_all &! output
        & ,prox,ph_error,ph_iter &! output
        & ) 
        ! solving charge balance with specific primary variables input; 
        ! here maqx is assumed to be concs. of free cations or H4SiO4 or SO42- or NO3-  
        ! gases are already treated with specific gas form.  
        implicit none
        integer,intent(in)::nz
        real(kind=8),intent(in)::kw,tc
        real(kind=8) so4th
        real(kind=8),dimension(nz)::so4x,prox_save,error_save,prox_save_newton,prox_save_bisec,prox_init
        real(kind=8),dimension(nz)::iosx_save,ios_new
        real(kind=8),dimension(nz),intent(in)::z,poro,sat
        real(kind=8),dimension(nz),intent(inout)::prox
        logical,intent(out)::ph_error

        real(kind=8),dimension(nz)::prox_max,prox_min,ph_add_order,prox_tmp1,prox_tmp2
        real(kind=8),dimension(nz)::f1_max,f1_min
        real(kind=8),dimension(nz)::df1,f1,f2,df2,df21,df12,d2f1
        real(kind=8),dimension(nz),intent(inout)::iosx
        real(kind=8) k_order,ph_inflex,a_order,c_order
        real(kind=8) error,tol,dconc 
        integer iter,iz,ispa,ispg

        integer,intent(in)::nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_nh3
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_no3
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_oxa
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_cl
        real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all

        real(kind=8),dimension(nsp_aq_all)::base_charge
        real(kind=8),dimension(nsp_aq_all,nz)::maqx_loc,maqf_loc
        real(kind=8),dimension(nsp_aq_all,nz)::dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2
        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc
        real(kind=8),dimension(nsp_aq_all,nz)::df1dmaq,df2dmaq,df1dmaqf,d2f1dmaqf
        real(kind=8),dimension(nsp_gas_all,nz)::df1dmgas,df2dmgas,d2f1dmgas

        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::dprodmaq_all
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::dprodmgas_all

        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::diosdmaq_all
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::diosdmgas_all

        real(kind=8),dimension(nsp_aq_all,nz)::dmaq,maqtmp_loc
        real(kind=8),dimension(nsp_gas_all,nz)::dmgas,mgastmp_loc
        real(kind=8),dimension(nz)::df1_dum,f1_dum,d2f1_dum,fact,f1_tmp,df1_tmp,d2f1_tmp
        real(kind=8),dimension(nz)::f1_tmp1,df1_tmp1,d2f1_tmp1
        real(kind=8),dimension(nz)::f1_tmp2,df1_tmp2,d2f1_tmp2
        real(kind=8),dimension(nsp_aq_all,nz)::df1dmaqf_dum,d2f1dmaqf_dum,df1dmaqf_tmp,d2f1dmaqf_tmp
        real(kind=8),dimension(nsp_gas_all,nz)::df1dmgas_dum,d2f1dmgas_dum,df1dmgas_tmp,d2f1dmgas_tmp
        real(kind=8),dimension(nz)::d2f2,df1df2,df2df1
        real(kind=8),dimension(nsp_aq_all,nz)::df2dmaqf,d2f2dmaqf
        real(kind=8),dimension(nsp_gas_all,nz)::d2f2dmgas

        real(kind=8),dimension(nsp_aq_all,nz)::dmaqft_dpro_loc,maqft_loc,dmaqft_dios_loc
        real(kind=8),dimension(nsp_aq_all,nsp_aq_all,nz)::dmaqft_dmaqf_loc
        real(kind=8),dimension(nsp_aq_all,nsp_gas_all,nz)::dmaqft_dmgas_loc

        integer iso4,iph,iph2,iph3
        integer :: nph = 3000
        ! integer :: nph2 = 1000
        integer :: nph2 = 100
        integer :: nph3 = 15

        integer,intent(out)::ph_iter

        logical,intent(in)::print_cb,act_ON
        character(500),intent(in)::print_loc
        logical so4_error,print_res
        logical bisec_chk,bisec_chk_ON,bisec_only,mod_ph_order,calc_simple,halley,first_chk_done

        real(kind=8),allocatable::amx(:,:),ymx(:)
        integer,allocatable::ipiv(:)
        integer info,nmx

        real(kind=8),parameter :: threshold = 10d0
        real(kind=8),parameter :: corr = exp(threshold)

        real(kind=8) ph_tmp,ph_fact,err1,err2,slp,slplog,ph_tmp_min,ph_tmp_max,slp_save
        real(kind=8) ph_max,ph_min 
        real(kind=8) u
        integer judge

        real(kind=8) f1_min_save,ph_f1min_save
        real(kind=8),parameter :: ph_init_min = 1d-20
        real(kind=8),parameter :: ph_init_max = 1d4 

        ! bisec_chk_ON = .false.
        bisec_chk_ON = .true.

        bisec_only = .false.
        ! bisec_only = .true.

        mod_ph_order = .false.
        ! mod_ph_order = .true.

        ! calc_simple = .false.
        calc_simple = .true.

        error = 1d4
        tol = 1d-6
        dconc = 1d0
        ph_add_order = 0d0
        ph_add_order = 2d0

        k_order = 0.5d0
        ph_inflex = 7d0
        a_order = 2d0
        c_order = 2d0

        ! where(-log10(prox)<3d0)
            ! ph_add_order=0d0
        ! endwhere 

        ! ph_add_order = a_order/(1d0+EXP(-2d0*k_order*(-log10(prox)-ph_inflex))) + c_order

        ! prox = 1d0 
        iter = 0

        if (any(isnan(maqx)) .or. any(isnan(maqc))) then 
            print*,'nan in input aqueosu species'
            stop
        endif 

        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maqx,mgasx,maqc,mgasc &
            & ,maqx_loc,mgasx_loc  &! output
            & )
            
        call get_base_charge( &
            & nsp_aq_all & 
            & ,chraq_all & 
            & ,base_charge &! output 
            & )
            
        call get_maqt_all( &
        ! call get_maqt_all_v2( &
            & nz,nsp_aq_all,nsp_gas_all &
            & ,chraq_all,chrgas_all &
            & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
            & ,mgasx_loc,maqx_loc,prox,iosx,tc &
            & ,dmaqft_dpro_loc,dmaqft_dmaqf_loc,dmaqft_dmgas_loc,dmaqft_dios_loc &! output
            & ,maqft_loc  &! output
            & )
            
        iso4 = findloc(chraq_all,'so4',dim=1)
        so4x = maqx_loc(iso4,:)*maqft_loc(iso4,:)
            
        maqf_loc = maqx_loc ! fixed free concs. 

        if (.not.act_ON) iosx = 0d0

        nmx = nz*2
        nmx = nz

        if (allocated(amx)) deallocate(amx)
        if (allocated(ymx)) deallocate(ymx)
        if (allocated(ipiv)) deallocate(ipiv)
        allocate(amx(nmx,nmx),ymx(nmx),ipiv(nmx))

        ph_error = .false.

        print_res = .false.

        prox_init = prox

        ! print*,'calc_pH'
        if (.not. print_cb) then
        ! if (.true.) then
            ! obtaining ph and so4f from scratch
        
            ! prox = 1d0 
            do while (error > tol)
            ! do while (error > tol*1d-4)

                prox_save = prox
                iosx_save = iosx
                
                call calc_charge_balance( &
                    & nz,nsp_aq_all,nsp_gas_all &
                    & ,chraq_all,chrgas_all &
                    & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                    & ,base_charge &
                    & ,mgasx_loc,maqf_loc &
                    & ,z,prox,iosx,tc &
                    & ,print_loc,print_res,ph_add_order &
                    & ,f1,df1,df1dmaqf,df1dmgas &!output
                    & ,d2f1,d2f1dmaqf,d2f1dmgas &!output
                    & ,f2,df2,df2dmaqf,df2dmgas &!output
                    & ,df1df2,df2df1,ios_new &!output
                    & )
                
                ! df1 = df1*prox
                
                if (any(isnan(f1)).or.any(isnan(df1))) then 
                    print*,'found nan during the course of ph calc: newton'
                    print *,any(isnan(f1)),any(isnan(df1))
                    print *,prox
                    print *
                    print *,f1
                    print *
                    ! if (any(isnan(f1))) print *, f1
                    if (any(isnan(df1))) print *, df1
                    if (act_ON) then
                        print *,iosx
                        print *
                        print *,f2
                        print *
                        if (any(isnan(df2))) print *, df2
                    endif 
                    ph_error = .true.
                    ! stop
                    return
                    exit
                    ! pause 
                endif 
                
                if (act_ON) then
                    if (any(isnan(f2)).or.any(isnan(df2)) ) then 
                        print*,'found nan during the course of ios calc: newton'
                        print *,any(isnan(f2)),any(isnan(df2))
                        print *,iosx
                        print *
                        print *,f2
                        print *
                        ! if (any(isnan(f2))) print *, f2
                        if (any(isnan(df2))) print *, df2
                        ph_error = .true.
                        ! stop
                        return
                        exit
                        ! pause 
                    endif 
                endif 
                
                if (nmx==nz) then 
                
                    where (prox -f1/df1>0d0)
                        prox = prox -f1/df1
                    elsewhere 
                        prox = prox*dexp( -f1/df1/prox )
                    endwhere
                    error = maxval(dabs(dexp( -f1/df1/prox )-1d0))
                    
                    if (act_ON) then 
                    
                        iosx = ios_new
                        
                        error = max( error, maxval(dabs(dexp( -f2/df2/iosx )-1d0)) )
                    else
                        iosx = 0d0
                    endif 
                endif 
                
                df1 = df1*prox
                df2df1 = df2df1*prox
                df2 = df2*iosx
                df1df2 = df1df2*iosx
                
                if (any(isnan(f1)).or.any(isnan(df1))) then 
                    print*,'found nan during the course of ph calc'
                    print *,any(isnan(f1)),any(isnan(df1))
                    print *,prox
                    ph_error = .true.
                    exit
                    ! pause 
                endif 
                
                if (act_ON) then
                    if (any(isnan(f2)).or.any(isnan(df2)) &
                        & .or.any(isnan(df1df2)).or.any(isnan(df2df1))) then 
                        print*,'found nan during the course of ios calc'
                        print *,any(isnan(f2)),any(isnan(df2)) &
                            & ,any(isnan(df1df2)),any(isnan(df2df1))
                        print *,iosx
                        ph_error = .true.
                        exit
                        ! pause 
                    endif 
                endif 
                
                
                if (nmx/=nz) then 
                    amx = 0d0
                    ymx = 0d0
                    
                    ymx(1:nz) = f1(:)
                    ymx(nz+1:nmx) = f2(:)
                    
                    do iz=1,nz
                        amx(iz,iz)=df1(iz)
                        amx(nz+iz,nz+iz)=df2(iz)
                        amx(iz,nz+iz)=df1df2(iz)
                        amx(nz+iz,iz)=df2df1(iz)
                    enddo 
                    ymx = -ymx
                    
                    call DGESV(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 
                    
                    prox = prox*exp( ymx(1:nz) )
                    iosx = iosx*exp( ymx(nz+1:nmx) )
                    
                    error = maxval(abs(exp( ymx )-1d0))
                    if (isnan(error) .or. info/=0) then 
                        print *,'error in error or dgesv'
                        error = 1d4
                        ph_error = .true.
                        exit 
                    endif 
                endif 

                
                ! error = maxval(dabs((prox_save-prox)/prox))
                ! error_save = dabs((prox_save-prox)/prox)
                
                ! if (any(prox == 0d0)) then 
                    ! error = 1d4
                    ! where (prox == 0d0)
                        ! prox = 1d-12
                        ! error_save = 1d4
                    ! endwhere
                ! endif 
                
                iter = iter + 1
                
                ! print*,iter,error
                
                if (iter > 3000) then 
                    print *,'iteration exceeds 3000 with newton method: error = ',error, ' tol = ',tol,halley
                    do iz=1,nz
                        print*,iz,-log10(prox_save(iz)),-log10(prox(iz)),dabs((prox_save(iz)-prox(iz))/prox(iz)),abs(f1(iz))
                        error_save(iz) = dabs((prox_save(iz)-prox(iz))/prox(iz))
                    enddo
                    print *
                    do iz=1,nz
                        print*,iz,-log10(iosx_save(iz)),-log10(iosx(iz)),dabs((iosx_save(iz)-iosx(iz))/iosx(iz)),abs(f2(iz))
                        error_save(iz) = dabs((iosx_save(iz)-iosx(iz))/iosx(iz))
                    enddo
                    ! print*,error
                    ! print*,prox
                    ph_error = .true.
                    if (.not.bisec_chk_ON) return
                endif 
                
                if (ph_error) exit 
            enddo  
            
            bisec_chk = .false.
            if ( bisec_chk_ON .and. ph_error ) bisec_chk = .true.
            
            ! tring brutal forcing 
            if (bisec_chk) then 
                prox_save_newton = prox
                do iz=1,nz
                    if (error_save(iz)<tol) cycle
                    
                    first_chk_done = .false.
                    ! check to where a root likely exists
                    ph_min = -2d0
                    ph_max = 16d0
                    ! if (error_save(iz)<1d-5) then
                        ! ph_min = -log10(prox_save_newton(iz))-2d0
                        ! ph_max = -log10(prox_save_newton(iz))+2d0                
                    ! endif 
                    prox_tmp1 = prox
                    prox_tmp2 = prox
                    
                    print *, 'not converged @ ',iz,ph_min,ph_max,error_save(iz),-log10(prox_save_newton(iz))
                    print *, maqf_loc(findloc(chraq_all,'ca',dim=1),iz) &
                        & ,maqf_loc(findloc(chraq_all,'no3',dim=1),iz) &
                        & ,maqf_loc(findloc(chraq_all,'oxa',dim=1),iz) 

                    f1_min_save = 1d100
                    ph_f1min_save = 1d100
                    
                    do iph3 = 1,nph3
                        ph_tmp_min = 1d-100
                        ph_tmp_max = 1d100
                        ! ph_tmp_min = ph_init_min
                        ! ph_tmp_max = ph_init_max
                        ph_tmp_min = 10d0**-ph_max
                        ph_tmp_max = 10d0**-ph_min

                        ! initially give slp a random negative value to be saved to slp_save
                        slp = -100d0
                        ! print *,'start from alkaline pH'
                        do iph2=1,nph2 ! start from alkaline pH
                            ph_tmp = ph_max + (ph_min - ph_max) &
                                & * (real(iph2,kind=8)-1d0)/(real(nph2,kind=8)-1d0) 
                            ph_tmp = 10d0**-ph_tmp
                            
                            dconc = ph_tmp*1d-6
                            
                            prox_tmp1(iz) = ph_tmp + dconc
                            prox_tmp2(iz) = ph_tmp - dconc
                            
                            call calc_charge_balance_point( &
                                & nz,nsp_aq_all,nsp_gas_all &
                                & ,chraq_all,chrgas_all &
                                & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                                & ,base_charge &
                                & ,mgasx_loc,maqf_loc &
                                & ,z,prox_tmp1,iz,iosx,tc &
                                & ,print_loc,print_res,ph_add_order &
                                & ,f1_tmp1,df1_tmp1,df1dmaqf_dum,df1dmgas_dum &!output
                                & ,d2f1_tmp1,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                                & )
                            
                            call calc_charge_balance_point( &
                                & nz,nsp_aq_all,nsp_gas_all &
                                & ,chraq_all,chrgas_all &
                                & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                                & ,base_charge &
                                & ,mgasx_loc,maqf_loc &
                                & ,z,prox_tmp2,iz,iosx,tc &
                                & ,print_loc,print_res,ph_add_order &
                                & ,f1_tmp2,df1_tmp2,df1dmaqf_dum,df1dmgas_dum &!output
                                & ,d2f1_tmp2,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                                & )
                            
                            err1 = dabs(dexp(-f1_tmp1(iz)/df1_tmp1(iz)/prox_tmp1(iz))-1d0)
                            err2 = dabs(dexp(-f1_tmp2(iz)/df1_tmp2(iz)/prox_tmp2(iz))-1d0)
                            
                            if (isnan(err1)) then
                                print *,'err1 is nan',f1_tmp1(iz),df1_tmp1(iz),prox_tmp1(iz)
                                stop
                            endif 
                            if (isnan(err2)) then 
                                print *,'err2 is nan',f1_tmp2(iz),df1_tmp2(iz),prox_tmp2(iz)
                                stop
                            endif 
                            
                            slp_save = slp
                            slp = ( err1 - err2) / (2d0*dconc)
                            slplog = ( err1 - err2) / (-dlog10(prox_tmp1(iz)) - (-dlog10(prox_tmp1(iz)))  )
                            
                            prox(iz) = ph_tmp

                            call calc_charge_balance_point( &
                                & nz,nsp_aq_all,nsp_gas_all &
                                & ,chraq_all,chrgas_all &
                                & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                                & ,base_charge &
                                & ,mgasx_loc,maqf_loc &
                                & ,z,prox,iz,iosx,tc &
                                & ,print_loc,print_res,ph_add_order &
                                & ,f1_dum,df1_dum,df1dmaqf_dum,df1dmgas_dum &!output
                                & ,d2f1_dum,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                                & )
                                
                            if (isnan(f1_dum(iz))) then 
                                print *,'f1_dum(iz) is nan',f1_dum(iz),ph_tmp
                                stop
                            endif 
                            
                            ! print *,-log10(ph_tmp),slp,f1_dum(iz),-log10(ph_tmp_min)
                            
                            ! if (slp <= 0d0 .and. f1_dum(iz) <= 0d0) ph_tmp_min = max(ph_tmp_min,ph_tmp)
                            if (slp_save <= 0d0 .and. slp <= 0d0 .and. f1_dum(iz) <= 0d0) ph_tmp_min = max(ph_tmp_min,ph_tmp)
                            
                            if (abs(f1_dum(iz)) < f1_min_save) then
                                f1_min_save = abs(f1_dum(iz))
                                ph_f1min_save = ph_tmp
                            endif 
                        
                        enddo 
                        
                        ! initially give slp a random positive value to be saved to slp_save
                        slp = 100d0
                        ! print *,'start from acidic pH'
                        do iph2=nph2,1,-1 ! start from acidic pH
                            ph_tmp = ph_max + (ph_min - ph_max) &
                                & * (real(iph2,kind=8)-1d0)/(real(nph2,kind=8)-1d0) 
                            ph_tmp = 10d0**-ph_tmp
                            
                            dconc = ph_tmp*1d-6
                            
                            prox_tmp1(iz) = ph_tmp + dconc
                            prox_tmp2(iz) = ph_tmp - dconc
                            
                            call calc_charge_balance_point( &
                                & nz,nsp_aq_all,nsp_gas_all &
                                & ,chraq_all,chrgas_all &
                                & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                                & ,base_charge &
                                & ,mgasx_loc,maqf_loc &
                                & ,z,prox_tmp1,iz,iosx,tc &
                                & ,print_loc,print_res,ph_add_order &
                                & ,f1_tmp1,df1_tmp1,df1dmaqf_dum,df1dmgas_dum &!output
                                & ,d2f1_tmp1,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                                & )
                            
                            call calc_charge_balance_point( &
                                & nz,nsp_aq_all,nsp_gas_all &
                                & ,chraq_all,chrgas_all &
                                & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                                & ,base_charge &
                                & ,mgasx_loc,maqf_loc &
                                & ,z,prox_tmp2,iz,iosx,tc &
                                & ,print_loc,print_res,ph_add_order &
                                & ,f1_tmp2,df1_tmp2,df1dmaqf_dum,df1dmgas_dum &!output
                                & ,d2f1_tmp2,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                                & )
                            
                            err1 = dabs(dexp(-f1_tmp1(iz)/df1_tmp1(iz)/prox_tmp1(iz))-1d0)
                            err2 = dabs(dexp(-f1_tmp2(iz)/df1_tmp2(iz)/prox_tmp2(iz))-1d0)
                            
                            if (isnan(err1)) then
                                print *,'err1 is nan',f1_tmp1(iz),df1_tmp1(iz),prox_tmp1(iz)
                                stop
                            endif 
                            if (isnan(err2)) then 
                                print *,'err2 is nan',f1_tmp2(iz),df1_tmp2(iz),prox_tmp2(iz)
                                stop
                            endif 
                            
                            slp_save = slp
                            slp = ( err1 - err2) / (2d0*dconc)
                            slplog = ( err1 - err2) / (-dlog10(prox_tmp1(iz)) - (-dlog10(prox_tmp1(iz)))  )
                            
                            prox(iz) = ph_tmp

                            print_res = .true.
                            print_res = .false.
                            call calc_charge_balance_point( &
                                & nz,nsp_aq_all,nsp_gas_all &
                                & ,chraq_all,chrgas_all &
                                & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                                & ,base_charge &
                                & ,mgasx_loc,maqf_loc &
                                & ,z,prox,iz,iosx,tc &
                                & ,print_loc,print_res,ph_add_order &
                                & ,f1_dum,df1_dum,df1dmaqf_dum,df1dmgas_dum &!output
                                & ,d2f1_dum,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                                & )
                            print_res = .false.
                            
                            if (isnan(f1_dum(iz))) then 
                                print *,'f1_dum(iz) is nan',f1_dum(iz),ph_tmp
                                stop
                            endif 
                                
                            ! print *,-log10(ph_tmp),slp,f1_dum(iz),-log10(ph_tmp_max)
                            
                            ! if (slp >= 0d0 .and. f1_dum(iz) >= 0d0) ph_tmp_max = min(ph_tmp_max,ph_tmp)
                            if (slp_save >= 0d0 .and. slp >= 0d0 .and. f1_dum(iz) >= 0d0) ph_tmp_max = min(ph_tmp_max,ph_tmp)
                            
                            ! if (abs(f1_dum(iz)) < f1_min_save) then
                                ! f1_min_save = abs(f1_dum(iz))
                                ! ph_f1min_save = ph_tmp
                            ! endif 
                        
                        enddo 
                        ph_max = -log10(ph_tmp_min)
                        ph_min = -log10(ph_tmp_max)
                        error = abs( 10d0**-ph_min -  10d0**-ph_max)/10d0**-ph_max
                        
                        ! if (iph3 /= nph3) then 
                            ! if ( abs((ph_tmp_min - ph_init_min)/ph_init_min) < 1d-6 &
                                ! & .and. abs((ph_tmp_max - ph_init_max)/ph_init_max) > 1d-6  &
                                ! & ) then
                                ! ph_tmp_min = ph_tmp_max * 1d-4
                                ! error =1d4
                            ! endif 
                            
                            ! if ( abs((ph_tmp_min - ph_init_min)/ph_init_min) > 1d-6 &
                                ! & .and. abs((ph_tmp_max - ph_init_max)/ph_init_max) < 1d-6  &
                                ! & ) then
                                ! ph_tmp_max = ph_tmp_min * 1d4
                                ! error =1d4
                            ! endif 
                        ! endif 
                            
                        print*, iph3, 'a root likely between'& 
                            & , ph_min , 'and', ph_max, '>>> error=', error 
                        
                        if (ph_min > ph_max) then 
                            print *, 'ph_min > ph_max detected: something is wrong in bracketing root'
                            print *, 'Possibility: there could be 2 solutions to charge balance equation'
                            print *, '--> discard ph_min or ph_max randomly and get a new ph_min or ph_max as a midpoint'
                            ! stop
                            if (error >= tol) then 
                                ! stop
                                
                                ph_tmp = ph_min
                                ph_min = ph_max
                                ph_max = ph_tmp
                                
                                call random_number(u)
                                judge = 0 + FLOOR(2*u)
                                
                                if (judge ==0) then
                                    ph_min = 0.5d0*(ph_min + ph_max)
                                else
                                    ph_max = 0.5d0*(ph_min + ph_max)
                                endif 
                                
                                ! ph_error = .true.
                                ! return
                            else
                                print *, ' error is small so do not care the above message'
                                prox(iz) = 10d0**(-0.5d0*(ph_max + ph_min))
                                first_chk_done = .true.
                                exit
                            endif 
                        endif 
                        
                        ! prox(iz) = prox_save_newton(iz)
                        
                        if (error < tol*1d-6) then 
                            print *, ' *** root found *** ',iz
                            prox(iz) = 10d0**(-0.5d0*(ph_max + ph_min))
                            first_chk_done = .true.
                            exit
                        endif 
                        
                        ! if (iph3 == nph3) then
                            ! print *,' *** too large error *** '
                            ! ph_tmp = ph_f1min_save
                            ! print *,' ... so adopt where error can be minimum? pH = ',-log10(ph_tmp)
                            ! prox(iz) = ph_tmp
                        ! endif 
                        
                    enddo
                    ! pause
                    cycle
                    if (first_chk_done) cycle
                    
                    ! trying to find solution where error = 0d0 instead of f1
                    prox_tmp1 = prox
                    prox_tmp2 = prox
                    
                    dconc = 1d-7
                    
                    first_chk_done = .false.
                    
                    do iph=1,nph
                        
                        ! dconc = prox(iz)*1d-6
                        
                        prox_tmp1(iz) = prox(iz) + dconc
                        prox_tmp2(iz) = prox(iz) - dconc
                        
                        call calc_charge_balance_point( &
                            & nz,nsp_aq_all,nsp_gas_all &
                            & ,chraq_all,chrgas_all &
                            & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                            & ,base_charge &
                            & ,mgasx_loc,maqf_loc &
                            & ,z,prox_tmp1,iz,iosx,tc &
                            & ,print_loc,print_res,ph_add_order &
                            & ,f1_tmp1,df1_tmp1,df1dmaqf_dum,df1dmgas_dum &!output
                            & ,d2f1_tmp1,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                            & )
                        
                        call calc_charge_balance_point( &
                            & nz,nsp_aq_all,nsp_gas_all &
                            & ,chraq_all,chrgas_all &
                            & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                            & ,base_charge &
                            & ,mgasx_loc,maqf_loc &
                            & ,z,prox_tmp2,iz,iosx,tc &
                            & ,print_loc,print_res,ph_add_order &
                            & ,f1_tmp2,df1_tmp2,df1dmaqf_dum,df1dmgas_dum &!output
                            & ,d2f1_tmp2,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                            & )
                        
                        err1 = dabs(dexp(-f1_tmp1(iz)/df1_tmp1(iz)/prox_tmp1(iz))-1d0)
                        err2 = dabs(dexp(-f1_tmp2(iz)/df1_tmp2(iz)/prox_tmp2(iz))-1d0)
                        
                        if (isnan(err1) .or. isnan(err2)) then 
                            print*,prox_tmp1(iz),f1_tmp1(iz),df1_tmp1(iz)
                            print*,prox_tmp2(iz),f1_tmp2(iz),df1_tmp2(iz)
                            ! stop
                            ! exit
                            err1 = f1_tmp1(iz)
                            err2 = f1_tmp2(iz)
                        endif 
                        
                        slp = ( err1 - err2) / (2d0*dconc)
                        slplog = ( err1 - err2) / (-dlog10(prox_tmp1(iz)) - (-dlog10(prox_tmp1(iz)))  )
                        
                        if (prox_tmp1(iz) - err1/slp > 0d0) then 
                            ph_tmp = prox_tmp1(iz) - err1/slp
                        else 
                            ph_tmp = -dlog10(prox_tmp1(iz)) - err1/slplog
                            ph_tmp = 10d0**(-ph_tmp)
                        endif 
                        
                        prox(iz) = ph_tmp

                        call calc_charge_balance_point( &
                            & nz,nsp_aq_all,nsp_gas_all &
                            & ,chraq_all,chrgas_all &
                            & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                            & ,base_charge &
                            & ,mgasx_loc,maqf_loc &
                            & ,z,prox,iz,iosx,tc &
                            & ,print_loc,print_res,ph_add_order &
                            & ,f1_dum,df1_dum,df1dmaqf_dum,df1dmgas_dum &!output
                            & ,d2f1_dum,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                            & )

                        ! print *, iph, -log10(ph_tmp),dabs(dexp(-f1_dum(iz)/df1_dum(iz)/prox(iz))-1d0),dabs(f1_dum(iz))

                        if ( dabs(dexp(-f1_dum(iz)/df1_dum(iz)/prox(iz))-1d0) < tol ) then 
                            first_chk_done = .true.
                            exit 
                        endif 
                    
                    enddo 
                    print *, 'new ph ', -dlog10(ph_tmp), 'old ph ',  -dlog10(prox_save_newton(iz))
                    
                    if (first_chk_done) cycle
                    
                    ph_tmp = prox_save_newton(iz)
                    ph_fact = 3d0
                    f1_tmp = 1d100
                    do iph=1,nph
                        prox(iz) = ph_fact*prox_save_newton(iz) &
                            & + (prox_save_newton(iz)/ph_fact - ph_fact*prox_save_newton(iz)) &
                            & * (real(iph,kind=8)-1d0)/(real(nph,kind=8)-1d0) 
                        call calc_charge_balance_point( &
                            & nz,nsp_aq_all,nsp_gas_all &
                            & ,chraq_all,chrgas_all &
                            & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                            & ,base_charge &
                            & ,mgasx_loc,maqf_loc &
                            & ,z,prox,iz,iosx,tc &
                            & ,print_loc,print_res,ph_add_order &
                            & ,f1_dum,df1_dum,df1dmaqf_dum,df1dmgas_dum &!output
                            & ,d2f1_dum,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                            & )
                        ! if (ph_tmp -f1_dum(iz)/df1_dum(iz) > 0d0) then 
                            ! ph_tmp = ph_tmp -f1_dum(iz)/df1_dum(iz) 
                        ! else 
                            ! ph_tmp = ph_tmp*exp(-f1_dum(iz)/df1_dum(iz)/ph_tmp)
                        ! endif 
                        ! ph_tmp = ph_tmp*exp(-f1_dum(iz)/df1_dum(iz)/prox(iz))
                        ! ph_tmp = ph_tmp &
                            ! & *dexp( -2d0*f1_dum(iz)*df1_dum(iz)/(2d0*df1_dum(iz)**2d0 - f1_dum(Iz)*d2f1_dum(iz) )/prox(iz))
                        ! prox(iz) = ph_tmp
                        ! print *, iph, -log10(ph_tmp),dabs(exp(-f1_dum(iz)/df1_dum(iz)/prox(iz))-1d0),dabs(f1_dum(iz))
                        ! if ( dabs(exp(-f1_dum(iz)/df1_dum(iz)/prox(iz))-1d0) < tol ) exit 
                        if ( abs(f1_dum(iz)) < abs(f1_tmp(iz)) ) then 
                            ph_tmp = prox(iz)
                            f1_tmp(iz) = f1_dum(iz)
                        endif 
                    enddo 
                    prox(iz) = ph_tmp
                    print *, 'new ph ', -dlog10(ph_tmp), 'old ph ',  -dlog10(prox_save_newton(iz))
                enddo 
            endif 
            
            
        endif 

        ph_iter = iter

        ph_error = .false.

        if (any(isnan(prox)) .or. any(prox<=0d0)) then     
            print *, (-log10(prox(iz)),iz=1,nz,nz/5)
            print*,'ph is nan or <= zero'
            ! prox = prox_init
            ph_error = .true.
            ! stop
        endif 

        if (print_cb) print_res = .true.

                
        call calc_charge_balance( &
            & nz,nsp_aq_all,nsp_gas_all &
            & ,chraq_all,chrgas_all &
            & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
            & ,base_charge &
            & ,mgasx_loc,maqf_loc &
            & ,z,prox,iosx,tc &
            & ,print_loc,print_res,ph_add_order &
            & ,f1,df1,df1dmaqf,df1dmgas &!output
            & ,d2f1,d2f1dmaqf,d2f1dmgas &!output
            & ,f2,df2,df2dmaqf,df2dmgas &!output
            & ,df1df2,df2df1,ios_new &!output
            & )

        ! if ( maxval(abs((ios_new - iosx)/iosx)) > tol ) then
            ! print *, 'Ionic strength calculation check failure'
            ! print *, maxval(abs((ios_new - iosx)/iosx))
            ! stop
        ! endif 

        ! stop

        ! ### CHECKING WIHT POINT SUBROUTINE WORKS ###

        ! f1_tmp = 0d0
        ! df1_tmp = 0d0 
        ! d2f1_tmp = 0d0
        ! df1dmaqf_tmp = 0d0
        ! df1dmgas_tmp = 0d0
        ! d2f1dmaqf_tmp = 0d0
        ! d2f1dmgas_tmp = 0d0 
        ! print_res = .false.
        ! do iz=1,nz
            ! call calc_charge_balance_point( &
                ! & nz,nsp_aq_all,nsp_gas_all &
                ! & ,chraq_all,chrgas_all &
                ! & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
                ! & ,base_charge &
                ! & ,mgasx_loc,maqf_loc &
                ! & ,z,prox,iz &
                ! & ,print_loc,print_res,ph_add_order &
                ! & ,f1_dum,df1_dum,df1dmaqf_dum,df1dmgas_dum &!output
                ! & ,d2f1_dum,d2f1dmaqf_dum,d2f1dmgas_dum &!output
                ! & )
            ! f1_tmp(iz) = f1_dum(iz)
            ! df1_tmp(iz) = df1_dum(iz) 
            ! d2f1_tmp(iz) = d2f1_dum(iz)
            ! df1dmaqf_tmp(:,iz) = df1dmaqf_dum(:,iz)
            ! df1dmgas_tmp(:,iz) = df1dmgas_dum(:,iz)
            ! d2f1dmaqf_tmp(:,iz) = d2f1dmaqf_dum(:,iz)
            ! d2f1dmgas_tmp(:,iz) = d2f1dmgas_dum(:,iz)

        ! enddo 

        ! if (maxval(abs((f1_tmp -f1)/f1)) > tol) then 
            ! print *, 'point wrong for f1?'
            ! stop
        ! endif 

        ! if (maxval(abs((df1_tmp - df1)/df1)) > tol) then 
            ! print *, 'point wrong for df1?'
            ! stop
        ! endif 

        ! if (maxval(abs((df1dmaqf_tmp - df1dmaqf)/df1dmaqf)) > tol) then 
            ! print *, 'point wrong for df1dmaqf?'
            ! stop
        ! endif 

        ! if (maxval(abs((df1dmgas_tmp - df1dmgas)/df1dmgas)) > tol) then 
            ! print *, 'point wrong for df1dmaqf?'
            ! stop
        ! endif 


        ! ### END CHECKING ###

        do ispa = 1, nsp_aq_all
            dprodmaq_all(ispa,:) = - df1dmaqf(ispa,:) / df1   
        enddo 

        do ispg = 1, nsp_gas_all
            dprodmgas_all(ispg,:) = - df1dmgas(ispg,:) /df1
        enddo 

        diosdmaq_all = 0d0
        diosdmgas_all = 0d0

        if (act_ON) then 
            do ispa = 1, nsp_aq_all
                diosdmaq_all(ispa,:) = - df2dmaqf(ispa,:) / df2   
            enddo 

            do ispg = 1, nsp_gas_all
                diosdmgas_all(ispg,:) = - df2dmgas(ispg,:) /df2
            enddo 
        endif 

        ! solving two equations analytically:
        ! df1/dmsp + df1/dph * dph/dmsp  = 0  
        ! df1/dmsp * Dmsp + df1/dph *Dph = 0
        ! df1/dmsp * Dmsp + df1/dph *Dph + (d2f1/dmsp2) * (Dmsp)^2 + (d2f1/d2ph) *(Dph)^2 = 0
        ! df1/dmsp * Dmsp + (d2f1/dmsp2) * (Dmsp)^2 + df1/dph *Dph +  (d2f1/d2ph) *(Dph)^2 = 0

        ! do ispa = 1, nsp_aq_all
            ! dprodmaq_all(ispa,:) = - (df2*df1dmaqf(ispa,:) - df1df2*df2dmaqf(ispa,:))/(df2*df1 - df1df2*df2df1)   
            ! diosdmaq_all(ispa,:) = - ( df2df1*df1dmaqf(ispa,:) - df1*df2dmaqf(ispa,:) )/(df2df1*df1df2 - df1*df2 ) 
        ! enddo 

        ! do ispg = 1, nsp_gas_all
            ! dprodmgas_all(ispg,:) = - (df2*df1dmgas(ispg,:) - df1df2*df2dmgas(ispg,:) )/(df2*df1 -df12*df2df1)
            ! diosdmgas_all(ispg,:) = - ( df2df1*df1dmgas(ispg,:) - df1*df2dmgas(ispg,:) )/(df2df1*df1df2 - df1*df2 )  
        ! enddo 

        return

    endsubroutine calc_pH_v7_4

    ! Charge balance calculation
    subroutine calc_charge_balance( &
        & nz,nsp_aq_all,nsp_gas_all &
        & ,chraq_all,chrgas_all &
        & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
        & ,base_charge &
        & ,mgasx_loc,maqf_loc &
        & ,z,prox,iosx,tc &
        & ,print_loc,print_res,ph_add_order &
        & ,f1,df1,df1dmaqf,df1dmgas &!output
        & ,d2f1,d2f1dmaqf,d2f1dmgas &!output
        & ,f2,df2,df2dmaqf,df2dmgas &!output
        & ,df1df2,df2df1,ios_new &!output
        & )
        implicit none

        integer,intent(in)::nz,nsp_aq_all,nsp_gas_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        real(kind=8),intent(in)::kw,tc
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl
        real(kind=8),dimension(nsp_gas_all,nz),intent(in)::mgasx_loc
        real(kind=8),dimension(nsp_aq_all,nz),intent(in)::maqf_loc
        real(kind=8),dimension(nsp_aq_all),intent(in)::base_charge
        real(kind=8),dimension(nz),intent(in)::z,prox,ph_add_order
        real(kind=8),dimension(nz),intent(in)::iosx
        real(kind=8),dimension(nz),intent(out)::f1,df1,d2f1
        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::df1dmaqf,d2f1dmaqf
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::df1dmgas,d2f1dmgas
        real(kind=8),dimension(nz),intent(out)::f2,df2,df1df2,df2df1,ios_new
        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::df2dmaqf
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::df2dmgas
        real(kind=8),dimension(nz)::d2f2
        real(kind=8),dimension(nsp_aq_all,nz)::d2f2dmaqf
        real(kind=8),dimension(nsp_gas_all,nz)::d2f2dmgas

        logical,intent(in)::print_res
        character(500),intent(in)::print_loc
        character(500)::path_tmp,index_tmp

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ispa,ispa_h,ispa_c,ispa_s,iz,ipco2,ipnh3,iso4,ioxa,ispa_no3,ino3,ispa_nh3,ispa_oxa,ispa_cl &
            & ,icl,icharge,ic1,ic2,ic3

        real(kind=8) kco2,k1,k2,knh3,k1nh3,rspa_h,rspa_s,rspa_no3,rspa_nh3,rspa_oxa,rspa_oxa_2,rspa_oxa_3 &
            & ,rspa_cl,rcharge
        ! real(kind=8) tc
        real(kind=8),dimension(nz)::pco2x,pnh3x,so4f,no3f,oxaf,clf
        real(kind=8),dimension(nz)::isf,fkw,fkeq,dfkw_dios,dfkeq_dios
        real(kind=8),dimension(nz)::gamma_tmp,dgamma_dios_tmp
        real(kind=8),dimension(4,nz)::gamma,dgamma_dios
        real(kind=8),dimension(nz)::f1_chk,ss_add,back

        character(1) chrint

        real(kind=8),dimension(nsp_gas_all,3,nz)::fkeqgas_h
        real(kind=8),dimension(nsp_aq_all,4,nz)::fkeqaq_h
        real(kind=8),dimension(nsp_aq_all,2,nz)::fkeqaq_c,fkeqaq_s,fkeqaq_no3,fkeqaq_nh3,fkeqaq_oxa,fkeqaq_cl
        #ifdef debug_phcalc
        logical::debug = .true. 
        #else
        logical::debug = .false. 
        #endif

        path_tmp = print_loc(:index(print_loc,'.txt')-5)
        index_tmp = print_loc(index(print_loc,'.txt')-4:)

        if (print_res) open(88,file = trim(adjustl(print_loc)),status='replace')
        if (print_res) then 
            if (print_loc == './ph.txt') then 
                open(99,file = './ph(eq).txt',status='replace')
            else
                open(99,file = trim(adjustl(path_tmp))//'(eq)'//trim(adjustl(index_tmp)),status='replace')
            endif 
        endif 

        ipco2   = findloc(chrgas_all,'pco2',dim=1)
        ipnh3   = findloc(chrgas_all,'pnh3',dim=1)
        iso4    = findloc(chraq_all,'so4',dim=1)
        ino3    = findloc(chraq_all,'no3',dim=1)
        ioxa    = findloc(chraq_all,'oxa',dim=1)
        icl     = findloc(chraq_all,'cl',dim=1)

        kco2    = keqgas_h(ipco2,ieqgas_h0)
        k1      = keqgas_h(ipco2,ieqgas_h1)
        k2      = keqgas_h(ipco2,ieqgas_h2)

        pco2x   = mgasx_loc(ipco2,:)


        knh3    = keqgas_h(ipnh3,ieqgas_h0)
        k1nh3   = keqgas_h(ipnh3,ieqgas_h1)

        pnh3x   = mgasx_loc(ipnh3,:)

        so4f    = maqf_loc(iso4,:)
        no3f    = maqf_loc(ino3,:)
        oxaf    = maqf_loc(ioxa,:)
        clf     = maqf_loc(icl,:)

        ss_add = ph_add_order

        f1 = 0d0
        df1 = 0d0
        d2f1 = 0d0
        df1dmaqf = 0d0
        df1dmgas = 0d0
        d2f1dmaqf = 0d0
        d2f1dmgas = 0d0

        f2 = 0d0
        df2 = 0d0
        d2f2 = 0d0
        df2dmaqf = 0d0
        df2dmgas = 0d0
        d2f2dmaqf = 0d0
        d2f2dmgas = 0d0

        df1df2=0d0
        df2df1=0d0

        back = 1d0
        back = 0d0

        fkeqaq_c=0d0;fkeqaq_s=0d0;fkeqaq_no3=0d0;fkeqaq_nh3=0d0;fkeqaq_oxa=0d0;fkeqaq_cl=0d0
        fkeqaq_h=0d0
        fkeqgas_h=0d0

        do icharge=1,4
            rcharge = 1d0*icharge
            call calc_gamma_davies(  &
                & nz,iosx,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma(icharge,:)=gamma_tmp(:)
            dgamma_dios(icharge,:)=dgamma_dios_tmp(:)
        enddo
            
        fkw = 1d0/gamma(1,:)/gamma(1,:) ! H2O = H+ + OH- <--> Kw = {H+}{OH-} <--> Kw/gamma/gamma = [H+][OH-]
        dfkw_dios = 1d0*(-2d0)*gamma(1,:)**(-3d0)*dgamma_dios(1,:)

        ! print *,fkw 
        ! print *,dfkw_dios 
        ! stop

        f1 = f1 + prox**(ss_add+1d0) - fkw*kw*prox**(ss_add-1d0) 
        df1 = df1 + (ss_add+1d0)*prox**ss_add - fkw*kw*(ss_add-1d0)*prox**(ss_add-2d0) 
        d2f1 = d2f1 + (ss_add+1d0)*ss_add*prox**(ss_add-1d0) &
            & - fkw*kw*(ss_add-1d0)*(ss_add-2d0)*prox**(ss_add-3d0) 
        df1df2 = df1df2 - dfkw_dios*kw*prox**(ss_add-1d0)
        f2 = f2 - 2d0*iosx*prox**(ss_add) + prox**(ss_add+1d0) + fkw*kw*prox**(ss_add-1d0) 
        df2 = df2 - 2d0*prox**(ss_add) + dfkw_dios*kw*prox**(ss_add-1d0)
        df2df1 = df2df1 - 2d0*iosx*ss_add*prox**(ss_add-1d0) &
            & + (ss_add+1d0)*prox**(ss_add) + fkw*kw*(ss_add-1d0)*prox**(ss_add-2d0) 
        if (print_res) write(88,'(3A11)', advance='no') 'z','h', 'oh'
        if (print_res) write(99,'(3A11)', advance='no') 'z','h', 'oh'
        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) print*,'nan found f1 and df1: point 1'

        ! adding charges coming from aq species in eq with gases
        ! pCO2 
        ! Kco2: CO2(g) = CO2(a) assume no correction for activity/fugacity 
        ! K1  : CO2(a) + H2O = HCO3- + H+ <--> K1 = {HCO3-}{H+}/{CO2(a)} <--> K1/gamma/gamma = [HCO3-][H+]/[CO2(a)]
        ! K2  : HCO3- = CO32- + H+ <--> K2 = {CO32-}{H+}/{HCO3-} <--> K2*gamma/gamma/gamma2 = [CO32-][H+]/[HCO3-]
        fkeq = 1d0/gamma(2,:)
        dfkeq_dios = -1d0/gamma(2,:)**2d0*dgamma_dios(2,:)

        f1 = f1  -  fkw*k1*kco2*pco2x*prox**(ss_add-1d0)  -  2d0*fkeq*fkw*k2*k1*kco2*pco2x*prox**(ss_add-2d0)
        df1 = df1  -  fkw*k1*kco2*pco2x*(ss_add-1d0)*prox**(ss_add-2d0)  -  2d0*fkeq*fkw*k2*k1*kco2*pco2x*(ss_add-2d0)*prox**(ss_add-3d0)
        d2f1 = d2f1  -  fkw*k1*kco2*pco2x*(ss_add-1d0)*(ss_add-2d0)*prox**(ss_add-3d0)  &
            & -  2d0*fkeq*fkw*k2*k1*kco2*pco2x*(ss_add-2d0)*(ss_add-3d0)*prox**(ss_add-4d0)
        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) -  fkw*k1*kco2*1d0*prox**(ss_add-1d0)  -  2d0*fkeq*fkw*k2*k1*kco2*1d0*prox**(ss_add-2d0)
        df1df2 = df1df2  + ( &
            & -  dfkw_dios*k1*kco2*pco2x*prox**(ss_add-1d0)  &
            & -  2d0*dfkeq_dios*fkw*k2*k1*kco2*pco2x*prox**(ss_add-2d0) &
            & -  2d0*fkeq*dfkw_dios*k2*k1*kco2*pco2x*prox**(ss_add-2d0) &
            & )
        f2 = f2  +  fkw*k1*kco2*pco2x*prox**(ss_add-1d0)  +  4d0*fkeq*fkw*k2*k1*kco2*pco2x*prox**(ss_add-2d0)
        df2 = df2  + ( &
            & +  dfkw_dios*k1*kco2*pco2x*prox**(ss_add-1d0)  &
            & +  4d0*dfkeq_dios*fkw*k2*k1*kco2*pco2x*prox**(ss_add-2d0) &
            & +  4d0*fkeq*dfkw_dios*k2*k1*kco2*pco2x*prox**(ss_add-2d0) &
            & )
        df2df1 = df2df1  &
            & +  fkw*k1*kco2*pco2x*(ss_add-1d0)*prox**(ss_add-2d0)  +  4d0*fkeq*fkw*k2*k1*kco2*pco2x*(ss_add-2d0)*prox**(ss_add-3d0)
        df2dmgas(ipco2,:) = df2dmgas(ipco2,:) +  fkw*k1*kco2*1d0*prox**(ss_add-1d0)  +  4d0*fkeq*fkw*k2*k1*kco2*1d0*prox**(ss_add-2d0)
        if (print_res) write(88,'(2A11)', advance='no') 'hco3','co3'
        if (print_res) write(99,'(2A11)', advance='no') 'hco3','co3'
        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) print*,'nan found f1 and df1: point 2'
        ! pNH3 
        ! k1nh3: NH4+ = NH3 + H+ (no change in thermodynamic const is necessary?)
        f1 = f1  +  pnh3x*knh3/k1nh3*prox**(ss_add+1d0)
        df1 = df1  +  pnh3x*knh3/k1nh3*(ss_add+1d0)*prox**ss_add
        d2f1 = d2f1  +  pnh3x*knh3/k1nh3*(ss_add+1d0)*ss_add*prox**(ss_add-1d0)
        df1dmgas(ipnh3,:) = df1dmgas(ipnh3,:)  +  1d0*knh3/k1nh3*prox**(ss_add+1d0)
        f2 = f2  +  pnh3x*knh3/k1nh3*prox**(ss_add+1d0)
        df2df1 = df2df1  +  pnh3x*knh3/k1nh3*(ss_add+1d0)*prox**ss_add
        df2dmgas(ipnh3,:) = df2dmgas(ipnh3,:)  +  1d0*knh3/k1nh3*prox**(ss_add+1d0)
        if (print_res) write(88,'(A11)', advance='no') 'nh4'
        if (print_res) write(99,'(A11)', advance='no') 'nh4'
        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) print*,'nan found f1 and df1: point 3'

        do ispa = 1, nsp_aq_all
            
            f1 = f1 + base_charge(ispa)*maqf_loc(ispa,:)*prox**(ss_add)
            df1 = df1 + ( &
                & + base_charge(ispa)*maqf_loc(ispa,:)*(ss_add)*prox**(ss_add-1d0)  &
                & )
            d2f1 = d2f1 + ( &
                & + base_charge(ispa)*maqf_loc(ispa,:)*(ss_add)*(ss_add-1d0)*prox**(ss_add-2d0)  &
                & )
            df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + base_charge(ispa)*1d0*prox**(ss_add)
            f2 = f2 + base_charge(ispa)**2d0*maqf_loc(ispa,:)*prox**(ss_add)
            df2df1 = df2df1 + ( &
                & + base_charge(ispa)**2d0*maqf_loc(ispa,:)*(ss_add)*prox**(ss_add-1d0)  &
                & )
            df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + base_charge(ispa)**2d0*1d0*prox**(ss_add)
            if (print_res) write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))
            if (print_res) write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))
            if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) print*,'nan found f1 and df1: point 4 | '//trim(adjustl(chraq_all(ispa)))
            
            ! account for speces associated with NH4+ (both anions and cations: X + NH4+ = XNH4+)
            do ispa_nh3 = 1,2
                if ( keqaq_nh3(ispa,ispa_nh3) > 0d0) then 
                    rspa_nh3 = real(ispa_nh3,kind=8)
                    ic1 = nint(abs(base_charge(ispa)))
                    ic2 = nint(abs(base_charge(ispa)+rspa_nh3))
                    if ( ic1>0 .and. ic2 > 0) then  
                        fkeq = gamma(ic1,:)*gamma(1,:)**rspa_nh3/gamma(ic2,:)
                        dfkeq_dios = ( &
                            & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_nh3/gamma(ic2,:) &
                            & + gamma(ic1,:)*rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                            &   /gamma(ic2,:) &
                            & + gamma(ic1,:)*gamma(1,:)**rspa_nh3*(-1d0) &
                            &   /gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                            & )
                    elseif ( ic1==0 .and. ic2 > 0) then  
                        fkeq = gamma(1,:)**rspa_nh3/gamma(ic2,:)
                        dfkeq_dios = ( &
                            & + rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                            &   /gamma(ic2,:) &
                            & + gamma(1,:)**rspa_nh3*(-1d0) &
                            &   /gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                            & )
                    elseif ( ic1>0 .and. ic2 == 0) then  
                        fkeq = gamma(ic1,:)*gamma(1,:)**rspa_nh3
                        dfkeq_dios = ( &
                            & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_nh3 &
                            & + gamma(ic1,:)*rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                            & )
                    elseif ( ic1==0 .and. ic2 == 0) then  
                        fkeq = gamma(1,:)**rspa_nh3
                        dfkeq_dios = ( &
                            & + rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                            & )
                    else    
                        print *, 'something is wrong'
                        stop
                    endif 
                    fkeqaq_nh3(ispa,ispa_nh3,:) = fkeq
                    f1 = f1 + (base_charge(ispa) + rspa_nh3)*fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:) &
                        & *(pnh3x*knh3/k1nh3)**rspa_nh3*prox**(rspa_nh3+ss_add)
                    df1 = df1 + ( & 
                        & + (base_charge(ispa) + rspa_nh3)*fkeq &
                        & *keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(pnh3x*knh3/k1nh3)**rspa_nh3 &
                        & *(rspa_nh3+ss_add)*prox**(rspa_nh3+ss_add-1d0) &
                        & )
                    d2f1 = d2f1 + ( & 
                        & + (base_charge(ispa) + rspa_nh3)*fkeq &
                        & *keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(rspa_nh3+ss_add)*(rspa_nh3+ss_add-1d0)*(pnh3x*knh3/k1nh3)**rspa_nh3 &
                        & *prox**(rspa_nh3+ss_add-2d0) &
                        & )
                    df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + (& 
                        & + (base_charge(ispa) + rspa_nh3)*fkeq*keqaq_nh3(ispa,ispa_nh3) &
                        & *1d0*(pnh3x*knh3/k1nh3)**rspa_nh3*prox**(rspa_nh3+ss_add) &
                        & )
                    df1dmgas(ipnh3,:) = df1dmgas(ipnh3,:) + (& 
                        & + (base_charge(ispa) + rspa_nh3)*fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:) &
                        & *(knh3/k1nh3)**rspa_nh3*rspa_nh3*rspa_nh3**(rspa_nh3-1d0)*prox**(rspa_nh3+ss_add) &
                        & )
                    df1df2 = df1df2 + (base_charge(ispa) + rspa_nh3)*dfkeq_dios*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:) &
                        & *(pnh3x*knh3/k1nh3)**rspa_nh3*prox**(rspa_nh3+ss_add)
                    f2 = f2 + (base_charge(ispa) + rspa_nh3)**2d0*fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:) &
                        & *(pnh3x*knh3/k1nh3)**rspa_nh3*prox**(rspa_nh3+ss_add)
                    df2df1 = df2df1 + ( & 
                        & + (base_charge(ispa) + rspa_nh3)**2d0*fkeq &
                        & *keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(pnh3x*knh3/k1nh3)**rspa_nh3 &
                        & *(rspa_nh3+ss_add)*prox**(rspa_nh3+ss_add-1d0) &
                        & )
                    df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + (& 
                        & + (base_charge(ispa) + rspa_nh3)**2d0*fkeq*keqaq_nh3(ispa,ispa_nh3) &
                        & *1d0*(pnh3x*knh3/k1nh3)**rspa_nh3*prox**(rspa_nh3+ss_add) &
                        & )
                    df2dmgas(ipnh3,:) = df2dmgas(ipnh3,:) + (& 
                        & + (base_charge(ispa) + rspa_nh3)**2d0*fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:) &
                        & *(knh3/k1nh3)**rspa_nh3*rspa_nh3*rspa_nh3**(rspa_nh3-1d0)*prox**(rspa_nh3+ss_add) &
                        & )
                    df2 = df2 + (base_charge(ispa) + rspa_nh3)**2d0*dfkeq_dios*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:) &
                        & *(pnh3x*knh3/k1nh3)**rspa_nh3*prox**(rspa_nh3+ss_add)
                    if (print_res) then 
                        write(chrint,'(I1)') ispa_nh3
                        write(88,'(A11)', advance='no') '(nh4)'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                        write(99,'(A11)', advance='no') '(nh4)'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                    endif 
                    if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                        write(chrint,'(I1)') ispa_nh3
                        print'("nan found f1 and df1: point 5 | ",A11)', '(nh4)'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                    endif 
                    
                endif 
            enddo 
            
            ! anions
            if ( &
                & trim(adjustl(chraq_all(ispa)))=='no3' &
                & .or. trim(adjustl(chraq_all(ispa)))=='so4' &
                & .or. trim(adjustl(chraq_all(ispa)))=='cl' &
                & .or. trim(adjustl(chraq_all(ispa)))=='ac' &
                & .or. trim(adjustl(chraq_all(ispa)))=='mes' &
                & .or. trim(adjustl(chraq_all(ispa)))=='im' &
                & .or. trim(adjustl(chraq_all(ispa)))=='tea' &
                ! & .or. trim(adjustl(chraq_all(ispa)))=='oxa' &
                & ) then 
                
                ! account for speces associated with H+ (X + H+ = XH+)
                do ispa_h = 1,2
                    if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                        rspa_h = real(ispa_h,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)+rspa_h))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**rspa_h/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_h/gamma(ic2,:) &
                                & + gamma(ic1,:)*rspa_h*gamma(1,:)**(rspa_h-1d0)*dgamma_dios(1,:) &
                                &   /gamma(ic2,:) &
                                & + gamma(ic1,:)*gamma(1,:)**rspa_h*(-1d0) &
                                &   /gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq = gamma(1,:)**rspa_h/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + rspa_h*gamma(1,:)**(rspa_h-1d0)*dgamma_dios(1,:) &
                                &   /gamma(ic2,:) &
                                & + gamma(1,:)**rspa_h*(-1d0) &
                                &   /gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**rspa_h
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_h &
                                & + gamma(ic1,:)*rspa_h*gamma(1,:)**(rspa_h-1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq = gamma(1,:)**rspa_h
                            dfkeq_dios = ( &
                                & + rspa_h*gamma(1,:)**(rspa_h-1d0)*dgamma_dios(1,:) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        fkeqaq_h(ispa,ispa_h,:) = fkeq
                        f1 = f1 + (base_charge(ispa) + rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add)
                        df1 = df1 + ( & 
                            & + (base_charge(ispa) + rspa_h)*fkeq &
                            &        *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(rspa_h+ss_add)*prox**(rspa_h+ss_add-1d0) &
                            & )
                        d2f1 = d2f1 + ( & 
                            & + (base_charge(ispa) + rspa_h)*fkeq &
                            &        *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(rspa_h+ss_add)*(rspa_h+ss_add-1d0)*prox**(rspa_h+ss_add-2d0) &
                            & )
                        df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + (& 
                            & + (base_charge(ispa) + rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(rspa_h+ss_add) &
                            & )
                        df1df2 = df1df2 + ( &
                            & + (base_charge(ispa) + rspa_h)*dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add) &
                            & )
                        f2 = f2 + (base_charge(ispa) + rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add)
                        df2df1 = df2df1 + ( & 
                            & + (base_charge(ispa) + rspa_h)**2d0*fkeq &
                            &        *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(rspa_h+ss_add)*prox**(rspa_h+ss_add-1d0) &
                            & )
                        df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + (& 
                            & + (base_charge(ispa) + rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(rspa_h+ss_add) &
                            & )
                        df2 = df2 + (base_charge(ispa) + rspa_h)**2d0*dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add)
                        if (print_res) then 
                            write(chrint,'(I1)') ispa_h
                            write(88,'(A11)', advance='no') 'h'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                            write(99,'(A11)', advance='no') 'h'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                        endif 
                        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                            write(chrint,'(I1)') ispa_h
                            print'("nan found f1 and df1: point 6 | ",A11)', 'h'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                        endif 
                        
                    endif 
                enddo
            ! oxalic acid
            elseif ( &
                & trim(adjustl(chraq_all(ispa)))=='oxa' &
                & .or. trim(adjustl(chraq_all(ispa)))=='glp' &
                & ) then
                do ispa_h = 1,2
                    if (ispa_h==1) then  ! OxaH- = Oxa= + H+ 
                        if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                            rspa_h = real(ispa_h,kind=8)
                            fkeq = 1d0/gamma(2,:)
                            dfkeq_dios = -1d0/gamma(2,:)**2d0*dgamma_dios(2,:)
                            fkeqaq_h(ispa,ispa_h,:) = fkeq
                            f1 = f1 + (base_charge(ispa) - rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h)
                            df1 = df1 + ( &
                                & + (base_charge(ispa) - rspa_h)*fkeq &
                                &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*prox**(ss_add-rspa_h-1d0) &
                                & )
                            d2f1 = d2f1 + ( &
                                & + (base_charge(ispa) - rspa_h)*fkeq &
                                &   *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*(ss_add-rspa_h-1d0)*prox**(ss_add-rspa_h-2d0) &
                                & )
                            df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( &
                                & + (base_charge(ispa) - rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(ss_add-rspa_h) &
                                & )
                            df1df2 = df1df2 + ( &
                                & + (base_charge(ispa) - rspa_h)*dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h) &
                                & )
                            f2 = f2 + (base_charge(ispa) - rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h)
                            df2df1 = df2df1 + ( &
                                & + (base_charge(ispa) - rspa_h)**2d0*fkeq &
                                &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*prox**(ss_add-rspa_h-1d0) &
                                & )
                            df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( &
                                & + (base_charge(ispa) - rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(ss_add-rspa_h) &
                                & )
                            df2 = df2 + ( &
                                & + (base_charge(ispa) - rspa_h)**2d0*dfkeq_dios &
                                &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h) &
                                & )
                            if (print_res) then 
                                write(chrint,'(I1)') ispa_h
                                write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(oh)'//trim(adjustl(chrint))
                                write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(oh)'//trim(adjustl(chrint))
                            endif 
                            if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                                write(chrint,'(I1)') ispa_h
                                print'("nan found f1 and df1: point 7 | ",A11)' &
                                    & ,trim(adjustl(chraq_all(ispa)))//'(oh)'//trim(adjustl(chrint))
                            endif 
                        endif 
                    elseif (ispa_h==2) then  ! OxaH- + H+ = OxaH2  
                        if ( keqaq_h(ispa,ispa_h) > 0d0) then  
                            rspa_h = real(ispa_h-1,kind=8)
                            fkeq = gamma(1,:)**2d0
                            dfkeq_dios = 2d0*gamma(1,:)*dgamma_dios(1,:)
                            fkeqaq_h(ispa,ispa_h,:)=fkeq
                            f1 = f1 + (base_charge(ispa) + rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add)
                            df1 = df1 + ( & 
                                & + (base_charge(ispa) + rspa_h)*fkeq &
                                &        *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(rspa_h+ss_add)*prox**(rspa_h+ss_add-1d0) &
                                & )
                            d2f1 = d2f1 + ( & 
                                & + (base_charge(ispa) + rspa_h)*fkeq &
                                &   *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(rspa_h+ss_add)*(rspa_h+ss_add-1d0)*prox**(rspa_h+ss_add-2d0) &
                                & )
                            df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + (& 
                                & + (base_charge(ispa) + rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(rspa_h+ss_add) &
                                & )
                            df1df2 = df1df2 + ( &
                                & + (base_charge(ispa) + rspa_h)*dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add) &
                                & )
                            f2 = f2 + (base_charge(ispa) + rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add)
                            df2df1 = df2df1 + ( & 
                                & + (base_charge(ispa) + rspa_h)**2d0*fkeq &
                                &        *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(rspa_h+ss_add)*prox**(rspa_h+ss_add-1d0) &
                                & )
                            df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + (& 
                                & + (base_charge(ispa) + rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(rspa_h+ss_add) &
                                & )
                            df2 = df2 + ( &
                                & + (base_charge(ispa) + rspa_h)**2d0*dfkeq_dios &
                                &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(rspa_h+ss_add) &
                                & )
                            if (print_res) then 
                                write(chrint,'(I1)') ispa_h-1
                                write(88,'(A11)', advance='no') 'h'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                                write(99,'(A11)', advance='no') 'h'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                            endif 
                            if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                                write(chrint,'(I1)') ispa_h-1
                                print'("nan found f1 and df1: point 8 | ",A11)', 'h'//trim(adjustl(chrint))//trim(adjustl(chraq_all(ispa)))
                            endif 
                            
                        endif 
                    endif 
                enddo
            ! cations
            else 
                ! account for hydrolysis speces (X + H2O = XOH- + H+)
                do ispa_h = 1,4
                    if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                        rspa_h = real(ispa_h,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-rspa_h))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq = gamma(ic1,:)/gamma(ic2,:)/gamma(1,:)**rspa_h
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)/gamma(ic2,:)/gamma(1,:)**rspa_h &
                                & + gamma(ic1,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:)/gamma(1,:)**rspa_h &
                                & + gamma(ic1,:)/gamma(ic2,:)*(-rspa_h)/gamma(1,:)**(rspa_h+1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq = 1d0/gamma(ic2,:)/gamma(1,:)**rspa_h
                            dfkeq_dios = ( &
                                & + 1d0*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:)/gamma(1,:)**rspa_h &
                                & + 1d0/gamma(ic2,:)*(-rspa_h)/gamma(1,:)**(rspa_h+1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq = gamma(ic1,:)/gamma(1,:)**rspa_h
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)/gamma(1,:)**rspa_h &
                                & + gamma(ic1,:)*(-rspa_h)/gamma(1,:)**(rspa_h+1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq = 1d0/gamma(1,:)**rspa_h
                            dfkeq_dios = ( &
                                & + 1d0*(-rspa_h)/gamma(1,:)**(rspa_h+1d0)*dgamma_dios(1,:) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        fkeqaq_h(ispa,ispa_h,:)=fkeq
                        f1 = f1 + (base_charge(ispa) - rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h)
                        df1 = df1 + ( &
                            & + (base_charge(ispa) - rspa_h)*fkeq &
                            &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*prox**(ss_add-rspa_h-1d0) &
                            & )
                        d2f1 = d2f1 + ( &
                            & + (base_charge(ispa) - rspa_h)*fkeq &
                            &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*(ss_add-rspa_h-1d0)*prox**(ss_add-rspa_h-2d0) &
                            & )
                        df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( &
                            & + (base_charge(ispa) - rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(ss_add-rspa_h) &
                            & )
                        df1df2 = df1df2 + ( & 
                            & + (base_charge(ispa) - rspa_h)*dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h) &
                            & )
                        f2 = f2 + (base_charge(ispa) - rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h)
                        df2df1 = df2df1 + ( &
                            & + (base_charge(ispa) - rspa_h)**2d0*fkeq &
                            &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*prox**(ss_add-rspa_h-1d0) &
                            & )
                        df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( &
                            & + (base_charge(ispa) - rspa_h)**2d0*fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**(ss_add-rspa_h) &
                            & )
                        df2 = df2 + (base_charge(ispa) - rspa_h)**2d0*dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h)
                        if (print_res) then 
                            write(chrint,'(I1)') ispa_h
                            write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(oh)'//trim(adjustl(chrint))
                            write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(oh)'//trim(adjustl(chrint))
                        endif 
                        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                            write(chrint,'(I1)') ispa_h
                            print'("nan found f1 and df1: point 9 | ",A11)', trim(adjustl(chraq_all(ispa)))//'(oh)'//trim(adjustl(chrint))
                            print* &
                                & , (base_charge(ispa) - rspa_h)*fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**(ss_add-rspa_h) &
                                & , (base_charge(ispa) - rspa_h)*fkeq &
                                &   *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(ss_add-rspa_h)*(ss_add-rspa_h-1d0)*prox**(ss_add-rspa_h-2d0) 
                        endif 
                    endif 
                enddo 
                ! account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                do ispa_c = 1,2
                    if ( keqaq_c(ispa,ispa_c) > 0d0) then 
                        if (ispa_c == 1) then ! with CO3-- (e.g., Mg2+ + CO32- = MgCO3 )
                            ic1 = nint(abs(base_charge(ispa)))
                            ic2 = nint(abs(base_charge(ispa)-2d0))
                            if ( ic1>0 .and. ic2 > 0) then  
                                fkeq = gamma(ic1,:)*gamma(2,:)/gamma(ic2,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(ic1,:)*gamma(2,:)/gamma(ic2,:) &
                                    & + gamma(ic1,:)*dgamma_dios(2,:)/gamma(ic2,:) &
                                    & + gamma(ic1,:)*gamma(2,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                    & )
                            elseif ( ic1==0 .and. ic2 > 0) then  
                                fkeq = gamma(2,:)/gamma(ic2,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(2,:)/gamma(ic2,:) &
                                    & + gamma(2,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                    & )
                            elseif ( ic1>0 .and. ic2 == 0) then  
                                fkeq = gamma(ic1,:)*gamma(2,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(ic1,:)*gamma(2,:) &
                                    & + gamma(ic1,:)*dgamma_dios(2,:) &
                                    & )
                            elseif ( ic1==0 .and. ic2 == 0) then  
                                fkeq = gamma(2,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(2,:) &
                                    & )
                            else    
                                print *, 'something is wrong'
                                stop
                            endif 
                            fkeqaq_c(ispa,ispa_c,:) = fkeq
                            f1 = f1 + (base_charge(ispa)-2d0)*fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-2d0)
                            df1 = df1 + ( & 
                                & + (base_charge(ispa)-2d0)*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(ss_add-2d0)*prox**(ss_add-3d0) &
                                & )
                            d2f1 = d2f1 + ( & 
                                & + (base_charge(ispa)-2d0)*fkeq &
                                & *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(ss_add-2d0)*(ss_add-3d0)*prox**(ss_add-4d0) &
                                & )
                            df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( & 
                                & + (base_charge(ispa)-2d0)*fkeq*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x*prox**(ss_add-2d0) &
                                & )
                            df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( & 
                                & + (base_charge(ispa)-2d0)*fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0*prox**(ss_add-2d0) &
                                & )
                            df1df2 = df1df2 + ( &
                                & + (base_charge(ispa)-2d0)*dfkeq_dios &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-2d0) &
                                & )
                            f2 = f2 + ( &
                                & + (base_charge(ispa)-2d0)**2d0*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-2d0) &
                                & )
                            df2df1 = df2df1 + ( & 
                                & + (base_charge(ispa)-2d0)**2d0*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(ss_add-2d0)*prox**(ss_add-3d0) &
                                & )
                            df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( & 
                                & + (base_charge(ispa)-2d0)**2d0*fkeq*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x*prox**(ss_add-2d0) &
                                & )
                            df2dmgas(ipco2,:) = df2dmgas(ipco2,:) + ( & 
                                & + (base_charge(ispa)-2d0)**2d0*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0*prox**(ss_add-2d0) &
                                & )
                            df2 = df2 + ( &
                                & + (base_charge(ispa)-2d0)**2d0*dfkeq_dios &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-2d0) &
                                & )
                            if (print_res) write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(co3)'
                            if (print_res) write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(co3)'
                            if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                                print*,'nan found f1 and df1: point 10 | '//trim(adjustl(chraq_all(ispa)))//'(co3)'
                            endif 
                        elseif (ispa_c == 2) then ! with HCO3- (e.g., Mg2+ + H+ + CO32- = MgHCO3+ )
                            ic1 = nint(abs(base_charge(ispa)))
                            ic2 = nint(abs(base_charge(ispa)-1d0))
                            if ( ic1>0 .and. ic2 > 0) then  
                                fkeq = gamma(ic1,:)*gamma(2,:)*gamma(1,:)/gamma(ic2,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(ic1,:)*gamma(2,:)*gamma(1,:)/gamma(ic2,:) &
                                    & + gamma(ic1,:)*dgamma_dios(2,:)*gamma(1,:)/gamma(ic2,:) &
                                    & + gamma(ic1,:)*gamma(2,:)*dgamma_dios(1,:)/gamma(ic2,:) &
                                    & + gamma(ic1,:)*gamma(2,:)*gamma(1,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                    & )
                            elseif ( ic1==0 .and. ic2 > 0) then  
                                fkeq = gamma(2,:)*gamma(1,:)/gamma(ic2,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(2,:)*gamma(1,:)/gamma(ic2,:) &
                                    & + gamma(2,:)*dgamma_dios(1,:)/gamma(ic2,:) &
                                    & + gamma(2,:)*gamma(1,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                    & )
                            elseif ( ic1>0 .and. ic2 == 0) then  
                                fkeq = gamma(ic1,:)*gamma(2,:)*gamma(1,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(ic1,:)*gamma(2,:)*gamma(1,:) &
                                    & + gamma(ic1,:)*dgamma_dios(2,:)*gamma(1,:) &
                                    & + gamma(ic1,:)*gamma(2,:)*dgamma_dios(1,:) &
                                    & )
                            elseif ( ic1==0 .and. ic2 == 0) then  
                                fkeq = gamma(2,:)*gamma(1,:)
                                dfkeq_dios = ( &
                                    & + dgamma_dios(2,:)*gamma(1,:) &
                                    & + gamma(2,:)*dgamma_dios(1,:) &
                                    & )
                            else    
                                print *, 'something is wrong'
                                stop
                            endif 
                            fkeqaq_c(ispa,ispa_c,:) = fkeq
                            f1 = f1 + (base_charge(ispa)-1d0)*fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-1d0)
                            df1 = df1 + ( & 
                                & + (base_charge(ispa)-1d0)*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(ss_add-1d0)*prox**(ss_add-2d0) &
                                & )
                            d2f1 = d2f1 + ( & 
                                & + (base_charge(ispa)-1d0)*fkeq &
                                & *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(ss_add-1d0)*(ss_add-2d0)*prox**(ss_add-3d0) &
                                & )
                            df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( & 
                                & + (base_charge(ispa)-1d0)*fkeq*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x*prox**(ss_add-1d0) &
                                & )
                            df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( & 
                                & + (base_charge(ispa)-1d0)*fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0*prox**(ss_add-1d0) &
                                & )
                            df1df2 = df1df2 + ( &
                                & + (base_charge(ispa)-1d0)*dfkeq_dios &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-1d0) &
                                & )
                            f2 = f2 + ( &
                                & + (base_charge(ispa)-1d0)**2d0*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-1d0) &
                                & )
                            df2df1 = df2df1 + ( & 
                                & + (base_charge(ispa)-1d0)**2d0*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(ss_add-1d0)*prox**(ss_add-2d0) &
                                & )
                            df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( & 
                                & + (base_charge(ispa)-1d0)**2d0*fkeq*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x*prox**(ss_add-1d0) &
                                & )
                            df2dmgas(ipco2,:) = df2dmgas(ipco2,:) + ( & 
                                & + (base_charge(ispa)-1d0)**2d0*fkeq &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0*prox**(ss_add-1d0) &
                                & )
                            df2 = df2 + ( &
                                & + (base_charge(ispa)-1d0)**2d0*dfkeq_dios &
                                &       *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*prox**(ss_add-1d0) &
                                & )
                            if (print_res) write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(hco3)'
                            if (print_res) write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(hco3)'
                            if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                                print*,'nan found f1 and df1: point 11 | '//trim(adjustl(chraq_all(ispa)))//'(hco3)'
                            endif 
                        endif 
                    endif 
                enddo 
                ! account for complexation with free SO4 (e.g., X + SO42- = XSO42-)
                do ispa_s = 1,2
                    if ( keqaq_s(ispa,ispa_s) > 0d0) then 
                        rspa_s = real(ispa_s,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-2d0*rspa_s))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq = gamma(ic1,:)*gamma(2,:)**rspa_s/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(2,:)**rspa_s/gamma(ic2,:) &
                                & + gamma(ic1,:)*rspa_s*gamma(2,:)**(rspa_s-1d0)*dgamma_dios(2,:)/gamma(ic2,:) &
                                & + gamma(ic1,:)*gamma(2,:)**rspa_s*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq = gamma(2,:)**rspa_s/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + rspa_s*gamma(2,:)**(rspa_s-1d0)*dgamma_dios(2,:)/gamma(ic2,:) &
                                & + gamma(2,:)**rspa_s*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq = gamma(ic1,:)*gamma(2,:)**rspa_s
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(2,:)**rspa_s &
                                & + gamma(ic1,:)*rspa_s*gamma(2,:)**(rspa_s-1d0)*dgamma_dios(2,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq = gamma(2,:)**rspa_s
                            dfkeq_dios = ( &
                                & + rspa_s*gamma(2,:)**(rspa_s-1d0)*dgamma_dios(2,:) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        fkeqaq_s(ispa,ispa_s,:) = fkeq
                        f1 = f1 + (base_charge(ispa)-2d0*rspa_s)*fkeq*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*prox**ss_add
                        df1 = df1 + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)*fkeq &
                            &       *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*ss_add*prox**(ss_add-1d0) & 
                            & )
                        d2f1 = d2f1 + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)*fkeq &
                            &       *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*ss_add*(ss_add-1d0)*prox**(ss_add-2d0) & 
                            & )
                        df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)*fkeq*keqaq_s(ispa,ispa_s)*1d0*so4f**rspa_s*prox**ss_add & 
                            & )
                        df1dmaqf(iso4,:) = df1dmaqf(iso4,:) + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)*fkeq*keqaq_s(ispa,ispa_s) &
                            & *maqf_loc(ispa,:)*rspa_s*so4f**(rspa_s-1d0)*prox**ss_add & 
                            & )
                        df1df2 = df1df2 + ( &
                            & + (base_charge(ispa)-2d0*rspa_s)*dfkeq_dios*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*prox**ss_add &
                            & )
                        f2 = f2 + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)**2d0*fkeq*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*prox**ss_add &
                            & )
                        df2df1 = df2df1 + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)**2d0*fkeq &
                            &       *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*ss_add*prox**(ss_add-1d0) & 
                            & )
                        df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)**2d0*fkeq*keqaq_s(ispa,ispa_s)*1d0*so4f**rspa_s*prox**ss_add & 
                            & )
                        df2dmaqf(iso4,:) = df2dmaqf(iso4,:) + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)**2d0*fkeq*keqaq_s(ispa,ispa_s) &
                            & *maqf_loc(ispa,:)*rspa_s*so4f**(rspa_s-1d0)*prox**ss_add & 
                            & )
                        df2 = df2 + ( & 
                            & + (base_charge(ispa)-2d0*rspa_s)**2d0*dfkeq_dios &
                            &       *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s*prox**ss_add &
                            & )
                        if (print_res) then 
                            write(chrint,'(I1)') ispa_s
                            write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(so4)'//trim(adjustl(chrint))
                            write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(so4)'//trim(adjustl(chrint))
                        endif 
                        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                            write(chrint,'(I1)') ispa_s
                            print*,'nan found f1 and df1: point 12 | '//trim(adjustl(chraq_all(ispa)))//'(so4)'//trim(adjustl(chrint))
                        endif 
                            
                    endif 
                enddo 
                ! accounting for complexation with free NO3 (e.g., X + NO3- = XNO3-)
                do ispa_no3 = 1,2
                    if ( keqaq_no3(ispa,ispa_no3) > 0d0) then 
                        rspa_no3 = real(ispa_no3,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-1d0*rspa_no3))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**rspa_no3/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_no3/gamma(ic2,:) &
                                & + gamma(ic1,:)*rspa_no3*gamma(1,:)**(rspa_no3-1d0)*dgamma_dios(1,:)/gamma(ic2,:) &
                                & + gamma(ic1,:)*gamma(1,:)**rspa_no3*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq = gamma(1,:)**rspa_no3/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + rspa_no3*gamma(1,:)**(rspa_no3-1d0)*dgamma_dios(1,:)/gamma(ic2,:) &
                                & + gamma(1,:)**rspa_no3*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**rspa_no3
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_no3 &
                                & + gamma(ic1,:)*rspa_no3*gamma(1,:)**(rspa_no3-1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq = gamma(1,:)**rspa_no3
                            dfkeq_dios = ( &
                                & + rspa_no3*gamma(1,:)**(rspa_no3-1d0)*dgamma_dios(1,:) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        fkeqaq_no3(ispa,ispa_no3,:) = fkeq
                        f1 = f1 + ( &
                            & + (base_charge(ispa)-1d0*rspa_no3)*fkeq &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*prox**ss_add &
                            & )
                        df1 = df1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)*fkeq &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*ss_add*prox**(ss_add-1d0) & 
                            & )
                        d2f1 = d2f1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)*fkeq &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*ss_add*(ss_add-1d0)*prox**(ss_add-2d0) & 
                            & )
                        df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)*fkeq*keqaq_no3(ispa,ispa_no3)*1d0*no3f**rspa_no3*prox**ss_add & 
                            & )
                        df1dmaqf(ino3,:) = df1dmaqf(ino3,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)*fkeq*keqaq_no3(ispa,ispa_no3) &
                            &   *maqf_loc(ispa,:)*rspa_no3*no3f**(rspa_no3-1d0)*prox**ss_add & 
                            & )
                        df1df2 = df1df2 + ( &
                            & + (base_charge(ispa)-1d0*rspa_no3)*dfkeq_dios &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*prox**ss_add &
                            & )
                        f2 = f2 + ( &
                            & + (base_charge(ispa)-1d0*rspa_no3)**2d0*fkeq &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*prox**ss_add &
                            & )
                        df2df1 = df2df1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)**2d0*fkeq &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*ss_add*prox**(ss_add-1d0) & 
                            & )
                        df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)**2d0*fkeq*keqaq_no3(ispa,ispa_no3)*1d0*no3f**rspa_no3*prox**ss_add & 
                            & )
                        df2dmaqf(ino3,:) = df2dmaqf(ino3,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_no3)**2d0*fkeq*keqaq_no3(ispa,ispa_no3) &
                            &   *maqf_loc(ispa,:)*rspa_no3*no3f**(rspa_no3-1d0)*prox**ss_add & 
                            & )
                        df2 = df2 + ( &
                            & + (base_charge(ispa)-1d0*rspa_no3)**2d0*dfkeq_dios &
                            &       *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3*prox**ss_add &
                            & )
                        if (print_res) then 
                            write(chrint,'(I1)') ispa_no3
                            write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(no3)'//trim(adjustl(chrint))
                            write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(no3)'//trim(adjustl(chrint))
                        endif 
                        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                            write(chrint,'(I1)') ispa_no3
                            print*,'nan found f1 and df1: point 13 | '//trim(adjustl(chraq_all(ispa)))//'(no3)'//trim(adjustl(chrint))
                        endif 
                            
                    endif 
                enddo 
                ! accounting for complexation with free Cl
                do ispa_cl = 1,2
                    if ( keqaq_cl(ispa,ispa_cl) > 0d0) then 
                        rspa_cl = real(ispa_cl,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-1d0*rspa_cl))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**rspa_cl/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_cl/gamma(ic2,:) &
                                & + gamma(ic1,:)*rspa_cl*gamma(1,:)**(rspa_cl-1d0)*dgamma_dios(1,:)/gamma(ic2,:) &
                                & + gamma(ic1,:)*gamma(1,:)**rspa_cl*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq = gamma(1,:)**rspa_cl/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + rspa_cl*gamma(1,:)**(rspa_cl-1d0)*dgamma_dios(1,:)/gamma(ic2,:) &
                                & + gamma(1,:)**rspa_cl*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**rspa_cl
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_cl &
                                & + gamma(ic1,:)*rspa_cl*gamma(1,:)**(rspa_cl-1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq = gamma(1,:)**rspa_cl
                            dfkeq_dios = ( &
                                & + rspa_cl*gamma(1,:)**(rspa_cl-1d0)*dgamma_dios(1,:) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        fkeqaq_cl(ispa,ispa_cl,:) = fkeq
                        f1 = f1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)*fkeq &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*prox**ss_add & 
                            & ) 
                        df1 = df1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)*fkeq &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*ss_add*prox**(ss_add-1d0) & 
                            & )
                        d2f1 = d2f1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)*fkeq &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*ss_add*(ss_add-1d0)*prox**(ss_add-2d0) & 
                            & )
                        df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)*fkeq*keqaq_cl(ispa,ispa_cl)*1d0*clf**rspa_cl*prox**ss_add & 
                            & )
                        df1dmaqf(icl,:) = df1dmaqf(icl,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)*fkeq*keqaq_cl(ispa,ispa_cl) &
                            &   *maqf_loc(ispa,:)*rspa_cl*clf**(rspa_cl-1d0)*prox**ss_add & 
                            & )
                        df1df2 = df1df2 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)*dfkeq_dios &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*prox**ss_add &
                            & )
                        f2 = f2 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)**2d0*fkeq &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*prox**ss_add &
                            & )
                        df2df1 = df2df1 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)**2d0*fkeq &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*ss_add*prox**(ss_add-1d0) & 
                            & )
                        df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)**2d0*fkeq*keqaq_cl(ispa,ispa_cl)*1d0*clf**rspa_cl*prox**ss_add & 
                            & )
                        df2dmaqf(icl,:) = df2dmaqf(icl,:) + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)**2d0*fkeq*keqaq_cl(ispa,ispa_cl) &
                            &   *maqf_loc(ispa,:)*rspa_cl*clf**(rspa_cl-1d0)*prox**ss_add & 
                            & )
                        df2 = df2 + ( & 
                            & + (base_charge(ispa)-1d0*rspa_cl)**2d0*dfkeq_dios &
                            &       *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl*prox**ss_add &
                            & )
                        if (print_res) then 
                            write(chrint,'(I1)') ispa_cl
                            write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(cl)'//trim(adjustl(chrint))
                            write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(cl)'//trim(adjustl(chrint))
                        endif 
                        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                            write(chrint,'(I1)') ispa_cl
                            print*,'nan found f1 and df1: point 14 | '//trim(adjustl(chraq_all(ispa)))//'(cl)'//trim(adjustl(chrint))
                        endif 
                            
                    endif 
                enddo 
                ! accounting for complexation with HOxa- 
                ! (e.g., Fe+3 + 2 OxaH- = Fe(Oxa)2- + 2 H+ )
                ! (Al3+ + H2O + HOxa- = Al(OH)Oxa + 2H+    | ispa_oxa=1)
                ! (Al3+ + 2H2O + HOxa- = Al(OH)2Oxa- + 3H+ | ispa_oxa=2)
                do ispa_oxa = 1,2
                    rspa_oxa   = real(ispa_oxa,kind=8)
                    rspa_oxa_2 = real(ispa_oxa,kind=8) * 2d0
                    rspa_oxa_3 = real(ispa_oxa,kind=8)
                    if (trim(adjustl(chraq_all(ispa)))=='al') then
                        rspa_oxa   = real(ispa_oxa,kind=8) + 1d0
                        rspa_oxa_2 = real(ispa_oxa,kind=8) + 2d0
                        rspa_oxa_3 = 1d0
                    endif 
                    if ( keqaq_oxa(ispa,ispa_oxa) > 0d0) then 
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-rspa_oxa_2))
                        if ( ic1>0 .and. ic2 > 0) then  
                            ! fkeq = gamma(ic1,:)*gamma(1,:)**rspa_oxa_3/gamma(ic2,:)/gamma(1,:)**rspa_oxa
                            fkeq = gamma(ic1,:)*gamma(1,:)**(rspa_oxa_3-rspa_oxa)/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**(rspa_oxa_3-rspa_oxa)/gamma(ic2,:) &
                                & + gamma(ic1,:)*(rspa_oxa_3-rspa_oxa)*gamma(1,:)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,:) &
                                &       /gamma(ic2,:) &
                                & + gamma(ic1,:)*gamma(1,:)**(rspa_oxa_3-rspa_oxa)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq = gamma(1,:)**(rspa_oxa_3-rspa_oxa)/gamma(ic2,:)
                            dfkeq_dios = ( &
                                & + (rspa_oxa_3-rspa_oxa)*gamma(1,:)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,:)/gamma(ic2,:) &
                                & + gamma(1,:)**(rspa_oxa_3-rspa_oxa)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq = gamma(ic1,:)*gamma(1,:)**(rspa_oxa_3-rspa_oxa)
                            dfkeq_dios = ( &
                                & + dgamma_dios(ic1,:)*gamma(1,:)**(rspa_oxa_3-rspa_oxa) &
                                & + gamma(ic1,:)*(rspa_oxa_3-rspa_oxa)*gamma(1,:)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,:) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then
                            fkeq = gamma(1,:)**(rspa_oxa_3-rspa_oxa)
                            dfkeq_dios = ( &
                                & + (rspa_oxa_3-rspa_oxa)*gamma(1,:)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,:) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        fkeqaq_oxa(ispa,ispa_oxa,:) = fkeq
                        f1 = f1 + (base_charge(ispa)-rspa_oxa_2)*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            ! & *maqf_loc(ispa,:)*oxaf**rspa_oxa*prox**ss_add
                            & *maqf_loc(ispa,:)*oxaf**rspa_oxa_3*prox**(ss_add-rspa_oxa)
                        df1 = df1 + ( & 
                            & + (base_charge(ispa)-rspa_oxa_2)*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            ! &   *maqf_loc(ispa,:)*oxaf**rspa_oxa*ss_add*prox**(ss_add-1d0) & 
                            &   *maqf_loc(ispa,:)*oxaf**rspa_oxa_3*(ss_add-rspa_oxa)*prox**(ss_add-rspa_oxa-1d0) & 
                            & )
                        d2f1 = d2f1 + ( & 
                            & + (base_charge(ispa)-rspa_oxa_2)*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            ! &   *maqf_loc(ispa,:)*oxaf**rspa_oxa*ss_add*(ss_add-1d0)*prox**(ss_add-2d0) & 
                            &   *maqf_loc(ispa,:)*oxaf**rspa_oxa_3 &
                            &   *(ss_add-rspa_oxa)*(ss_add-rspa_oxa-1d0)*prox**(ss_add-rspa_oxa-2d0) & 
                            & )
                        df1dmaqf(ispa,:) = df1dmaqf(ispa,:) + ( & 
                            ! & + (base_charge(ispa)-2d0*rspa_oxa)*keqaq_oxa(ispa,ispa_oxa)*1d0*oxaf**rspa_oxa*prox**ss_add & 
                            & + (base_charge(ispa)-rspa_oxa_2)*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            &   *1d0*oxaf**rspa_oxa_3*prox**(ss_add-rspa_oxa) & 
                            & )
                        df1dmaqf(ioxa,:) = df1dmaqf(ioxa,:) + ( & 
                            ! & + (base_charge(ispa)-2d0*rspa_oxa)*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*1d0**rspa_oxa*prox**ss_add & 
                            & + (base_charge(ispa)-rspa_oxa_2)*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            &   *maqf_loc(ispa,:)*rspa_oxa_3*oxaf**(rspa_oxa_3-1d0)*prox**(ss_add-rspa_oxa) & 
                            & )
                        df1df2 = df1df2 + (base_charge(ispa)-rspa_oxa_2)*dfkeq_dios*keqaq_oxa(ispa,ispa_oxa) &
                            ! & *maqf_loc(ispa,:)*oxaf**rspa_oxa*prox**ss_add
                            & *maqf_loc(ispa,:)*oxaf**rspa_oxa_3*prox**(ss_add-rspa_oxa)
                        f2 = f2 + (base_charge(ispa)-rspa_oxa_2)**2d0*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            ! & *maqf_loc(ispa,:)*oxaf**rspa_oxa*prox**ss_add
                            & *maqf_loc(ispa,:)*oxaf**rspa_oxa_3*prox**(ss_add-rspa_oxa)
                        df2df1 = df2df1 + ( & 
                            & + (base_charge(ispa)-rspa_oxa_2)**2d0*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            ! &   *maqf_loc(ispa,:)*oxaf**rspa_oxa*ss_add*prox**(ss_add-1d0) & 
                            &   *maqf_loc(ispa,:)*oxaf**rspa_oxa_3*(ss_add-rspa_oxa)*prox**(ss_add-rspa_oxa-1d0) & 
                            & )
                        df2dmaqf(ispa,:) = df2dmaqf(ispa,:) + ( & 
                            ! & + (base_charge(ispa)-2d0*rspa_oxa)*keqaq_oxa(ispa,ispa_oxa)*1d0*oxaf**rspa_oxa*prox**ss_add & 
                            & + (base_charge(ispa)-rspa_oxa_2)**2d0*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            &   *1d0*oxaf**rspa_oxa_3*prox**(ss_add-rspa_oxa) & 
                            & )
                        df2dmaqf(ioxa,:) = df2dmaqf(ioxa,:) + ( & 
                            ! & + (base_charge(ispa)-2d0*rspa_oxa)*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*1d0**rspa_oxa*prox**ss_add & 
                            & + (base_charge(ispa)-rspa_oxa_2)**2d0*fkeq*keqaq_oxa(ispa,ispa_oxa) &
                            &   *maqf_loc(ispa,:)*rspa_oxa_3*oxaf**(rspa_oxa_3-1d0)*prox**(ss_add-rspa_oxa) & 
                            & )
                        df2 = df2 + (base_charge(ispa)-rspa_oxa_2)**2d0*dfkeq_dios*keqaq_oxa(ispa,ispa_oxa) &
                            ! & *maqf_loc(ispa,:)*oxaf**rspa_oxa*prox**ss_add
                            & *maqf_loc(ispa,:)*oxaf**rspa_oxa_3*prox**(ss_add-rspa_oxa)
                        if (print_res) then 
                            write(chrint,'(I1)') ispa_oxa
                            write(88,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(oxa)'//trim(adjustl(chrint))
                            write(99,'(A11)', advance='no') trim(adjustl(chraq_all(ispa)))//'(oxa)'//trim(adjustl(chrint))
                        endif 
                        if ( debug.and.(any(isnan(f1)).or.any(isnan(df1))) ) then 
                            write(chrint,'(I1)') ispa_oxa
                            print*,'nan found f1 and df1: point 15 | '//trim(adjustl(chraq_all(ispa)))//'(oxa)'//trim(adjustl(chrint))
                        endif 
                            
                    endif 
                enddo 
            endif 
        enddo     

        ios_new = f2 + 2d0*iosx*prox**(ss_add)
        ios_new = 0.5d0*ios_new/prox**(ss_add)

        ! Note (3/31/2023): ios_new should be independent of iosx (input) because the term 2d0*iosx*prox**(ss_add) was subtracted initially

        if (print_res) write(88,'(A11)', advance='no') 'I'
        if (print_res) write(99,'(A11)', advance='no') 'I'
        if (print_res) write(88,'(A11)') 'tot_charge'
        if (print_res) write(99,'(A11)') 'tot_charge'

        f1_chk = 0d0
        ss_add = 0d0

        fkeq = 1d0/gamma(2,:)
        fkw = 1d0/gamma(1,:)/gamma(1,:) ! H2O = H+ + OH- <--> Kw = {H+}{OH-} <--> Kw/gamma/gamma = [H+][OH-]
        if (print_res) then
            do iz = 1, nz
                f1_chk(iz) = f1_chk(iz) + prox(iz)**(ss_add(iz)+1d0) - fkw(iz)*kw*prox(iz)**(ss_add(iz)-1d0)
                write(88,'(3E25.16)', advance='no') z(iz),prox(iz), fkw(iz)*kw/prox(iz)
                write(99,'(3E25.16)', advance='no') z(iz),prox(iz), -fkw(iz)*kw/prox(iz)

                ! adding charges coming from aq species in eq with gases
                ! pCO2
                f1_chk(iz) = f1_chk(iz)  -  fkw(iz)*k1*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-1d0) &
                    & -  2d0*fkw(iz)*fkeq(iz)*k2*k1*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-2d0)
                write(88,'(2E25.16)', advance='no')    fkw(iz)*k1*kco2*pco2x(iz)/prox(iz) &
                    & ,  fkw(iz)*fkeq(iz)*k2*k1*kco2*pco2x(iz)/prox(iz)**2d0
                write(99,'(2E25.16)', advance='no')  - fkw(iz)*k1*kco2*pco2x(iz)/prox(iz)  &
                    & , -2d0*fkw(iz)*fkeq(iz)*k2*k1*kco2*pco2x(iz)/prox(iz)**2d0
                ! pNH3
                f1_chk(iz) = f1_chk(iz)  +  pnh3x(iz)*knh3/k1nh3*prox(iz)**(ss_add(iz)+1d0)
                write(88,'(E25.16)', advance='no')    pnh3x(iz)*knh3/k1nh3*prox(iz)
                write(99,'(E25.16)', advance='no')    pnh3x(iz)*knh3/k1nh3*prox(iz)

                do ispa = 1, nsp_aq_all
                    
                    f1_chk(iz) = f1_chk(iz) + base_charge(ispa)*maqf_loc(ispa,iz)*prox(iz)**(ss_add(iz))
                    write(88,'(E25.16)', advance='no') maqf_loc(ispa,iz) 
                    write(99,'(E25.16)', advance='no') base_charge(ispa)*maqf_loc(ispa,iz) 
                    
                    ! account for speces associated with NH4+ (both anions and cations)
                    do ispa_nh3 = 1,2
                        if ( keqaq_nh3(ispa,ispa_nh3) > 0d0) then 
                            rspa_nh3 = real(ispa_nh3,kind=8)
                            f1_chk(iz) = f1_chk(iz) &
                                & + (base_charge(ispa) + rspa_nh3)*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                                & *(pnh3x(iz)*knh3/k1nh3)**rspa_nh3*prox(iz)**(rspa_nh3+ss_add(iz)) &
                                & *fkeqaq_nh3(ispa,ispa_nh3,iz)
                            write(88,'(E25.16)', advance='no') keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                                & *(pnh3x(iz)*knh3/k1nh3)**rspa_nh3*prox(iz)**rspa_nh3&
                                & *fkeqaq_nh3(ispa,ispa_nh3,iz)
                            write(99,'(E25.16)', advance='no') (base_charge(ispa) + rspa_nh3) &
                                & *keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                                & *(pnh3x(iz)*knh3/k1nh3)**rspa_nh3*prox(iz)**rspa_nh3&
                                & *fkeqaq_nh3(ispa,ispa_nh3,iz)
                        endif 
                    enddo 
                    
                    ! annions
                    if ( &
                        & trim(adjustl(chraq_all(ispa)))=='no3' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='so4' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='cl' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='ac' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='mes' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='im' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='tea' &
                        ! & .or. trim(adjustl(chraq_all(ispa)))=='oxa' &
                        & ) then 
                        
                        ! account for speces associated with H+
                        do ispa_h = 1,2
                            if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                                rspa_h = real(ispa_h,kind=8)
                                f1_chk(iz) = f1_chk(iz) &
                                    & + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**(rspa_h+ss_add(iz)) &
                                    & *fkeqaq_h(ispa,ispa_h,iz)
                                write(88,'(E25.16)', advance='no') keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**rspa_h &
                                    & *fkeqaq_h(ispa,ispa_h,iz)
                                write(99,'(E25.16)', advance='no') (base_charge(ispa) + rspa_h) &
                                    & *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**rspa_h &
                                    & *fkeqaq_h(ispa,ispa_h,iz)
                            endif 
                        enddo 
                    ! oxalic acid 
                    elseif ( &
                        & trim(adjustl(chraq_all(ispa)))=='oxa' &
                        & .or. trim(adjustl(chraq_all(ispa)))=='glp' &
                        & ) then 
                        do ispa_h = 1,2
                            if (ispa_h==1) then 
                                if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                                    rspa_h = real(ispa_h,kind=8)
                                    f1_chk(iz) = f1_chk(iz) &
                                        & + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h) &
                                        &       *maqf_loc(ispa,iz)*prox(iz)**(ss_add(iz)-rspa_h)  &
                                        &       *fkeqaq_h(ispa,ispa_h,iz)
                                    write(88,'(E25.16)', advance='no') keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)/prox(iz)**rspa_h &
                                        & *fkeqaq_h(ispa,ispa_h,iz)
                                    write(99,'(E25.16)', advance='no') (base_charge(ispa) - rspa_h) &
                                        & *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)/prox(iz)**rspa_h &
                                        & *fkeqaq_h(ispa,ispa_h,iz)
                                endif 
                            elseif (ispa_h==2)then
                                if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                                    rspa_h = real(ispa_h-1,kind=8)
                                    f1_chk(iz) = f1_chk(iz) &
                                        & + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h) &
                                        &       *maqf_loc(ispa,iz)*prox(iz)**(rspa_h+ss_add(iz)) &
                                        &       *fkeqaq_h(ispa,ispa_h,iz)
                                    write(88,'(E25.16)', advance='no') keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**rspa_h &
                                        & *fkeqaq_h(ispa,ispa_h,iz)
                                    write(99,'(E25.16)', advance='no') (base_charge(ispa) + rspa_h) &
                                        & *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**rspa_h &
                                        & *fkeqaq_h(ispa,ispa_h,iz)
                                endif 
                            endif 
                        enddo 
                    ! cations
                    else 
                        ! account for hydrolysis speces
                        do ispa_h = 1,4
                            if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                                rspa_h = real(ispa_h,kind=8)
                                f1_chk(iz) = f1_chk(iz) &
                                    & + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**(ss_add(iz)-rspa_h) &
                                    & *fkeqaq_h(ispa,ispa_h,iz)
                                write(88,'(E25.16)', advance='no') keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)/prox(iz)**rspa_h &
                                    & *fkeqaq_h(ispa,ispa_h,iz)
                                write(99,'(E25.16)', advance='no') (base_charge(ispa) - rspa_h) &
                                    & *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)/prox(iz)**rspa_h &
                                    & *fkeqaq_h(ispa,ispa_h,iz)
                            endif 
                        enddo 
                        ! account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                        do ispa_c = 1,2
                            if ( keqaq_c(ispa,ispa_c) > 0d0) then 
                                if (ispa_c == 1) then ! with CO3--
                                    f1_chk(iz) = f1_chk(iz) + (base_charge(ispa)-2d0) &
                                        & *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-2d0) &
                                        & *fkeqaq_c(ispa,ispa_c,iz)
                                    write(88,'(E25.16)', advance='no') &
                                        & keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
                                        & *fkeqaq_c(ispa,ispa_c,iz)
                                    write(99,'(E25.16)', advance='no') (base_charge(ispa)-2d0) &
                                        & *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
                                        & *fkeqaq_c(ispa,ispa_c,iz)
                                elseif (ispa_c == 2) then ! with HCO3-
                                    f1_chk(iz) = f1_chk(iz) + (base_charge(ispa)-1d0) &
                                        & *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-1d0) &
                                        & *fkeqaq_c(ispa,ispa_c,iz)
                                    write(88,'(E25.16)', advance='no') &
                                        & keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz) &
                                        & *fkeqaq_c(ispa,ispa_c,iz)
                                    write(99,'(E25.16)', advance='no') (base_charge(ispa)-1d0) &
                                        & *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz) &
                                        & *fkeqaq_c(ispa,ispa_c,iz)
                                endif 
                            endif 
                        enddo 
                        ! account for complexation with free SO4
                        do ispa_s = 1,2
                            if ( keqaq_s(ispa,ispa_s) > 0d0) then 
                                rspa_s = real(ispa_s,kind=8)
                                f1_chk(iz) = f1_chk(iz)  + (base_charge(ispa)-2d0*rspa_s) &
                                    & *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)**rspa_s*prox(iz)**ss_add(iz) &
                                    & *fkeqaq_s(ispa,ispa_s,iz)
                                write(88,'(E25.16)', advance='no') keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)**rspa_s &
                                    & *fkeqaq_s(ispa,ispa_s,iz)
                                write(99,'(E25.16)', advance='no') (base_charge(ispa)-2d0*rspa_s) &
                                    & *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)**rspa_s &
                                    & *fkeqaq_s(ispa,ispa_s,iz)
                            endif 
                        enddo 
                        ! account for complexation with free NO3
                        do ispa_no3 = 1,2
                            if ( keqaq_no3(ispa,ispa_no3) > 0d0) then 
                                rspa_no3 = real(ispa_no3,kind=8)
                                f1_chk(iz) = f1_chk(iz)  + (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                                    & *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*no3f(iz)**rspa_no3*prox(iz)**ss_add(iz) &
                                    & *fkeqaq_no3(ispa,ispa_no3,iz)
                                write(88,'(E25.16)', advance='no') keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*no3f(iz)**rspa_no3 &
                                    & *fkeqaq_no3(ispa,ispa_no3,iz)
                                write(99,'(E25.16)', advance='no') (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                                    & *keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*no3f(iz)**rspa_no3 &
                                    & *fkeqaq_no3(ispa,ispa_no3,iz)
                            endif 
                        enddo 
                        ! account for complexation with free Cl
                        do ispa_cl = 1,2
                            if ( keqaq_cl(ispa,ispa_cl) > 0d0) then 
                                rspa_cl = real(ispa_cl,kind=8)
                                f1_chk(iz) = f1_chk(iz)  + (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                                    & *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*clf(iz)**rspa_cl*prox(iz)**ss_add(iz) &
                                    & *fkeqaq_cl(ispa,ispa_cl,iz)
                                write(88,'(E25.16)', advance='no') keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*clf(iz)**rspa_cl &
                                    & *fkeqaq_cl(ispa,ispa_cl,iz)
                                write(99,'(E25.16)', advance='no') (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                                    & *keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*clf(iz)**rspa_cl &
                                    & *fkeqaq_cl(ispa,ispa_cl,iz)
                            endif 
                        enddo 
                        ! account for complexation with Hoxa-
                        do ispa_oxa = 1,2
                            rspa_oxa   = real(ispa_oxa,kind=8)
                            rspa_oxa_2 = real(ispa_oxa,kind=8) * 2d0
                            rspa_oxa_3 = real(ispa_oxa,kind=8)
                            if (trim(adjustl(chraq_all(ispa)))=='al') then
                                rspa_oxa   = real(ispa_oxa,kind=8) + 1d0
                                rspa_oxa_2 = real(ispa_oxa,kind=8) + 2d0
                                rspa_oxa_3 = 1d0
                            endif 
                            if ( keqaq_oxa(ispa,ispa_oxa) > 0d0) then 
                                f1_chk(iz) = f1_chk(iz)  + (base_charge(ispa)-rspa_oxa_2) &
                                    & *keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,iz)*oxaf(iz)**rspa_oxa_3*prox(iz)**(ss_add(iz)-rspa_oxa) &
                                    & *fkeqaq_oxa(ispa,ispa_oxa,iz)
                                write(88,'(E25.16)', advance='no') &
                                    & keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,iz)*oxaf(iz)**rspa_oxa_3/prox(iz)**rspa_oxa &
                                    & *fkeqaq_oxa(ispa,ispa_oxa,iz)
                                write(99,'(E25.16)', advance='no') (base_charge(ispa)-rspa_oxa_2) &
                                    & *keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,iz)*oxaf(iz)**rspa_oxa_3/prox(iz)**rspa_oxa &
                                    & *fkeqaq_oxa(ispa,ispa_oxa,iz)
                            endif 
                        enddo 
                    endif 
                enddo     
                ! case to save the input ionic strength, which is zero when act_ON = .false.
                ! write(88,'(E25.16)', advance='no') iosx(iz)
                ! write(99,'(E25.16)', advance='no') iosx(iz)
                ! case to save the newly calculated input ionic strength regardless of act_ON = .true. or .false.
                write(88,'(E25.16)', advance='no') ios_new(iz)
                write(99,'(E25.16)', advance='no') ios_new(iz)
                write(88,'(E25.16)') f1_chk(iz)
                write(99,'(E25.16)') f1_chk(iz)
            enddo 
        endif 

        if (print_res) close(88)
        if (print_res) close(99)

    endsubroutine calc_charge_balance

    ! Charge balance calculation at a point
    subroutine calc_charge_balance_point( &
        & nz,nsp_aq_all,nsp_gas_all &
        & ,chraq_all,chrgas_all &
        & ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &
        & ,base_charge &
        & ,mgasx_loc,maqf_loc &
        & ,z,prox,iz,iosx,tc &
        & ,print_loc,print_res,ph_add_order &
        & ,f1,df1,df1dmaqf,df1dmgas &!output
        & ,d2f1,d2f1dmaqf,d2f1dmgas &!output
        & )
        implicit none

        integer,intent(in)::nz,nsp_aq_all,nsp_gas_all,iz
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        real(kind=8),intent(in)::kw,tc
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl
        real(kind=8),dimension(nsp_gas_all,nz),intent(in)::mgasx_loc
        real(kind=8),dimension(nsp_aq_all,nz),intent(in)::maqf_loc
        real(kind=8),dimension(nsp_aq_all),intent(in)::base_charge
        real(kind=8),dimension(nz),intent(in)::z,prox,ph_add_order,iosx
        real(kind=8),dimension(nz),intent(out)::f1,df1,d2f1
        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::df1dmaqf,d2f1dmaqf
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::df1dmgas,d2f1dmgas

        logical,intent(in)::print_res
        character(500),intent(in)::print_loc

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ispa,ispa_h,ispa_c,ispa_s,ipco2,ipnh3,iso4,ioxa,ispa_no3,ino3,ispa_nh3,ispa_oxa,icl,ispa_cl

        real(kind=8) kco2,k1,k2,knh3,k1nh3,rspa_h,rspa_s,rspa_no3,rspa_nh3,rspa_oxa,rspa_oxa_2,rspa_oxa_3 &
            & ,rspa_cl
        real(kind=8),dimension(nz)::pco2x,pnh3x,so4f,no3f,oxaf,clf
        real(kind=8),dimension(nz)::f1_chk,ss_add,back

        integer icharge,ic1,ic2
        ! real(kind=8) tc
        real(kind=8) rcharge
        real(kind=8),dimension(nz)::gamma_tmp,dgamma_dios_tmp,fkw,fkeq,dfkw_dios,dfkeq_dios 
        real(kind=8),dimension(4,nz)::gamma,dgamma_dios 

        character(1) chrint

        ipco2   = findloc(chrgas_all,'pco2',dim=1)
        ipnh3   = findloc(chrgas_all,'pnh3',dim=1)
        iso4    = findloc(chraq_all,'so4',dim=1)
        ino3    = findloc(chraq_all,'no3',dim=1)
        ioxa    = findloc(chraq_all,'oxa',dim=1)
        icl     = findloc(chraq_all,'cl',dim=1)

        kco2    = keqgas_h(ipco2,ieqgas_h0)
        k1      = keqgas_h(ipco2,ieqgas_h1)
        k2      = keqgas_h(ipco2,ieqgas_h2)

        pco2x   = mgasx_loc(ipco2,:)


        knh3    = keqgas_h(ipnh3,ieqgas_h0)
        k1nh3   = keqgas_h(ipnh3,ieqgas_h1)

        pnh3x   = mgasx_loc(ipnh3,:)

        so4f    = maqf_loc(iso4,:)
        no3f    = maqf_loc(ino3,:)
        oxaf    = maqf_loc(ioxa,:)
        clf     = maqf_loc(icl,:)

        ss_add = ph_add_order

        f1 = 0d0
        df1 = 0d0
        d2f1 = 0d0
        df1dmaqf = 0d0
        df1dmgas = 0d0
        d2f1dmaqf = 0d0
        d2f1dmgas = 0d0

        back = 1d0
        back = 0d0

        do icharge=1,4
            rcharge = 1d0*icharge
            call calc_gamma_davies(  &
                & nz,iosx,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma(icharge,:)=gamma_tmp(:)
            dgamma_dios(icharge,:)=dgamma_dios_tmp(:)
        enddo

        fkw(iz) = 1d0/gamma(1,iz)/gamma(1,iz) ! H2O = H+ + OH- <--> Kw = {H+}{OH-} <--> Kw/gamma/gamma = [H+][OH-]
        dfkw_dios(iz) = 1d0*(-2d0)*gamma(1,iz)**(-3d0)*dgamma_dios(1,iz)

        f1(iz) = f1(iz) + prox(iz)**(ss_add(iz)+1d0) - fkw(iz)*kw*prox(iz)**(ss_add(iz)-1d0) &
            & + back(iz)*prox(iz)**(ss_add(iz))- back(iz)*prox(iz)**(ss_add(iz))
        df1(iz) = df1(iz) + (ss_add(iz)+1d0)*prox(iz)**(ss_add(iz)) &
            & - fkw(iz)*kw*(ss_add(iz)-1d0)*prox(iz)**(ss_add(iz)-2d0) &
            & + ss_add(iz)*back(iz)*prox(iz)**(ss_add(iz)-1d0) &
            & - ss_add(iz)*back(iz)*prox(iz)**(ss_add(iz)-1d0)
        d2f1(iz) = d2f1(iz) + (ss_add(iz)+1d0)*(ss_add(iz))*prox(iz)**(ss_add(iz)-1d0) &
            & - fkw(iz)*kw*(ss_add(iz)-1d0)*(ss_add(iz)-2d0)*prox(iz)**(ss_add(iz)-3d0) &
            & + ss_add(iz)*(ss_add(iz)-1d0)*back(iz)*prox(iz)**(ss_add(iz)-2d0) &
            & - ss_add(iz)*(ss_add(iz)-1d0)*back(iz)*prox(iz)**(ss_add(iz)-2d0)

        ! adding charges coming from aq species in eq with gases
        ! pCO2
        fkeq(iz) = 1d0/gamma(2,iz)
        dfkeq_dios(iz) = -1d0/gamma(2,iz)**2d0*dgamma_dios(2,iz)

        f1(iz) = f1(iz)  -  fkw(iz)*k1*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-1d0) &
            & -  2d0*fkeq(iz)*k2*fkw(iz)*k1*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-2d0)
        df1(iz) = df1(iz)  &
            & -  fkw(iz)*k1*kco2*pco2x(iz)*(ss_add(iz)-1d0)*prox(iz)**(ss_add(iz)-2d0) &
            & -  2d0*fkeq(iz)*k2*fkw(iz)*k1*kco2*pco2x(iz)*(ss_add(iz)-2d0)*prox(iz)**(ss_add(iz)-3d0)
        d2f1(iz) = d2f1(iz)  &
            & -  fkw(iz)*k1*kco2*pco2x(iz)*(ss_add(iz)-1d0)*(ss_add(iz)-2d0)*prox(iz)**(ss_add(iz)-3d0) &
            & -  2d0*fkeq(iz)*k2*fkw(iz)*k1*kco2*pco2x(iz)*(ss_add(iz)-2d0)*(ss_add(iz)-3d0)*prox(iz)**(ss_add(iz)-4d0)
        df1dmgas(ipco2,iz) = df1dmgas(ipco2,iz)  -  fkw(iz)*k1*kco2*1d0*prox(iz)**(ss_add(iz)-1d0) &
            & -  2d0*fkeq(iz)*k2*fkw(iz)*k1*kco2*1d0*prox(iz)**(ss_add(iz)-2d0)
        ! pNH3
        f1(iz) = f1(iz)  +  pnh3x(iz)*knh3/k1nh3*prox(iz)**(ss_add(iz)+1d0)
        df1(iz) = df1(iz)  +  pnh3x(iz)*knh3/k1nh3*(ss_add(iz)+1d0)*prox(iz)**(ss_add(iz))
        d2f1(iz) = d2f1(iz)  +  pnh3x(iz)*knh3/k1nh3*(ss_add(iz)+1d0)*(ss_add(iz))*prox(iz)**(ss_add(iz)-1d0)
        df1dmgas(ipnh3,iz) = df1dmgas(ipnh3,iz)  +  1d0*knh3/k1nh3*prox(iz)**(ss_add(iz)+1d0)

        if (print_res) write(*,fmt='(a,L)', advance='no') 'after gas',f1(iz)>0d0

        do ispa = 1, nsp_aq_all
            
            f1(iz) = f1(iz) + base_charge(ispa)*maqf_loc(ispa,iz)*prox(iz)**(ss_add(iz))
            df1(iz) = df1(iz) + base_charge(ispa)*maqf_loc(ispa,iz)*(ss_add(iz))*prox(iz)**(ss_add(iz)-1d0)
            d2f1(iz) = d2f1(iz) + base_charge(ispa)*maqf_loc(ispa,iz)*(ss_add(iz))*(ss_add(iz)-1d0)*prox(iz)**(ss_add(iz)-2d0)
            df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) + base_charge(ispa)*1d0*prox(iz)**(ss_add(iz))
            
            if (print_res .and. f1(iz)<0d0) then 
                write(*,fmt='(1x,a,1x,a)', advance='no') 'base-aq',chraq_all(ispa)
            endif 
            
            ! account for speces associated with NH4+ (both anions and cations)
            do ispa_nh3 = 1,2
                if ( keqaq_nh3(ispa,ispa_nh3) > 0d0) then 
                    rspa_nh3 = real(ispa_nh3,kind=8)
                    ic1 = nint(abs(base_charge(ispa)))
                    ic2 = nint(abs(base_charge(ispa)+rspa_nh3))
                    if ( ic1>0 .and. ic2 > 0) then  
                        fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_nh3/gamma(ic2,iz)
                        dfkeq_dios(iz) = ( &
                            & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_nh3/gamma(ic2,iz) &
                            & + gamma(ic1,iz)*rspa_nh3*gamma(1,iz)**(rspa_nh3-1d0)*dgamma_dios(1,iz) &
                            &   /gamma(ic2,iz) &
                            & + gamma(ic1,iz)*gamma(1,iz)**rspa_nh3*(-1d0) &
                            &   /gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                            & )
                    elseif ( ic1==0 .and. ic2 > 0) then  
                        fkeq(iz) = gamma(1,iz)**rspa_nh3/gamma(ic2,iz)
                        dfkeq_dios(iz) = ( &
                            & + rspa_nh3*gamma(1,iz)**(rspa_nh3-1d0)*dgamma_dios(1,iz) &
                            &   /gamma(ic2,iz) &
                            & + gamma(1,iz)**rspa_nh3*(-1d0) &
                            &   /gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                            & )
                    elseif ( ic1>0 .and. ic2 == 0) then  
                        fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_nh3
                        dfkeq_dios(iz) = ( &
                            & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_nh3 &
                            & + gamma(ic1,iz)*rspa_nh3*gamma(1,iz)**(rspa_nh3-1d0)*dgamma_dios(1,iz) &
                            & )
                    elseif ( ic1==0 .and. ic2 == 0) then  
                        fkeq(iz) = gamma(1,iz)**rspa_nh3
                        dfkeq_dios(iz) = ( &
                            & + rspa_nh3*gamma(1,iz)**(rspa_nh3-1d0)*dgamma_dios(1,iz) &
                            & )
                    else    
                        print *, 'something is wrong'
                        stop
                    endif 
                    f1(iz) = f1(iz) &
                        & + (base_charge(ispa) + rspa_nh3)*fkeq(iz)*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                        & * (pnh3x(iz)*knh3/k1nh3)**rspa_nh3*prox(iz)**(rspa_nh3+ss_add(iz))
                    df1(iz) = df1(iz) &
                        & + (base_charge(ispa) + rspa_nh3)*fkeq(iz)*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                        & * (pnh3x(iz)*knh3/k1nh3)**rspa_nh3*(rspa_nh3+ss_add(iz))*prox(iz)**(rspa_nh3+ss_add(iz)-1d0)
                    d2f1(iz) = d2f1(iz) &
                        & + (base_charge(ispa) + rspa_nh3)*fkeq(iz)*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                        & * (pnh3x(iz)*knh3/k1nh3)**rspa_nh3*(rspa_nh3+ss_add(iz))*(rspa_nh3+ss_add(iz)-1d0) &
                        & * prox(iz)**(rspa_nh3+ss_add(iz)-2d0)
                    df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) &
                        & + (base_charge(ispa) + rspa_nh3)*fkeq(iz)*keqaq_nh3(ispa,ispa_nh3)*1d0 &
                        & * (pnh3x(iz)*knh3/k1nh3)**rspa_nh3*prox(iz)**(rspa_nh3+ss_add(iz))
                    df1dmgas(ipnh3,iz) = df1dmgas(ipnh3,iz) &
                        & + (base_charge(ispa) + rspa_nh3)*fkeq(iz)*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,iz) &
                        & * (knh3/k1nh3)**rspa_nh3*prox(iz)**(rspa_nh3+ss_add(iz)) &
                        & * rspa_nh3*pnh3x(iz)**(rspa_nh3-1d0)
                endif 
            enddo 
            
            if (print_res .and. f1(iz)<0d0) then 
                write(*,fmt='(1x,a,1x,a)', advance='no') 'nh4-aq',chraq_all(ispa)
            endif 
            
            ! annions
            if ( &
                & trim(adjustl(chraq_all(ispa)))=='no3' &
                & .or. trim(adjustl(chraq_all(ispa)))=='so4' &
                & .or. trim(adjustl(chraq_all(ispa)))=='cl' &
                & .or. trim(adjustl(chraq_all(ispa)))=='ac' &
                & .or. trim(adjustl(chraq_all(ispa)))=='mes' &
                & .or. trim(adjustl(chraq_all(ispa)))=='im' &
                & .or. trim(adjustl(chraq_all(ispa)))=='tea' &
                ! & .or. trim(adjustl(chraq_all(ispa)))=='oxa' &
                & ) then 
                
                ! account for speces associated with H+
                do ispa_h = 1,2
                    if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                        rspa_h = real(ispa_h,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)+rspa_h))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_h/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_h/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*rspa_h*gamma(1,iz)**(rspa_h-1d0)*dgamma_dios(1,iz) &
                                &   /gamma(ic2,iz) &
                                & + gamma(ic1,iz)*gamma(1,iz)**rspa_h*(-1d0) &
                                &   /gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(1,iz)**rspa_h/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + rspa_h*gamma(1,iz)**(rspa_h-1d0)*dgamma_dios(1,iz) &
                                &   /gamma(ic2,iz) &
                                & + gamma(1,iz)**rspa_h*(-1d0) &
                                &   /gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_h
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_h &
                                & + gamma(ic1,iz)*rspa_h*gamma(1,iz)**(rspa_h-1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(1,iz)**rspa_h
                            dfkeq_dios(iz) = ( &
                                & + rspa_h*gamma(1,iz)**(rspa_h-1d0)*dgamma_dios(1,iz) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        f1(iz) = f1(iz) &
                            & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**(rspa_h+ss_add(iz))
                        df1(iz) = df1(iz) &
                            & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(rspa_h+ss_add(iz)) &
                            & * prox(iz)**(rspa_h+ss_add(iz)-1d0)
                        d2f1(iz) = d2f1(iz) &
                            & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(rspa_h+ss_add(iz)) &
                            & * (rspa_h+ss_add(iz)-1d0)*prox(iz)**(rspa_h+ss_add(iz)-2d0)
                        df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) &
                            & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*1d0*prox(iz)**(rspa_h+ss_add(iz))
                    endif 
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'h-aq',chraq_all(ispa)
                endif 
            ! oxalic acid 
            elseif ( &
                & trim(adjustl(chraq_all(ispa)))=='oxa' &
                & .or. trim(adjustl(chraq_all(ispa)))=='glp' &
                & ) then 
                do ispa_h = 1,2
                    if (ispa_h==1) then 
                        if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                            rspa_h = real(ispa_h,kind=8)
                            fkeq(iz) = 1d0/gamma(2,iz)
                            dfkeq_dios(iz) = -1d0/gamma(2,iz)**2d0*dgamma_dios(2,iz)
                            f1(iz) = f1(iz) &
                                & + (base_charge(ispa) - rspa_h)*fkeq(iz) &
                                &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**(ss_add(iz)-rspa_h)
                            df1(iz) = df1(iz) &
                                & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(ss_add(iz)-rspa_h) &
                                & * prox(iz)**(ss_add(iz)-rspa_h-1d0)
                            d2f1(iz) = d2f1(iz) &
                                & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(ss_add(iz)-rspa_h) &
                                & * (ss_add(iz)-rspa_h-1d0)* prox(iz)**(ss_add(iz)-rspa_h-2d0)
                            df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) &
                                & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*1d0*prox(iz)**(ss_add(iz)-rspa_h)
                        endif 
                    elseif(ispa_h==2)then
                        if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                            rspa_h = real(ispa_h-1,kind=8)
                            fkeq(iz) = gamma(1,iz)**2d0
                            dfkeq_dios(iz) = 2d0*gamma(1,iz)*dgamma_dios(1,iz)
                            f1(iz) = f1(iz) &
                                & + (base_charge(ispa) + rspa_h)*fkeq(iz) &
                                &       *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**(rspa_h+ss_add(iz))
                            df1(iz) = df1(iz) &
                                & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(rspa_h+ss_add(iz)) &
                                & * prox(iz)**(rspa_h+ss_add(iz)-1d0)
                            d2f1(iz) = d2f1(iz) &
                                & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(rspa_h+ss_add(iz)) &
                                & * (rspa_h+ss_add(iz)-1d0)*prox(iz)**(rspa_h+ss_add(iz)-2d0)
                            df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) &
                                & + (base_charge(ispa) + rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*1d0*prox(iz)**(rspa_h+ss_add(iz))
                        endif 
                    endif
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'h-aq',chraq_all(ispa)
                endif 
            ! cations
            else 
                ! account for hydrolysis speces
                do ispa_h = 1,4
                    if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                        rspa_h = real(ispa_h,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-rspa_h))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(ic1,iz)/gamma(ic2,iz)/gamma(1,iz)**rspa_h
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)/gamma(ic2,iz)/gamma(1,iz)**rspa_h &
                                & + gamma(ic1,iz)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz)/gamma(1,iz)**rspa_h &
                                & + gamma(ic1,iz)/gamma(ic2,iz)*(-rspa_h)/gamma(1,iz)**(rspa_h+1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq(iz) = 1d0/gamma(ic2,iz)/gamma(1,iz)**rspa_h
                            dfkeq_dios(iz) = ( &
                                & + 1d0*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz)/gamma(1,iz)**rspa_h &
                                & + 1d0/gamma(ic2,iz)*(-rspa_h)/gamma(1,iz)**(rspa_h+1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(ic1,iz)/gamma(1,iz)**rspa_h
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)/gamma(1,iz)**rspa_h &
                                & + gamma(ic1,iz)*(-rspa_h)/gamma(1,iz)**(rspa_h+1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq(iz) = 1d0/gamma(1,iz)**rspa_h
                            dfkeq_dios(iz) = ( &
                                & + 1d0*(-rspa_h)/gamma(1,iz)**(rspa_h+1d0)*dgamma_dios(1,iz) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        f1(iz) = f1(iz) &
                            & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)**(ss_add(iz)-rspa_h)
                        df1(iz) = df1(iz) &
                            & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(ss_add(iz)-rspa_h) &
                            & * prox(iz)**(ss_add(iz)-rspa_h-1d0)
                        d2f1(iz) = d2f1(iz) &
                            & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*(ss_add(iz)-rspa_h) &
                            & * (ss_add(iz)-rspa_h-1d0)* prox(iz)**(ss_add(iz)-rspa_h-2d0)
                        df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) &
                            & + (base_charge(ispa) - rspa_h)*fkeq(iz)*keqaq_h(ispa,ispa_h)*1d0*prox(iz)**(ss_add(iz)-rspa_h)
                    endif 
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'h-aq',chraq_all(ispa)
                endif 
                ! account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                do ispa_c = 1,2
                    if ( keqaq_c(ispa,ispa_c) > 0d0) then 
                        if (ispa_c == 1) then ! with CO3--
                            ic1 = nint(abs(base_charge(ispa)))
                            ic2 = nint(abs(base_charge(ispa)-2d0))
                            if ( ic1>0 .and. ic2 > 0) then  
                                fkeq(iz) = gamma(ic1,iz)*gamma(2,iz)/gamma(ic2,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(ic1,iz)*gamma(2,iz)/gamma(ic2,iz) &
                                    & + gamma(ic1,iz)*dgamma_dios(2,iz)/gamma(ic2,iz) &
                                    & + gamma(ic1,iz)*gamma(2,iz)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                    & )
                            elseif ( ic1==0 .and. ic2 > 0) then  
                                fkeq(iz) = gamma(2,iz)/gamma(ic2,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(2,iz)/gamma(ic2,iz) &
                                    & + gamma(2,iz)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                    & )
                            elseif ( ic1>0 .and. ic2 == 0) then  
                                fkeq(iz) = gamma(ic1,iz)*gamma(2,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(ic1,iz)*gamma(2,iz) &
                                    & + gamma(ic1,iz)*dgamma_dios(2,iz) &
                                    & )
                            elseif ( ic1==0 .and. ic2 == 0) then  
                                fkeq(iz) = gamma(2,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(2,iz) &
                                    & )
                            else    
                                print *, 'something is wrong'
                                stop
                            endif 
                            f1(iz) = f1(iz) + (base_charge(ispa)-2d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-2d0)
                            df1(iz) = df1(iz) + (base_charge(ispa)-2d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*(ss_add(iz)-2d0) &
                                & *prox(iz)**(ss_add(iz)-3d0)
                            d2f1(iz) = d2f1(iz) + (base_charge(ispa)-2d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*(ss_add(iz)-2d0) &
                                & *(ss_add(iz)-3d0)*prox(iz)**(ss_add(iz)-4d0)
                            df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) + (base_charge(ispa)-2d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-2d0)
                            df1dmgas(ipco2,iz) = df1dmgas(ipco2,iz) + (base_charge(ispa)-2d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*1d0*prox(iz)**(ss_add(iz)-2d0)
                        elseif (ispa_c == 2) then ! with HCO3-
                            ic1 = nint(abs(base_charge(ispa)))
                            ic2 = nint(abs(base_charge(ispa)-1d0))
                            if ( ic1>0 .and. ic2 > 0) then  
                                fkeq(iz) = gamma(ic1,iz)*gamma(2,iz)*gamma(1,iz)/gamma(ic2,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(ic1,iz)*gamma(2,iz)*gamma(1,iz)/gamma(ic2,iz) &
                                    & + gamma(ic1,iz)*dgamma_dios(2,iz)*gamma(1,iz)/gamma(ic2,iz) &
                                    & + gamma(ic1,iz)*gamma(2,iz)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                    & + gamma(ic1,iz)*gamma(2,iz)*gamma(1,iz)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                    & )
                            elseif ( ic1==0 .and. ic2 > 0) then  
                                fkeq(iz) = gamma(2,iz)*gamma(1,iz)/gamma(ic2,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(2,iz)*gamma(1,iz)/gamma(ic2,iz) &
                                    & + gamma(2,iz)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                    & + gamma(2,iz)*gamma(1,iz)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                    & )
                            elseif ( ic1>0 .and. ic2 == 0) then  
                                fkeq(iz) = gamma(ic1,iz)*gamma(2,iz)*gamma(1,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(ic1,iz)*gamma(2,iz)*gamma(1,iz) &
                                    & + gamma(ic1,iz)*dgamma_dios(2,iz)*gamma(1,iz) &
                                    & + gamma(ic1,iz)*gamma(2,iz)*dgamma_dios(1,iz) &
                                    & )
                            elseif ( ic1==0 .and. ic2 == 0) then  
                                fkeq(iz) = gamma(2,iz)*gamma(1,iz)
                                dfkeq_dios(iz) = ( &
                                    & + dgamma_dios(2,iz)*gamma(1,iz) &
                                    & + gamma(2,iz)*dgamma_dios(1,iz) &
                                    & )
                            else    
                                print *, 'something is wrong'
                                stop
                            endif 
                            f1(iz) = f1(iz) + (base_charge(ispa)-1d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-1d0)
                            df1(iz) = df1(iz) + (base_charge(ispa)-1d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*(ss_add(iz)-1d0) &
                                & *prox(iz)**(ss_add(iz)-2d0)
                            d2f1(iz) = d2f1(iz) + (base_charge(ispa)-1d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*(ss_add(iz)-1d0) &
                                & *(ss_add(iz)-2d0)*prox(iz)**(ss_add(iz)-3d0)
                            df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz) + (base_charge(ispa)-1d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x(iz)*prox(iz)**(ss_add(iz)-1d0)
                            df1dmgas(ipco2,iz) = df1dmgas(ipco2,iz) + (base_charge(ispa)-1d0) &
                                & *fkeq(iz)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*1d0*prox(iz)**(ss_add(iz)-1d0)
                        endif 
                    endif 
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'co2-aq',chraq_all(ispa)
                endif 
                ! account for complexation with free SO4
                do ispa_s = 1,2
                    if ( keqaq_s(ispa,ispa_s) > 0d0) then 
                        rspa_s = real(ispa_s,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-2d0*rspa_s))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(2,iz)**rspa_s/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(2,iz)**rspa_s/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*rspa_s*gamma(2,iz)**(rspa_s-1d0)*dgamma_dios(2,iz)/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*gamma(2,iz)**rspa_s*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(2,iz)**rspa_s/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + rspa_s*gamma(2,iz)**(rspa_s-1d0)*dgamma_dios(2,iz)/gamma(ic2,iz) &
                                & + gamma(2,iz)**rspa_s*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(2,iz)**rspa_s
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(2,iz)**rspa_s &
                                & + gamma(ic1,iz)*rspa_s*gamma(2,iz)**(rspa_s-1d0)*dgamma_dios(2,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(2,iz)**rspa_s
                            dfkeq_dios(iz) = ( &
                                & + rspa_s*gamma(2,iz)**(rspa_s-1d0)*dgamma_dios(2,iz) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        f1(iz) = f1(iz)  + (base_charge(ispa)-2d0*rspa_s) &
                            & *fkeq(iz)*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)**rspa_s*prox(iz)**ss_add(iz)
                        df1(iz) = df1(iz)  + (base_charge(ispa)-2d0*rspa_s) &
                            & *fkeq(iz)*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)**rspa_s*ss_add(iz)*prox(iz)**(ss_add(iz)-1d0)
                        d2f1(iz) = d2f1(iz)  + (base_charge(ispa)-2d0*rspa_s) &
                            & *fkeq(iz)*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)**rspa_s*ss_add(iz) &
                            & *(ss_add(iz)-1d0)*prox(iz)**(ss_add(iz)-2d0)
                        df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz)  + (base_charge(ispa)-2d0*rspa_s) &
                            & *fkeq(iz)*keqaq_s(ispa,ispa_s)*1d0*so4f(iz)**rspa_s*prox(iz)**ss_add(iz)
                        df1dmaqf(iso4,iz) = df1dmaqf(iso4,iz)  + (base_charge(ispa)-2d0*rspa_s) &
                            & *fkeq(iz)*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*rspa_s*so4f(iz)**(rspa_s-1d0)*prox(iz)**ss_add(iz)
                    endif 
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'so4-aq',chraq_all(ispa)
                endif 
                ! account for complexation with free NO3
                do ispa_no3 = 1,2
                    if ( keqaq_no3(ispa,ispa_no3) > 0d0) then 
                        rspa_no3 = real(ispa_no3,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-1d0*rspa_no3))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_no3/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_no3/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*rspa_no3*gamma(1,iz)**(rspa_no3-1d0)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*gamma(1,iz)**rspa_no3*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(1,iz)**rspa_no3/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + rspa_no3*gamma(1,iz)**(rspa_no3-1d0)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                & + gamma(1,iz)**rspa_no3*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_no3
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_no3 &
                                & + gamma(ic1,iz)*rspa_no3*gamma(1,iz)**(rspa_no3-1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(1,iz)**rspa_no3
                            dfkeq_dios(iz) = ( &
                                & + rspa_no3*gamma(1,iz)**(rspa_no3-1d0)*dgamma_dios(1,iz) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        f1(iz) = f1(iz)  + (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                            & *fkeq(iz)*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*no3f(iz)**rspa_no3*prox(iz)**ss_add(iz)
                        df1(iz) = df1(iz)  + (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                            & *fkeq(iz)*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*no3f(iz)**rspa_no3*ss_add(iz)*prox(iz)**(ss_add(iz)-1d0)
                        d2f1(iz) = d2f1(iz)  + (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                            & *fkeq(iz)*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*no3f(iz)**rspa_no3*ss_add(iz) &
                            & *(ss_add(iz)-1d0)*prox(iz)**(ss_add(iz)-2d0)
                        df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz)  + (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                            & *fkeq(iz)*keqaq_no3(ispa,ispa_no3)*1d0*no3f(iz)**rspa_no3*prox(iz)**ss_add(iz)
                        df1dmaqf(ino3,iz) = df1dmaqf(ino3,iz)  + (base_charge(ispa)+base_charge(ino3)*rspa_no3) &
                            & *fkeq(iz)*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,iz)*rspa_no3*no3f(iz)**(rspa_no3-1d0)*prox(iz)**ss_add(iz)
                    endif 
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'no3-aq',chraq_all(ispa)
                endif 
                ! account for complexation with free Cl
                do ispa_cl = 1,2
                    if ( keqaq_cl(ispa,ispa_cl) > 0d0) then 
                        rspa_cl = real(ispa_cl,kind=8)
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-1d0*rspa_cl))
                        if ( ic1>0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_cl/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_cl/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*rspa_cl*gamma(1,iz)**(rspa_cl-1d0)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*gamma(1,iz)**rspa_cl*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(1,iz)**rspa_cl/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + rspa_cl*gamma(1,iz)**(rspa_cl-1d0)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                & + gamma(1,iz)**rspa_cl*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_cl
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**rspa_cl &
                                & + gamma(ic1,iz)*rspa_cl*gamma(1,iz)**(rspa_cl-1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(1,iz)**rspa_cl
                            dfkeq_dios(iz) = ( &
                                & + rspa_cl*gamma(1,iz)**(rspa_cl-1d0)*dgamma_dios(1,iz) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        f1(iz) = f1(iz)  + (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                            & *fkeq(iz)*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*clf(iz)**rspa_cl*prox(iz)**ss_add(iz)
                        df1(iz) = df1(iz)  + (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                            & *fkeq(iz)*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*clf(iz)**rspa_cl*ss_add(iz)*prox(iz)**(ss_add(iz)-1d0)
                        d2f1(iz) = d2f1(iz)  + (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                            & *fkeq(iz)*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*clf(iz)**rspa_cl*ss_add(iz) &
                            & *(ss_add(iz)-1d0)*prox(iz)**(ss_add(iz)-2d0)
                        df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz)  + (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                            & *fkeq(iz)*keqaq_cl(ispa,ispa_cl)*1d0*clf(iz)**rspa_cl*prox(iz)**ss_add(iz)
                        df1dmaqf(icl,iz) = df1dmaqf(icl,iz)  + (base_charge(ispa)+base_charge(icl)*rspa_cl) &
                            & *fkeq(iz)*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,iz)*rspa_cl*clf(iz)**(rspa_cl-1d0)*prox(iz)**ss_add(iz)
                    endif 
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'cl-aq',chraq_all(ispa)
                endif 
                ! account for complexation with HOxa-
                do ispa_oxa = 1,2
                    rspa_oxa   = real(ispa_oxa,kind=8)
                    rspa_oxa_2 = real(ispa_oxa,kind=8) * 2d0
                    rspa_oxa_3 = real(ispa_oxa,kind=8)
                    if (trim(adjustl(chraq_all(ispa)))=='al') then
                        rspa_oxa   = real(ispa_oxa,kind=8) + 1d0
                        rspa_oxa_2 = real(ispa_oxa,kind=8) + 2d0
                        rspa_oxa_3 = 1d0
                    endif 
                        
                    if ( keqaq_oxa(ispa,ispa_oxa) > 0d0) then 
                        ic1 = nint(abs(base_charge(ispa)))
                        ic2 = nint(abs(base_charge(ispa)-rspa_oxa_2))
                        if ( ic1>0 .and. ic2 > 0) then  
                            ! fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**rspa_oxa_3/gamma(ic2,iz)/gamma(1,iz)**rspa_oxa
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa)/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa)/gamma(ic2,iz) &
                                & + gamma(ic1,iz)*(rspa_oxa_3-rspa_oxa)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,iz) &
                                &       /gamma(ic2,iz) &
                                & + gamma(ic1,iz)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 > 0) then  
                            fkeq(iz) = gamma(1,iz)**(rspa_oxa_3-rspa_oxa)/gamma(ic2,iz)
                            dfkeq_dios(iz) = ( &
                                & + (rspa_oxa_3-rspa_oxa)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,iz)/gamma(ic2,iz) &
                                & + gamma(1,iz)**(rspa_oxa_3-rspa_oxa)*(-1d0)/gamma(ic2,iz)**2d0*dgamma_dios(ic2,iz) &
                                & )
                        elseif ( ic1>0 .and. ic2 == 0) then  
                            fkeq(iz) = gamma(ic1,iz)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa)
                            dfkeq_dios(iz) = ( &
                                & + dgamma_dios(ic1,iz)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa) &
                                & + gamma(ic1,iz)*(rspa_oxa_3-rspa_oxa)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,iz) &
                                & )
                        elseif ( ic1==0 .and. ic2 == 0) then
                            fkeq(iz) = gamma(1,iz)**(rspa_oxa_3-rspa_oxa)
                            dfkeq_dios(iz) = ( &
                                & + (rspa_oxa_3-rspa_oxa)*gamma(1,iz)**(rspa_oxa_3-rspa_oxa-1d0)*dgamma_dios(1,iz) &
                                & )
                        else    
                            print *, 'something is wrong'
                            stop
                        endif 
                        f1(iz) = f1(iz)  + (base_charge(ispa)-rspa_oxa_2) &
                            & *fkeq(iz)*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,iz)*oxaf(iz)**rspa_oxa_3*prox(iz)**(ss_add(iz)-rspa_oxa)
                        df1(iz) = df1(iz)  + (base_charge(ispa)-rspa_oxa_2)*fkeq(iz)*keqaq_oxa(ispa,ispa_oxa) &
                            & *maqf_loc(ispa,iz)*oxaf(iz)**rspa_oxa_3*(ss_add(iz)-rspa_oxa)*prox(iz)**(ss_add(iz)-rspa_oxa-1d0)
                        d2f1(iz) = d2f1(iz)  + (base_charge(ispa)-rspa_oxa_2) &
                            & *fkeq(iz)*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,iz)*oxaf(iz)**rspa_oxa_3*(ss_add(iz)-rspa_oxa) &
                            & *(ss_add(iz)-rspa_oxa-1d0)*prox(iz)**(ss_add(iz)-rspa_oxa-2d0)
                        df1dmaqf(ispa,iz) = df1dmaqf(ispa,iz)  + (base_charge(ispa)-rspa_oxa_2) &
                            & *fkeq(iz)*keqaq_oxa(ispa,ispa_oxa)*1d0*oxaf(iz)**rspa_oxa_3*prox(iz)**(ss_add(iz)-rspa_oxa)
                        df1dmaqf(ioxa,iz) = df1dmaqf(ioxa,iz)  + (base_charge(ispa)-rspa_oxa_2)*fkeq(iz)*keqaq_oxa(ispa,ispa_oxa) &
                            & *maqf_loc(ispa,iz)*rspa_oxa_3*oxaf(iz)**(rspa_oxa_3-1d0)*prox(iz)**(ss_add(iz)-rspa_oxa)
                    endif  
                enddo 
                if (print_res .and. f1(iz)<0d0) then 
                    write(*,fmt='(1x,a,1x,a)', advance='no') 'oxa-aq',chraq_all(ispa)
                endif 
            endif 
            
            ! if (print_res) then 
                ! if (ispa==nsp_aq_all) then
                    ! write(*,fmt='(a,a,L)') 'after aq',chraq_all(ispa),f1(iz)>0d0
                ! else
                    ! write(*,fmt='(a,a,L)', advance='no') 'after aq',chraq_all(ispa),f1(iz)>0d0
                ! endif 
            ! endif 
            if (print_res .and. ispa==nsp_aq_all) then 
                write(*,fmt='(1x,a,1x,a,1x,L)') 'end-aq',chraq_all(ispa),f1(iz)>0d0
            endif 
        enddo


    endsubroutine calc_charge_balance_point

end module scepter_equilibrium 