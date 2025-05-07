program weathering
    use scepter_constants
    use scepter_variables
    use scepter_weathering_main
    use scepter_IO
    use scepter_input
    use scepter_physics
    use scepter_psd
    use scepter_concentration
    implicit none 

    CALL getcwd(cwd)
    WRITE(*,*) TRIM(cwd)

    call getarg(0,path)
    WRITE(*,*) TRIM(path)
    ! call get_command_argument(0, cmd)
    ! WRITE(*,*) TRIM(cmd)
    ! path2 = path(:index(path,'weathering')-2)
    path2 = path(:len(trim(path))-len('scepter')-1)
    WRITE(*,*) TRIM(path2)

    CALL chdir(TRIM(path2))
    CALL getcwd(path)
    WRITE(*,*) TRIM(path)

    call get_variables_num( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! output
        & )

    print *,nsp_sld,nsp_aq,nsp_gas,nrxn_ext

    allocate(chraq(nsp_aq),chrsld(nsp_sld),chrgas(nsp_gas),chrrxn_ext(nrxn_ext))
    allocate(chrsld_kinspc(nsld_kinspc),kin_sld_spc(nsld_kinspc))
        
    call get_variables( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! input
        & ,chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc &! output
        & ) 
        
    print *,chraq
    print *,chrsld 
    print *,chrgas 
    print *,chrrxn_ext 
    print *,chrsld_kinspc 
    print *,kin_sld_spc 
    
    ! pause
    sim_name = 'chkchk'

    call get_bsdvalues( &
        & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
        & ,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau &
        & )
        
    call weathering_main( &
        & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,ttot,plant_rain,rainpowder_2nd  &! input
        & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save &! input
        & ,count_dtunchanged_Max,tc,step_tau &! input 
        & ,nsld_kinspc,chrsld_kinspc,kin_sld_spc &! input
        & )

    contains 

    subroutine makegrid(beta,nz,ztot,dz,z,regular_grid)  !  making grid, after Hoffmann & Chiang, 2000
        implicit none
        integer(kind=4),intent(in) :: nz
        logical,intent(in)::regular_grid
        real(kind=8),intent(in)::beta,ztot
        real(kind=8),intent(out)::dz(nz),z(nz)
        integer(kind=4) iz

        do iz = 1, nz 
            z(iz) = iz*ztot/nz  ! regular grid 
            if (iz==1) then
                dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
            endif
            if (iz/=1) then 
                dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(:iz-1))
            endif
        enddo

        if (regular_grid) then 
            dz = ztot/nz  ! when implementing regular grid
        endif 

        do iz=1,nz  ! depth is defined at the middle of individual layers 
            if (iz==1) z(iz)=dz(iz)*0.5d0  
            if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
        enddo

    endsubroutine makegrid

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    subroutine sld_kin( &
        & nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv_tmp &! input
        & ,nsp_gas_all,chrgas_all,mgas_loc &! input
        & ,nsp_aq_all,chraq_all,maqf_loc &! input
        & ,mineral,dev_sp &! input 
        & ,kin,dkin_dmsp &! output
        & ) 
        implicit none

        integer,intent(in)::nz
        real(kind=8),intent(in)::rg,tc,sec2yr,tempk_0,mv_tmp,kw,kho
        real(kind=8),dimension(nz),intent(in)::prox

        real(kind=8) :: cal2j = 4.184d0 

        character(5),intent(in)::mineral,dev_sp
        real(kind=8),dimension(nz),intent(out)::kin,dkin_dmsp
        real(kind=8) mh,moh,kinn_ref,kinh_ref,kinoh_ref,ean,eah,eaoh,tc_ref

        integer,intent(in)::nsp_gas_all,nsp_aq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        real(kind=8),dimension(nsp_gas_all,nz),intent(in)::mgas_loc
        real(kind=8),dimension(nsp_aq_all,nz),intent(in)::maqf_loc

        real(kind=8),dimension(nz) :: pco2,alf
        real(kind=8) mco2,kinco2_ref,eaco2,q10,kref

        ! real(kind=8) k_arrhenius


        kin = 0d0
        dkin_dmsp = 0d0

        select case(trim(adjustl(mineral)))
            ! case('ka','al2o3') ! corundum dissolution rate is assumed to be the same as kaolinite (cf., Carroll-Webb and Walther, 1988)
            case('ka') ! corundum dissolution rate is assumed to be the same as kaolinite (cf., Carroll-Webb and Walther, 1988)
                mh = 0.777d0
                moh = -0.472d0
                kinn_ref = 10d0**(-13.18d0)*sec2yr
                kinh_ref = 10d0**(-11.31d0)*sec2yr
                kinoh_ref = 10d0**(-17.05d0)*sec2yr
                ean = 22.2d0
                eah = 65.9d0
                eaoh = 17.9d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('ab')
                mh = 0.457d0
                moh = -0.572d0
                kinn_ref = 10d0**(-12.56d0)*sec2yr
                kinh_ref = 10d0**(-10.16d0)*sec2yr
                kinoh_ref = 10d0**(-15.6d0)*sec2yr
                ean = 69.8d0
                eah = 65d0
                eaoh = 71d0
                ! above is from table 13 ( this could be a confused mixture of linear and non-linear regression parameters in Table 1)
                ! following is linear regression result of Table 1
                mh = 0.457d0
                moh = -0.572d0
                kinn_ref = 10d0**(-12.04d0)*sec2yr
                kinh_ref = 10d0**(-9.87d0)*sec2yr
                kinoh_ref = 10d0**(-16.98d0)*sec2yr
                ean = 69.8d0
                eah = 65d0
                eaoh = 71d0
                ! then non-linear regression
                ! mh = 0.317d0
                ! moh = -0.471d0
                ! kinn_ref = 10d0**(-12.56d0)*sec2yr
                ! kinh_ref = 10d0**(-10.16d0)*sec2yr
                ! kinoh_ref = 10d0**(-15.6d0)*sec2yr
                ! ean = 65d0
                ! eah = 65d0
                ! eaoh = 66.5d0
                !!! 
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('kfs','sdn') ! sanidine assumed to follow rate law of K-feldspar
                mh = 0.5d0
                moh = -0.823d0
                kinn_ref = 10d0**(-12.41d0)*sec2yr
                kinh_ref = 10d0**(-10.06d0)*sec2yr
                ! kinoh_ref = 10d0**(-9.68d0)*sec2yr*kw**(-moh)
                kinoh_ref = 10d0**(-21.2d0)*sec2yr
                ! ean = 9.08*cal2j
                ! eah = 12.4d0*cal2j
                ! eaoh = 22.5d0*cal2j
                ean = 38d0
                eah = 51.7d0
                eaoh = 94.1d0
                tc_ref = 25d0
                ! Brantley et al 2008 used data from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('leu') 
                mh = 0.7d0
                moh = -0.2d0
                kinn_ref = 10d0**(-9.20d0)*sec2yr
                kinh_ref = 10d0**(-6.00d0)*sec2yr
                kinoh_ref = 10d0**(-10.66d0)*sec2yr
                ean = 75.5d0
                eah = 132.2d0
                eaoh = 56.6d0
                tc_ref = 25d0
                ! Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('fo')
                mh = 0.47d0
                moh = 0d0
                kinn_ref = 10d0**(-10.64d0)*sec2yr
                kinh_ref = 10d0**(-6.85d0)*sec2yr
                kinoh_ref = 0d0
                ean = 79d0
                eah = 67.2d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('fa')
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-12.80d0)*sec2yr
                kinh_ref = 10d0**(-4.80d0)*sec2yr
                kinoh_ref = 0d0
                ean = 94.4d0
                eah = 94.4d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('an')
                mh = 1.411d0
                moh = 0d0
                kinn_ref = 10d0**(-9.12d0)*sec2yr
                kinh_ref = 10d0**(-3.5d0)*sec2yr
                kinoh_ref = 0d0
                ean = 17.8d0
                eah = 16.6d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('la')
                mh = 0.626d0
                moh = 0d0
                moh = -0.57d0
                kinn_ref = 10d0**(-10.91d0)*sec2yr
                kinh_ref = 10d0**(-7.87d0)*sec2yr
                kinoh_ref = 0d0 ! original data 
                ! kinoh_ref = 10d0**(-15.57d0)*sec2yr ! added by Beerling et al. 2020 from albite data of Palandri and Kharaka 2004
                ean = 45.2d0
                eah = 42.1d0
                eaoh = 0d0 ! original data
                ! eaoh = 71d0 ! Beerling
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('and')
                mh = 0.541d0
                moh = 0d0
                kinn_ref = 10d0**(-11.47d0)*sec2yr
                kinh_ref = 10d0**(-8.88d0)*sec2yr
                kinoh_ref = 0d0
                ean = 57.4d0
                eah = 53.5d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('olg')
                mh = 0.457d0
                moh = 0d0
                kinn_ref = 10d0**(-11.84d0)*sec2yr
                kinh_ref = 10d0**(-9.67d0)*sec2yr
                kinoh_ref = 0d0
                ean = 69.8d0
                eah = 65.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('by')
                mh = 1.018d0
                moh = 0d0
                kinn_ref = 10d0**(-9.82d0)*sec2yr
                kinh_ref = 10d0**(-5.85d0)*sec2yr
                kinoh_ref = 0d0
                ean = 31.5d0
                eah = 29.3d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('cc')
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-5.81d0)*sec2yr
                kinh_ref = 10d0**(-0.3d0)*sec2yr
                kinoh_ref = 0d0
                ean = 23.5d0
                eah = 14.4d0
                eaoh = 0d0
                tc_ref = 25d0
                ! adding co2 mechanism
                mco2 = 1d0
                kinco2_ref = 10d0**(-3.48d0)*sec2yr
                eaco2 = 35.4d0
                pco2 = mgas_loc(findloc(chrgas_all,'pco2',dim=1),:)
                ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & + pco2**mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case('pco2')
                        dkin_dmsp = ( & 
                            & + mco2*pco2**(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('arg')
                ! assumed to be the same as those for cc
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-5.81d0)*sec2yr
                kinh_ref = 10d0**(-0.3d0)*sec2yr
                kinoh_ref = 0d0
                ean = 23.5d0
                eah = 14.4d0
                eaoh = 0d0
                tc_ref = 25d0
                ! adding co2 mechanism
                mco2 = 1d0
                kinco2_ref = 10d0**(-3.48d0)*sec2yr
                eaco2 = 35.4d0
                pco2 = mgas_loc(findloc(chrgas_all,'pco2',dim=1),:)
                ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & + pco2**mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case('pco2')
                        dkin_dmsp = ( & 
                            & + mco2*pco2**(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('dlm') ! for disordered dolomite
                mh = 0.500d0
                moh = 0d0
                kinn_ref = 10d0**(-7.53d0)*sec2yr
                kinh_ref = 10d0**(-3.19d0)*sec2yr
                kinoh_ref = 0d0
                ean = 52.2d0
                eah = 36.1d0
                eaoh = 0d0
                tc_ref = 25d0
                ! adding co2 mechanism
                mco2 = 0.5d0
                kinco2_ref = 10d0**(-5.11d0)*sec2yr
                eaco2 = 34.8d0
                pco2 = mgas_loc(findloc(chrgas_all,'pco2',dim=1),:)
                ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & + pco2**mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case('pco2')
                        dkin_dmsp = ( & 
                            & + mco2*pco2**(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('gb','amal')
                mh = 0.992d0
                moh = -0.784d0
                kinn_ref = 10d0**(-11.50d0)*sec2yr
                kinh_ref = 10d0**(-7.65d0)*sec2yr
                kinoh_ref = 10d0**(-16.65d0)*sec2yr
                ean = 61.2d0
                eah = 47.5d0
                eaoh = 80.1d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('amsi')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-12.23d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 74.5d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('phsi')
                mh = 1d0
                moh = -0.33d0
                kinn_ref = 5d-18*1d4*sec2yr ! given as mol/cm2/sec; converted to mol/m2/yr
                kinh_ref = 6d-16*1d4*sec2yr ! given as mol/cm2/sec; converted to mol/m2/yr
                kinoh_ref = 3.5d-13*1d4*sec2yr *kw**(-moh) ! given as mol/cm2/sec; converted to mol/m2/yr
                ean = 74.5d0 
                eah = 74.5d0 
                eaoh = 74.5d0 
                tc_ref = 25d0
                ! from Fraysse et al. 2009 except for Ea which is taken from value for Palandri and Kharaka 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('qtz')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-13.40d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 90.9d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                mh = 0.309d0 ! brantley 2008
                moh = -0.411d0 ! brantley 2008
                kinn_ref = 10d0**(-13.40d0)*sec2yr
                kinh_ref = 4.34d-12 *sec2yr ! brantley 2008
                kinoh_ref = 6.06d-10 * kw**(-moh) *sec2yr ! brantley 2008
                ean = 90.9d0
                eah = 90.9d0 ! assumed to be the same as ean
                eaoh = 90.9d0 ! assumed to be the same as ean
                tc_ref = 25d0
                ! pH neutral range from Palandri and Kharaka, 2004 pH dependence from Brantley et al 2008
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('gt','amfe3')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-7.94d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 86.5d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('hm')
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-14.60d0)*sec2yr
                kinh_ref = 10d0**(-9.39d0)*sec2yr
                kinoh_ref = 0d0
                ean = 66.2d0
                eah = 66.2d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('ct')
                mh = 0d0
                moh = -0.23d0
                kinn_ref = 10d0**(-12d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 10d0**(-13.58d0)*sec2yr
                ean = 73.5d0
                eah = 0d0
                eaoh = 73.5d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('mscv')
                mh = 0.370d0
                moh = -0.22d0
                kinn_ref = 10d0**(-13.55d0)*sec2yr
                kinh_ref = 10d0**(-11.85d0)*sec2yr
                kinoh_ref = 10d0**(-13.55d0)*sec2yr
                ean = 22d0
                eah = 22d0
                eaoh = 22d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('plgp')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-12.4d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 29d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('cabd','ill','kbd','nabd','mgbd','casp','ksp','nasp','mgsp') ! illite kinetics is assumed to be the same as smectite (Bibi et al., 2011)
                mh = 0.34d0
                moh = -0.4d0
                kinn_ref = 10d0**(-12.78d0)*sec2yr
                kinh_ref = 10d0**(-10.98d0)*sec2yr
                kinoh_ref = 10d0**(-16.52d0)*sec2yr
                ean = 35d0
                eah = 23.6d0
                eaoh = 58.9d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('nph','anl') ! analcime kinetics is assumed to be the same as nepherine (cf. Ragnarsdottir, GCA, 1993)
                mh = 1.130d0
                moh = -0.200d0
                kinn_ref = 10d0**(-8.56d0)*sec2yr
                kinh_ref = 10d0**(-2.73d0)*sec2yr
                kinoh_ref = 10d0**(-10.76d0)*sec2yr
                ean = 65.4d0
                eah = 62.9d0
                eaoh = 37.8d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('dp')
                mh = 0.71d0
                moh = 0d0
                kinn_ref = 10d0**(-11.11d0)*sec2yr
                kinh_ref = 10d0**(-6.36d0)*sec2yr
                kinoh_ref = 0d0
                ean = 50.6d0
                eah = 96.1d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('hb','cpx','agt')
                mh = 0.70d0
                moh = 0d0
                kinn_ref = 10d0**(-11.97d0)*sec2yr
                kinh_ref = 10d0**(-6.82d0)*sec2yr
                kinoh_ref = 0d0
                ean = 78.0d0
                eah = 78.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for augite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('en','opx','fer')
                mh = 0.60d0
                moh = 0d0
                kinn_ref = 10d0**(-12.72d0)*sec2yr
                kinh_ref = 10d0**(-9.02d0)*sec2yr
                kinoh_ref = 0d0
                ean = 80.0d0
                eah = 80.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for enstatite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('jd')
                mh = 0.70d0
                moh = 0d0
                kinn_ref = 10d0**(-9.50d0)*sec2yr
                kinh_ref = 10d0**(-6.00d0)*sec2yr
                kinoh_ref = 0d0
                ean = 94.4d0
                eah = 132.2d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for jadeite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('wls')
                mh = 0.40d0
                moh = 0d0
                kinn_ref = 10d0**(-8.88d0)*sec2yr
                kinh_ref = 10d0**(-5.37d0)*sec2yr
                kinoh_ref = 0d0
                ean = 54.7d0
                eah = 54.7d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for wollastonite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('tm')
                mh = 0.70d0
                moh = 0d0
                kinn_ref = 10d0**(-10.60d0)*sec2yr
                kinh_ref = 10d0**(-8.40d0)*sec2yr
                kinoh_ref = 0d0
                ean = 94.4d0
                eah = 18.9d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('ep')
                mh = 0.338d0
                moh = -0.556d0
                kinn_ref = 10d0**(-11.99d0)*sec2yr
                kinh_ref = 10d0**(-10.60d0)*sec2yr
                kinoh_ref = 10d0**(-17.33d0)*sec2yr
                ean = 70.7d0
                eah = 71.1d0
                eaoh = 79.1d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('clch')
                mh = 0.5d0
                moh = 00
                kinn_ref = 10d0**(-12.52d0)*sec2yr
                kinh_ref = 10d0**(-11.11d0)*sec2yr
                kinoh_ref = 0d0
                ean = 88.0d0
                eah = 88.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('cdr')
                mh = 1.0d0
                moh = 00
                kinn_ref = 10d0**(-11.20d0)*sec2yr
                kinh_ref = 10d0**(-3.8d0)*sec2yr
                kinoh_ref = 0d0
                ean = 28.3d0
                eah = 113.3d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('antp','splt')
                mh = 0.440d0
                moh = 0d0
                kinn_ref = 10d0**(-14.24d0)*sec2yr
                kinh_ref = 10d0**(-11.94d0)*sec2yr
                kinoh_ref = 0d0
                ean = 51.0d0
                eah = 51.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for anthophyllite from Palandri and Kharaka, 2004
                ! sepiolite dissolution is assumed to be close to anthophyllite dissolution (Mulders and Oelkers 2020 GCA)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('gps')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-2.79d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 0d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            ! case('caso4') ! from Palandri and Kharaka, 2004
                ! mh = 0d0
                ! moh = 0d0
                ! kinn_ref = 10d0**(-3.19d0)*sec2yr
                ! kinh_ref = 0d0
                ! kinoh_ref = 0d0
                ! ean = 14.3d0
                ! eah = 0d0
                ! eaoh = 0d0
                ! tc_ref = 25d0
                ! kin = ( & 
                    ! & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    ! & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    ! & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    ! & ) 
                ! select case(trim(adjustl(dev_sp)))
                    ! case('pro')
                        ! dkin_dmsp = ( & 
                            ! & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            ! & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            ! & ) 
                    ! case default 
                        ! dkin_dmsp = 0d0
                ! endselect 
                
            case('fe2o','mgo','k2o','cao','na2o','al2o3','sio2','caso4')
                kin = ( &
                    & 1d0/0.01d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 1 year as in Chen et al. (2010, AFM) 
                    & )
                ! assuming randomly fast dissolution rate for oxides
                dkin_dmsp = 0d0
            
            case('gbas')
                mh = 1.013d0
                moh = -0.258d0
                kinn_ref = 0d0
                kinh_ref = 10d0**(-4.27d0)*sec2yr
                kinoh_ref = 10d0**(-11.00d0)*sec2yr
                ean = 0d0
                eah = 39.7d0
                eaoh = 38.4d0
                tc_ref = 25d0
                ! from Pollyea and Rimstidt, 2017 (this is normalized to geometric surface area)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                ! from Brantley et al. (2008)
                mh = 1.16d0
                moh = -0.16d0
                eah = 47.5d0
                eaoh = 47.5d0
                kinh_ref = 588d0*sec2yr
                kinoh_ref = 0.0822d0*kw**(-moh)*sec2yr
                kin = ( &
                    & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**mh &
                    & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**moh &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                        ! from Brantley et al. (2008)
                        dkin_dmsp = ( &
                            & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**(mh-1d0)*mh &
                            & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**(moh-1d0)*moh &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('cbas')
                mh = 0.68d0
                moh = -0.286d0
                kinn_ref = 0d0
                kinh_ref = 10d0**(-6.15d0)*sec2yr
                kinoh_ref = 10d0**(-11.83d0)*sec2yr
                ean = 0d0
                eah = 40.1d0
                eaoh = 32.9d0
                tc_ref = 25d0
                ! from Pollyea and Rimstidt, 2017 (this is normalized to geometric surface area)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                ! ----- from Brantley et al. (2008)
                mh = 1.16d0
                moh = -0.16d0
                eah = 47.5d0
                eaoh = 47.5d0
                kinh_ref = 588d0*sec2yr
                kinoh_ref = 0.0822d0*kw**(-moh)*sec2yr
                kin = ( &
                    & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**mh &
                    & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**moh &
                    & ) 
                ! ----- from Oelkers and Gislason (2001) & Gislason & Oelkers (2003)
                alf = maqf_loc(findloc(chraq_all,'al',dim=1),:)
                eah = 25.5d0 ! from Gislason and Oelkers (2003) (but for geometric surface based equation)
                tc_ref = 25d0
                mh = 0.35d0 ! OG01
                ! mh = 1d0/3d0 ! GO03
                kinh_ref = 10d0**(-11.64d0)*1d4*sec2yr ! units converted from mol/cm2/sec to mol/m2/yr (OG01)
                ! kinh_ref = 10d0**(-5.6d0)*1d4/92d0*sec2yr ! units converted from mol/cm2/sec to mol/m2/yr (GO03); 92 is surface roughness
                kin = (&
                    & k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg)*(prox**3d0/alf)**mh &
                    ! & kinh_ref*exp(-eah/(rg*(tc+tempk_0)))*(prox**3d0/alf)**mh &
                    & )
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                        ! ----- from Brantley et al. (2008)
                        dkin_dmsp = ( &
                            & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**(mh-1d0)*mh &
                            & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**(moh-1d0)*moh &
                            & ) 
                        ! ----- from Oelkers and Gislason (2001)
                        dkin_dmsp = ( &
                            & k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg)*prox**(3d0*mh-1d0)*3d0*mh *(1d0/alf)**mh &
                            ! & kinh_ref*exp(-eah/(rg*(tc+tempk_0)))*prox**(3d0*mh-1d0)*3d0*mh *(1d0/alf)**mh &
                            & ) 
                    case('al')
                        ! ----- from Oelkers and Gislason (2001)
                        dkin_dmsp = ( &
                            & k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg)*prox**(3d0*mh)*(-mh) *(1d0/alf)**(mh+1d0) &
                            ! & kinh_ref*exp(-eah/(rg*(tc+tempk_0)))*prox**(3d0*mh)*(-mh) *(1d0/alf)**(mh+1d0) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
                
            case('py')
                mh = 0d0
                moh = -0.11d0
                kinn_ref = 0d0
                kinh_ref = 0d0
                kinoh_ref = 10.0d0**(-8.19d0)*sec2yr*kho**0.5d0
                ean = 0d0
                eah = 0d0
                eaoh = 57d0
                tc_ref = 15d0
                ! from Williamson and Rimstidt (1994)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
                
            case('g1')
                kin = ( &
                    & 1d0/1d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 1 year as in Chen et al. (2010, AFM) 
                    & )
                ! adding temperature dependence in the form of Q10
                kref = 1d0/1d0
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                dkin_dmsp = 0d0
                
            case('g2')
                kin = ( &
                    & 1d0/8d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 8 year as in Chen et al. (2010, AFM) 
                    & )
                ! adding temperature dependence in the form of Q10
                kref = 1d0/8d0
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                dkin_dmsp = 0d0
                
            case('g3')
                kin = ( &
                    & 1d0/1d3 &! mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 1k year  
                    & )
                ! adding temperature dependence in the form of Q10
                kref = 1d0/1d3
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                dkin_dmsp = 0d0
                
            case('amnt')
                kin = ( &
                    & 1d0/0.01d0 &! just a value assumed; turnover time of 0.1 year for NH4NO3 
                    & )
                ! adding temperature dependence in the form of Q10
                ! kref = 0.01d0/1d0
                kref = 1d0/0.1d0
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                ! kin = kref
                dkin_dmsp = 0d0
                
            case('kcl','gac','mesmh','ims','teas','naoh','naglp','cacl2','nacl')
                kin = ( &
                    & 1d0/0.01d0 &! just a value assumed; turnover time of 0.1 year for NH4NO3 
                    & )
                ! adding temperature dependence in the form of Q10
                ! kref = 0.01d0/1d0
                kref = 1d0/0.0001d0
                tc_ref = 15d0
                q10 = 1d0
                kin = k_q10(kref,tc,tc_ref,q10)
                ! kin = kref
                dkin_dmsp = 0d0
                
            case default 
                kin =0d0
                dkin_dmsp = 0d0

        endselect  


    endsubroutine sld_kin

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    subroutine sld_therm( &
        & rg,tc,tempk_0,ss_x,ss_y,ss_z &! input
        & ,mineral &! input
        & ,therm &! output
        & ) 
        implicit none

        real(kind=8),intent(in)::rg,tc,tempk_0,ss_x,ss_y,ss_z
        real(kind=8) :: cal2j = 4.184d0 
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

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    subroutine calc_omega_v5( &
        & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst & 
        & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &
        & ,maqx,maqc,mgasx,mgasc,mgasth_all,prox,iosx,tc &
        & ,keqsld_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 &
        & ,staq_all,stgas_all &
        & ,mineral &
        & ,domega_dmaq_all,domega_dmgas_all,domega_dpro_loc,domega_dios_loc &! output
        & ,omega,omega_error &! output
        & )
        ! this subroutine assumes to receive free ions
        implicit none
        integer,intent(in)::nz
        real(kind=8):: k1,k2,kco2,po2th,mo2g1,mo2g2,mo2g3,keq_tmp,ss_x,ss_pro,ss_pco2,mo2_tmp,tc
        real(kind=8),dimension(nz),intent(in):: prox,iosx
        real(kind=8),dimension(nz):: pco2x,po2x
        real(kind=8),dimension(nz),intent(out)::omega
        logical,intent(out)::omega_error
        character(5),intent(in):: mineral

        integer,intent(in)::nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3
        real(kind=8),dimension(nsp_sld_all),intent(in)::keqsld_all
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::staq_all
        real(kind=8),dimension(nsp_sld_all,nsp_gas_all),intent(in)::stgas_all

        real(kind=8),dimension(nz),intent(out)::domega_dpro_loc,domega_dios_loc
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::domega_dmgas_all
        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::domega_dmaq_all

        real(kind=8),dimension(nsp_aq_all,nz)::maqx_loc,maqf_loc
        real(kind=8),dimension(nsp_aq_all,nz)::dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2
        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

        integer ieqaq_co3,ieqaq_hco3
        data ieqaq_co3,ieqaq_hco3/1,2/

        integer ieqaq_so4,ieqaq_so42
        data ieqaq_so4,ieqaq_so42/1,2/

        integer ispa,ipco2,ipo2
        ! real(kind=8)::thon = 1d0
        real(kind=8)::thon = -1d100

        integer icharge
        real(kind=8),dimension(nz)::fkeq,gamma_tmp,dgamma_dios_tmp
        real(kind=8),dimension(4,nz)::gamma,dgamma_dios
        real(kind=8) rcharge

        logical::act_ON = .true.
        ! logical::act_ON = .false.

        mo2g1 = keqsld_all(findloc(chrsld_all,'g1',dim=1))
        mo2g2 = keqsld_all(findloc(chrsld_all,'g2',dim=1))
        mo2g3 = keqsld_all(findloc(chrsld_all,'g3',dim=1))

        po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

        kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
        k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
        k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)

        ! assuming maqx = maqf, what is obtained is maqf_loc instead of maqx_loc
        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maqx,mgasx,maqc,mgasc &
            & ,maqf_loc,mgasx_loc  &! output
            & )


        pco2x = mgasx_loc(findloc(chrgas_all,'pco2',dim=1),:)
        po2x = mgasx_loc(findloc(chrgas_all,'po2',dim=1),:)

        ipco2 = findloc(chrgas_all,'pco2',dim=1)
        ipo2 = findloc(chrgas_all,'po2',dim=1)

        do icharge=1,4
            rcharge = 1d0*icharge
            call calc_gamma_davies(  &
                & nz,iosx,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma(icharge,:)=gamma_tmp(:)
            dgamma_dios(icharge,:)=dgamma_dios_tmp(:)
        enddo

        domega_dmaq_all  =0d0
        domega_dmgas_all =0d0
        domega_dpro_loc  =0d0
        domega_dios_loc  =0d0

        select case(trim(adjustl(mineral)))

            ! case default ! (almino)silicates & oxides
            case ( &
                & 'fo','ab','an','ka','gb','ct','fa','gt','cabd','dp','hb','kfs','amsi','hm','ill','anl','nph' &
                & ,'qtz','tm','la','by','olg','and','cpx','en','fer','opx','mgbd','kbd','nabd','mscv','plgp','antp' &
                & ,'agt','jd','wls','phsi','splt','casp','ksp','nasp','mgsp','fe2o','mgo','k2o','cao','na2o','al2o3' &
                & ,'gbas','cbas','ep','clch','sdn','cdr','leu','amal','amfe3','sio2' &
                & )  ! (almino)silicates & oxides
                keq_tmp = keqsld_all(findloc(chrsld_all,mineral,dim=1))
                omega = 1d0
                ss_pro = 0d0
                fkeq = 1d0
                do ispa = 1,nsp_aq_all
                    if (staq_all(findloc(chrsld_all,mineral,dim=1),ispa) > 0d0) then 

                        omega = omega*maqf_loc(ispa,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                        
                        ! derivatives are first given as d(log omega)/dc 
                        domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:) + ( &
                            & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)/maqf_loc(ispa,:)*1d0 &
                            & )

                        selectcase(trim(adjustl(chraq_all(ispa)))) 
                            case('na','k')
                                ss_pro = ss_pro + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                fkeq  = fkeq * gamma(1,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                ! derivatives are first given as d(log gamma)/dios 
                                domega_dios_loc = domega_dios_loc &
                                    & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)*dgamma_dios(1,:)/gamma(1,:)
                            case('fe2','ca','mg')
                                ss_pro = ss_pro + 2d0*staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                fkeq  = fkeq * gamma(2,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                domega_dios_loc = domega_dios_loc &
                                    & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)*dgamma_dios(2,:)/gamma(2,:)
                            case('fe3','al')
                                ss_pro = ss_pro + 3d0*staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                fkeq  = fkeq * gamma(3,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                domega_dios_loc = domega_dios_loc &
                                    & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)*dgamma_dios(3,:)/gamma(3,:)
                        endselect
                    endif 
                enddo 
                
                if (ss_pro > 0d0) then 
                    omega = omega / prox**ss_pro
                    
                    fkeq  = fkeq / gamma(1,:)**ss_pro

                    ! derivatives are first given as d(log omega)/dc 
                    domega_dpro_loc = domega_dpro_loc - ss_pro/prox 
                    ! derivatives are first given as d(log gamma)/dios 
                    domega_dios_loc = domega_dios_loc - ss_pro*dgamma_dios(1,:)/gamma(1,:)
                endif 
                
                if (keq_tmp > 0d0) then 
                    if (.not.act_ON) omega = omega / keq_tmp
                    if (act_ON)      omega = omega / keq_tmp * fkeq
                endif         
                
                ! derivatives are now d(omega)/dc ( = d(omega)/d(log omega) * d(log omega)/dc = omega * d(log omega)/dc)
                do ispa=1,nsp_aq_all
                    domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:)*omega(:)
                enddo 
                domega_dpro_loc = domega_dpro_loc*omega
                
                if (.not.act_ON) domega_dios_loc = 0d0
                if (act_ON)      domega_dios_loc = domega_dios_loc*omega
                
            case('cc','arg','dlm') ! carbonates
                keq_tmp = keqsld_all(findloc(chrsld_all,mineral,dim=1))
                ss_pco2 = stgas_all(findloc(chrsld_all,mineral,dim=1),findloc(chrgas_all,'pco2',dim=1))
                omega = 1d0
                fkeq = 1d0
                
                do ispa = 1,nsp_aq_all
                    if (staq_all(findloc(chrsld_all,mineral,dim=1),ispa) > 0d0) then 
                        omega = omega*maqf_loc(ispa,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                        
                        ! derivatives are first given as d(log omega)/dc 
                        domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:) + ( &
                            & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)/maqf_loc(ispa,:)*1d0 &
                            & )
                        
                        fkeq  = fkeq * gamma(2,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                        domega_dios_loc = domega_dios_loc &
                            & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)*dgamma_dios(2,:)/gamma(2,:)
                    endif 
                enddo 
                
                if (ss_pco2 > 0d0) then
                    omega = omega*(k1*k2*kco2*pco2x/(prox**2d0))**ss_pco2
                    
                    ! derivatives are first given as d(log omega)/dc 
                    domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:) + ss_pco2/pco2x 
                    domega_dpro_loc = domega_dpro_loc - 2d0*ss_pco2/prox 
                    
                    fkeq  = fkeq * gamma(2,:)**ss_pco2
                    domega_dios_loc = domega_dios_loc + ss_pco2*dgamma_dios(2,:)/gamma(2,:)
                endif 
                
                if (keq_tmp > 0d0) then 
                    if (.not.act_ON) omega = omega / keq_tmp
                    if (act_ON)      omega = omega / keq_tmp * fkeq
                endif     
                
                ! derivatives are now d(omega)/dc ( = d(omega)/d(log omega) * d(log omega)/dc = omega * d(log omega)/dc)
                do ispa=1,nsp_aq_all
                    domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:)*omega(:)
                enddo 
                domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:)*omega(:)
                domega_dpro_loc = domega_dpro_loc*omega
                
                if (.not.act_ON) domega_dios_loc = 0d0
                if (act_ON)      domega_dios_loc = domega_dios_loc*omega
                
            case('gps','nacl','caso4','cacl2','kcl') ! salts (sulfates/chlorides)
            ! CaSO4*2H2O = Ca+2 + SO4-2 + 2H2O   
            ! NaCl = Na+ + Cl-
            ! KCl = K+ + Cl-
            ! CaCl2  =  Ca2+ + 2 Cl-
                keq_tmp = keqsld_all(findloc(chrsld_all,mineral,dim=1))
                omega = 1d0
                fkeq = 1d0
                
                do ispa = 1,nsp_aq_all
                    if (staq_all(findloc(chrsld_all,mineral,dim=1),ispa) > 0d0) then 
                        omega = omega*maqf_loc(ispa,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                        
                        ! derivatives are first given as d(log omega)/dc 
                        domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:) + ( &
                            & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)/maqf_loc(ispa,:)*1d0 &
                            & )
                            
                        selectcase(trim(adjustl(chraq_all(ispa)))) 
                            case('na','k','cl')
                                gamma_tmp = gamma(1,:)
                                fkeq  = fkeq * gamma(1,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                domega_dios_loc = domega_dios_loc &
                                    & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)*dgamma_dios(1,:)/gamma(1,:)
                            case('ca','mg','so4')
                                fkeq  = fkeq * gamma(2,:)**staq_all(findloc(chrsld_all,mineral,dim=1),ispa)
                                domega_dios_loc = domega_dios_loc &
                                    & + staq_all(findloc(chrsld_all,mineral,dim=1),ispa)*dgamma_dios(2,:)/gamma(2,:)
                        endselect
                    endif 
                enddo 
                
                if (keq_tmp > 0d0) then 
                    if (.not.act_ON) omega = omega / keq_tmp
                    if (act_ON)      omega = omega / keq_tmp * fkeq
                endif     
                
                ! derivatives are now d(omega)/dc ( = d(omega)/d(log omega) * d(log omega)/dc = omega * d(log omega)/dc)
                do ispa=1,nsp_aq_all
                    domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:)*omega(:)
                enddo 
                
                if (.not.act_ON) domega_dios_loc = 0d0
                if (act_ON)      domega_dios_loc = domega_dios_loc*omega
                
            !!! other minerals that are assumed not to be controlled by distance from equilibrium i.e. omega
            
            case('py') ! sulfides (assumed to be totally controlled by kinetics)
            ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
            ! i.e., 1.0 - omega_py = po2x**0.5 
                ! omega = 1d0 - po2x**0.5d0
                omega = 1d0 - po2x**0.5d0*merge(0d0,1d0,po2x<po2th*thon)
                domega_dmgas_all(ipo2,:) = - 0.5d0*po2x**(-0.5d0)*merge(0d0,1d0,po2x<po2th*thon)
                
            case('om','omb')
                omega = 1d0 ! these are not used  
                
            case('g1','g2','g3','amnt')
            ! omega is defined so that kg1*poro*hr*mvg1*1d-6*mg1x*(1d0-omega_g1) = kg1*poro*hr*mvg1*1d-6*mg1x*po2x/(po2x+mo2)
            ! i.e., 1.0 - omega_g1 = po2x/(po2x+mo2) 
                if (trim(adjustl(mineral)) == 'g1') mo2_tmp = mo2g1
                if (trim(adjustl(mineral)) == 'g2') mo2_tmp = mo2g2
                if (trim(adjustl(mineral)) == 'g3') mo2_tmp = mo2g3
                omega = 1d0 - po2x/(po2x+mo2_tmp)*merge(0d0,1d0,po2x < po2th*thon)
                domega_dmgas_all(ipo2,:) = ( &
                    & - 1d0/(po2x+mo2_tmp)*merge(0d0,1d0,po2x < po2th*thon) &
                    & - po2x*(-1d0)/(po2x+mo2_tmp)**2d0*merge(0d0,1d0,po2x < po2th*thon) &
                    & )
            
            case('gac','mesmh','ims','teas','naoh','naglp') ! reacting in any case
                omega = 0d0
            
            case('inrt') ! not reacting in any case
                omega = 1d0
                
            case default 
                ! this should not be selected
                omega = 1d0
                print *, '*** CAUTION: mineral (',mineral,') saturation state is not defined --- > pause'
                pause
                
        endselect

        omega_error = .false.
        if (any(isnan(omega))) then 
            print *,'nan in calc_omega_v4',any(isnan(omega)),mineral
            omega_error = .true.
            ! stop
        endif 

    endsubroutine calc_omega_v5

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 
    subroutine calc_rxn_ext_dev_3( &
        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
        & ,rg,tempk_0,tc,iosx &!input
        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
        & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
        & ,rxn_name,sp_name &! input 
        & ,rxn_ext,drxnext_dmsp,rxnext_error &! output
        & )
        implicit none
        integer,intent(in)::nz
        real(kind=8),dimension(nz),intent(in):: poro,sat,prox,iosx
        real(kind=8),dimension(nz),intent(out):: drxnext_dmsp
        real(kind=8),dimension(nz),intent(out):: rxn_ext
        character(5),intent(in)::rxn_name,sp_name
        logical,intent(out)::rxnext_error

        integer,intent(in)::nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst 

        character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst

        real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl
        real(kind=8),dimension(nrxn_ext_all,nz),intent(in)::krxn1_ext_all,krxn2_ext_all

        integer,intent(in)::nsp_sld,nsp_sld_cnst,nsp_sld_all

        character(5),dimension(nsp_sld),intent(in)::chrsld
        character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all

        real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx
        real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc
        real(kind=8),dimension(nsp_sld_all),intent(in)::msldth_all,mv_all

        real(kind=8),intent(in)::rho_grain,kw,rg,tempk_0,tc

        real(kind=8),dimension(nsp_sld,nz),intent(in):: hr

        ! local variables
        real(kind=8):: po2th,fe2th,mwtom,g1th,g2th,g3th,mvpy,fe3th,knh3,k1nh3,ko2,v_tmp,km_tmp1,km_tmp2,km_tmp3  &
            & ,kn2o,k1fe2,k1fe2co3,k1fe2hco3,k1fe2so4,kco2,k1,k2
        real(kind=8),dimension(nz):: po2x,vmax,mo2,fe2x,koxa,vmax2,mom2,komb,beta,omx,ombx &
            & ,mo2g1,mo2g2,mo2g3,kg1,kg2,kg3,g1x,g2x,g3x,pyx,fe3x,koxpy,pnh3x,nh4x,dnh4_dpro,dnh4_dpnh3 &
            & ,no3x,pn2ox,dv_dph_tmp,fe2f,dfe2f_dfe2,dfe2f_dpco2,dfe2f_dpro,dfe2f_dso4f,pco2x,hrpy,oxax &
            & ,vmax_tmp
        real(kind=8),dimension(nsp_aq_all,nz)::maqx_loc,maqft_loc,dmaqft_dpro_loc,dmaqft_dios_loc
        real(kind=8),dimension(nsp_aq_all,nsp_aq_all,nz)::dmaqft_dmaqf_loc
        real(kind=8),dimension(nsp_aq_all,nsp_gas_all,nz)::dmaqft_dmgas_loc
        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc
        real(kind=8),dimension(nsp_sld_all,nz)::msldx_loc

        ! real(kind=8):: thon = 1d0
        real(kind=8):: thon = -1d100

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

        integer ieqaq_co3,ieqaq_hco3
        data ieqaq_co3,ieqaq_hco3/1,2/

        integer ieqaq_so4,ieqaq_so42
        data ieqaq_so4,ieqaq_so42/1,2/

        character(5) sp_tmp
        character(25) scheme

        ! ... need to clean up the following mess at some day ...

        vmax = krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)
        mo2 = krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)

        po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

        koxa = krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) 

        fe2th = maqth_all(findloc(chraq_all,'fe2',dim=1))

        fe3th = maqth_all(findloc(chraq_all,'fe3',dim=1))

        g1th = msldth_all(findloc(chrsld_all,'g1',dim=1))
        g2th = msldth_all(findloc(chrsld_all,'g2',dim=1))
        g3th = msldth_all(findloc(chrsld_all,'g3',dim=1))

        ko2 = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)

        kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
        k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
        k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)

        knh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h0)
        k1nh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h1)

        kn2o = keqgas_h(findloc(chrgas_all,'pn2o',dim=1),ieqgas_h0)

        k1fe2 = keqaq_h(findloc(chraq_all,'fe2',dim=1),ieqaq_h1)
        k1fe2co3 = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_co3)
        k1fe2hco3  = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_hco3)
        k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)

        vmax2 = krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)
        mom2 = krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)

        komb = krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)
        beta = krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)

        koxpy = krxn1_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1),:)

        kg2 = krxn1_ext_all(findloc(chrrxn_ext_all,'g2k',dim=1),:)

        mvpy = mv_all(findloc(chrsld_all,'py',dim=1))

        hrpy = 0d0
        if (any(chrsld=='py')) then 
            hrpy = hr(findloc(chrsld,'py',dim=1),:)
        endif 


        ! getting mgasx_loc & maqx_loc
        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maqx,mgasx,maqc,mgasc &
            & ,maqx_loc,mgasx_loc  &! output
            & )

        ! getting maqft_loc and its derivatives
        call get_maqt_all( &
        ! call get_maqt_all_v2( &
            & nz,nsp_aq_all,nsp_gas_all &
            & ,chraq_all,chrgas_all &
            & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
            & ,mgasx_loc,maqx_loc,prox,iosx,tc &
            & ,dmaqft_dpro_loc,dmaqft_dmaqf_loc,dmaqft_dmgas_loc,dmaqft_dios_loc &! output
            & ,maqft_loc  &! output
            & )

        call get_msldx_all( &
            & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
            & ,chrsld,chrsld_all,chrsld_cnst &
            & ,msldx,msldc &
            & ,msldx_loc  &! output
            & )
            
        po2x    = mgasx_loc(findloc(chrgas_all,'po2',dim=1),:)
        pco2x   = mgasx_loc(findloc(chrgas_all,'pco2',dim=1),:)
        pnh3x   = mgasx_loc(findloc(chrgas_all,'pnh3',dim=1),:)
        pn2ox   = mgasx_loc(findloc(chrgas_all,'pn2o',dim=1),:)

        oxax    = maqx_loc(findloc(chraq_all,'oxa',dim=1),:)    *maqft_loc(findloc(chraq_all,'oxa',dim=1),:)
        fe2x    = maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    *maqft_loc(findloc(chraq_all,'fe2',dim=1),:)
        fe3x    = maqx_loc(findloc(chraq_all,'fe3',dim=1),:)    *maqft_loc(findloc(chraq_all,'fe3',dim=1),:)
        no3x    = maqx_loc(findloc(chraq_all,'no3',dim=1),:)    *maqft_loc(findloc(chraq_all,'no3',dim=1),:)
        fe2f    = maqx_loc(findloc(chraq_all,'fe2',dim=1),:)

        omx     = msldx_loc(findloc(chrsld_all,'om',dim=1),:)
        ombx    = msldx_loc(findloc(chrsld_all,'omb',dim=1),:)
        g1x     = msldx_loc(findloc(chrsld_all,'g1',dim=1),:)
        g2x     = msldx_loc(findloc(chrsld_all,'g2',dim=1),:)
        g3x     = msldx_loc(findloc(chrsld_all,'g3',dim=1),:)
        pyx     = msldx_loc(findloc(chrsld_all,'py',dim=1),:)


        nh4x = pnh3x*knh3*prox/k1nh3
        dnh4_dpro = pnh3x*knh3*1d0/k1nh3
        dnh4_dpnh3 = 1d0*knh3*prox/k1nh3

        select case(trim(adjustl(rxn_name)))

            case('resp')
                rxn_ext = vmax*po2x/(po2x+mo2)
                ! rxn_ext = vmax*merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
                
                select case(trim(adjustl(sp_name)))
                    case('po2')
                        drxnext_dmsp = (&
                            & vmax*1d0/(po2x+mo2) &
                            & +vmax*po2x*(-1d0)/(po2x+mo2)**2d0 &
                            & )
                    case default
                        drxnext_dmsp = 0d0
                endselect 
            
            case('oxao2')
                ! parameterization adopted by Lawrence et al. 2014; Perez-Fodich and Derry, 2019
                ! decay const of 0.5 yr-1
                rxn_ext = ( &
                    & + poro*sat*1d3*0.5d0 &
                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                    & *maqft_loc(findloc(chraq_all,'oxa',dim=1),:) &
                    & )
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & *dmaqft_dpro_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & )
                            case('ios')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & *dmaqft_dios_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & ) 
                            case('oxa')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & *dmaqft_dmaqf_loc( &
                                    &       findloc(chraq_all,'oxa',dim=1) &
                                    &       ,findloc(chraq_all,'oxa',dim=1),:) &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *1d0 &
                                    & *maqft_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & ) 
                            case default
                                if (any( trim(adjustl(sp_name)) == chraq ) ) then
                                    drxnext_dmsp = ( &
                                        & + poro*sat*1d3*0.5d0 &
                                        & * maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                        & * dmaqft_dmaqf_loc( &
                                        &       findloc(chraq_all,'oxa',dim=1)  &
                                        &       ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                        & ) 
                                elseif (any( trim(adjustl(sp_name)) == chrgas ) ) then
                                    drxnext_dmsp = ( &
                                        & + poro*sat*1d3*0.5d0 &
                                        & * maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                        & * dmaqft_dmgas_loc( &
                                        &       findloc(chraq_all,'oxa',dim=1)  &
                                        &       ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:) &
                                        & ) 
                                else
                                    drxnext_dmsp = 0d0
                                endif 
                        endselect 
                
            case('fe2o2')
                ! scheme = 'full' ! reflecting individual rate consts for different Fe2+ species (after Kanzaki and Murakami 2016)
                scheme = 'default' ! as a function of pH and pO2
                
                selectcase(trim(adjustl(scheme)))
                    case('full')
                        rxn_ext = ( &
                            & + poro*sat*1d3*fe2f*( &
                            & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                            & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                            & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                            & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x/prox &
                            & )*po2x &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2*(-1d0)/prox**2d0 &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0 &
                                    & )*po2x &
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x/prox &
                                    & )*1d0 &
                                    & )
                            case('pco2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*1d0/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*1d0/prox &
                                    & )*po2x &
                                    & )
                            case('fe2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*1d0*( &
                                    & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x/prox &
                                    & )*po2x &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect
                        
                    case default
                        rxn_ext = ( &
                            & poro*sat*1d3*fe2x*po2x &
                            & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                            & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & poro*sat*1d3*fe2x*po2x &
                                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)*2.0d0*(kw*(-1d0)/prox**2d0)) &
                                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & poro*sat*1d3*fe2x*1d0 &
                                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                    & )
                            case('fe2')
                                drxnext_dmsp = ( &
                                    & poro*sat*1d3*po2x &
                                    & * ( & 
                                    & + 1d0    & 
                                    & * maqft_loc(findloc(chraq_all,'fe2',dim=1),:)   &
                                    & + maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    & 
                                    & * dmaqft_dmaqf_loc( &
                                    &       findloc(chraq_all,'fe2',dim=1) &
                                    &       ,findloc(chraq_all,'fe2',dim=1),:)   &
                                    & ) & 
                                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                    & )
                            case default
                                ! drxnext_dmsp = 0d0
                                if (any( trim(adjustl(sp_name)) == chraq ) ) then
                                    drxnext_dmsp = ( &
                                        & poro*sat*1d3*po2x &
                                        & * ( & 
                                        & + maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    & 
                                        & * dmaqft_dmaqf_loc( &
                                        &       findloc(chraq_all,'fe2',dim=1) &
                                        &       ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:)   &
                                        & ) & 
                                        & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                        & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                        & )
                                elseif (any( trim(adjustl(sp_name)) == chrgas ) ) then
                                    drxnext_dmsp = ( &
                                        & poro*sat*1d3*po2x &
                                        & * ( & 
                                        & + maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    & 
                                        & * dmaqft_dmgas_loc( &
                                        &       findloc(chraq_all,'fe2',dim=1) &
                                        &       ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:)   &
                                        & ) & 
                                        & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                        & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                        & )
                                else
                                    drxnext_dmsp = 0d0
                                endif 
                        endselect
                endselect 
            
            case('omomb')
                rxn_ext = vmax2 & ! mg C / soil g /yr
                    & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                    & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                    & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
                
                select case(trim(adjustl(sp_name)))
                    case('om')
                        drxnext_dmsp = ( &
                            & vmax2 & ! mg C / soil g /yr
                            & *1d0*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                            & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & + vmax2 & ! mg C / soil g /yr
                            & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                            & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & *(-1d0)/(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3)**2d0 &
                            & *(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & ) 
                    case('omb')
                        drxnext_dmsp = ( &
                            & vmax2 & ! mg C / soil g /yr
                            & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                            & *1d0*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & ) 
                    case default
                        drxnext_dmsp = 0d0
                endselect
            
            case('ombto')
                rxn_ext = komb*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**beta &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
                
                select case(trim(adjustl(sp_name)))
                    case('omb')
                        drxnext_dmsp = ( &
                            & komb*beta*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**(beta-1d0) &
                            & *(1d0-poro)*rho_grain*rho_grain*1d6 &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & ) 
            
                    case default
                        drxnext_dmsp = 0d0
                        
                endselect 
                        
            
            case('pyfe3') 
                rxn_ext = ( &
                    & koxpy*poro*hrpy*mvpy*pyx*fe3x**0.93d0*fe2x**(-0.40d0) &
                    & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                    & /(1d0 - poro) &
                    & )
                
                select case(trim(adjustl(sp_name)))
                    case('py')
                        drxnext_dmsp = (&
                            & koxpy*poro*hrpy*mvpy*1d0*fe3x**0.93d0*fe2x**(-0.40d0) &
                            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                            & )
                    case('fe3')
                        drxnext_dmsp = (&
                            & koxpy*poro*hrpy*mvpy*pyx*(0.93d0)*fe3x**(0.93d0-1d0)*fe2x**(-0.40d0) &
                            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                            & )
                    case('fe2')
                        drxnext_dmsp = (&
                            & koxpy*poro*hrpy*mvpy*pyx*fe3x**0.93d0*(-0.4d0)*fe2x**(-0.40d0-1d0) &
                            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                            & )
                    case default
                        drxnext_dmsp = 0d0
                endselect 
                
            case('amo2o')
                scheme = 'maggi08' ! Maggi et al. (2008) wihtout baterial, pH and water saturation functions
                ! scheme = 'Fennel' ! from biogem_box_geochem.f90 in GENIE model referring to Fennel et al. 2005 with a correction 
                ! scheme = 'FennelOLD' ! from biogem_box_geochem.f90 in GENIE model referring to Fennel et al. 2005 without a correction 
                ! scheme = 'Ozaki' ! from biogem_box_geochem.f90 in GENIE model referring to Ozaki et al. [EPSL ... ?]
                
                select case(trim(adjustl(scheme)))
                    case('maggi08')
                        v_tmp = 9.53d-6*60d0*60d0*24d0*365d0 ! (~300 /yr)
                        ! v_tmp = v_tmp/100d0 ! (~3 /yr; default value produces too much nitrate (pH goes down to ~1)
                        km_tmp1 = 14d-5
                        km_tmp2 = 2.41d-5
                        rxn_ext = ( &
                            & v_tmp &
                            & *nh4x/(nh4x + km_tmp1 ) &
                            & *po2x*ko2/(po2x*ko2 + km_tmp2 ) &
                            & *min(2d0*sat,1d0) &
                            ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                            & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                            & )
                        
                        dv_dph_tmp = 0d0
                        where (3d0 < -log10(prox) .and. -log10(prox) < 7d0)
                            dv_dph_tmp = 0.25d0
                        elsewhere (7d0 < -log10(prox) .and. -log10(prox) < 11d0)
                            dv_dph_tmp = -0.25d0
                        elsewhere (-log10(prox) == 7d0)
                            dv_dph_tmp = 0d0
                        elsewhere (-log10(prox) == 3d0)
                            dv_dph_tmp = 0.125d0
                        elsewhere (-log10(prox) == 11d0)
                            dv_dph_tmp = -0.125d0
                        elsewhere 
                            dv_dph_tmp = 0d0
                        endwhere 
                        
                        ! when using modified version using normal distribution with sigma = 1 
                        dv_dph_tmp = exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &
                            & *-0.5d0*2d0*((-log10(prox)-7d0)/1d0)  &
                            & *(-1d0) &
                            & *1d0/log(10d0)/prox
                        
                        select case(trim(adjustl(sp_name)))
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *nh4x/(nh4x + km_tmp1 ) &
                                    & *( &
                                    & 1d0*ko2/(po2x*ko2 + km_tmp2) &
                                    & + po2x*ko2*(-1d0)/(po2x*ko2 + km_tmp2)**2d0 * ko2 &
                                    & ) &
                                    & *min(2d0*sat,1d0) &
                                    ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pnh3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & dnh4_dpnh3/(nh4x + km_tmp1 ) &
                                    & + nh4x*(-1d0)/(nh4x + km_tmp1 )**2d0 * dnh4_dpnh3 &
                                    & ) &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & *min(2d0*sat,1d0) &
                                    ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & dnh4_dpro/(nh4x + km_tmp1 ) &
                                    & + nh4x*(-1d0)/(nh4x + km_tmp1 )**2d0 * dnh4_dpro &
                                    & ) &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & *min(2d0*sat,1d0) &
                                    ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & + &
                                    & v_tmp &
                                    & *nh4x/(nh4x + km_tmp1 ) &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *dv_dph_tmp &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                    case('Fennel','FennelOLD')
                        if (trim(adjustl(scheme))== 'Fennel') v_tmp = 6.0d0 ! /yr
                        if (trim(adjustl(scheme))== 'FennelOLD') v_tmp = 0.16667d0 ! /yr
                        km_tmp2 = 2.0D-05
                        rxn_ext = ( &
                            & v_tmp &
                            & *nh4x &
                            & *po2x*ko2/(po2x*ko2 + km_tmp2 ) &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *nh4x &
                                    & *( &
                                    & 1d0*ko2/(po2x*ko2 + km_tmp2) &
                                    & + po2x*ko2*(-1d0)/(po2x*ko2 + km_tmp2)**2d0 * ko2 &
                                    & ) &
                                    & )
                            case('pnh3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *dnh4_dpnh3 &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * dnh4_dpro &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                    case('Ozaki')
                        v_tmp = 18250.0d0/1027.649d0 ! /yr
                        rxn_ext = ( &
                            & v_tmp &
                            & *nh4x &
                            & *po2x*ko2 &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *nh4x &
                                    & * 1d0*ko2 &
                                    & )
                            case('pnh3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *dnh4_dpnh3 &
                                    & *po2x*ko2 &
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * dnh4_dpro &
                                    & *po2x*ko2 &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                    endselect
                        
            case('g2n0','g2n21') 
                ! overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
                ! first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O)   
                ! (assuming that oxidation by N2O governs overall denitrification)
                scheme = 'maggi08' ! Maggi et al. (2008) wihtout baterial, pH and water saturation functions; vmax from oxidation by N2O (rate-limiting)
                
                select case(trim(adjustl(scheme)))
                    case('maggi08')
                        v_tmp = 1.23d-7*60d0*60d0*24d0*365d0
                        km_tmp1 = 10d-5 * 1d6 ! mol L-1 converted to mol m-3
                        km_tmp2 = 11.3d-5
                        km_tmp3 = 2.52d-5
                        rxn_ext = ( &
                            & v_tmp &
                            & *g2x/(g2x + km_tmp1 ) &
                            & *no3x/(no3x + km_tmp2 ) &
                            & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                            & *min(2d0*sat,1d0) &
                            & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                            & )
                        
                        ! when using modified version using normal distribution with sigma = 1 
                        dv_dph_tmp = exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &
                            & *-0.5d0*2d0*((-log10(prox)-7d0)/1d0)  &
                            & *(-1d0) &
                            & *1d0/log(10d0)/prox
                            
                        select case(trim(adjustl(sp_name)))
                            case('g2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & 1d0/(g2x + km_tmp1 ) &
                                    & + g2x*(-1d0)/(g2x + km_tmp1 )**2d0 * 1d0 &
                                    & ) &
                                    & *no3x/(no3x + km_tmp2 ) &
                                    & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('no3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & * ( & 
                                    & 1d0/(no3x + km_tmp2 ) &
                                    & + no3x*(-1d0)/(no3x + km_tmp2 )**2d0 * 1d0 &
                                    & ) &
                                    & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *no3x/(no3x + km_tmp2 ) &
                                    & *km_tmp3*(-1d0)/(po2x*ko2 + km_tmp3 )**2d0 * ko2 &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *no3x/(no3x + km_tmp2 ) &
                                    & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *dv_dph_tmp &! modified version using normal distribution with sigma = 1 
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                endselect 
                        
            case('g2n22') ! 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O)  
                scheme = 'maggi08' ! Maggi et al. (2008) wihtout baterial, pH and water saturation functions
                
                select case(trim(adjustl(scheme)))
                    case('maggi08')
                        v_tmp = 1.23d-7*60d0*60d0*24d0*365d0
                        km_tmp1 = 10d-5 * 1d6 ! mol L-1 converted to mol m-3
                        km_tmp2 = 11.3d-5
                        km_tmp3 = 2.52d-5
                        rxn_ext = ( &
                            & v_tmp &
                            & *g2x/(g2x + km_tmp1 ) &
                            & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                            ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                            & *km_tmp3/(no3x + km_tmp3 ) &
                            & *min(2d0*sat,1d0) &
                            & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                            & )
                        
                        ! when using modified version using normal distribution with sigma = 1 
                        dv_dph_tmp = exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &
                            & *-0.5d0*2d0*((-log10(prox)-7d0)/1d0)  &
                            & *(-1d0) &
                            & *1d0/log(10d0)/prox
                            
                        select case(trim(adjustl(sp_name)))
                            case('g2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & 1d0/(g2x + km_tmp1 ) &
                                    & + g2x*(-1d0)/(g2x + km_tmp1 )**2d0 * 1d0 &
                                    & ) &
                                    & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *km_tmp3/(no3x + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pn2o')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & * ( & 
                                    & kn2o/(kn2o*pn2ox + km_tmp2 ) &
                                    & + kn2o*pn2ox*(-1d0)/(kn2o*pn2ox + km_tmp2 )**2d0 * kn2o &
                                    & ) &
                                    ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *km_tmp3/(no3x + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            ! case('po2')
                                ! drxnext_dmsp = ( &
                                    ! & v_tmp &
                                    ! & *g2x/(g2x + km_tmp1 ) &
                                    ! & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3*(-1d0)/(po2x*ko2 + km_tmp3 )**2d0 * ko2 &
                                    ! & *min(2d0*sat,1d0) &
                                    ! & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    ! & )
                            case('no3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3*(-1d0)/(po2x*ko2 + km_tmp3 )**2d0 * ko2 &
                                    & *km_tmp3*(-1d0)/(no3x + km_tmp3 )**2d0 * 1d0 &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *km_tmp3/(no3x + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *dv_dph_tmp &! modified version using normal distribution with sigma = 1 
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                endselect 
            
            case('g2k','g2ca','g2mg')
                ! cation uptake by g2 
                ! 0.0395 Eq mol-1 C harvested from Kantzas et al. 2022 
                ! for now OM is parameterized as NPP x 1.5 
                ! harvest = NPP * HI * RS /DF / C (Monfreda et al. 2008)
                ! where 
                !   HI = harvest index (0.85 as assumed by Kantzas et al. 2022)
                !   RS = root:shoot ratio (0.25 assumed) 
                !   DF = dry proportion of economic yield (0.9 assumed, cf. Kroodsma and Field, 2006)  
                !   C  = carbon content (0.45 g C / g dry matter) 
                !   
                ! So, vmax = 0.0395 * kresp * g2 / 1.5 * 0.85 * 0.25 / 0.9 / 0.45 
                ! And assume ion uptake kinetics by Fageria et al. 2010 for wheat (data from Barber 1995)
                ! Michaelis-Menten formulation 
                !   maq/(Kmaq + maq) 
                !       where 
                !           maq is ion conc. (M)
                !           Kmaq is Michaeris constants (M)
                ! for now, only consider K, Ca and Mg uptake
                ! Maximum uptake is distributed between these ions by K:Ca:Mg = 7:1.6:0.4
                ! and Kmaq for K, Ca and Mg are 7e-6, 5e-6 and 1e-6 M, respectively.
                
                select case(trim(adjustl(rxn_name)))
                    case('g2k')
                        sp_tmp = 'k    '
                        km_tmp1 = 7d0/11d0
                        km_tmp2 = 7d-6
                    case('g2ca')
                        sp_tmp = 'ca   '
                        km_tmp1 = 1.6d0/11d0
                        km_tmp2 = 5d-6
                    case('g2mg')
                        sp_tmp = 'mg   '
                        km_tmp1 = 0.4d0/11d0
                        km_tmp2 = 1d-6
                    case default
                        print *, '*** FUNDAMENTAL ERROR in extra reaction: cation uptake'
                        stop
                endselect  
                
                vmax_tmp = 0.0395d0 * kg2  / 1.5d0 * 0.85d0 * 0.25d0 / 0.9d0 / 0.45d0
                
                rxn_ext = ( &
                    & + vmax_tmp * g2x  &
                    & *po2x/(mo2 + po2x)  &
                    & *km_tmp1  &
                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & /( km_tmp2 &
                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & ) &
                    & )
                    
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dpro_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & +  ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *(-1d0)/( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & )**2d0 &
                                    & ) &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dpro_loc(findloc(chraq_all,sp_tmp,dim=1),:) 
                            case('ios')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dios_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & +  ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *(-1d0)/( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & )**2d0 &
                                    & ) &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dios_loc(findloc(chraq_all,sp_tmp,dim=1),:) 
                            case('g2')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * 1d0   &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * g2x  &
                                    & *1d0/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & + ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x*(-1d0)/(mo2 + po2x)**2d0  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & + ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dmgas_loc( &
                                        & findloc(chraq_all,sp_tmp,dim=1) &
                                        & ,findloc(chrgas_all,'po2',dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) & 
                                    & + ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *(-1d0)/( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & )**2d0 &
                                    & ) &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dmgas_loc( &
                                    & findloc(chraq_all,sp_tmp,dim=1) &
                                    & ,findloc(chrgas_all,'po2',dim=1),:) 
                            case default
                                if (any( trim(adjustl(sp_name)) == chraq ) ) then
                                    if ( sp_name == sp_tmp ) then 
                                        drxnext_dmsp = ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                                & findloc(chraq_all,sp_tmp,dim=1) &
                                                & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                            & /( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & ) &
                                            & ) & 
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *(-1d0)/( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )**2d0 &
                                            & ) &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                            & findloc(chraq_all,sp_tmp,dim=1) &
                                            & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *1d0 &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & /( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & ) &
                                            & ) &
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *(-1d0)/( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )**2d0 &
                                            & *1d0 &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )
                                    else
                                        drxnext_dmsp = ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                                & findloc(chraq_all,sp_tmp,dim=1) &
                                                & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                            & /( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & ) &
                                            & ) & 
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *(-1d0)/( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )**2d0 &
                                            & ) &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                            & findloc(chraq_all,sp_tmp,dim=1) &
                                            & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:)
                                    endif 
                                elseif (any( trim(adjustl(sp_name)) == chrgas ) ) then
                                    drxnext_dmsp = ( &
                                        & + vmax_tmp * g2x  &
                                        & *po2x/(mo2 + po2x)  &
                                        & *km_tmp1  &
                                        & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *dmaqft_dmgas_loc( &
                                            & findloc(chraq_all,sp_tmp,dim=1) &
                                            & ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:) &
                                        & /( km_tmp2 &
                                        & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & ) &
                                        & ) & 
                                        & + ( &
                                        & + vmax_tmp * g2x  &
                                        & *po2x/(mo2 + po2x)  &
                                        & *km_tmp1  &
                                        & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *(-1d0)/( km_tmp2 &
                                        & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & )**2d0 &
                                        & ) &
                                        & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *dmaqft_dmgas_loc( &
                                        & findloc(chraq_all,sp_tmp,dim=1) &
                                        & ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:) 
                                else
                                    drxnext_dmsp = 0d0
                                endif 
                        endselect 
                
                
            case default 
                rxn_ext = 0d0
                drxnext_dmsp = 0d0
                
        endselect

        rxnext_error = .false.
        if (any(isnan(rxn_ext)) .or. any(isnan(drxnext_dmsp))) then 
            print *,'nan in calc_rxn_ext_dev_3'
            print *,'any(isnan(rxn_ext)) | any(isnan(drxnext_dmsp))'
            print *,rxn_name,sp_name,any(isnan(rxn_ext)),any(isnan(drxnext_dmsp)) &
                & ,any(isnan(dmaqft_dios_loc(findloc(chraq_all,sp_tmp,dim=1),:) ))
            rxnext_error = .true.
            ! stop
        endif 

    endsubroutine calc_rxn_ext_dev_3

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    subroutine alsilicate_aq_gas_1D_v3_2( &
        ! new input 
        & nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext &
        & ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  &
        & ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi &
        & ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp &
        & ,stgas_ext,stgas_dext,staq_ext,stsld_ext,staq_dext,stsld_dext &
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst &
        & ,chraq_cnst,chraq_all,chrgas_cnst,chrgas_all,chrsld_all &
        & ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
        & ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &
        & ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain,msldth_all,mv_all,staq_all,stgas_all &
        & ,trans,display,chrflx,sld_enforce &! input
        & ,nsld_kinspc,chrsld_kinspc,kin_sld_spc &! input
        & ,precstyle,solmod,fkin &! input
        !  old inputs
        & ,hr,poro,z,dz,w_btm,sat,pro,poroprev,tora,v,tol,it,nflx,kw,maqft_prev,disp & 
        & ,ucv,torg,cplprec,rg,tc,sec2yr,tempk_0,proi,poroi,up,dwn,cnr,adf,msldunit  &
        & ,ads_ON,maqfads_prev,keqcec_all,keqiex_all,cec_pH_depend,aq_close,ios,act_ON,beta_all & 
        ! old inout
        & ,dt,flgback,w &    
        ! output 
        & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec,rxnsld,flx_co2sp,maqft & 
        & ,maqfads,msldf_loc,beta_loc,iosx &
        & )
        ! this is an attempt to calculate mass balance based on specific primary variables for aq. species.  
        implicit none 

        integer,intent(in)::nz,nflx
        real(kind=8),intent(in)::w_btm,tol,kw,ucv,rho_grain,rg,tc,sec2yr,tempk_0,proi,poroi
        real(kind=8),dimension(nz),intent(in)::poro,z,sat,tora,v,poroprev,dz,torg,pro,up,dwn,cnr,adf,disp,ios
        real(kind=8),dimension(nz),intent(out)::prox
        real(kind=8),dimension(nz),intent(out)::iosx
        real(kind=8),dimension(nz),intent(inout)::w
        integer,intent(inout)::it
        integer iter
        logical,intent(in)::cplprec,display
        logical,intent(inout)::flgback
        character(3),intent(in)::msldunit
        real(kind=8),intent(in)::dt
        real(kind=8) error

        integer,intent(in)::nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext,nsld_kinspc
        character(5),dimension(nsp_sld),intent(in)::chrsld
        character(5),dimension(nsp_sld_2),intent(in)::chrsld_2
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_ph),intent(in)::chraq_ph
        character(5),dimension(nsp_gas_ph),intent(in)::chrgas_ph
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nrxn_ext),intent(in)::chrrxn_ext
        character(5),dimension(nsld_kinspc),intent(in)::chrsld_kinspc
        real(kind=8),dimension(nsp_sld),intent(in)::msldi,msldth,mv
        real(kind=8),dimension(nsp_aq),intent(in)::maqi,maqth,daq 
        real(kind=8),dimension(nsp_gas),intent(in)::mgasi,mgasth,dgasa,dgasg,khgasi
        real(kind=8),dimension(nsp_gas)::dgasi,dgasn
        real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
        real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas
        real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp 
        real(kind=8),dimension(nsp_sld,nz),intent(inout)::ksld
        real(kind=8),dimension(nz,nz,nsp_sld),intent(in)::trans
        real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx,omega,nonprec,rxnsld
        real(kind=8),dimension(nsp_sld,nz)::domega_dpro,dmsld,dksld_dpro,drxnsld_dmsld,dksld_dso4f,domega_dso4f,dksld_dios,domega_dios
        real(kind=8),dimension(nsp_sld,nsp_aq,nz)::domega_dmaq,dksld_dmaq,drxnsld_dmaq
        real(kind=8),dimension(nsp_sld,nsp_gas,nz)::domega_dmgas,dksld_dmgas,drxnsld_dmgas
        real(kind=8),dimension(nsp_sld,nflx,nz),intent(out)::flx_sld
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maq,maqsupp,maqft_prev,maqfads_prev
        real(kind=8),dimension(nsp_aq,nz),intent(inout)::maqx,maqft,maqfads 
        real(kind=8),dimension(nsp_aq,nz)::dprodmaq,dmaq,maqf,dmaqft_dpro,dmaqfads_dpro,maqx_save,dmaqx,dmaqft_dios,dmaqfads_dios
        real(kind=8),dimension(nsp_aq,nz)::diosdmaq
        real(kind=8),dimension(nsp_aq,nsp_aq,nz)::dmaqft_dmaqf,dmaqfads_dmaqf
        real(kind=8),dimension(nsp_aq,nsp_gas,nz)::dmaqft_dmgas,dmaqfads_dmgas
        real(kind=8),dimension(nsp_aq,nsp_sld,nz)::dmaqfads_dmsld
        real(kind=8),dimension(nsp_aq,nflx,nz),intent(out)::flx_aq
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,mgassupp
        real(kind=8),dimension(nsp_gas,nz),intent(inout)::mgasx 
        real(kind=8),dimension(nsp_gas,nz)::khgasx,khgas,dgas,agasx,agas,rxngas,dkhgas_dpro,dprodmgas,dmgas,dso4fdmgas,dkhgas_dso4f &
            & ,mgasx_save,dmgasx,dkhgas_dios
        real(kind=8),dimension(nsp_gas,nz)::diosdmgas
        real(kind=8),dimension(nsp_gas,nsp_aq,nz)::dkhgas_dmaq,ddgas_dmaq,dagas_dmaq,drxngas_dmaq 
        real(kind=8),dimension(nsp_gas,nsp_sld,nz)::drxngas_dmsld 
        real(kind=8),dimension(nsp_gas,nsp_gas,nz)::dkhgas_dmgas,ddgas_dmgas,dagas_dmgas,drxngas_dmgas 
        real(kind=8),dimension(nsp_gas,nflx,nz),intent(out)::flx_gas 
        real(kind=8),dimension(nrxn_ext,nz),intent(inout)::rxnext
        real(kind=8),dimension(nrxn_ext,nz)::drxnext_dpro,drxnext_dso4f,drxnext_dios
        real(kind=8),dimension(nrxn_ext,nsp_gas),intent(in)::stgas_ext,stgas_dext
        real(kind=8),dimension(nrxn_ext,nsp_aq),intent(in)::staq_ext,staq_dext
        real(kind=8),dimension(nrxn_ext,nsp_sld),intent(in)::stsld_ext,stsld_dext
        real(kind=8),dimension(nrxn_ext,nsp_gas,nz)::drxnext_dmgas
        real(kind=8),dimension(nrxn_ext,nsp_aq,nz)::drxnext_dmaq
        real(kind=8),dimension(nrxn_ext,nsp_sld,nz)::drxnext_dmsld
        real(kind=8),dimension(nsld_kinspc),intent(in)::kin_sld_spc

        integer,intent(in)::nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst,nsp_sld_cnst
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::staq_all
        real(kind=8),dimension(nsp_sld_all,nsp_gas_all),intent(in)::stgas_all
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_no3
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_nh3
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_oxa
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_cl
        real(kind=8),dimension(nsp_sld_all),intent(in)::keqsld_all,msldth_all,mv_all

        real(kind=8),dimension(nsp_sld_all),intent(in)::keqcec_all,beta_all
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::keqiex_all
        logical,dimension(nsp_sld_all),intent(in)::cec_pH_depend
        real(kind=8),dimension(nsp_sld_all,nz),intent(out)::msldf_loc,beta_loc

        real(kind=8),dimension(nsp_aq_all,nz)::dprodmaq_all,dso4fdmaq_all,diosdmaq_all
        real(kind=8),dimension(nsp_gas_all,nz)::dprodmgas_all,dso4fdmgas_all,diosdmgas_all

        real(kind=8),dimension(nz)::domega_dpro_loc,domega_dso4f_loc,domega_dios_loc
        real(kind=8),dimension(nsp_gas_all,nz)::domega_dmgas_all
        real(kind=8),dimension(nsp_aq_all,nz)::domega_dmaq_all
        real(kind=8),dimension(nsp_aq_all,nz)::maqft_loc,dmaqft_dpro_loc,maqf_loc,maqx_loc,maqx_loc_tmp,dmaqft_dios_loc
        real(kind=8),dimension(nsp_aq_all,nsp_aq_all,nz)::dmaqft_dmaqf_loc
        real(kind=8),dimension(nsp_aq_all,nsp_gas_all,nz)::dmaqft_dmgas_loc
        real(kind=8),dimension(nsp_aq_all,nz)::dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2
        real(kind=8),dimension(nsp_aq_cnst,nz)::maqcx

        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc
        real(kind=8),dimension(nsp_gas_all,nz)::khgas_all,khgasx_all,dkhgas_dpro_all,dkhgas_dso4f_all,dkhgas_dios_all
        real(kind=8),dimension(nsp_gas_all,nsp_aq_all,nz)::dkhgas_dmaq_all
        real(kind=8),dimension(nsp_gas_all,nsp_gas_all,nz)::dkhgas_dmgas_all

        real(kind=8),dimension(nsp_sld_all,nz)::msldx_loc
        real(kind=8),dimension(nsp_aq_all,nz)::maqfads_loc,dmaqfads_dpro_loc
        real(kind=8),dimension(nsp_aq_all,nsp_aq_all,nz)::dmaqfads_dmaqf_loc
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz)::dmaqfads_dmsld_loc

        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz)::maqfads_sld_loc
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nsp_aq_all,nz)::dmaqfads_sld_dmaqf_loc
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz)::dmaqfads_sld_dmsld_loc
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz)::dmaqfads_sld_dpro_loc

        real(kind=8),dimension(nsp_aq,nsp_sld,nz)::maqfads_sld
        real(kind=8),dimension(nsp_aq,nsp_sld,nsp_aq,nz)::dmaqfads_sld_dmaqf
        real(kind=8),dimension(nsp_aq,nsp_sld,nsp_gas,nz)::dmaqfads_sld_dmgas
        real(kind=8),dimension(nsp_aq,nsp_sld,nz)::dmaqfads_sld_dmsld
        real(kind=8),dimension(nsp_aq,nsp_sld,nz)::dmaqfads_sld_dpro,dmaqfads_sld_dios

        character(5),dimension(nflx),intent(in)::chrflx

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

        integer ieqaq_co3,ieqaq_hco3
        data ieqaq_co3,ieqaq_hco3/1,2/

        integer ieqaq_so4,ieqaq_so42
        data ieqaq_so4,ieqaq_so42/1,2/

        integer,intent(in)::nrxn_ext_all

        character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all

        real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
        real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all
        real(kind=8),dimension(nrxn_ext_all,nz),intent(in)::krxn1_ext_all,krxn2_ext_all

        real(kind=8),dimension(4,nflx,nz),intent(out)::flx_co2sp

        integer iz,row,ie,ie2,iflx,isps,ispa,ispg,ispa2,ispg2,col,irxn,isps2,iiz,isps_kinspc,row_w,col_w
        integer izp,izn
        integer::itflx,iadv,idif,irain,ires
        integer::ph_iter,ph_iter2
        data itflx,iadv,idif,irain/1,2,3,4/

        integer,dimension(nsp_sld)::irxn_sld 
        integer,dimension(nrxn_ext)::irxn_ext 

        real(kind=8),dimension(nsp_sld,nz),intent(in)::hr

        real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp &
            & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,omega_tmp_th,rxn_ext_tmp &
            & ,edif_tmp,edif_tmp_n,edif_tmp_p,khco2n_tmp,pco2n_tmp,edifn_tmp,caqsupp_tmp,kco2,k1,k2,kho,sw_red &
            & ,flx_max,flx_max_max,proi_tmp,knh3,k1nh3,kn2o,wp_tmp,w_tmp,sporo_tmp,sporop_tmp,sporoprev_tmp  &
            & ,mn_tmp,wn_tmp,sporon_tmp,caqdif_tmp_n

        real(kind=8),parameter::infinity = huge(0d0)
        real(kind=8),parameter::fact = 1d-3
        real(kind=8),parameter::dconc = 1d-14
        real(kind=8),parameter::maxfact = 1d200
        ! real(kind=8),parameter::threshold = log(maxfact)
        real(kind=8),parameter::threshold = 10d0
        ! real(kind=8),parameter::threshold = 3d0
        ! real(kind=8),parameter::corr = 1.5d0
        real(kind=8),parameter::corr = exp(threshold)

        real(kind=8),dimension(nz)::dummy,dummy2,dummy3,kin,dkin_dmsp,dumtest,sporo,prox_save,iosx_save

        logical print_cb,ph_error,omega_error,rxnext_error,ads_error
        character(500) print_loc
        character(20) chrfmt

        integer,parameter :: iter_max = 50
        ! integer,parameter :: iter_max = 300

        integer :: nz_disp = 10

        real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz),emx3(nsp3*nz)
        integer ipiv3(nsp3*nz)
        integer info 

        external DGESV

        logical::chkflx = .true.
        logical::dt_norm = .true.
        logical::kin_iter = .true.
        logical::new_gassol = .true.
        ! logical::new_gassol = .false.
        logical,intent(in)::ads_ON != .true.
        ! logical::ads_ON = .false.

        logical::ph_precalc = .true.
        ! logical::ph_precalc = .false.

        !*** Previously aq species could diffuse out of the soil surface as in marine sediment; now shut this down
        !*** For gas species, it must assume no flux of dissolved gaseuous species via diffusion
        !       Previously surface diffusion considers both gas + aq species: 
        !           0.5*[dgasi + dgas(1)], where dgasi is gas only but dgas(1) is mixture of aq + gas diffusions 
        !       To eliminate aq diffusion at the surface more completely, a change has made to calculate the surface diffusion coefficient as 
        !           0.5*[dgasi + dgasn(1)], where dgasn(1) considers only gas diffusion at the topmost soil layer excluding aqueous diffusion
        logical::aq_diff_close = .true. 
        ! logical::aq_diff_close = .false.

        !*** An attempt to implement closed system for gas species;
        !       Previously, gas exchange is allowed at the topmost soil layer via diffusion
        !       The switch below is to enable shutting down this exchange.
        ! logical::gas_close = .true.
        logical::gas_close = .false.

        ! logical::sld_enforce = .false.
        logical,intent(in)::sld_enforce != .true.

        ! logical::aq_close = .false.
        logical,intent(in)::aq_close != .true.
        logical,intent(in)::act_ON 

        character(10),dimension(nsp_sld),intent(in):: precstyle 
        real(kind=8),dimension(nsp_sld,nz),intent(in):: solmod,fkin ! factor to modify solubility used only to implement rxn rate law as defined by Emmanuel and Ague, 2011
        real(kind=8) msld_seed ,fact2
        ! real(kind=8):: fact_tol = 1d-3
        real(kind=8):: fact_tol = 1d-4
        real(kind=8):: dt_th = 1d-6
        real(kind=8):: flx_tol = 1d-4 != tol*fact_tol*(z(nz)+0.5d0*dz(nz))
        ! real(kind=8):: flx_tol = 1d-3 ! desparate to make things converge 
        ! real(kind=8):: flx_max_tol = 1d-9 != tol*fact_tol*(z(nz)+0.5d0*dz(nz)) ! working for most cases but not when spinup with N cycles
        real(kind=8):: flx_max_tol = 1d-6 != tol*fact_tol*(z(nz)+0.5d0*dz(nz)) 
        real(kind=8):: flx_max_max_tol = 1d-6 != tol*fact_tol*(z(nz)+0.5d0*dz(nz)) 
        integer solve_sld 

        real(kind=8):: sat_lim_prec = 1d50 ! maximum value of saturation state for minerals that can precipitate 
        real(kind=8):: sat_lim_noprec = 2d0 ! maximum value of saturation state for minerals that cannot precipitate 

        !-----------------------------------------------

        if (aq_close) chkflx = .false.

        !! added to enable colosed gas system 
        gas_close = .false.
        if (aq_close) gas_close = .true.

        msld_seed = 1d-20

        ! flx_tol = tol*fact_tol*(z(nz)+0.5d0*dz(nz))
        ! flx_tol = 1d-4

        if (sld_enforce) then 
            solve_sld = 0
        else
            solve_sld = 1
        endif 

        sw_red = 1d0
        sw_red = -1d100

        do isps=1,nsp_sld
            irxn_sld(isps) = 4+isps
        enddo 

        do irxn=1,nrxn_ext
            irxn_ext(irxn) = 4+nsp_sld+irxn
        enddo 

        ires = nflx

        print_cb = .false. 
        print_loc = './ph.txt'

        kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
        k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
        k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)

        kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)

        knh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h0)
        k1nh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h1)

        kn2o = keqgas_h(findloc(chrgas_all,'pn2o',dim=1),ieqgas_h0)

        sporo = 1d0 - poro
        if (msldunit=='blk') sporo = 1d0

        ! so4fprev = so4f

        ! w = win
            
        nonprec = 1d0 ! primary minerals only dissolve
        if (cplprec)then
            do isps = 1, nsp_sld
                if (any(chrsld_2 == chrsld(isps))) then  
                    nonprec(isps,:) = 0d0 ! allowing precipitation for secondary phases
                endif 
            enddo
        endif 

        prox = pro
        iosx = ios

        dummy = 0d0
        dummy2 = 0d0

        error = 1d4
        iter = 0

        ! print *, 'starting silciate calculation'
        ! stop

        do while ((.not.isnan(error)).and.(error > tol*fact_tol))

            amx3=0.0d0
            ymx3=0.0d0 
            emx3=0.0d0 
            
            flx_sld = 0d0
            flx_aq = 0d0
            flx_gas = 0d0
            
            ! precalculation of pH when iteration is not first time
            if (ph_precalc .and. iter/=0) then 
                dmaqx = maqx - maqx_save
                dmgasx = mgasx - mgasx_save
                do iz=1,nz
                    prox(iz) = prox(iz) * exp( &
                        & sum(dprodmaq(:,iz)*dmaqx(:,iz))/prox_save(iz) &
                        & + sum(dprodmgas(:,iz)*dmgasx(:,iz))/prox_save(iz) &
                        & )
                    iosx(iz) = iosx(iz) * exp( &
                        & sum(diosdmaq(:,iz)*dmaqx(:,iz))/iosx_save(iz) &
                        & + sum(diosdmgas(:,iz)*dmgasx(:,iz))/iosx_save(iz) &
                        & )
                enddo 
            endif 
            
            if (.not.act_ON) iosx = 0d0
            
            ! pH calculation and its derivative wrt aq and gas species
            
            ! print_cb = .true. 
            
            call calc_pH_v7_4( &
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
            
            ! print *
            ! print *, -log10(prox)
            ! print *, iosx
            ! print *,diosdmaq_all
            ! print *
            ! print *,diosdmgas_all
            ! stop

            if (ph_error) then 
                print *, 'error issued from ph calculation: raising flag and return to main' 
                flgback = .true.
                return
            endif 
            
            ! *** sanity check 
            if (any(isnan(prox)) .or. any(prox<=0d0)) then    
                print *, ' NAN or <=0 H+ conc.',any(isnan(prox)),any(prox<=0d0)
                print *,prox
                stop
            endif 
            
            dprodmaq = 0d0
            diosdmaq = 0d0
            do ispa=1,nsp_aq
                if (any (chraq_ph == chraq(ispa))) then 
                    dprodmaq(ispa,:)=dprodmaq_all(findloc(chraq_all,chraq(ispa),dim=1),:)
                    diosdmaq(ispa,:)=diosdmaq_all(findloc(chraq_all,chraq(ispa),dim=1),:)
                endif 
            enddo 
            
            dprodmgas = 0d0
            diosdmgas = 0d0
            do ispg=1,nsp_gas
                if (any (chrgas_ph == chrgas(ispg))) then 
                    dprodmgas(ispg,:)=dprodmgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                    diosdmgas(ispg,:)=diosdmgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                endif 
            enddo 
            
            ! saving maqx and mgasx
            maqx_save = maqx
            mgasx_save = mgasx
            prox_save = prox
            iosx_save = iosx
            
            ! getting mgasx_loc & maqx_loc
            call get_maqgasx_all( &
                & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
                & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
                & ,maqx,mgasx,maqc,mgasc &
                & ,maqx_loc,mgasx_loc  &! output
                & )
            
            ! getting maqft_loc and its derivatives
            call get_maqt_all( &
                & nz,nsp_aq_all,nsp_gas_all &
                & ,chraq_all,chrgas_all &
                & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
                & ,mgasx_loc,maqx_loc,prox,iosx,tc &
                & ,dmaqft_dpro_loc,dmaqft_dmaqf_loc,dmaqft_dmgas_loc,dmaqft_dios_loc &! output
                & ,maqft_loc  &! output
                & )
            ! if (any(isnan(maqft_loc))) then 
                ! print *,'nan in maqft_loc'
                ! print *,maqft_loc
                ! stop
            ! endif 
            maqft = 0d0
            
            dmaqft_dpro = 0d0
            dmaqft_dios = 0d0
            dmaqft_dmaqf = 0d0
            dmaqft_dmgas = 0d0
            do ispa=1,nsp_aq
                maqft(ispa,:)=maqft_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
                
                dmaqft_dpro(ispa,:)=dmaqft_dpro_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
                if (act_ON) dmaqft_dios(ispa,:)=dmaqft_dios_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
                
                do ispa2=1,nsp_aq
                    dmaqft_dmaqf(ispa,ispa2,:) = ( &
                        & + dmaqft_dmaqf_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chraq_all,chraq(ispa2),dim=1),:) &
                        & + dmaqft_dpro(ispa,:)*dprodmaq(ispa2,:) &
                        & + dmaqft_dios(ispa,:)*diosdmaq(ispa2,:) & 
                        & )
                enddo 
                do ispg=1,nsp_gas
                    dmaqft_dmgas(ispa,ispg,:) = ( &
                        & + dmaqft_dmgas_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrgas_all,chrgas(ispg),dim=1),:) &
                        & + dmaqft_dpro(ispa,:)*dprodmgas(ispg,:) &
                        & + dmaqft_dios(ispa,:)*diosdmgas(ispg,:) &
                        & )
                enddo 
            enddo
            
            ! print *,dprodmaq(findloc(chraq,'na',dim=1),:)
            ! print *
            ! print *,dmaqft_dmaqf(findloc(chraq,'na',dim=1),findloc(chraq,'na',dim=1),:)
            ! print *
            ! print *,dprodmaq(findloc(chraq,'na',dim=1),:) &
                ! & /(maqx(findloc(chraq,'na',dim=1),:)*dmaqft_dmaqf(findloc(chraq,'na',dim=1),findloc(chraq,'na',dim=1),:) & 
                ! & + 1d0*maqft(findloc(chraq,'na',dim=1),:) )
            
            
            
            !!!  for adsorption 
            if (ads_ON) then 
                call get_msldx_all( &
                    & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
                    & ,chrsld,chrsld_all,chrsld_cnst &
                    & ,msldx,msldc &
                    & ,msldx_loc  &! output
                    & )

                ! call get_maqads_all_v3( &
                call get_maqads_all_v4( &
                ! call get_maqads_all_v4a( &
                    & nz,nsp_aq_all,nsp_sld_all &
                    & ,chraq_all,chrsld_all &
                    & ,keqcec_all,keqiex_all,cec_pH_depend,beta_all &
                    & ,msldx_loc,maqx_loc,prox &
                    & ,dmaqfads_sld_dpro_loc,dmaqfads_sld_dmaqf_loc,dmaqfads_sld_dmsld_loc &! output
                    & ,msldf_loc,maqfads_sld_loc,beta_loc,ads_error  &! output
                    & )

                if (ads_error) then 
                    print *, 'error issued from adsorption calculation: raising flag and return to main' 
                    flgback = .true.
                    return
                endif 
                
                ! print *,maqfads_loc(findloc(chraq_all,'na',dim=1),:)
                ! print *,dmaqfads_dpro_loc(findloc(chraq_all,'na',dim=1),:)
                ! print *,dmaqfads_dmaqf_loc(findloc(chraq_all,'na',dim=1),findloc(chraq_all,'na',dim=1),:)
                ! print *,dmaqfads_dmsld_loc(findloc(chraq_all,'na',dim=1),findloc(chrsld_all,'ka',dim=1),:) 
                
                maqfads_sld = 0d0
                dmaqfads_sld_dpro = 0d0
                dmaqfads_sld_dios = 0d0
                dmaqfads_sld_dmaqf = 0d0
                dmaqfads_sld_dmgas = 0d0
                dmaqfads_sld_dmsld = 0d0
                do ispa=1,nsp_aq
                    do isps=1,nsp_sld
                        maqfads_sld(ispa,isps,:) &
                            & =maqfads_sld_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrsld_all,chrsld(isps),dim=1),:)
                        
                        dmaqfads_sld_dpro(ispa,isps,:) &
                            & =dmaqfads_sld_dpro_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrsld_all,chrsld(isps),dim=1),:)

                        do ispa2=1,nsp_aq
                            dmaqfads_sld_dmaqf(ispa,isps,ispa2,:) &
                                & = dmaqfads_sld_dmaqf_loc( &
                                &       findloc(chraq_all,chraq(ispa),dim=1) &
                                &       ,findloc(chrsld_all,chrsld(isps),dim=1) &
                                &       ,findloc(chraq_all,chraq(ispa2),dim=1) &
                                &       ,:) &
                                &   +  dmaqfads_sld_dpro(ispa,isps,:)*dprodmaq(ispa2,:) 
                        enddo 
                        
                        do ispg=1,nsp_gas
                            dmaqfads_sld_dmgas(ispa,isps,ispg,:) &
                                & = dmaqfads_sld_dpro(ispa,isps,:)*dprodmgas(ispg,:) 
                        enddo
                        
                        dmaqfads_sld_dmsld(ispa,isps,:) &
                            & = dmaqfads_sld_dmsld_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrsld_all,chrsld(isps),dim=1),:) 
                        
                    enddo
                enddo 
            else
                maqfads_sld = 0d0
                dmaqfads_sld_dpro = 0d0
                dmaqfads_sld_dios = 0d0
                dmaqfads_sld_dmaqf = 0d0
                dmaqfads_sld_dmgas = 0d0
                dmaqfads_sld_dmsld = 0d0
            endif 
                
            maqfads = 0d0
            dmaqfads_dpro = 0d0
            dmaqfads_dios = 0d0
            dmaqfads_dmaqf = 0d0
            dmaqfads_dmgas = 0d0
            dmaqfads_dmsld = 0d0
            do ispa=1,nsp_aq
                do iz=1,nz
                    maqfads(ispa,iz) = sum(maqfads_sld(ispa,:,iz))
                    dmaqfads_dpro(ispa,iz) = sum(dmaqfads_sld_dpro(ispa,:,iz))
                    if (act_ON) dmaqfads_dios(ispa,iz) = sum(dmaqfads_sld_dios(ispa,:,iz))
                    do ispa2=1,nsp_aq
                        dmaqfads_dmaqf(ispa,ispa2,iz) = sum(dmaqfads_sld_dmaqf(ispa,:,ispa2,iz))
                    enddo 
                    do ispg=1,nsp_gas
                        dmaqfads_dmgas(ispa,ispg,iz) = sum(dmaqfads_sld_dmgas(ispa,:,ispg,iz))
                    enddo 
                enddo
                dmaqfads_dmsld(ispa,:,:) = dmaqfads_sld_dmsld(ispa,:,:)
            enddo 
            ! stop
            
            ! recalculation of rate constants for mineral reactions
            if (kin_iter) then 
                ksld = 0d0
                dksld_dpro = 0d0
                dksld_dios = 0d0
                dksld_dmaq = 0d0
                dksld_dmgas = 0d0
                
                do isps =1,nsp_sld 
                    call sld_kin( &
                        & nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) &! input
                        & ,nsp_gas_all,chrgas_all,mgasx_loc &! input
                        & ,nsp_aq_all,chraq_all,maqx_loc &! input
                        & ,chrsld(isps),'pro  ' &! input 
                        & ,kin,dkin_dmsp &! output
                        & ) 
                    ksld(isps,:) = kin              *fkin(isps,:)
                    dksld_dpro(isps,:) = dkin_dmsp  *fkin(isps,:)
                    
                    do ispa = 1,nsp_aq
                        if (any (chraq_ph == chraq(ispa)) .or. staq(isps,ispa)/=0d0 ) then 
                            call sld_kin( &
                                & nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) &! input
                                & ,nsp_gas_all,chrgas_all,mgasx_loc &! input
                                & ,nsp_aq_all,chraq_all,maqx_loc &! input
                                & ,chrsld(isps),chraq(ispa) &! input 
                                & ,kin,dkin_dmsp &! output
                                & ) 
                            dksld_dmaq(isps,ispa,:) = dkin_dmsp *fkin(isps,:) + ( &
                                & dksld_dpro(isps,:)*dprodmaq(ispa,:) &
                                & )
                        endif 
                    enddo 
                    
                    do ispg = 1,nsp_gas
                        if (any (chrgas_ph == chrgas(ispg)) .or. stgas(isps,ispg)/=0d0) then 
                            call sld_kin( &
                                & nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) &! input
                                & ,nsp_gas_all,chrgas_all,mgasx_loc &! input
                                & ,nsp_aq_all,chraq_all,maqx_loc &! input
                                & ,chrsld(isps),chrgas(ispg) &! input 
                                & ,kin,dkin_dmsp &! output
                                & ) 
                            dksld_dmgas(isps,ispg,:) = dkin_dmsp *fkin(isps,:) + ( &
                                & dksld_dpro(isps,:)*dprodmgas(ispg,:) &
                                & )
                        endif 
                    enddo 
                
                enddo 
            else 
                dksld_dpro = 0d0
                dksld_dios = 0d0
                dksld_dmaq = 0d0
                dksld_dmgas = 0d0
            endif 
            
            ! if kin const. is specified in input file 
            if (nsld_kinspc > 0) then 
                do isps_kinspc=1,nsld_kinspc    
                    if ( any( chrsld == chrsld_kinspc(isps_kinspc))) then 
                        select case (trim(adjustl(chrsld_kinspc(isps_kinspc))))
                            case('g1','g2','g3') ! for OMs, turn over year needs to be provided [yr]
                                if (kin_sld_spc(isps_kinspc)/=0d0) then  
                                    ksld(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = ( &                   
                                        & 1d0/kin_sld_spc(isps_kinspc) &
                                        & ) 
                                else
                                    ksld(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = kin_sld_spc(isps_kinspc)
                                endif 
                                dksld_dpro(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = 0d0
                                dksld_dios(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = 0d0
                                dksld_dmaq(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:,:) = 0d0
                                dksld_dmgas(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:,:) = 0d0
                            case default ! otherwise, usual rate constant [mol/m2/yr]
                                ksld(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = ( &                            
                                    & kin_sld_spc(isps_kinspc) &
                                    & ) 
                                dksld_dpro(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = 0d0
                                dksld_dios(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = 0d0
                                dksld_dmaq(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:,:) = 0d0
                                dksld_dmgas(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:,:) = 0d0
                        end select 
                    endif 
                enddo 
            endif 
            
            ! print *,dksld_dmaq(findloc(chrsld,'ab',dim=1),findloc(chraq,'na',dim=1),:)
            ! print *
            ! print *,dksld_dmaq(findloc(chrsld,'ab',dim=1),findloc(chraq,'no3',dim=1),:)
            ! print *
                    
            ! *** sanity check *** 
            if (any(isnan(ksld)) .or. any(ksld>infinity)) then 
                print *,' *** found insanity in ksld: listing below -- '
                do isps=1,nsp_sld
                    do iz=1,nz
                        if (isnan(ksld(isps,iz)) .or. ksld(isps,iz)>infinity) print*,chrsld(isps),iz,ksld(isps,iz)
                    enddo
                enddo 
                stop
            ! else 
                ! print *,' *** found sanity in ksld -- '
            endif 
            ! print *, 'main loop'
            ! print *, ksld(findloc(chrsld,'kfs',dim=1),:)
            
            ! saturation state calc. and their derivatives wrt aq and gas species
            
            ! print *,'ksld',ksld(findloc(chrsld,'gt',dim=1),:)
            
            omega = 0d0
            domega_dpro = 0d0
            domega_dios = 0d0
            domega_dmaq = 0d0
            domega_dmgas = 0d0
            
            do isps =1, nsp_sld
                
                dummy = 0d0
                domega_dpro_loc = 0d0
                domega_dios_loc = 0d0
                call calc_omega_v5( &
                    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst & 
                    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &
                    & ,maqx,maqc,mgasx,mgasc,mgasth_all,prox,iosx,tc &
                    & ,keqsld_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 &
                    & ,staq_all,stgas_all &
                    & ,chrsld(isps) &
                    & ,domega_dmaq_all,domega_dmgas_all,domega_dpro_loc,domega_dios_loc &! output
                    & ,dummy,omega_error &! output
                    & )
                if (omega_error) then
                    flgback = .true.
                    return 
                endif 
                omega(isps,:) = dummy
                domega_dpro(isps,:) = domega_dpro_loc
                if (act_ON) domega_dios(isps,:) = domega_dios_loc
                
                do ispa = 1, nsp_aq
                    if (any (chraq_ph == chraq(ispa)) .or. staq(isps,ispa)/=0d0 ) then 
                    
                        domega_dmaq(isps,ispa,:) = domega_dmaq_all(findloc(chraq_all,chraq(ispa),dim=1),:)+ ( &
                            & + domega_dpro(isps,:)*dprodmaq(ispa,:) &
                            & + domega_dios(isps,:)*diosdmaq(ispa,:) &
                            & )
                        
                        
                    endif 
                enddo
                do ispg = 1, nsp_gas
                    if (any (chrgas_ph == chrgas(ispg)) .or. stgas(isps,ispg)/=0d0) then 
                    
                        domega_dmgas(isps,ispg,:) = domega_dmgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)+ ( &
                            & + domega_dpro(isps,:)*dprodmgas(ispg,:) &
                            & + domega_dios(isps,:)*diosdmgas(ispg,:) &
                            & )
                    endif 
                enddo
            enddo 
            
            
            ! *** reducing saturation ***    
            do isps=1,nsp_sld
                dummy = 0d0
                if (any(chrsld_2 == chrsld(isps))) then  ! chrsld(isps) is included in secondary minerals
                    ! cycle
                    do iz=1,nz
                        if (omega(isps,iz)>=sat_lim_prec) omega(isps,iz) = sat_lim_prec
                    enddo 
                else
                    ! omega(isps,:) = dummy
                    do iz=1,nz
                        if (omega(isps,iz)>=sat_lim_noprec) omega(isps,iz) = sat_lim_noprec
                    enddo 
                endif 
            enddo 
                        
            ! *** sanity check ***     
            if (any(isnan(omega))) then 
                print *,' *** found NAN in omega: listing below -- '
                do isps=1,nsp_sld
                    do iz=1,nz
                        if (isnan(omega(isps,iz))) print*,chrsld(isps),iz,omega(isps,iz)
                    enddo
                enddo 
                stop
            endif 
            if (any(omega>infinity)) then 
                print *,' *** found INF in omega  '
                stop
                print *,' *** proceed maximum saturation 1d+100 if precipitating while 1d1 if not'
                do isps=1,nsp_sld
                    dummy = 0d0
                    if (any(omega(isps,:)>infinity)) then 
                        dummy = omega(isps,:)
                        if (any(chrsld_2 == chrsld(isps))) then  ! chrsld(isps) is included in secondary minerals
                            print *,chrsld(isps),' (precipitation allowed)'
                            where(dummy>infinity)
                                dummy = sat_lim_prec
                            endwhere
                        else
                            print *,chrsld(isps),' (precipitation not allowed)'
                            where(dummy>infinity)
                                dummy = sat_lim_noprec
                            endwhere
                        endif 
                        if (any(dummy>infinity)) then 
                            print *, 'somthing is wrong'
                            stop
                        endif 
                        omega(isps,:) = dummy
                    endif 
                enddo 
            endif 
            
            
            ! adding reactions that are not based on dis/prec of minerals
            rxnext = 0d0
            drxnext_dpro = 0d0
            drxnext_dios = 0d0
            drxnext_dmaq = 0d0
            drxnext_dmgas = 0d0
            drxnext_dmsld = 0d0
            
            do irxn=1,nrxn_ext
                dummy = 0d0
                dummy2 = 0d0
                call calc_rxn_ext_dev_3( &
                    & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                    & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                    & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                    & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                    & ,rg,tempk_0,tc,iosx &!input
                    & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
                    & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
                    & ,chrrxn_ext(irxn),'pro  ' &! input 
                    & ,dummy,dummy2,rxnext_error &! output
                    & )
                if (rxnext_error) then
                    flgback = .true.
                    return 
                endif 
                rxnext(irxn,:) = dummy
                drxnext_dpro(irxn,:) = dummy2
                
                if (act_ON) then 
                    dummy = 0d0
                    dummy2 = 0d0
                    call calc_rxn_ext_dev_3( &
                        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                        & ,rg,tempk_0,tc,iosx &!input
                        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
                        & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
                        & ,chrrxn_ext(irxn),'ios  ' &! input 
                        & ,dummy,dummy2,rxnext_error &! output
                        & )
                    if (rxnext_error) then
                        flgback = .true.
                        return 
                    endif 
                    drxnext_dios(irxn,:) = dummy2
                endif 
                
                do ispg=1,nsp_gas
                    ! if (stgas_dext(irxn,ispg)==0d0) cycle
                    
                    dummy = 0d0
                    dummy2 = 0d0
                    call calc_rxn_ext_dev_3( &
                        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                        & ,rg,tempk_0,tc,iosx &!input
                        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
                        & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
                        & ,chrrxn_ext(irxn),chrgas(ispg) &! input 
                        & ,dummy,dummy2,rxnext_error &! output
                        & )
                    if (rxnext_error) then
                        flgback = .true.
                        return 
                    endif 
                    drxnext_dmgas(irxn,ispg,:) = dummy2 + (&
                        & + drxnext_dpro(irxn,:)*dprodmgas(ispg,:) &
                        & + drxnext_dios(irxn,:)*diosdmgas(ispg,:) &
                        & )
                enddo 
                
                do ispa=1,nsp_aq
                    ! if (staq_dext(irxn,ispa)==0d0) cycle
                    
                    dummy = 0d0
                    dummy2 = 0d0
                    call calc_rxn_ext_dev_3( &
                        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                        & ,rg,tempk_0,tc,iosx &!input
                        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
                        & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
                        & ,chrrxn_ext(irxn),chraq(ispa) &! input 
                        & ,dummy,dummy2,rxnext_error &! output
                        & )
                    if (rxnext_error) then
                        flgback = .true.
                        return 
                    endif 
                    drxnext_dmaq(irxn,ispa,:) = dummy2 + ( &
                        & + drxnext_dpro(irxn,:)*dprodmaq(ispa,:) &
                        & + drxnext_dios(irxn,:)*diosdmaq(ispa,:) &
                        & )
                enddo 
                
                do isps=1,nsp_sld
                    ! if (stsld_dext(irxn,isps)==0d0) cycle
                    
                    dummy = 0d0
                    dummy2 = 0d0
                    call calc_rxn_ext_dev_3( &
                        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                        & ,rg,tempk_0,tc,iosx &!input
                        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
                        & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
                        & ,chrrxn_ext(irxn),chrsld(isps) &! input 
                        & ,dummy,dummy2,rxnext_error &! output
                        & )
                    if (rxnext_error) then
                        flgback = .true.
                        return 
                    endif 
                    drxnext_dmsld(irxn,isps,:) = dummy2
                enddo 
            enddo 
            
            ! gas tansport
            khgas = 0d0
            khgasx = 0d0
            dkhgas_dmaq = 0d0
            dkhgas_dmgas = 0d0
            ! added
            dkhgas_dpro = 0d0
            dkhgas_dios = 0d0
            
            if (new_gassol) then 
                call calc_khgas_all_v2( &
                    & nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst &
                    & ,chraq_all,chrgas_all,chraq_cnst,chrgas_cnst,chraq,chrgas &
                    & ,maq,mgas,maqx,mgasx,maqc,mgasc &
                    & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3  &
                    & ,pro,prox,ios,iosx,tc &
                    & ,khgas_all,khgasx_all,dkhgas_dpro_all,dkhgas_dmaq_all,dkhgas_dmgas_all,dkhgas_dios_all &!output
                    & )
                    
                do ispg=1,nsp_gas
                    khgas(ispg,:)=khgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                    khgasx(ispg,:)=khgasx_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                    dkhgas_dpro(ispg,:)=dkhgas_dpro_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                    if (act_ON) dkhgas_dios(ispg,:)=dkhgas_dios_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                    do ispa=1,nsp_aq
                        dkhgas_dmaq(ispg,ispa,:)= ( &
                            & + dkhgas_dmaq_all(findloc(chrgas_all,chrgas(ispg),dim=1),findloc(chraq_all,chraq(ispa),dim=1),:) &
                            & + dkhgas_dpro(ispg,:)*dprodmaq(ispa,:) &
                            & + dkhgas_dios(ispg,:)*diosdmaq(ispa,:) &
                            & )
                    enddo 
                    do ispg2=1,nsp_gas
                        dkhgas_dmgas(ispg,ispg2,:)= ( &
                            & + dkhgas_dmgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),findloc(chrgas_all,chrgas(ispg2),dim=1),:) &
                            & + dkhgas_dpro(ispg,:)*dprodmgas(ispg2,:) & 
                            & + dkhgas_dios(ispg,:)*diosdmgas(ispg2,:) &
                            & )
                    enddo 
                enddo 
            endif
            
            dgas = 0d0
            ddgas_dmaq = 0d0
            ddgas_dmgas = 0d0
            
            agas = 0d0
            agasx = 0d0
            dagas_dmaq = 0d0
            dagas_dmgas = 0d0
            
            do ispg = 1, nsp_gas
                
                if (.not. new_gassol) then ! old way to calc solubility (to be removed?)
                    select case (trim(adjustl(chrgas(ispg))))
                        case('pco2')
                            khgas(ispg,:) = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
                            khgasx(ispg,:) = kco2*(1d0+k1/prox + k1*k2/prox/prox)
                    
                            dkhgas_dpro(ispg,:) = kco2*(k1*(-1d0)/prox**2d0 + k1*k2*(-2d0)/prox**3d0)
                        case('po2')
                            khgas(ispg,:) = kho ! previous value; should not change through iterations 
                            khgasx(ispg,:) = kho
                    
                            dkhgas_dpro(ispg,:) = 0d0
                        case('pnh3')
                            khgas(ispg,:) = knh3*(1d0+pro/k1nh3) ! previous value; should not change through iterations 
                            khgasx(ispg,:) = knh3*(1d0+prox/k1nh3)
                    
                            dkhgas_dpro(ispg,:) = knh3*(1d0/k1nh3)
                        case('pn2o')
                            khgas(ispg,:) = kn2o ! previous value; should not change through iterations 
                            khgasx(ispg,:) = kn2o
                    
                            dkhgas_dpro(ispg,:) = 0d0
                    endselect 
                endif 
                
                dgas(ispg,:) = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg)+poro*sat*khgasx(ispg,:)*1d3*(tora*dgasa(ispg)+disp)  !! effective gas + aq diffusion
                dgasi(ispg) = ucv*1d3*dgasg(ispg)   !! gas diffusion alone in air 
                dgasn(ispg) = ucv*poro(1)*(1.0d0-sat(1))*1d3*torg(1)*dgasg(ispg)  ! gas diffusion alone in soil air at the upper most layer
                
                agas(ispg,:)= ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khgas(ispg,:)*1d3
                agasx(ispg,:)= ucv*poro*(1.0d0-sat)*1d3+poro*sat*khgasx(ispg,:)*1d3
                
                do ispa = 1,nsp_aq 
                    if (.not. new_gassol) dkhgas_dmaq(ispg,ispa,:) = dkhgas_dpro(ispg,:)*dprodmaq(ispa,:) ! old way to calc solubility (to be removed?)
                    ddgas_dmaq(ispg,ispa,:) = poro*sat*dkhgas_dmaq(ispg,ispa,:)*1d3*(tora*dgasa(ispg)+disp)
                    dagas_dmaq(ispg,ispa,:) =  poro*sat*dkhgas_dmaq(ispg,ispa,:)*1d3
                enddo 
                
                do ispg2 = 1,nsp_gas 
                    if (.not. new_gassol) dkhgas_dmgas(ispg,ispg2,:) = dkhgas_dpro(ispg,:)*dprodmgas(ispg2,:) ! old way to calc solubility (to be removed?)
                    ddgas_dmgas(ispg,ispg2,:) = poro*sat*dkhgas_dmgas(ispg,ispg2,:)*1d3*(tora*dgasa(ispg)+disp)
                    dagas_dmgas(ispg,ispg2,:) =  poro*sat*dkhgas_dmgas(ispg,ispg2,:)*1d3
                enddo 
            enddo 
            
            ! sld phase reactions
            
            rxnsld = 0d0
            drxnsld_dmsld = 0d0
            drxnsld_dmaq = 0d0
            drxnsld_dmgas = 0d0
            
            call sld_rxn( &
                & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
                & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle,solmod &! input
                & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas,chrsld &! input
                & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
                & ) 
                        
            ! *** sanity check ***     
            if (any(isnan(rxnsld)) .or. any(rxnsld>infinity)) then 
                print *,' *** found insanity in rxnsld: listing below -- '
                do isps=1,nsp_sld
                    do iz=1,nz
                        if (isnan(rxnsld(isps,iz)) .or. rxnsld(isps,iz)>infinity) print*,chrsld(isps),iz,rxnsld(isps,iz)
                    enddo
                enddo 
                stop
            endif 
            
            ! gas reactions 
            
            rxngas = 0d0
            drxngas_dmaq = 0d0
            drxngas_dmsld = 0d0
            drxngas_dmgas = 0d0
                
            do ispg = 1, nsp_gas
                do isps = 1, nsp_sld
                    rxngas(ispg,:) =  rxngas(ispg,:) + (&
                        ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                        ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & + stgas(isps,ispg)*rxnsld(isps,:) &
                        & )
                    drxngas_dmsld(ispg,isps,:) =  drxngas_dmsld(ispg,isps,:) + (&
                        ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & + stgas(isps,ispg)*drxnsld_dmsld(isps,:) &
                        & )
                    do ispg2 = 1,nsp_gas
                        drxngas_dmgas(ispg,ispg2,:) =  drxngas_dmgas(ispg,ispg2,:) + (&
                            ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg2,:)) &
                            ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            ! & + stgas(isps,ispg)*dksld_dmgas(isps,ispg2,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + stgas(isps,ispg)*drxnsld_dmgas(isps,ispg2,:) &
                            & )
                    enddo 
                    do ispa = 1,nsp_aq
                        drxngas_dmaq(ispg,ispa,:) =  drxngas_dmaq(ispg,ispa,:) + ( &
                            ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            ! & + stgas(isps,ispg)*dksld_dmaq(isps,ispa,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + stgas(isps,ispg)*drxnsld_dmaq(isps,ispa,:) &
                            & )
                    enddo 
                enddo 
            enddo 
                    
            if (.not.sld_enforce) then 

                do iz = 1, nz  !================================
                        
                    izp = iz+1
                    izn = iz-1
                    
                    if (iz==1)  izn = iz
                    if (iz==nz) izp = iz
                    
                    do isps = 1, nsp_sld
                    
                        row = nsp3*(iz-1)+isps
                        
                        m_tmp           = msldx(isps,iz) 
                        mth_tmp         = msldth(isps) 
                        mi_tmp          = msldi(isps)
                        mp_tmp          = msldx(isps,izp)
                        msupp_tmp       = msldsupp(isps,iz) 
                        rxn_ext_tmp     = sum(stsld_ext(:,isps)*rxnext(:,iz))
                        mprev_tmp       = msld(isps,iz)  
                        w_tmp           = w(iz) 
                        wp_tmp          = w(izp) 
                        sporo_tmp       = 1d0-poro(iz)
                        sporop_tmp      = 1d0-poro(izp) 
                        sporoprev_tmp   = 1d0-poroprev(iz)
                        
                        if (iz==nz) then 
                            mp_tmp      = mi_tmp
                            wp_tmp      = w_btm 
                            sporop_tmp  = 1d0- poroi
                        endif 
                        
                        if (msldunit == 'blk') then 
                            sporo_tmp       = 1d0
                            sporop_tmp      = 1d0
                            sporoprev_tmp   = 1d0
                        endif 

                        amx3(row,row) = ( &
                            & 1d0 *  sporo_tmp /merge(1d0,dt,dt_norm)     &
                            ! & + adf(iz)*up(iz)*sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    &
                            ! & - adf(iz)*dwn(iz)*sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    &
                            & + sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    &
                            & + drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) &
                            & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & * merge(1.0d0,m_tmp,m_tmp<mth_tmp*sw_red)

                        ymx3(row) = ( &
                            & ( sporo_tmp*m_tmp - sporoprev_tmp*mprev_tmp )/merge(1d0,dt,dt_norm) &
                            & - ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  &
                            ! & - adf(iz)*up(iz)*( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  &
                            ! & - adf(iz)*dwn(iz)*( sporo_tmp*w_tmp* m_tmp - sporon_tmp*wn_tmp*mn_tmp )/dz(iz)*merge(dt,1d0,dt_norm)  &
                            ! & - adf(iz)*cnr(iz)*( sporop_tmp*wp_tmp*mp_tmp - sporon_tmp*wn_tmp*mn_tmp )/dz(iz)*merge(dt,1d0,dt_norm)  &
                            & + rxnsld(isps,iz)*merge(dt,1d0,dt_norm) &
                            & -msupp_tmp*merge(dt,1d0,dt_norm)  &
                            & -rxn_ext_tmp*merge(dt,1d0,dt_norm)  &
                            & ) &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                            
                        if (iz/=nz) amx3(row,row+nsp3) = ( &
                            & (- sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & (- adf(iz)*up(iz)* sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & +(- adf(iz)*cnr(iz)* sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(1.0d0,mp_tmp,m_tmp<mth_tmp*sw_red)
                            
                        ! if (iz/=1) amx3(row,row-nsp3) = ( &
                            ! & (+ adf(iz)*dwn(iz)* sporon_tmp*wn_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & +(+ adf(iz)*cnr(iz)* sporon_tmp*wn_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & ) &
                            ! & *merge(1.0d0,mn_tmp,m_tmp<mth_tmp*sw_red)
                        
                        do ispa = 1, nsp_aq
                            col = nsp3*(iz-1) + nsp_sld + ispa
                            
                            amx3(row,col ) = ( &
                                & + drxnsld_dmaq(isps,ispa,iz)*merge(dt,1d0,dt_norm) &
                                & - sum(stsld_ext(:,isps)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *maqx(ispa,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        enddo 
                        
                        do ispg = 1, nsp_gas 
                            col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                            amx3(row,col) = ( &
                                & + drxnsld_dmgas(isps,ispg,iz)*merge(dt,1d0,dt_norm) &
                                & - sum(stsld_ext(:,isps)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *mgasx(ispg,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        enddo 
                        
                        do isps2 = 1,nsp_sld 
                            if (isps2 == isps) cycle
                            col = nsp3*(iz-1)+ isps2

                            amx3(row,col) = ( &
                                & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps2,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *msldx(isps2,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        enddo 

        #ifdef calcw_full
                        col =  nsp3*(iz-1)+ nsp3
                        amx3(row,col) = ( &
                            & - ( - sporo_tmp* m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  &
                            & ) &
                            ! & * w_tmp &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                            
                        if (iz/=nz) amx3(row,col+nsp3) = ( &
                            & (- sporop_tmp*mp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & (- adf(iz)*up(iz)* sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & +(- adf(iz)*cnr(iz)* sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            ! & *wp_tmp  &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
        #endif 
                        ! diffusion terms are filled with transition matrices 
                        ! if (turbo2(isps).or.labs(isps)) then
                            ! do iiz = 1, nz
                                ! col = nsp3*(iiz-1)+isps
                                ! if (trans(iiz,iz,isps)==0d0) cycle
                                ! amx3(row,col) = amx3(row,col) &
                                    ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz)
                                ! ymx3(row) = ymx3(row) &
                                    ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz)
                                    
                                ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                                    ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz) &
                                    ! & )
                            ! enddo
                        ! else
                            ! do iiz = 1, nz
                                ! col = nsp3*(iiz-1)+isps
                                ! if (trans(iiz,iz,isps)==0d0) cycle
                                    
                                ! amx3(row,col) = amx3(row,col) -trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                                    ! & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                                ! ymx3(row) = ymx3(row) - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                                    ! & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                                    
                                ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                                    ! & - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                                    ! & )
                            ! enddo
                        ! endif
                        
                        ! modifications with porosity and dz are made in make_trans subroutine
                        do iiz = 1, nz
                            col = nsp3*(iiz-1)+isps
                            if (trans(iiz,iz,isps)==0d0) cycle
                                
                            amx3(row,col) = amx3(row,col) - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz)* merge(dt,1d0,dt_norm) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                            ymx3(row) = ymx3(row) - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz)* merge(dt,1d0,dt_norm) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                                
                            flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                                & - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz) &
                                & )
                        enddo
                        
                        flx_sld(isps,itflx,iz) = ( &
                            & ( sporo_tmp*m_tmp- sporoprev_tmp*mprev_tmp)/dt &
                            & )
                        flx_sld(isps,iadv,iz) = ( &
                            & - ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)  &
                            ! & - adf(iz)*up(iz)*( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)  &
                            ! & - adf(iz)*dwn(iz)*( sporo_tmp*w_tmp* m_tmp - sporon_tmp*wn_tmp*mn_tmp )/dz(iz)  &
                            ! & - adf(iz)*cnr(iz)*( sporop_tmp*wp_tmp*mp_tmp - sporon_tmp*wn_tmp*mn_tmp )/dz(iz)  &
                            & )
                        flx_sld(isps,irxn_sld(isps),iz) = ( &
                            & + rxnsld(isps,iz) &
                            & )
                        flx_sld(isps,irain,iz) = (&
                            & - msupp_tmp  &
                            & )
                        flx_sld(isps,irxn_ext(:),iz) = (&
                                & - stsld_ext(:,isps)*rxnext(:,iz)  &
                                & )
                        flx_sld(isps,ires,iz) = sum(flx_sld(isps,:,iz))
                        if (isnan(flx_sld(isps,ires,iz))) then 
                            print *,chrsld(isps),iz,(flx_sld(isps,iflx,iz),iflx=1,nflx)
                        endif 
                    enddo 
                end do  !================================
            
            endif 
            
        #ifdef calcw_full
            do iz=1,nz
                row = nsp3*(iz-1) + nsp3
                        
                w_tmp = w(iz) 
                wp_tmp = w(min(nz,iz+1)) 
                sporo_tmp = 1d0-poro(iz)
                sporop_tmp = 1d0-poro(min(nz,iz+1)) 
                sporoprev_tmp = 1d0-poroprev(iz)
                wn_tmp = w(max(1,iz-1))
                sporon_tmp = 1d0-poro(max(1,iz-1))
                
                if (iz==1) then 
                    wn_tmp = 0d0
                    sporon_tmp = 0d0
                endif 
                
                if (iz==nz) then 
                    wp_tmp = w_btm 
                    sporop_tmp = 1d0- poroi
                endif 
                        
                ymx3(row) = ( &
                    & ( sporo_tmp - sporoprev_tmp )/merge(1d0,dt,dt_norm) &
                    & - ( sporop_tmp*wp_tmp - sporo_tmp*w_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  &
                    & ) 
                    
                amx3(row,row) = amx3(row,row) + ( &
                    & - (- sporo_tmp*1d0)/dz(iz)*merge(dt,1d0,dt_norm)  &
                    & ) &
                    ! & *w_tmp &
                    & *1d0
                    
                if (iz/=nz) amx3(row,row+nsp3) = amx3(row,row) + ( &
                    & - ( sporop_tmp*1d0 )/dz(iz)*merge(dt,1d0,dt_norm)  &
                    & ) &
                    ! & *wp_tmp &
                    & *1d0
                    
                do isps = 1, nsp_sld
                    
                    col = nsp3*(iz-1)+isps
                    
                    k_tmp = ksld(isps,iz)
                    mv_tmp = mv(isps)
                    omega_tmp = omega(isps,iz)
                    omega_tmp_th = omega_tmp*nonprec(isps,iz)
                    m_tmp = msldx(isps,iz) 
                    mth_tmp = msldth(isps) 
                    mi_tmp = msldi(isps)
                    mp_tmp = msldx(isps,min(nz,iz+1))
                    msupp_tmp = msldsupp(isps,iz) 
                    rxn_ext_tmp = sum(stsld_ext(:,isps)*rxnext(:,iz))
                    mprev_tmp = msld(isps,iz)  

                    ymx3(row) = ymx3(row) + ( &
                        & + rxnsld(isps,iz)*merge(dt,1d0,dt_norm) &
                        & -msupp_tmp*merge(dt,1d0,dt_norm)  &
                        & -rxn_ext_tmp*merge(dt,1d0,dt_norm)  &
                        & ) &
                        & * mv(isps) * 1d-6 &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)

                    amx3(row,col) = amx3(row,col) + ( &  
                        & + drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) &
                        & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & * mv(isps) * 1d-6 &
                        & * merge(1.0d0,m_tmp,m_tmp<mth_tmp*sw_red)
                        
                    do iiz = 1, nz
                        col = nsp3*(iiz-1)+isps
                        ymx3(row) = ymx3(row) + ( &
                            & - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz)* merge(dt,1d0,dt_norm) &
                            & ) &
                            & * mv(isps) * 1d-6 &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        
                        amx3(row,col) = amx3(row,col) + ( &
                            & - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz)* merge(dt,1d0,dt_norm) &
                            & ) &
                            & * mv(isps) * 1d-6 &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    enddo 
                        
                    do ispa = 1, nsp_aq
                        col = nsp3*(iz-1) + nsp_sld + ispa
                        
                        amx3(row,col ) = amx3(row,col ) + ( &
                            & + drxnsld_dmaq(isps,ispa,iz)*merge(dt,1d0,dt_norm) &
                            & - sum(stsld_ext(:,isps)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *maqx(ispa,iz) &
                            & * mv(isps) * 1d-6 &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    enddo 
                    
                    do ispg = 1, nsp_gas 
                        col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                        amx3(row,col) = amx3(row,col ) + ( &
                            & + drxnsld_dmgas(isps,ispg,iz)*merge(dt,1d0,dt_norm) &
                            & - sum(stsld_ext(:,isps)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *mgasx(ispg,iz) &
                            & * mv(isps) * 1d-6 &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    enddo 
                    
                    do isps2 = 1,nsp_sld 
                        if (isps2 == isps) cycle
                        col = nsp3*(iz-1)+ isps2

                        amx3(row,col) = amx3(row,col ) + ( &
                            & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps2,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *msldx(isps2,iz) &
                            & * mv(isps) * 1d-6 &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    enddo 
                    
                enddo 
            
            
            enddo
        #endif 
            

            do iz = 1, nz
                        
                izp = iz+1
                izn = iz-1
                
                if (iz==1)  izn = iz
                if (iz==nz) izp = iz
                
                do ispa = 1, nsp_aq

                    row = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa
                    
                    caq_tmp         = maqx(ispa,iz) * maqft(ispa,iz)
                    caq_tmp_prev    = maq(ispa,iz) * maqft_prev(ispa,iz)
                    caq_tmp_p       = maqx(ispa,izp) * maqft(ispa,izp)
                    caq_tmp_n       = maqx(ispa,izn) * maqft(ispa,izn)
                    
                    
                    d_tmp           = daq(ispa)
                    caqdif_tmp_n    = maqx(ispa,izn) * maqft(ispa,izn)
                    caqth_tmp       = maqth(ispa)
                    caqi_tmp        = maqi(ispa)
                    caqsupp_tmp     = maqsupp(ispa,iz) 
                    rxn_ext_tmp     = sum(staq_ext(:,ispa)*rxnext(:,iz))
                    rxn_tmp         = sum(staq(:,ispa)*rxnsld(:,iz))
                    drxndisp_tmp    = sum(staq(:,ispa)*drxnsld_dmaq(:,ispa,iz))
                    
                    if (iz==1 .and. (.not. aq_close) ) caq_tmp_n = caqi_tmp
                    if (iz==1 .and. (.not. aq_diff_close) ) caqdif_tmp_n = caqi_tmp
                        
                    edif_tmp    = 1d3*poro(iz )*sat(iz )*( tora(iz )*d_tmp + disp(iz ) )
                    edif_tmp_p  = 1d3*poro(izp)*sat(izp)*( tora(izp)*d_tmp + disp(izp) )
                    edif_tmp_n  = 1d3*poro(izn)*sat(izn)*( tora(izn)*d_tmp + disp(izn) )

                    amx3(row,row) = ( &
                        & + (poro(iz)*sat(iz)*1d3*1d0*maqft(ispa,iz))/merge(1d0,dt,dt_norm)  &
                        & + (poro(iz)*sat(iz)*1d3*maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa,iz))/merge(1d0,dt,dt_norm)  &
                        & -(0.5d0*(edif_tmp +edif_tmp_p) &
                        &   *merge(0d0,-1d0*maqft(ispa,iz),iz==nz)/( 0.5d0*(dz(iz)+dz(izp)) ) &
                        & -0.5d0*(edif_tmp +edif_tmp_n) &
                        &   *merge(0d0,1d0*maqft(ispa,iz),iz==1 .and. aq_diff_close)/( 0.5d0*(dz(iz)+dz(izn)) ))/dz(iz) &
                        & *merge(dt,1d0,dt_norm) &
                        & -(0.5d0*(edif_tmp +edif_tmp_p) &
                        &   *merge(0d0,-maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa,iz),iz==nz)/( 0.5d0*(dz(iz)+dz(izp)) ) &
                        & -0.5d0*(edif_tmp +edif_tmp_n) &
                        &   *merge(0d0,maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa,iz),iz==1 .and.aq_diff_close) &
                        &   /(0.5d0*(dz(iz)+dz(izn)))  )/dz(iz) &
                        & *merge(dt,1d0,dt_norm) &
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(1d0*maqft(ispa,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                        & -drxndisp_tmp*merge(dt,1d0,dt_norm) &
                        & - sum(staq_ext(:,ispa)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(1.0d0,maqx(ispa,iz),caq_tmp<caqth_tmp*sw_red)

                    ymx3(row) = ( &
                        & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/merge(1d0,dt,dt_norm)  &
                        & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(izp))) &
                        & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caqdif_tmp_n)/(0.5d0*(dz(iz)+dz(izn))))/dz(iz) &
                        & *merge(dt,1d0,dt_norm) &
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*merge(dt,1d0,dt_norm) &
                        & - rxn_tmp*merge(dt,1d0,dt_norm) &
                        & - caqsupp_tmp*merge(dt,1d0,dt_norm) &
                        & - rxn_ext_tmp*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)

                    if (iz/=1) then 
                        amx3(row,row-nsp3) = ( &
                            & -(-0.5d0*(edif_tmp +edif_tmp_n) &
                            &   *(-1d0*maqft(ispa,izn))/(0.5d0*(dz(iz)+dz(izn))))/dz(iz) &
                            & *merge(dt,1d0,dt_norm) &
                            & -(-0.5d0*(edif_tmp +edif_tmp_n) &
                            &   *(-maqx(ispa,izn)*dmaqft_dmaqf(ispa,ispa,izn))/(0.5d0*(dz(iz)+dz(izn))))/dz(iz) &
                            & *merge(dt,1d0,dt_norm) &
                            & + poro(iz)*sat(iz)*1d3*v(iz)*(-1d0*maqft(ispa,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & + poro(iz)*sat(iz)*1d3*v(iz) &
                            &   *(-maqx(ispa,izn)*dmaqft_dmaqf(ispa,ispa,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *maqx(ispa,izn)  &
                            & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                    endif 
                    
                    if (iz/=nz) then 
                        amx3(row,row+nsp3) = ( &
                            & -(0.5d0*(edif_tmp +edif_tmp_p) &
                            &   *(1d0*maqft(ispa,izp))/(0.5d0*(dz(iz)+dz(izp))))/dz(iz) &
                            & *merge(dt,1d0,dt_norm) &
                            & -(0.5d0*(edif_tmp +edif_tmp_p) &
                            &   *(maqx(ispa,izp)*dmaqft_dmaqf(ispa,ispa,izp))/(0.5d0*(dz(iz)+dz(izp))))/dz(iz) &
                            & *merge(dt,1d0,dt_norm) &
                            & ) &
                            & *maqx(ispa,izp) &
                            & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                    endif 
                    
                    if (.not.sld_enforce) then 
                        do isps = 1, nsp_sld
                            col = nsp3*(iz-1)+ isps
                            
                            amx3(row, col) = (     & 
                                ! & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz)) &
                                ! & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)*merge(dt,1d0,dt_norm)  &
                                & - staq(isps,ispa)*drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) &
                                & - sum(staq_ext(:,ispa)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *msldx(isps,iz) &
                                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                        enddo 
                    endif  
                    
                    do ispa2 = 1, nsp_aq
                        col = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa2
                        
                        if (ispa2 == ispa) cycle
                        
                        amx3(row,col) = amx3(row,col) + (     & 
                            & (poro(iz)*sat(iz)*1d3*maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa2,iz))/merge(1d0,dt,dt_norm)  &
                            & -(0.5d0*(edif_tmp +edif_tmp_p) &
                            &   *merge(0d0,-maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa2,iz),iz==nz)/( 0.5d0*(dz(iz)+dz(izp)) ) &
                            & -0.5d0*(edif_tmp +edif_tmp_n) & 
                            &   * merge(0d0,maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa2,iz),iz==1 .and. aq_diff_close) &
                            &   /( 0.5d0*(dz(iz)+dz(izn)) ))/dz(iz) &
                            & *merge(dt,1d0,dt_norm) &
                            & + poro(iz)*sat(iz)*1d3*v(iz)*(maqx(ispa,iz)*dmaqft_dmaqf(ispa,ispa2,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & - sum(staq(:,ispa)*drxnsld_dmaq(:,ispa2,iz))*merge(dt,1d0,dt_norm) &
                            & - sum(staq_ext(:,ispa)*drxnext_dmaq(:,ispa2,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *maqx(ispa2,iz) &
                            & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                            
                        if (iz/=1) then 
                            amx3(row,col-nsp3) = amx3(row,col-nsp3) + ( &
                                & -(-0.5d0*(edif_tmp +edif_tmp_n) &
                                &   *(-maqx(ispa,izn)*dmaqft_dmaqf(ispa,ispa2,izn)) &
                                &   /(0.5d0*(dz(iz)+dz(izn))))/dz(iz) &
                                &   *merge(dt,1d0,dt_norm) &
                                & + poro(iz)*sat(iz)*1d3*v(iz) &
                                &   *(-maqx(ispa,izn)*dmaqft_dmaqf(ispa,ispa2,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *maqx(ispa2,izn)  &
                                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                        endif 
                        
                        if (iz/=nz) then 
                            amx3(row,col+nsp3) = amx3(row,col+nsp3) + ( &
                                & -(0.5d0*(edif_tmp +edif_tmp_p) &
                                &   *(maqx(ispa,izp)*dmaqft_dmaqf(ispa,ispa2,izp)) &
                                &   /(0.5d0*(dz(iz)+dz(izp))))/dz(iz) &
                                &   *merge(dt,1d0,dt_norm) &
                                & ) &
                                & *maqx(ispa2,izp) &
                                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                        endif 
                    
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg
                        
                        amx3(row,col) = amx3(row,col) + (     & 
                            & (poro(iz)*sat(iz)*1d3*maqx(ispa,iz)*dmaqft_dmgas(ispa,ispg,iz))/merge(1d0,dt,dt_norm)  &
                            & -(0.5d0*(edif_tmp +edif_tmp_p) &
                            &   *merge(0d0,-maqx(ispa,iz)*dmaqft_dmgas(ispa,ispg,iz),iz==nz)/( 0.5d0*(dz(iz)+dz(izp)) ) &
                            & -0.5d0*(edif_tmp +edif_tmp_n) &
                            &   * merge(0d0,maqx(ispa,iz)*dmaqft_dmgas(ispa,ispg,iz),iz==1 .and. aq_diff_close) &
                            &   /( 0.5d0*(dz(iz)+dz(izn)) ))/dz(iz) &
                            &   *merge(dt,1d0,dt_norm) &
                            & + poro(iz)*sat(iz)*1d3*v(iz)*(maqx(ispa,iz)*dmaqft_dmgas(ispa,ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & - sum(staq(:,ispa)*drxnsld_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                            & - sum(staq_ext(:,ispa)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *mgasx(ispg,iz) &
                            & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                            
                        if (iz/=1) then 
                            amx3(row,col-nsp3) = amx3(row,col-nsp3) + ( &
                                & -(-0.5d0*(edif_tmp +edif_tmp_n) &
                                &   *(-maqx(ispa,izn)*dmaqft_dmgas(ispa,ispg,izn)) &
                                &   /(0.5d0*(dz(iz)+dz(izn))))/dz(iz) &
                                & *merge(dt,1d0,dt_norm) &
                                & + poro(iz)*sat(iz)*1d3*v(iz) &
                                &   *(-maqx(ispa,izn)*dmaqft_dmgas(ispa,ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *mgasx(ispg,izn)  &
                                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                        endif 
                        
                        if (iz/=nz) then 
                            amx3(row,col+nsp3) = amx3(row,col+nsp3) + ( &
                                & -(0.5d0*(edif_tmp +edif_tmp_p) &
                                &   *(maqx(ispa,izp)*dmaqft_dmgas(ispa,ispg,izp)) &
                                &   /(0.5d0*(dz(iz)+dz(izp))))/dz(iz) &
                                & *merge(dt,1d0,dt_norm) &
                                & ) &
                                & *mgasx(ispg,izp) &
                                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                        endif 
                    enddo 
                    
                    ! attempt to include adsorption 
                    if (ads_ON) then 
                        ! assuming sold conc. is given in mol per bul m3 
                        m_tmp       = maqx(ispa,iz) * maqfads(ispa,iz)
                        mprev_tmp   = maq(ispa,iz) * maqfads_prev(ispa,iz)
                        mp_tmp      = maqx(ispa,izp) * maqfads(ispa,izp)
                        mth_tmp     = caqth_tmp 
                        mi_tmp      = caqi_tmp
                        w_tmp       = w(iz) 
                        wp_tmp      = w(izp)
                        
                        
                        if (iz==nz) then 
                            mp_tmp = maqx(ispa,nz) * maqfads(ispa,nz) ! no gradient  
                            wp_tmp = w_btm 
                        endif 

                        amx3(row,row) = amx3(row,row) + ( &
                            & + 1d0*maqfads(ispa,iz) /merge(1d0,dt,dt_norm)     &
                            & + maqx(ispa,iz)*dmaqfads_dmaqf(ispa,ispa,iz) /merge(1d0,dt,dt_norm)     &
                            & + w_tmp *1d0*maqfads(ispa,iz) /dz(iz)*merge(dt,1d0,dt_norm)    &
                            & + w_tmp *maqx(ispa,iz)*dmaqfads_dmaqf(ispa,ispa,iz) /dz(iz)*merge(dt,1d0,dt_norm)    &
                            & ) &
                            & * merge(1.0d0,maqx(ispa,iz),m_tmp<mth_tmp*sw_red)
                        
                        ymx3(row) = ymx3(row) + ( &
                            & + ( m_tmp - mprev_tmp )/merge(1d0,dt,dt_norm) &
                            & - ( wp_tmp*mp_tmp - w_tmp* m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  &
                            & ) &
                            & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)

                        flx_aq(ispa,itflx,iz) = flx_aq(ispa,itflx,iz) + ( &
                            & ( m_tmp - mprev_tmp )/dt  &
                            & ) 
                        flx_aq(ispa,iadv,iz) = flx_aq(ispa,iadv,iz) + ( &
                            & - ( wp_tmp*mp_tmp - w_tmp* m_tmp)/dz(iz) &
                            & ) 
                            
                        if (iz/=nz) amx3(row,row+nsp3) = amx3(row,row+nsp3) + ( &
                            & + (- 1d0*maqfads(ispa,izp)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & + (- maqx(ispa,izp)*dmaqfads_dmaqf(ispa,ispa,izp)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(1.0d0,maqx(ispa,izp),m_tmp<mth_tmp*sw_red)
                        
                        if (iz==nz) amx3(row,row) = amx3(row,row) + ( &
                            & + (- 1d0*maqfads(ispa,nz)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & + (- maqx(ispa,nz)*dmaqfads_dmaqf(ispa,ispa,nz)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(1.0d0,maqx(ispa,nz),m_tmp<mth_tmp*sw_red)
                        
                        do ispa2 = 1, nsp_aq
                            col = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa2
                        
                            if (ispa2 == ispa) cycle

                            amx3(row,col) = amx3(row,col) + ( &
                                & + maqx(ispa,iz)*dmaqfads_dmaqf(ispa,ispa2,iz) /merge(1d0,dt,dt_norm)     &
                                & + w_tmp *maqx(ispa,iz)*dmaqfads_dmaqf(ispa,ispa2,iz) /dz(iz)*merge(dt,1d0,dt_norm)    &
                                & ) &
                                & *maqx(ispa2,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                            
                            if (iz/=nz) amx3(row,col+nsp3) = amx3(row,col+nsp3) + ( &
                                & + (- maqx(ispa,izp)*dmaqfads_dmaqf(ispa,ispa2,izp)*wp_tmp/dz(iz)) &
                                &   *merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(1.0d0,maqx(ispa2,izp),m_tmp<mth_tmp*sw_red)
                        
                            if (iz==nz) amx3(row,col) = amx3(row,col) + ( &
                                & + (- maqx(ispa,nz)*dmaqfads_dmaqf(ispa,ispa2,nz)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(1.0d0,maqx(ispa2,nz),m_tmp<mth_tmp*sw_red)
                        enddo 
                        
                        do ispg = 1, nsp_gas
                            col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg

                            amx3(row,col) = amx3(row,col) + ( &
                                & + maqx(ispa,iz)*dmaqfads_dmgas(ispa,ispg,iz) /merge(1d0,dt,dt_norm)     &
                                & + w_tmp *maqx(ispa,iz)*dmaqfads_dmgas(ispa,ispg,iz) /dz(iz)*merge(dt,1d0,dt_norm)    &
                                & ) &
                                & *mgasx(ispg,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                            
                            if (iz/=nz) amx3(row,col+nsp3) = amx3(row,col+nsp3) + ( &
                                & + (- maqx(ispa,izp)*dmaqfads_dmgas(ispa,ispg,izp)*wp_tmp/dz(iz)) &
                                &   *merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(1.0d0,mgasx(ispg,izp),m_tmp<mth_tmp*sw_red)
                    
                            if (iz==nz) amx3(row,col) = amx3(row,col) + ( &
                                & + (- maqx(ispa,nz)*dmaqfads_dmgas(ispa,ispg,nz)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(1.0d0,mgasx(ispg,nz),m_tmp<mth_tmp*sw_red)
                        enddo 
                        
                        if (.not.sld_enforce) then 
                            do isps = 1,nsp_sld 
                                
                                ! if (all(maqfads_sld(ispa,isps,iz) == 0d0)) cycle
                                
                                col = nsp3*(iz-1)+ isps

                                amx3(row,col) = amx3(row,col) + ( &
                                    & + maqx(ispa,iz)*dmaqfads_dmsld(ispa,isps,iz) /merge(1d0,dt,dt_norm)     &
                                    & + w_tmp *maqx(ispa,iz)*dmaqfads_dmsld(ispa,isps,iz) /dz(iz)*merge(dt,1d0,dt_norm)    &
                                    & ) &
                                    & *msldx(isps,iz) &
                                    & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                                
                                if (iz/=nz) amx3(row,col+nsp3) = amx3(row,col+nsp3) + ( &
                                    & + (- maqx(ispa,izp)*dmaqfads_dmsld(ispa,isps,izp)*wp_tmp/dz(iz)) &
                                    &   *merge(dt,1d0,dt_norm) &
                                    & ) &
                                    & *merge(1.0d0,msldx(isps,izp),m_tmp<mth_tmp*sw_red)
                        
                                if (iz==nz) amx3(row,col) = amx3(row,col) + ( &
                                    & + (- maqx(ispa,nz)*dmaqfads_dmsld(ispa,isps,nz)*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                                    & ) &
                                    & *merge(1.0d0,msldx(isps,nz),m_tmp<mth_tmp*sw_red)
                                    

                                ! modifications with porosity and dz are made in make_trans subroutine
                                do iiz = 1, nz
                                    col = nsp3*(iiz-1)+isps
                                    if (trans(iiz,iz,isps)==0d0) cycle
                                        
                                    amx3(row,col) = amx3(row,col) &
                                        & - trans(iiz,iz,isps)*maqx(ispa,iiz)*dmaqfads_sld_dmsld(ispa,isps,iiz)* merge(dt,1d0,dt_norm) &
                                        & *merge(0.0d0,msldx(isps,iiz),m_tmp<mth_tmp*sw_red)
                                    ymx3(row) = ymx3(row) &
                                        & - trans(iiz,iz,isps)*maqx(ispa,iiz)*maqfads_sld(ispa,isps,iiz)* merge(dt,1d0,dt_norm) &
                                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                                        
                                    flx_aq(ispa,idif,iz) = flx_aq(ispa,idif,iz) + ( &
                                        & - trans(iiz,iz,isps)*maqx(ispa,iiz)*maqfads_sld(ispa,isps,iiz) &
                                        & )
                                        
                                    col = nsp3*(iiz-1)+ nsp_sld*solve_sld + ispa
                                    
                                    amx3(row,col) = amx3(row,col) &
                                        & - trans(iiz,iz,isps)*maqx(ispa,iiz)*dmaqfads_sld_dmaqf(ispa,isps,ispa,iiz) &
                                        &   * merge(dt,1d0,dt_norm) *merge(0.0d0,maqx(ispa,iiz),m_tmp<mth_tmp*sw_red)  &
                                        & - trans(iiz,iz,isps)*1d0*maqfads_sld(ispa,isps,iiz)* merge(dt,1d0,dt_norm) &
                                        &   *merge(0.0d0,maqx(ispa,iiz),m_tmp<mth_tmp*sw_red)
                                        
                                    do ispa2 = 1, nsp_aq
                                        col = nsp3*(iiz-1)+ nsp_sld*solve_sld + ispa2
                                        
                                        if (ispa2 ==ispa) cycle
                                        
                                        amx3(row,col) = amx3(row,col) &
                                            & -trans(iiz,iz,isps)*maqx(ispa,iiz)*dmaqfads_sld_dmaqf(ispa,isps,ispa2,iiz) &
                                            &   *merge(dt,1d0,dt_norm)*merge(0.0d0,maqx(ispa2,iiz),m_tmp<mth_tmp*sw_red)
                                        
                                    enddo 
                                        
                                    do ispg = 1, nsp_gas
                                        col = nsp3*(iiz-1) + nsp_sld*solve_sld + nsp_aq + ispg
                                        
                                        amx3(row,col) = amx3(row,col) &
                                            & -trans(iiz,iz,isps)*maqx(ispa,iiz)*dmaqfads_sld_dmgas(ispa,isps,ispg,iiz) &
                                            &   *merge(dt,1d0,dt_norm)*merge(0.0d0,mgasx(ispg,iiz),m_tmp<mth_tmp*sw_red)
                                        
                                    enddo 
                                    
                                enddo
                            enddo 
                        endif 
                        
                    endif ! End of if ads_ON 
                            
                    flx_aq(ispa,itflx,iz) = flx_aq(ispa,itflx,iz) + (&
                        & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                        & ) 
                    flx_aq(ispa,iadv,iz) = flx_aq(ispa,iadv,iz) + (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                        & ) 
                    flx_aq(ispa,idif,iz) = flx_aq(ispa,idif,iz) + (&
                        & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(izp))) &
                        & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caqdif_tmp_n)/(0.5d0*(dz(iz)+dz(izn))))/dz(iz) &
                        & ) 
                    flx_aq(ispa,irxn_sld(:),iz) = (& 
                        ! & -staq(:,ispa)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                        ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
                        & - staq(:,ispa)*rxnsld(:,iz) &
                        & ) 
                    flx_aq(ispa,irain,iz) = (&
                        & - caqsupp_tmp &
                        & ) 
                    flx_aq(ispa,irxn_ext(:),iz) = (&
                        & - staq_ext(:,ispa)*rxnext(:,iz) &
                        & ) 
                    flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz))
                    if (isnan(flx_aq(ispa,ires,iz))) then 
                        print *,chraq(ispa),iz,(flx_aq(ispa,iflx,iz),iflx=1,nflx)
                    endif 
                    
                    amx3(row,:) = amx3(row,:)*fact 
                    ymx3(row) = ymx3(row)*fact 
                
                enddo 
                
            end do  ! ==============================
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2 & pO2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! print *,drxngas_dmaq(findloc(chrgas,'pco2',dim=1),findloc(chraq,'ca',dim=1),:)
            
            do iz = 1, nz
                        
                izp = iz+1
                izn = iz-1
                
                if (iz==1)  izn = iz
                if (iz==nz) izp = iz    
                
                do ispg = 1, nsp_gas
                
                    row = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg            
                    
                    pco2n_tmp   = mgasx(ispg,izn)
                    khco2n_tmp  = khgasx(ispg,izn)
                    edifn_tmp   = dgas(ispg,izn)
                    if (iz == 1 .and. (.not. gas_close) ) then 
                        pco2n_tmp   = mgasi(ispg)
                        khco2n_tmp  = khgasi(ispg)
                        edifn_tmp   = dgasi(ispg)
                    endif 

                    amx3(row,row) = ( &
                        & (agasx(ispg,iz) + dagas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,izp))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(izp))) &
                        & +0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,izp)-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(izp))) &
                        & + merge( &
                        &   -0.5d0*(dgasi(ispg)+dgasn(ispg))*(merge(0d0,1d0,gas_close))/(0.5d0*(dz(iz)+dz(izn))) &
                        &   , &
                        & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(merge(0d0,1d0,iz==1 .and. gas_close))/(0.5d0*(dz(iz)+dz(izn))) &
                        & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                        &   ,iz==1 .and. aq_diff_close &
                        &       ) &
                        &       )/dz(iz)  &
                        & *merge(dt,1d0,dt_norm) &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*1d0)/dz(iz)*merge(dt,1d0,dt_norm) &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/dz(iz) *merge(dt,1d0,dt_norm) &
                        & -sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                        & -drxngas_dmgas(ispg,ispg,iz)*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                    
                    ymx3(row) = ( &
                        & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/merge(1d0,dt,dt_norm) &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,izp))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(izp))) &
                        & - merge( &
                        &   0.5d0*(dgasi(ispg)+dgasn(ispg))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                        &   ,0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn)))  &
                        &   ,iz==1 .and. aq_diff_close) &
                        &       )/dz(iz)  &
                        & *merge(dt,1d0,dt_norm) &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*merge(dt,1d0,dt_norm) &
                        & -sum(stgas_ext(:,ispg)*rxnext(:,iz))*merge(dt,1d0,dt_norm) &
                        & -rxngas(ispg,iz)*merge(dt,1d0,dt_norm) &
                        & -mgassupp(ispg,iz)*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                    
                    
                    if (iz/=nz) then 
                        amx3(row,row+nsp3) = ( &
                                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,izp))*(1d0)/(0.5d0*(dz(iz)+dz(izp))) &
                                & + 0.5d0*(ddgas_dmgas(ispg,ispg,izp))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                                &       /(0.5d0*(dz(iz)+dz(izp))))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(0.0d0,mgasx(ispg,izp),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                    endif 
                    
                    if (iz/=1) then 
                        amx3(row,row-nsp3) = ( &
                            & -(- 0.5d0*(dgas(ispg,iz)+dgas(ispg,izn))*(-1d0)/(0.5d0*(dz(iz)+dz(izn))) &
                            & - 0.5d0*(ddgas_dmgas(ispg,ispg,izn))*(mgasx(ispg,iz)-mgasx(ispg,izn)) &
                            &       /(0.5d0*(dz(iz)+dz(izn))))/dz(iz)*merge(dt,1d0,dt_norm)  &
                            & +poro(iz)*sat(iz)*v(iz)*1d3*(-khgasx(ispg,izn)*1d0)/dz(iz)*merge(dt,1d0,dt_norm) &
                            & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg,izn)*mgasx(ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(0.0d0,mgasx(ispg,izn),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                    endif 
                    
                    if (.not.sld_enforce) then 
                        do isps = 1,nsp_sld
                            col = nsp3*(iz-1) + isps 
                            amx3(row,col) = ( &
                                & -drxngas_dmsld(ispg,isps,iz)*merge(dt,1d0,dt_norm) &
                                & -sum(stgas_ext(:,ispg)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                        enddo 
                    endif 
                    
                    do ispa = 1, nsp_aq
                        col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa 
                        amx3(row,col) = ( &
                            & (dagas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) &
                            ! & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                            ! &       /(0.5d0*(dz(iz)+dz(izp))) &
                            ! & - 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) )/dz(iz)  &
                            ! & *merge(dt,1d0,dt_norm) &
                            & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,izp)-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(izp))) &
                            & + merge( &
                            &   0d0 &
                            &   ,-0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                            &   ,iz==1 .and. aq_diff_close &
                            &       ) &
                            &       )/dz(iz)  &
                            & *merge(dt,1d0,dt_norm) &
                            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & -drxngas_dmaq(ispg,ispa,iz)*merge(dt,1d0,dt_norm) &
                            & -sum(stgas_ext(:,ispg)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                        
                        
                        if (iz/=nz) then 
                            amx3(row,col+nsp3) = ( &
                                & -( 0.5d0*(ddgas_dmaq(ispg,ispa,izp))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                                &       /(0.5d0*(dz(iz)+dz(izp))))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(0.0d0,maqx(ispa,izp),mgasx(ispg,iz)<mgasth(ispg)*sw_red)            
                        endif 
                        
                        if (iz/=1) then 
                            amx3(row,col-nsp3) = ( &
                                & -(- 0.5d0*(ddgas_dmaq(ispg,ispa,izn))*(mgasx(ispg,iz)-mgasx(ispg,izn)) &
                                &       /(0.5d0*(dz(iz)+dz(izn))))/dz(iz)*merge(dt,1d0,dt_norm)  &
                                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,izn)*mgasx(ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(0.0d0,maqx(ispa,izn),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                        endif  
                    enddo 
                    
                    do ispg2 = 1, nsp_gas
                        if (ispg == ispg2) cycle
                        col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg2
                        amx3(row,col) = ( &
                            & (dagas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) &
                            & -( 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                            &       /(0.5d0*(dz(iz)+dz(izp))) &
                            ! & - 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,iz)-pco2n_tmp) &
                            ! &       /(0.5d0*(dz(iz)+dz(izn))) &
                            & + merge( &
                            &   0d0 &
                            &   ,-0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                            &   ,iz==1 .and. aq_diff_close &
                            &       ) &
                            & )/dz(iz)*merge(dt,1d0,dt_norm)  &
                            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                            & -drxngas_dmgas(ispg,ispg2,iz)*merge(dt,1d0,dt_norm) &
                            & -sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg2,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(1.0d0,mgasx(ispg2,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                        
                        if (iz/=nz) then 
                            amx3(row,col+nsp3) = ( &
                                & -( 0.5d0*(ddgas_dmgas(ispg,ispg2,izp))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                                &       /(0.5d0*(dz(iz)+dz(izp))))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(0.0d0,mgasx(ispg2,izp),mgasx(ispg,iz)<mgasth(ispg)*sw_red)            
                        endif 
                        
                        if (iz/=1) then 
                            amx3(row,col-nsp3) = ( &
                                & -(- 0.5d0*(ddgas_dmgas(ispg,ispg2,izn))*(mgasx(ispg,iz)-mgasx(ispg,izn)) &
                                &       /(0.5d0*(dz(iz)+dz(izn))))/dz(iz)*merge(dt,1d0,dt_norm)  &
                                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg2,izn)*mgasx(ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *merge(0.0d0,mgasx(ispg2,izn),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                        endif  
                    enddo 
                    
                    if (amx3(row,row)==0d0) then 
                        print *,amx3(row,row),mgasx(ispg,iz)<mgasth(ispg)*sw_red,mgasx(ispg,iz) 
                        print *, &
                        & (agasx(ispg,iz) + dagas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz)) &
                        & ,-( 0.5d0*(dgas(ispg,iz)+dgas(ispg,izp))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(izp))) &
                        & +0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,izp)-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(izp))) &
                        & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(izn))) &
                        & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) )/dz(iz)  &
                        & ,+poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*1d0)/dz(iz) &
                        & ,+poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/dz(iz) &
                        & ,-sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg,iz)) &
                        & ,-drxngas_dmgas(ispg,ispg,iz) 
                    endif 
                    
                    flx_gas(ispg,itflx,iz) = ( &
                        & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/dt &
                        & )         
                    flx_gas(ispg,idif,iz) = ( &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,izp))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(izp))) &
                        & - merge( &
                        &   0.5d0*(dgasi(ispg)+dgasn(ispg))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                        &   ,0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                        &   ,iz==1 .and. aq_diff_close &
                        &   ) &
                        & )/dz(iz)  &
                        & )
                    flx_gas(ispg,iadv,iz) = ( &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz) &
                        & )
                    flx_gas(ispg,irxn_ext(:),iz) = -stgas_ext(:,ispg)*rxnext(:,iz)
                    flx_gas(ispg,irain,iz) = - mgassupp(ispg,iz)
                    flx_gas(ispg,irxn_sld(:),iz) = ( &
                        ! & -stgas(:,ispg)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                        ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0)  &
                        & - stgas(:,ispg)*rxnsld(:,iz) &
                        & ) 
                    flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
                    
                    if (any(isnan(flx_gas(ispg,:,iz)))) then
                        ! print *,flx_gas(ispg,:,iz)
                        print *,'NAN detected in flx_gas'
                    endif 
                    
                    ! amx3(row,:) = amx3(row,:)/alpha(iz)
                    ! ymx3(row) = ymx3(row)/alpha(iz)
                enddo 

            end do 
            
            fact2= maxval(abs(amx3))
            
            amx3 = amx3/fact2
            ymx3 = ymx3/fact2
            
            ymx3=-1.0d0*ymx3

            if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
            ! if (.true.) then 
                print*,'error in mtx'
                print*,'any(isnan(amx3)),any(isnan(ymx3))'
                print*,any(isnan(amx3)),any(isnan(ymx3))

                if (any(isnan(ymx3))) then 
                    do ie = 1,nsp3*(nz)
                        if (isnan(ymx3(ie))) then 
                            print*,'NAN is here...',ie
                        endif
                    enddo
                endif


                if (any(isnan(amx3))) then 
                    do ie = 1,nsp3*(nz)
                        do ie2 = 1,nsp3*(nz)
                            if (isnan(amx3(ie,ie2))) then 
                                print*,'NAN is here...',ie,ie2
                            endif
                        enddo
                    enddo
                endif
                
        #ifdef errmtx_printout
                open(unit=11,file='amx.txt',status = 'replace')
                open(unit=12,file='ymx.txt',status = 'replace')
                do ie = 1,nsp3*(nz)
                    write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                    write(12,*) ymx3(ie)
                enddo 
                close(11)
                close(12) 
        #endif 
                
                flgback = .true.
                ! pause
                exit

                stop
            endif

            call DGESV(nsp3*(Nz),int(1),amx3,nsp3*(Nz),IPIV3,ymx3,nsp3*(Nz),INFO) 

            if (any(isnan(ymx3))) then
                print*,'error in soultion'
                
        #ifdef errmtx_printout
                open(unit=11,file='amx.txt',status = 'replace')
                open(unit=12,file='ymx.txt',status = 'replace')
                do ie = 1,nsp3*(nz)
                    write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                    write(12,*) ymx3(ie)
                enddo 
                close(11)
                close(12)   
        #endif     
                
                flgback = .true.
                ! pause
                exit
                
                
            endif

            do iz = 1, nz
                if (.not.sld_enforce) then 
                    do isps = 1, nsp_sld
                        row = isps + nsp3*(iz-1)

                        if (isnan(ymx3(row))) then 
                            print *,'nan at', iz,z(iz),chrsld(isps)
                            stop
                        endif
                        
                        ! emx3(row) = (1d0-poro(iz))*msldx(isps,iz)*exp(ymx3(row)) -(1d0-poro(iz))*msldx(isps,iz)
                        emx3(row) = msldx(isps,iz)*exp(ymx3(row)) - msldx(isps,iz)

                        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                            msldx(isps,iz) = msldx(isps,iz)*corr
                        else if (ymx3(row) < -threshold) then 
                            msldx(isps,iz) = msldx(isps,iz)/corr
                        else   
                            msldx(isps,iz) = msldx(isps,iz)*exp(ymx3(row))
                        endif
                        
                        if ( msldx(isps,iz)<msldth(isps)) then ! too small trancate value and not be accounted for error 
                            msldx(isps,iz)=msldth(isps)
                            ymx3(row) = 0d0
                        endif
                    enddo 
                endif 
                
                do ispa = 1, nsp_aq
                    row = ispa + nsp_sld*solve_sld + nsp3*(iz-1)

                    if (isnan(ymx3(row))) then 
                        print *,'nan at', iz,z(iz),chraq(ispa)
                        stop
                    endif
                    
                    if (ads_ON) then 
                        emx3(row) = poro(iz)*sat(iz)*1d3*maqft(ispa,iz)*maqx(ispa,iz)*exp(ymx3(row)) &
                            & + maqfads(ispa,iz)*maqx(ispa,iz)*exp(ymx3(row)) &
                            & - poro(iz)*sat(iz)*1d3*maqft(ispa,iz)*maqx(ispa,iz) &
                            & - maqfads(ispa,iz)*maqx(ispa,iz) 
                    else
                        emx3(row) = poro(iz)*sat(iz)*1d3*maqft(ispa,iz)*maqx(ispa,iz)*exp(ymx3(row)) &
                            & - poro(iz)*sat(iz)*1d3*maqft(ispa,iz)*maqx(ispa,iz)
                    endif 
                    
                    ! emx3(row) = emx3(row)*1d3

                    if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                        maqx(ispa,iz) = maqx(ispa,iz)*corr
                    else if (ymx3(row) < -threshold) then 
                        maqx(ispa,iz) = maqx(ispa,iz)/corr
                    else   
                        maqx(ispa,iz) = maqx(ispa,iz)*exp(ymx3(row))
                    endif
                    
                    if (maqx(ispa,iz)<maqth(ispa)) then ! too small trancate value and not be accounted for error 
                        maqx(ispa,iz)=maqth(ispa)
                        ymx3(row) = 0d0
                    endif
                enddo 
                
                do ispg = 1, nsp_gas
                    row = ispg + nsp_aq + nsp_sld*solve_sld + nsp3*(iz-1)

                    if (isnan(ymx3(row))) then 
                        print *,'nan at', iz,z(iz),chrgas(ispg)
                        stop
                    endif
                    
                    emx3(row) =agasx(ispg,iz)* mgasx(ispg,iz)*exp(ymx3(row)) - agasx(ispg,iz)*mgasx(ispg,iz) 

                    if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                        mgasx(ispg,iz) = mgasx(ispg,iz)*corr
                    else if (ymx3(row) < -threshold) then 
                        mgasx(ispg,iz) = mgasx(ispg,iz)/corr
                    else   
                        mgasx(ispg,iz) = mgasx(ispg,iz)*exp(ymx3(row))
                    endif
                    
                    if (mgasx(ispg,iz)<mgasth(ispg)) then ! too small trancate value and not be accounted for error 
                        mgasx(ispg,iz)=mgasth(ispg)
                        ymx3(row) = 0d0
                    endif
                enddo 
                
        #ifdef calcw_full
                row =  nsp3*(iz-1) + nsp3
                if (isnan(ymx3(row))) then 
                    print *,'nan at', iz,z(iz),'w'
                    stop
                endif
                
                emx3(row) = w(iz)*exp(ymx3(row)) - w(iz) 
                emx3(row) = abs(ymx3(row))  
                
                w(iz) = w(iz) + ymx3(row)
                
                ! if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                    ! w(iz) = w(iz)*corr
                ! else if (ymx3(row) < -threshold) then 
                    ! w(iz) = w(iz)/corr
                ! else   
                    ! w(iz) = w(iz)*exp(ymx3(row))
                ! endif
                
                ! if (mgasx(ispg,iz)<mgasth(ispg)) then ! too small trancate value and not be accounted for error 
                    ! mgasx(ispg,iz)=mgasth(ispg)
                    ! ymx3(row) = 0d0
                ! endif
        #endif 

            end do 

            if (fact_tol == 1d0) then 
                error = maxval(exp(abs(ymx3))) - 1.0d0
            else 
                error = maxval((abs(emx3)))
            endif 
            
            if (isnan(error)) error = 1d4

            if (isnan(error).or.info/=0 .or. any(isnan(msldx)) .or. any(isnan(maqx)).or. any(isnan(mgasx))) then 
                error = 1d3
                print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
                print*, 'isnan(error), info/=0,any(isnan(msldx)),any(isnan(maqx)),any(isnan(mgasx))'
                print*,isnan(error),info,any(isnan(msldx)),any(isnan(maqx)),any(isnan(mgasx))
                
        #ifdef errmtx_printout
                open(unit=11,file='amx.txt',status = 'replace')
                open(unit=12,file='ymx.txt',status = 'replace')
                do ie = 1,nsp3*(nz)
                    write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                    write(12,*) ymx3(ie)
                enddo 
                close(11)
                close(12)   
        #endif 
                
                ! dt = dt/10d0
                flgback = .true.
                ! pause
                exit
                
                
                ! stop
            endif

            if (display) then 
                print '(a,E11.3,a,i0,a,E11.3)', 'iteration error = ',error, ', iteration = ',iter,', time step [yr] = ',dt
            endif      
            iter = iter + 1 

            if (iter > iter_Max ) then
                ! dt = dt/1.01d0
                ! dt = dt/10d0
                if (dt==0d0) then 
                    print *, 'dt==0d0; stop'
                
        #ifdef errmtx_printout
                    open(unit=11,file='amx.txt',status = 'replace')
                    open(unit=12,file='ymx.txt',status = 'replace')
                    do ie = 1,nsp3*(nz)
                        write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                        write(12,*) ymx3(ie)
                    enddo 
                    close(11)
                    close(12)      
        #endif 
                    stop
                endif 
                flgback = .true.
                
                exit 
            end if
            
        #ifdef dispiter
                write(chrfmt,'(i0)') nz_disp
                chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
                
                print *
                print *,' [concs] '
                print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
                if (nsp_aq>0) then 
                    print *,' < aq species >'
                    do ispa = 1, nsp_aq
                        print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/nz_disp)
                    enddo 
                endif 
                if (nsp_sld>0) then 
                    print *,' < sld species >'
                    do isps = 1, nsp_sld
                        print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/nz_disp)
                    enddo 
                endif 
                if (nsp_gas>0) then 
                    print *,' < gas species >'
                    do ispg = 1, nsp_gas
                        print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/nz_disp)
                    enddo 
                endif 
                print *
        #endif     

        enddo

        ! just addint flx calculation at the end 

            
        flx_sld = 0d0
        flx_aq = 0d0
        flx_gas = 0d0

        flx_co2sp = 0d0

        ! pH calculation and its derivative wrt aq and gas species

        call calc_pH_v7_4( &
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
            
        if (ph_error) then 
            print *, 'error issued from ph calculation (after main iteration in alsilicate_aq_gas_1D_v3_2)'
            print *, '---> raising flag and return to main' 
            flgback = .true.
            return
        endif 

        ! *** sanity check 
        if (any(isnan(prox)) .or. any(prox<=0d0)) then    
            print *, ' NAN or <=0 H+ conc. (after main iteration)',any(isnan(prox)),any(prox<=0d0)
            print *,prox
            stop
        endif 
            
        ! getting mgasx_loc & maqx_loc
        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maqx,mgasx,maqc,mgasc &
            & ,maqx_loc,mgasx_loc  &! output
            & )

        ! getting maqft_loc and its derivatives
        call get_maqt_all( &
        ! call get_maqt_all_v2( &
            & nz,nsp_aq_all,nsp_gas_all &
            & ,chraq_all,chrgas_all &
            & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
            & ,mgasx_loc,maqx_loc,prox,iosx,tc &
            & ,dmaqft_dpro_loc,dmaqft_dmaqf_loc,dmaqft_dmgas_loc,dmaqft_dios_loc &! output
            & ,maqft_loc  &! output
            & )

        maqft = 0d0
        do ispa=1,nsp_aq
            maqft(ispa,:)=maqft_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
        enddo 

        !!!  for adsorption 
        if (ads_ON) then 
            call get_msldx_all( &
                & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
                & ,chrsld,chrsld_all,chrsld_cnst &
                & ,msldx,msldc &
                & ,msldx_loc  &! output
                & )

            ! call get_maqads_all_v3( &
            call get_maqads_all_v4( &
            ! call get_maqads_all_v4a( &
                & nz,nsp_aq_all,nsp_sld_all &
                & ,chraq_all,chrsld_all &
                & ,keqcec_all,keqiex_all,cec_pH_depend,beta_all &
                & ,msldx_loc,maqx_loc,prox &
                & ,dmaqfads_sld_dpro_loc,dmaqfads_sld_dmaqf_loc,dmaqfads_sld_dmsld_loc &! output
                & ,msldf_loc,maqfads_sld_loc,beta_loc,ads_error  &! output
                & )

            if (ads_error) then 
                print *, 'error issued from adsorption calculation: raising flag and return to main' 
                flgback = .true.
                return
            endif 

            maqfads_sld = 0d0
            do ispa=1,nsp_aq
                do isps=1,nsp_sld
                    maqfads_sld(ispa,isps,:) &
                        & =maqfads_sld_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrsld_all,chrsld(isps),dim=1),:)
                enddo
            enddo 
        else
            maqfads_sld = 0d0
        endif 

        do ispa=1,nsp_aq
            ! maqfads(ispa,:)=maqfads_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
            do iz=1,nz
                maqfads(ispa,iz) = sum(maqfads_sld(ispa,:,iz))
            enddo
        enddo
        ! recalculation of rate constants for mineral reactions

        if (kin_iter) then 

            ksld = 0d0

            do isps =1,nsp_sld 
                call sld_kin( &
                    & nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) &! input
                    & ,nsp_gas_all,chrgas_all,mgasx_loc &! input
                    & ,nsp_aq_all,chraq_all,maqx_loc &! input
                    & ,chrsld(isps),'pro  ' &! input 
                    & ,kin,dkin_dmsp &! output
                    & ) 
                ksld(isps,:) = kin * fkin(isps,:)
            enddo 

        endif 
            
        ! if kin const. is specified in input file 
        if (nsld_kinspc > 0) then 
            do isps_kinspc=1,nsld_kinspc    
                if ( any( chrsld == chrsld_kinspc(isps_kinspc))) then 
                    select case (trim(adjustl(chrsld_kinspc(isps_kinspc))))
                        case('g1','g2','g3') ! for OMs, turn over year needs to be provided [yr]
                            if (kin_sld_spc(isps_kinspc)/=0d0) then  
                                ksld(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = ( &                   
                                    & 1d0/kin_sld_spc(isps_kinspc) &
                                    & ) 
                            else
                                ksld(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = ( &                            
                                    & kin_sld_spc(isps_kinspc) &
                                    & ) 
                            endif 
                        case default ! otherwise, usual rate constant [mol/m2/yr]
                            ksld(findloc(chrsld,chrsld_kinspc(isps_kinspc),dim=1),:) = ( &                            
                                & kin_sld_spc(isps_kinspc) &
                                & ) 
                    end select 
                endif 
            enddo 
        endif 

        ! *** sanity check
        if (any(isnan(ksld)) .or. any(ksld>infinity)) then 
            print *,' *** found insanity in ksld (after main loop): listing below -- '
            do isps=1,nsp_sld
                do iz=1,nz
                    if (isnan(ksld(isps,iz)) .or. ksld(isps,iz)>infinity) print*,chrsld(isps),iz,ksld(isps,iz)
                enddo
            enddo 
            stop
        endif 

        ! saturation state calc. and their derivatives wrt aq and gas species

        omega = 0d0

        do isps =1, nsp_sld
            dummy = 0d0
            call calc_omega_v5( &
                & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst & 
                & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &
                & ,maqx,maqc,mgasx,mgasc,mgasth_all,prox,iosx,tc &
                & ,keqsld_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 &
                & ,staq_all,stgas_all &
                & ,chrsld(isps) &
                & ,domega_dmaq_all,domega_dmgas_all,domega_dpro_loc,domega_dios_loc &! output
                & ,dummy,omega_error &! output
                & )
            omega(isps,:) = dummy
        enddo 

        ! *** reducing saturation ***    
        do isps=1,nsp_sld
            dummy = 0d0
            if (any(chrsld_2 == chrsld(isps))) then  ! chrsld(isps) is included in secondary minerals
                ! cycle
                do iz=1,nz
                    if (omega(isps,iz)>=sat_lim_prec) omega(isps,iz) = sat_lim_prec
                enddo 
            else
                ! omega(isps,:) = dummy
                do iz=1,nz
                    if (omega(isps,iz)>=sat_lim_noprec) omega(isps,iz) = sat_lim_noprec
                enddo 
            endif 
        enddo 

        rxnsld = 0d0
            
        call sld_rxn( &
            & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
            & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle,solmod &! input
            & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas,chrsld &! input
            & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
            & ) 

        ! adding reactions that are not based on dis/prec of minerals
        rxnext = 0d0

        do irxn=1,nrxn_ext
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev_3( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                & ,rg,tempk_0,tc,iosx &!input
                & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
                & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
                & ,chrrxn_ext(irxn),'pro  ' &! input 
                & ,dummy,dummy2,rxnext_error &! output
                & )
            if (rxnext_error) then
                flgback = .true.
                exit 
            endif 
            rxnext(irxn,:) = dummy
        enddo 

        if (.not.sld_enforce)then 
            do iz = 1, nz  !================================
                        
                izp = iz+1
                izn = iz-1
                
                if (iz==1)  izn = iz
                if (iz==nz) izp = iz
                
                do isps = 1, nsp_sld
                    
                    m_tmp           = msldx(isps,iz) 
                    mth_tmp         = msldth(isps) 
                    mi_tmp          = msldi(isps)
                    mp_tmp          = msldx(isps,izp)
                    msupp_tmp       = msldsupp(isps,iz) 
                    rxn_ext_tmp     = sum(stsld_ext(:,isps)*rxnext(:,iz))
                    mprev_tmp       = msld(isps,iz)  
                    w_tmp           = w(iz) 
                    wp_tmp          = w(izp) 
                    sporo_tmp       = 1d0-poro(iz)
                    sporop_tmp      = 1d0-poro(izp) 
                    sporoprev_tmp   = 1d0-poroprev(iz)
                    
                    if (iz==nz) then 
                        mp_tmp      = mi_tmp
                        wp_tmp      = w_btm 
                        sporop_tmp  = 1d0- poroi
                    endif 
                    
                    if (msldunit == 'blk') then 
                        sporo_tmp       = 1d0
                        sporop_tmp      = 1d0
                        sporoprev_tmp   = 1d0
                    endif 
                    
                    ! diffusion terms are filled with transition matrices 
                    ! if (turbo2(isps).or.labs(isps)) then
                        ! do iiz = 1, nz
                            ! if (trans(iiz,iz,isps)==0d0) cycle
                                
                            ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                                ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz) &
                                ! & )
                        ! enddo
                    ! else
                        ! do iiz = 1, nz
                            ! if (trans(iiz,iz,isps)==0d0) cycle
                                
                            ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                                ! & - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                                ! & )
                        ! enddo
                    ! endif
                    
                    do iiz = 1, nz
                        if (trans(iiz,iz,isps)==0d0) cycle
                            
                        flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                            & - trans(iiz,iz,isps)*msldx(isps,iiz) * sporo(iiz) &
                            & )
                    enddo
                    
                    flx_sld(isps,itflx,iz) = ( &
                        & (sporo_tmp*m_tmp - sporoprev_tmp*mprev_tmp)/dt &
                        & )
                    flx_sld(isps,iadv,iz) = ( &
                        & - ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)  &
                        ! & - adf(iz)*up(iz)*( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)  &
                        ! & - adf(iz)*dwn(iz)*( sporo_tmp*w_tmp* m_tmp - sporon_tmp*wn_tmp*mn_tmp )/dz(iz)  &
                        ! & - adf(iz)*cnr(iz)*( sporop_tmp*wp_tmp*mp_tmp - sporon_tmp*wn_tmp*mn_tmp )/dz(iz)  &
                        & )
                    ! flx_sld(isps,irxn_sld(isps),iz) = ( &
                        ! & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                        ! & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                        ! & )
                    flx_sld(isps,irxn_sld(isps),iz) = ( &
                        & + rxnsld(isps,iz) &
                        & )
                    flx_sld(isps,irain,iz) = (&
                        & - msupp_tmp  &
                        & )
                    flx_sld(isps,irxn_ext(:),iz) = (&
                            & - stsld_ext(:,isps)*rxnext(:,iz)  &
                            & )
                    flx_sld(isps,ires,iz) = sum(flx_sld(isps,:,iz))
                    if (isnan(flx_sld(isps,ires,iz))) then 
                        print *,chrsld(isps),iz,(flx_sld(isps,iflx,iz),iflx=1,nflx)
                    endif   
                    
                enddo 
            end do  !================================
        endif 

        do iz = 1, nz
                        
            izp = iz+1
            izn = iz-1
            
            if (iz==1)  izn = iz
            if (iz==nz) izp = iz
            
            do ispa = 1, nsp_aq
                
                caq_tmp         = maqx(ispa,iz) * maqft(ispa,iz)
                caq_tmp_prev    = maq(ispa,iz) * maqft_prev(ispa,iz)
                caq_tmp_p       = maqx(ispa,izp) * maqft(ispa,izp)
                caq_tmp_n       = maqx(ispa,izn) * maqft(ispa,izn)
                
                
                d_tmp           = daq(ispa)
                caqdif_tmp_n    = maqx(ispa,izn) * maqft(ispa,izn)
                caqth_tmp       = maqth(ispa)
                caqi_tmp        = maqi(ispa)
                caqsupp_tmp     = maqsupp(ispa,iz) 
                rxn_ext_tmp     = sum(staq_ext(:,ispa)*rxnext(:,iz))
                rxn_tmp         = sum(staq(:,ispa)*rxnsld(:,iz))
                drxndisp_tmp    = sum(staq(:,ispa)*drxnsld_dmaq(:,ispa,iz))
                
                if (iz==1 .and. (.not. aq_close) ) caq_tmp_n = caqi_tmp
                if (iz==1 .and. (.not. aq_diff_close) ) caqdif_tmp_n = caqi_tmp
                    
                edif_tmp    = 1d3*poro(iz )*sat(iz )*( tora(iz )*d_tmp + disp(iz ) )
                edif_tmp_p  = 1d3*poro(izp)*sat(izp)*( tora(izp)*d_tmp + disp(izp) )
                edif_tmp_n  = 1d3*poro(izn)*sat(izn)*( tora(izn)*d_tmp + disp(izn) )
                
                ! attempt to include adsorption 
                if (ads_ON) then 
                    ! assuming sold conc. is given in mol per bul m3 
                    m_tmp       = maqx(ispa,iz) * maqfads(ispa,iz)
                    mprev_tmp   = maq(ispa,iz) * maqfads_prev(ispa,iz)
                    mp_tmp      = maqx(ispa,izp) * maqfads(ispa,izp)
                    mth_tmp     = caqth_tmp 
                    mi_tmp      = caqi_tmp
                    w_tmp       = w(iz) 
                    wp_tmp      = w(izp)
                    
                    
                    if (iz==nz) then 
                        mp_tmp = maqx(ispa,nz) * maqfads(ispa,nz) ! no gradient  
                        wp_tmp = w_btm 
                    endif 

                    flx_aq(ispa,itflx,iz) = flx_aq(ispa,itflx,iz) + ( &
                        & ( m_tmp - mprev_tmp )/dt  &
                        & ) 
                    flx_aq(ispa,iadv,iz) = flx_aq(ispa,iadv,iz) + ( &
                        & - ( wp_tmp*mp_tmp - w_tmp* m_tmp)/dz(iz) &
                        & ) 
                    
                    if (.not.sld_enforce) then 
                        do isps = 1,nsp_sld 
                            ! modifications with porosity and dz are made in make_trans subroutine
                            do iiz = 1, nz
                                if (trans(iiz,iz,isps)==0d0) cycle
                                    
                                flx_aq(ispa,idif,iz) = flx_aq(ispa,idif,iz) + ( &
                                    & - trans(iiz,iz,isps)*maqx(ispa,iiz)*maqfads_sld(ispa,isps,iiz) &
                                    & )
                            enddo
                        enddo 
                    endif 
                    
                endif 
                        
                flx_aq(ispa,itflx,iz) = flx_aq(ispa,itflx,iz) + (&
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                    & ) 
                flx_aq(ispa,iadv,iz) = flx_aq(ispa,iadv,iz) + (&
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                    & ) 
                flx_aq(ispa,idif,iz) = flx_aq(ispa,idif,iz) + (&
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(izp))) &
                    & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caqdif_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & ) 
                flx_aq(ispa,irxn_sld(:),iz) = (& 
                    ! & -staq(:,ispa)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                    ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
                    & - staq(:,ispa)*rxnsld(:,iz) &
                    & ) 
                flx_aq(ispa,irain,iz) = (&
                    & - caqsupp_tmp &
                    & ) 
                flx_aq(ispa,irxn_ext(:),iz) = (&
                    & - staq_ext(:,ispa)*rxnext(:,iz) &
                    & ) 
                flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz))
                if (isnan(flx_aq(ispa,ires,iz))) then 
                    print *,chraq(ispa),iz,(flx_aq(ispa,iflx,iz),iflx=1,nflx)
                endif 
            
            enddo 
            
        end do  ! ==============================

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2 & pO2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        khgas = 0d0
        khgasx = 0d0
        ! added
        if (new_gassol) then 
            call calc_khgas_all_v2( &
                & nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst &
                & ,chraq_all,chrgas_all,chraq_cnst,chrgas_cnst,chraq,chrgas &
                & ,maq,mgas,maqx,mgasx,maqc,mgasc &
                & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3  &
                & ,pro,prox,ios,iosx,tc &
                & ,khgas_all,khgasx_all,dkhgas_dpro_all,dkhgas_dmaq_all,dkhgas_dmgas_all,dkhgas_dios_all &!output
                & )
                
            do ispg=1,nsp_gas
                khgas(ispg,:)=khgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
                khgasx(ispg,:)=khgasx_all(findloc(chrgas_all,chrgas(ispg),dim=1),:)
            enddo 
        endif 

        dgas = 0d0

        agas = 0d0
        agasx = 0d0

        rxngas = 0d0

        do ispg = 1, nsp_gas
            
            if (.not.new_gassol) then ! to be removed?
                select case (trim(adjustl(chrgas(ispg))))
                    case('pco2')
                        khgas(ispg,:) = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
                        khgasx(ispg,:) = kco2*(1d0+k1/prox + k1*k2/prox/prox)
                    case('po2')
                        khgas(ispg,:) = kho ! previous value; should not change through iterations 
                        khgasx(ispg,:) = kho
                    case('pnh3')
                        khgas(ispg,:) = knh3*(1d0+pro/k1nh3) ! previous value; should not change through iterations 
                        khgasx(ispg,:) = knh3*(1d0+prox/k1nh3)
                    case('pn2o')
                        khgas(ispg,:) = kn2o ! previous value; should not change through iterations 
                        khgasx(ispg,:) = kn2o
                endselect 
            endif 
            
            dgas(ispg,:) = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg)+poro*sat*khgasx(ispg,:)*1d3*(tora*dgasa(ispg)+disp)
            dgasi(ispg) = ucv*1d3*dgasg(ispg) 
            
            agas(ispg,:)= ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khgas(ispg,:)*1d3
            agasx(ispg,:)= ucv*poro*(1.0d0-sat)*1d3+poro*sat*khgasx(ispg,:)*1d3
            
            do isps = 1, nsp_sld
                rxngas(ispg,:) =  rxngas(ispg,:) + (&
                    ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + stgas(isps,ispg)*rxnsld(isps,:) &
                    & )
            enddo 
        enddo 

        ! print *,drxngas_dmaq(findloc(chrgas,'pco2',dim=1),findloc(chraq,'ca',dim=1),:)

        do iz = 1, nz
                    
            izp = iz+1
            izn = iz-1
            
            if (iz==1)  izn = iz
            if (iz==nz) izp = iz    
            
            do ispg = 1, nsp_gas            
                
                pco2n_tmp   = mgasx(ispg,izn)
                khco2n_tmp  = khgasx(ispg,izn)
                edifn_tmp   = dgas(ispg,izn)
                if (iz == 1 .and. (.not. gas_close)) then 
                    pco2n_tmp   = mgasi(ispg)
                    khco2n_tmp  = khgasi(ispg)
                    edifn_tmp   = dgasi(ispg)
                endif 
                
                flx_gas(ispg,itflx,iz) = ( &
                    & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/dt &
                    & )         
                flx_gas(ispg,idif,iz) = ( &
                    & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,izp))*(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(izp))) &
                    & - merge( &
                    &   0.5d0*(dgasi(ispg)+dgasn(ispg))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                    &   ,0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                    &   ,iz==1 .and. aq_diff_close &
                    &   ) &
                    & )/dz(iz)  &
                    & )
                flx_gas(ispg,iadv,iz) = ( &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz) &
                    & )
                flx_gas(ispg,irxn_ext(:),iz) = -stgas_ext(:,ispg)*rxnext(:,iz)
                flx_gas(ispg,irain,iz) = - mgassupp(ispg,iz)
                flx_gas(ispg,irxn_sld(:),iz) = ( &
                    ! & -stgas(:,ispg)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                    ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
                    & - stgas(:,ispg)*rxnsld(:,iz) &
                    & )
                flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
                
                if (any(isnan(flx_gas(ispg,:,iz)))) then
                    ! print *,flx_gas(ispg,:,iz)
                    print *,'NAN detected in flx_gas'
                endif 
            enddo 
            
            if (any(chrgas=='pco2')) then 
                ispg = findloc(chrgas,'pco2',dim=1)
                
                pco2n_tmp = mgasx(ispg,izn)
                proi_tmp = prox(izn)
                if (iz == 1) then 
                    pco2n_tmp = mgasi(ispg)
                    proi_tmp = proi
                endif 
                
                ! gaseous CO2
                
                edifn_tmp = ucv*poro(izn)*(1.0d0-sat(izn))*1d3*torg(izn)*dgasg(ispg)
                if (iz==1) edifn_tmp = dgasi(ispg)
                
                flx_co2sp(1,itflx,iz) = ( &
                    & (ucv*poro(iz)*(1.0d0-sat(Iz))*1d3*mgasx(ispg,iz)-ucv*poroprev(iz)*(1.0d0-sat(Iz))*1d3*mgas(ispg,iz))/dt &
                    & )  
                flx_co2sp(1,idif,iz) = ( &
                    & -( 0.5d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(Iz)*dgasg(ispg) &
                    &       +ucv*poro(izp)*(1.0d0-sat(izp))*1d3*torg(izp)*dgasg(ispg)) &
                    &       *(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(izp))) &
                    & - 0.5d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(Iz)*dgasg(ispg) + edifn_tmp) &
                    &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                    &  )/dz(iz)  &
                    & ) 
                flx_co2sp(1,irxn_ext(:),iz) = -stgas_ext(:,ispg)*rxnext(:,iz)
                flx_co2sp(1,irain,iz) = - mgassupp(ispg,iz)
                flx_co2sp(1,irxn_sld(:),iz) = ( &
                    & - stgas(:,ispg)*rxnsld(:,iz) &
                    & )
                    
                ! dissolved CO2
                
                edifn_tmp = poro(izn)*sat(izn)*kco2*1d3*(tora(izn)*dgasa(ispg)+disp(izn))
                if (iz==1) edifn_tmp = 0d0
                
                flx_co2sp(2,itflx,iz) = ( &
                    & (poro(iz)*sat(iz)*kco2*1d3*mgasx(ispg,iz)-poroprev(iz)*sat(iz)*kco2*1d3*mgas(ispg,iz))/dt &
                    & )  
                flx_co2sp(2,idif,iz) = ( &
                    & -( 0.5d0*(poro(iz)*sat(iz)*kco2*1d3*(tora(iz)*dgasa(ispg)+disp(iz)) &
                    &       +poro(izp)*sat(izp)*kco2*1d3*(tora(izp)*dgasa(ispg)+disp(izp))) &
                    &       *(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(izp))) &
                    & - merge( &
                    &   0d0  &
                    &   ,0.5d0*(poro(iz)*sat(iz)*kco2*1d3*(tora(iz)*dgasa(ispg) +disp(iz)) + edifn_tmp) &
                    &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                    &   ,iz==1 .and.aq_diff_close &
                    &   ) &
                    & )/dz(iz)  &
                    & ) 
                flx_co2sp(2,iadv,iz) = ( &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(kco2*mgasx(ispg,iz)- kco2*pco2n_tmp)/dz(iz) &
                    & )
                    
                ! HCO3-
                
                edifn_tmp = poro(izn)*sat(izn)*kco2*k1/prox(izn)*1d3 &
                    & *(tora(izn)*dgasa(ispg)+disp(izn))
                if (iz==1) edifn_tmp = 0d0
                
                flx_co2sp(3,itflx,iz) = ( &
                    & (poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*mgasx(ispg,iz)-poroprev(iz)*sat(iz)*kco2*k1/pro(iz)*1d3*mgas(ispg,iz))/dt &
                    & )  
                flx_co2sp(3,idif,iz) = ( &
                    & -( 0.5d0*(poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*(tora(iz)*dgasa(ispg)+disp(iz)) &
                    &       +poro(izp)*sat(izp)*kco2*k1/prox(izp)*1d3 &
                    &       *(tora(izp)*dgasa(ispg)+disp(izp)) ) &
                    &       *(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(izp))) &
                    & - merge( &
                    &   0d0 &
                    &   ,0.5d0*(poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*(tora(iz)*dgasa(ispg)+disp(iz)) + edifn_tmp) &
                    &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                    &   ,iz==1 .and.aq_diff_close &
                    &   ) &
                    & )/dz(iz)  &
                    & ) 
                flx_co2sp(3,iadv,iz) = ( &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*( &
                    &       kco2*k1/prox(iz)*mgasx(ispg,iz) &
                    &       - kco2*k1/proi_tmp*pco2n_tmp)/dz(iz) &
                    & )
                    
                ! CO32-
                
                edifn_tmp = poro(izn)*sat(izn)*kco2*k1*k2/prox(izn)**2d0*1d3 &
                    & *(tora(izn)*dgasa(ispg)+disp(izn))
                if (iz==1) edifn_tmp = 0d0
                
                flx_co2sp(4,itflx,iz) = (  &
                    & (poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)**2d0*1d3*mgasx(ispg,iz) &
                    &       -poroprev(iz)*sat(iz)*kco2*k1*k2/pro(iz)**2d0*1d3*mgas(ispg,iz))/dt &
                    & )  
                flx_co2sp(4,idif,iz) = ( &
                    & -( 0.5d0*(poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)**2d0*1d3*(tora(iz)*dgasa(ispg)+disp(iz)) &
                    &       +poro(izp)*sat(izp)*kco2*k1*k2/prox(izp)**2d0*1d3 &
                    &       *(tora(izp)*dgasa(ispg)+disp(izp))) &
                    &       *(mgasx(ispg,izp)-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(izp))) &
                    & - merge( & 
                    &   0d0 &
                    &   ,0.5d0*(poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)**2d0*1d3*(tora(iz)*dgasa(ispg)+disp(iz)) + edifn_tmp) &
                    &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(izn))) &
                    &   ,iz==1 .and.aq_diff_close &
                    &   ) &
                    & )/dz(iz)  &
                    & ) 
                flx_co2sp(4,iadv,iz) = ( &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*( &
                    &       kco2*k1*k2/prox(iz)**2d0*mgasx(ispg,iz) &
                    &       - kco2*k1*k2/proi_tmp**2d0*pco2n_tmp)/dz(iz) &
                    & )
                    
                    
                flx_co2sp(1,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
                flx_co2sp(2,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
                flx_co2sp(3,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
                flx_co2sp(4,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
            endif 
            
        end do 

        ! because still total and so4f are traced in the main subroutine
        ! if (any(chraq == 'so4')) then 
            ! so4f = maqx(findloc(chraq,'so4',dim=1),:)
        ! elseif (any(chraq_cnst == 'so4')) then 
            ! so4f = maqc(findloc(chraq_cnst,'so4',dim=1),:)
        ! endif 
        ! returning maqx as total 
        ! maqx = maqx*maqft
            
        #ifdef dispiter
                
        print *
        print *,' [saturation & pH] '
        if (nsp_sld>0) then 
            print *,' < sld species omega >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        print *,' < pH >'
        print trim(adjustl(chrfmt)), 'ph', (-log10(prox(iz)),iz=1,nz, nz/nz_disp)
        print *

        write(chrfmt,'(i0)') nflx
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'

        print *
        print *,' [fluxes] '
        print trim(adjustl(chrfmt)),'time',(chrflx(iflx),iflx=1,nflx)

        write(chrfmt,'(i0)') nflx
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        if (nsp_aq>0) then 
            print *,' < aq species >'
            do ispa = 1, nsp_aq
                print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
            enddo 
        endif 
        if (nsp_sld>0) then 
            print *,' < sld species >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
            enddo 
        endif 
        if (nsp_gas>0) then 
            print *,' < gas species >'
            do ispg = 1, nsp_gas
                print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
            enddo 
        endif 
        print *
        #endif     

        if (chkflx .and. dt > dt_th) then 
            flx_max_max = 0d0
            do isps = 1, nsp_sld

                flx_max = 0d0
                do iflx = 1, nflx
                    flx_max = max(flx_max,abs(sum(flx_sld(isps,iflx,:)*dz)))
                enddo 
                
                flx_max_max = max(flx_max_max,flx_max)
            enddo 

            do ispa = 1, nsp_aq

                flx_max = 0d0
                do iflx = 1, nflx
                    flx_max = max(flx_max,abs(sum(flx_aq(ispa,iflx,:)*dz)))
                enddo 
                flx_max_max = max(flx_max_max,flx_max)
            enddo 

            do ispg = 1, nsp_gas

                flx_max = 0d0
                do iflx = 1, nflx
                    flx_max = max(flx_max,abs(sum(flx_gas(ispg,iflx,:)*dz)))
                enddo 
                flx_max_max = max(flx_max_max,flx_max)
            enddo 
            
            if (aq_close .and. flx_max_max < flx_max_max_tol) return
            
            if (.not.sld_enforce) then 
                do isps = 1, nsp_sld

                    flx_max = 0d0
                    do iflx = 1, nflx
                        flx_max = max(flx_max,abs(sum(flx_sld(isps,iflx,:)*dz)))
                    enddo 
                    
                    if (flx_max/flx_max_max > flx_max_tol .and.  abs(sum(flx_sld(isps,ires,:)*dz))/flx_max > flx_tol ) then 
                        print *
                        print *, '*** too large error in mass balance of sld phases'
                        print *,'sp          = ',chrsld(isps)
                        print *,'flx_max_tol = ',  flx_max_tol
                        print *,'flx_max_max = ',  flx_max_max
                        print *,'flx_max     = ',  flx_max
                        print *
                        do iflx=1,nflx
                            print *, chrflx(iflx)//'       = ', sum(flx_sld(isps,iflx,:)*dz(:))
                        enddo
                        print *
                        ! pause
                        flgback = .true.
                        return
                    
                        open(unit=11,file='amx.txt',status = 'replace')
                        open(unit=12,file='ymx.txt',status = 'replace')
                        do ie = 1,nsp3*(nz)
                            write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                            write(12,*) ymx3(ie)
                        enddo 
                        close(11)
                        close(12)     
                        
                        stop
                        ! dt = dt/10d0
                    endif 
                enddo 
            endif 

            do ispa = 1, nsp_aq

                flx_max = 0d0
                do iflx = 1, nflx
                    flx_max = max(flx_max,abs(sum(flx_aq(ispa,iflx,:)*dz)))
                enddo 
                
                if (flx_max/flx_max_max > flx_max_tol  .and. abs(sum(flx_aq(ispa,ires,:)*dz))/flx_max > flx_tol ) then 
                    print *
                    print *, '*** too large error in mass balance of aq phases'
                    print *,'sp          = ',chraq(ispa)
                    print *,'flx_max_tol = ',  flx_max_tol
                    print *,'flx_max_max = ',  flx_max_max
                    print *,'flx_max     = ',  flx_max
                    print *
                    do iflx=1,nflx
                        print *, chrflx(iflx)//'       = ', sum(flx_aq(ispa,iflx,:)*dz(:))
                    enddo
                    print *
                    ! pause
                    flgback = .true.
                    return
                
                    open(unit=11,file='amx.txt',status = 'replace')
                    open(unit=12,file='ymx.txt',status = 'replace')
                    do ie = 1,nsp3*(nz)
                        write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                        write(12,*) ymx3(ie)
                    enddo 
                    close(11)
                    close(12)     
                    
                    stop
                    ! dt = dt/10d0
                endif 
            enddo 

            do ispg = 1, nsp_gas

                flx_max = 0d0
                do iflx = 1, nflx
                    flx_max = max(flx_max,abs(sum(flx_gas(ispg,iflx,:)*dz)))
                enddo 
                
                if (flx_max/flx_max_max > flx_max_tol  .and. abs(sum(flx_gas(ispg,ires,:)*dz))/flx_max > flx_tol ) then 
                    print *
                    print *, '*** too large error in mass balance of gas phases'
                    print *,'sp          = ',chrgas(ispg)
                    print *,'flx_max_tol = ',  flx_max_tol
                    print *,'flx_max_max = ',  flx_max_max
                    print *,'flx_max     = ',  flx_max
                    print *
                    do iflx=1,nflx
                        print *, chrflx(iflx)//'       = ', sum(flx_gas(ispg,iflx,:)*dz(:))
                    enddo
                    print *
                    ! pause
                    flgback = .true.
                    return
                
                    open(unit=11,file='amx.txt',status = 'replace')
                    open(unit=12,file='ymx.txt',status = 'replace')
                    do ie = 1,nsp3*(nz)
                        write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                        write(12,*) ymx3(ie)
                    enddo 
                    close(11)
                    close(12)     
                    
                    ! dt = dt/10d0
                    stop
                endif 
            enddo 
        endif 

    endsubroutine alsilicate_aq_gas_1D_v3_2

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    subroutine sld_rxn( &
        & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
        & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle,solmod &! input
        & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas,chrsld &! input
        & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
        & ) 
        implicit none 

        integer,intent(in)::nz,nsp_sld,nsp_aq,nsp_gas
        real(kind=8),intent(in)::msld_seed,dt
        real(kind=8),dimension(nz),intent(in)::poro,sat,dz
        real(kind=8),dimension(nsp_sld,nz),intent(in)::hr
        real(kind=8),dimension(nsp_sld),intent(in)::mv,msldth
        real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
        real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas
        real(kind=8),dimension(nsp_aq),intent(in)::maqth
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maq
        real(kind=8),dimension(nsp_gas),intent(in)::mgasth
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,agas
        real(kind=8),dimension(nsp_sld,nz),intent(in)::ksld,omega,nonprec,msldx,msld,solmod
        real(kind=8),dimension(nsp_sld,nsp_aq,nz),intent(in)::dksld_dmaq,domega_dmaq
        real(kind=8),dimension(nsp_sld,nsp_gas,nz),intent(in)::dksld_dmgas,domega_dmgas
        character(10),dimension(nsp_sld),intent(in)::precstyle
        character(5),dimension(nsp_sld),intent(in)::chrsld
        real(kind=8),dimension(nsp_sld,nz),intent(out)::rxnsld,drxnsld_dmsld
        real(kind=8),dimension(nsp_sld,nsp_aq,nz),intent(out)::drxnsld_dmaq
        real(kind=8),dimension(nsp_sld,nsp_gas,nz),intent(out)::drxnsld_dmgas

        integer ispa,isps,ispg,iz
        real(kind=8),dimension(nsp_sld,nz)::maxdis,maxprec

        real(kind=8)::auth_th = 1d2
        real(kind=8),parameter::infinity = huge(0d0)
            
            

        ! *** sanity check
        if (any(isnan(ksld)) .or. any(ksld>infinity)) then 
            print *,' *** found insanity in ksld (in sld_rxn): listing below -- '
            do isps=1,nsp_sld
                do iz=1,nz
                    if (isnan(ksld(isps,iz)) .or. ksld(isps,iz)>infinity) print*,chrsld(isps),iz,ksld(isps,iz)
                enddo
            enddo 
            stop
        endif 
        ! print *, 'in sld_rxn'
        ! print *, ksld(findloc(chrsld,'kfs',dim=1),:)

            
        rxnsld = 0d0
        drxnsld_dmsld = 0d0
        drxnsld_dmaq = 0d0
        drxnsld_dmgas = 0d0

        maxdis = 1d200
        maxprec = -1d200
        do isps=1,nsp_sld
            maxdis(isps,:) = min(maxdis(isps,:),(msld(isps,:)-msldth(isps))/dt)
            do ispa = 1,nsp_aq
                if (staq(isps,ispa)<0d0) then 
                    maxdis(isps,:) = min(maxdis(isps,:), -1d0/staq(isps,ispa)*poro*sat*1d3*(maq(ispa,:)-maqth(ispa))/dt )
                endif 
            enddo 
            do ispg = 1,nsp_gas
                if (stgas(isps,ispg)<0d0) then 
                    maxdis(isps,:) = min(maxdis(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:)*(mgas(ispg,:)-mgasth(ispg))/dt )
                endif 
            enddo 
            do ispa = 1,nsp_aq
                if (staq(isps,ispa)>0d0) then 
                    maxprec(isps,:) = max(maxprec(isps,:), -1d0/staq(isps,ispa)*poro*sat*1d3*(maq(ispa,:)-maqth(ispa))/dt )
                endif 
            enddo 
            do ispg = 1,nsp_gas
                if (stgas(isps,ispg)>0d0) then 
                    maxprec(isps,:) = max(maxprec(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:)*(mgas(ispg,:)-mgasth(ispg))/dt )
                endif 
            enddo 
        enddo 

        do isps = 1,nsp_sld
            select case(trim(adjustl(precstyle(isps))))
            
                case ('full_lim') 
                    
                    do iz = 1,nz
                        if (1d0-omega(isps,iz) > 0d0) then 
                            rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) 
                            if (rxnsld(isps,iz)> maxdis(isps,iz)) then 
                                rxnsld(isps,iz) = maxdis(isps,iz)
                                drxnsld_dmsld(isps,iz) = 0d0
                                drxnsld_dmaq(isps,:,iz) = 0d0
                                drxnsld_dmgas(isps,:,iz) = 0d0
                            else 
                                drxnsld_dmsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz)) 
                                drxnsld_dmaq(isps,:,iz) = ( &
                                    & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmaq(isps,:,iz)) &
                                    & +dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                                    & )
                                drxnsld_dmgas(isps,:,iz) = ( &
                                    & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmgas(isps,:,iz)) &
                                    & +dksld_dmgas(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                                    & )
                            endif 
                        elseif (1d0-omega(isps,iz) < 0d0) then 
                            if (nonprec(isps,iz)==1d0) then 
                                rxnsld(isps,iz) = 0d0
                                drxnsld_dmsld(isps,iz) = 0d0
                                drxnsld_dmaq(isps,:,iz) = 0d0
                                drxnsld_dmgas(isps,:,iz) = 0d0
                            elseif (nonprec(isps,iz)==0d0) then 
                                rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz))
                                if (rxnsld(isps,iz) < maxprec(isps,iz)) then 
                                    rxnsld(isps,iz) = maxprec(isps,iz)
                                    drxnsld_dmsld(isps,iz) = 0d0
                                    drxnsld_dmaq(isps,:,iz) = 0d0
                                    drxnsld_dmgas(isps,:,iz) = 0d0
                                else
                                    rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz))
                                    drxnsld_dmsld(isps,iz) = 0d0
                                    drxnsld_dmaq(isps,:,iz) = ( &
                                        & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*(-domega_dmaq(isps,:,iz)) &
                                        & +dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz)) &
                                        & )
                                    drxnsld_dmgas(isps,:,iz) = ( &
                                        & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*(-domega_dmgas(isps,:,iz)) &
                                        & +dksld_dmgas(isps,:,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz)) &
                                        & )
                                endif 
                            endif 
                        endif 
                    enddo 
                    
            
                case ('full') 
                    rxnsld(isps,:) = ( &
                        & + min(ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:))/(1d0-poro),maxdis(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                        &  + max(ksld(isps,:)*poro*hr(isps,:)*(1d0-omega(isps,:)),maxprec(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                        & )
                    
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                        & )
                        
                    do ispa = 1,nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            &  + ksld(isps,:)*poro*hr(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            &  + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            & )
                    enddo 
                        
                    do ispg = 1,nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            & /(1d0-poro) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            & /(1d0-poro) &
                            &  + ksld(isps,:)*poro*hr(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            &  + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            & )
                    enddo 
                    
                    
                case('seed')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('decay')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*msldx(isps,:)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('2/3noporo')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6)**(2d0/3d0) &
                        &       *(2d0/3d0)*msldx(isps,:)**(-1d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('2/3')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:)*(mv(isps)*1d-6)**(2d0/3d0) &
                        &       *(2d0/3d0)*msldx(isps,:)**(-1d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('psd_full')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*hr(isps,:)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('emmanuel')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)*solmod(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case default
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    if (any(isnan(rxnsld(isps,:)))) then 
                        print *, 'NAN in rxnsld: ',chrsld(isps)
                        print *, 'ksld(isps,:)',ksld(isps,:)
                        print *, 'hr(isps,:)',hr(isps,:)
                        print *, 'msldx(isps,:)',msldx(isps,:)
                        print *, 'omega(isps,:)',omega(isps,:)
                        stop
                    endif 
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    ! correcting for solid fraction available to porewater 
                    
                    ! rxnsld(isps,:) = rxnsld(isps,:)/(1d0-poro)
                    ! drxnsld_dmsld(isps,:) = drxnsld_dmsld(isps,:)/(1d0-poro)
                    ! drxnsld_dmaq(isps,:,:) = drxnsld_dmaq(isps,:,:)/(1d0-poro)
                    ! drxnsld_dmgas(isps,:,:) = drxnsld_dmgas(isps,:,:)/(1d0-poro)
                    
                    
                    ! attempt to add authigenesis above some threshould for omega
                    ! do iz=1,nz 
                        ! if (nonprec(isps,iz)==0d0 .and. omega(isps,iz) > auth_th) then 
                            ! rxnsld(isps,iz) = rxnsld(isps,iz) + ( &
                                ! & + ksld(isps,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                ! & )
                            
                            ! do ispa = 1, nsp_aq
                                ! drxnsld_dmaq(isps,ispa,iz) = drxnsld_dmaq(isps,ispa,iz) + ( &
                                    ! & + ksld(isps,iz)*poro(iz)*hr(iz)*(-domega_dmaq(isps,ispa,iz)) &
                                    ! & + dksld_dmaq(isps,ispa,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                    ! & )
                            ! enddo 
                            
                            ! do ispg = 1, nsp_gas
                                ! drxnsld_dmgas(isps,ispg,iz) = drxnsld_dmgas(isps,ispg,iz) + ( &
                                    ! & + ksld(isps,iz)*poro(iz)*hr(iz)*(-domega_dmgas(isps,ispg,iz)) &
                                    ! & + dksld_dmgas(isps,ispg,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                    ! & )
                            ! enddo 
                        ! endif 
                    ! enddo 
                    
                    ! print *, 'max-rxnflx', isps, sum (ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*dz)
            
            endselect
        enddo   
    
    endsubroutine sld_rxn

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    subroutine calc_khgas_all_v2( &
        & nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst &
        & ,chraq_all,chrgas_all,chraq_cnst,chrgas_cnst,chraq,chrgas &
        & ,maq,mgas,maqx,mgasx,maqc,mgasc &
        & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3  &
        & ,pro,prox,ios,iosx,tc &
        & ,khgas,khgasx,dkhgas_dpro,dkhgas_dmaq,dkhgas_dmgas,dkhgas_dios &!output
        & )
        implicit none

        ! input 
        integer,intent(in)::nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_gas),intent(in)::chrgas
        real(kind=8),intent(in)::tc
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx,maq
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx,mgas
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nz),intent(in)::pro,prox,ios,iosx
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3
        ! output 
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::khgas,khgasx,dkhgas_dpro,dkhgas_dios
        real(kind=8),dimension(nsp_gas_all,nsp_gas_all,nz),intent(out)::dkhgas_dmgas
        real(kind=8),dimension(nsp_gas_all,nsp_aq_all,nz),intent(out)::dkhgas_dmaq

        ! local 
        real(kind=8),dimension(nsp_aq_all,nz)::maqx_loc,maq_loc
        real(kind=8),dimension(nsp_aq_all,nz)::maqf_loc,maqf_loc_prev
        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc,mgas_loc

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ispg,ispa,ispa_c,ipco2,ipnh3,io2,in2o,ispa_nh3

        real(kind=8) kco2,k1,k2,knh3,k1nh3,kho,kn2o,rspa_nh3
        real(kind=8),dimension(nz)::pnh3,pnh3x

        integer icharge,ic1,ic2
        real(kind=8) rcharge
        real(kind=8),dimension(nz)::gamma_tmp,dgamma_dios_tmp
        real(kind=8),dimension(nz)::fkw,fkeq,fkw_prev,fkeq_prev,dfkw_dios,dfkeq_dios
        real(kind=8),dimension(4,nz)::gamma,dgamma_dios,gamma_prev,dgamma_dios_prev
        real(kind=8),dimension(nsp_aq_all)::base_charge


        ipco2 = findloc(chrgas_all,'pco2',dim=1)
        ipnh3 = findloc(chrgas_all,'pnh3',dim=1)
        io2 = findloc(chrgas_all,'po2',dim=1)
        in2o = findloc(chrgas_all,'pn2o',dim=1)

        kco2 = keqgas_h(ipco2,ieqgas_h0)
        k1 = keqgas_h(ipco2,ieqgas_h1)
        k2 = keqgas_h(ipco2,ieqgas_h2)

        knh3 = keqgas_h(ipnh3,ieqgas_h0)
        k1nh3 = keqgas_h(ipnh3,ieqgas_h1)

        kho = keqgas_h(io2,ieqgas_h0)

        kn2o = keqgas_h(in2o,ieqgas_h0)

        khgas = 0d0
        khgasx = 0d0

        dkhgas_dpro = 0d0
        dkhgas_dios = 0d0
        dkhgas_dmgas = 0d0
        dkhgas_dmaq = 0d0



        do icharge=1,4
            rcharge = 1d0*icharge
            call calc_gamma_davies(  &
                & nz,iosx,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma(icharge,:)=gamma_tmp(:)
            dgamma_dios(icharge,:)=dgamma_dios_tmp(:)
            call calc_gamma_davies(  &
                & nz,ios,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma_prev(icharge,:)=gamma_tmp(:)
            dgamma_dios_prev(icharge,:)=dgamma_dios_tmp(:)
        enddo
            
        call get_base_charge( &
            & nsp_aq_all & 
            & ,chraq_all & 
            & ,base_charge &! output 
            & )

        do ispg = 1, nsp_gas_all
            select case (trim(adjustl(chrgas_all(ispg))))
                case('pco2')
                    ! Kco2: CO2(g) = CO2(a) assume no correction for activity/fugacity 
                    ! K1  : CO2(a) + H2O = HCO3- + H+ <--> K1 = {HCO3-}{H+}/{CO2(a)} <--> K1/gamma/gamma = [HCO3-][H+]/[CO2(a)]
                    ! K2  : HCO3- = CO32- + H+ <--> K2 = {CO32-}{H+}/{HCO3-} <--> K2*gamma/gamma/gamma2 = [CO32-][H+]/[HCO3-]    
                    fkw = 1d0/gamma(1,:)/gamma(1,:) ! H2O = H+ + OH- <--> Kw = {H+}{OH-} <--> Kw/gamma/gamma = [H+][OH-]
                    fkw_prev = 1d0/gamma_prev(1,:)/gamma_prev(1,:) ! H2O = H+ + OH- <--> Kw = {H+}{OH-} <--> Kw/gamma/gamma = [H+][OH-]
                    dfkw_dios = 1d0*(-2d0)*gamma(1,:)**(-3d0)*dgamma_dios(1,:)
                    fkeq = 1d0/gamma(2,:)
                    fkeq_prev = 1d0/gamma_prev(2,:)
                    dfkeq_dios = -1d0/gamma(2,:)**2d0*dgamma_dios(2,:)
                    khgas(ispg,:) = kco2*(1d0+fkw_prev*k1/pro + fkw_prev*fkeq_prev*k1*k2/pro/pro) ! previous value; should not change through iterations 
                    khgasx(ispg,:) = kco2*(1d0+fkw*k1/prox + fkw*fkeq*k1*k2/prox/prox)
                    
                    dkhgas_dpro(ispg,:) = kco2*(fkw*k1*(-1d0)/prox**2d0 + fkw*fkeq*k1*k2*(-2d0)/prox**3d0)
                    dkhgas_dios(ispg,:) = kco2*(dfkw_dios*k1/prox + dfkw_dios*fkeq*k1*k2/prox/prox + fkw*dfkeq_dios*k1*k2/prox/prox)
                    
                    ! obtain previous data 
                    call get_maqgasx_all( &
                        & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
                        & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
                        & ,maq,mgas,maqc,mgasc &
                        & ,maqf_loc_prev,mgas_loc  &! output
                        & )
                        
                    call get_maqgasx_all( &
                        & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
                        & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
                        & ,maqx,mgasx,maqc,mgasc &
                        & ,maqf_loc,mgasx_loc  &! output
                        & )
                        
                    ! account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                    do ispa = 1, nsp_aq_all
                        do ispa_c = 1,2
                            if ( keqaq_c(ispa,ispa_c) > 0d0) then 
                                if (ispa_c == 1) then ! with CO3--
                                    ic1 = nint(abs(base_charge(ispa)))
                                    ic2 = nint(abs(base_charge(ispa)-2d0))
                                    if ( ic1>0 .and. ic2 > 0) then  
                                        fkeq = gamma(ic1,:)*gamma(2,:)/gamma(ic2,:)
                                        fkeq_prev = gamma_prev(ic1,:)*gamma_prev(2,:)/gamma_prev(ic2,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(ic1,:)*gamma(2,:)/gamma(ic2,:) &
                                            & + gamma(ic1,:)*dgamma_dios(2,:)/gamma(ic2,:) &
                                            & + gamma(ic1,:)*gamma(2,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                            & )
                                    elseif ( ic1==0 .and. ic2 > 0) then  
                                        fkeq = gamma(2,:)/gamma(ic2,:)
                                        fkeq_prev = gamma_prev(2,:)/gamma_prev(ic2,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(2,:)/gamma(ic2,:) &
                                            & + gamma(2,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                            & )
                                    elseif ( ic1>0 .and. ic2 == 0) then  
                                        fkeq = gamma(ic1,:)*gamma(2,:)
                                        fkeq_prev = gamma_prev(ic1,:)*gamma_prev(2,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(ic1,:)*gamma(2,:) &
                                            & + gamma(ic1,:)*dgamma_dios(2,:) &
                                            & )
                                    elseif ( ic1==0 .and. ic2 == 0) then  
                                        fkeq = gamma(2,:)
                                        fkeq_prev = gamma_prev(2,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(2,:) &
                                            & )
                                    endif 
                                    khgas(ispg,:) = khgas(ispg,:) + ( &
                                        & + fkeq_prev*keqaq_c(ispa,ispa_c)*maqf_loc_prev(ispa,:)*k1*k2*kco2*pro**(-2d0) &
                                        & )
                                    khgasx(ispg,:) = khgasx(ispg,:) + ( &
                                        & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*prox**(-2d0) &
                                        & )
                                    dkhgas_dpro(ispg,:) = dkhgas_dpro(ispg,:) + ( &
                                        & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*(-2d0)*prox**(-3d0) &
                                        & )
                                    dkhgas_dmaq(ispg,ispa,:) = dkhgas_dmaq(ispg,ispa,:) + ( &
                                        & + fkeq*keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox**(-2d0) &
                                        & *1d0 &
                                        & )
                                    dkhgas_dios(ispg,:) = dkhgas_dios(ispg,:) + ( &
                                        & + dfkeq_dios*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*prox**(-2d0) &
                                        & )
                                elseif (ispa_c == 2) then ! with HCO3-
                                    ic1 = nint(abs(base_charge(ispa)))
                                    ic2 = nint(abs(base_charge(ispa)-1d0))
                                    if ( ic1>0 .and. ic2 > 0) then  
                                        fkeq = gamma(ic1,:)*gamma(2,:)*gamma(1,:)/gamma(ic2,:)
                                        fkeq_prev = gamma_prev(ic1,:)*gamma_prev(2,:)*gamma_prev(1,:)/gamma_prev(ic2,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(ic1,:)*gamma(2,:)*gamma(1,:)/gamma(ic2,:) &
                                            & + gamma(ic1,:)*dgamma_dios(2,:)*gamma(1,:)/gamma(ic2,:) &
                                            & + gamma(ic1,:)*gamma(2,:)*dgamma_dios(1,:)/gamma(ic2,:) &
                                            & + gamma(ic1,:)*gamma(2,:)*gamma(1,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                            & )
                                    elseif ( ic1==0 .and. ic2 > 0) then  
                                        fkeq = gamma(2,:)*gamma(1,:)/gamma(ic2,:)
                                        fkeq_prev = gamma_prev(2,:)*gamma_prev(1,:)/gamma_prev(ic2,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(2,:)*gamma(1,:)/gamma(ic2,:) &
                                            & + gamma(2,:)*dgamma_dios(1,:)/gamma(ic2,:) &
                                            & + gamma(2,:)*gamma(1,:)*(-1d0)/gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                            & )
                                    elseif ( ic1>0 .and. ic2 == 0) then  
                                        fkeq = gamma(ic1,:)*gamma(2,:)*gamma(1,:)
                                        fkeq_prev = gamma_prev(ic1,:)*gamma_prev(2,:)*gamma_prev(1,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(ic1,:)*gamma(2,:)*gamma(1,:) &
                                            & + gamma(ic1,:)*dgamma_dios(2,:)*gamma(1,:) &
                                            & + gamma(ic1,:)*gamma(2,:)*dgamma_dios(1,:) &
                                            & )
                                    elseif ( ic1==0 .and. ic2 == 0) then  
                                        fkeq = gamma(2,:)*gamma(1,:)
                                        fkeq_prev = gamma_prev(2,:)*gamma_prev(1,:)
                                        dfkeq_dios = ( &
                                            & + dgamma_dios(2,:)*gamma(1,:) &
                                            & + gamma(2,:)*dgamma_dios(1,:) &
                                            & )
                                    endif 
                                    khgas(ispg,:) = khgas(ispg,:) + ( &
                                        & + fkeq_prev*keqaq_c(ispa,ispa_c)*maqf_loc_prev(ispa,:)*k1*k2*kco2*pro**(-1d0) & 
                                        & )
                                    khgasx(ispg,:) = khgasx(ispg,:) + ( &
                                        & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*prox**(-1d0) & 
                                        & )
                                    dkhgas_dpro(ispg,:) = dkhgas_dpro(ispg,:) + ( &
                                        & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*(-1d0)*prox**(-2d0) &
                                        & )
                                    dkhgas_dmaq(ispg,ispa,:) = dkhgas_dmaq(ispg,ispa,:) + ( &
                                        & + fkeq*keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox**(-1d0) & 
                                        & *1d0 &
                                        & )
                                    dkhgas_dios(ispg,:) = dkhgas_dios(ispg,:) + ( &
                                        & + dfkeq_dios*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*prox**(-1d0) & 
                                        & )
                                endif 
                            endif 
                        enddo 
                    enddo 
                    
                case('po2')
                    khgas(ispg,:) = kho ! previous value; should not change through iterations 
                    khgasx(ispg,:) = kho

                case('pnh3')
                    khgas(ispg,:) = knh3*(1d0+pro/k1nh3) ! previous value; should not change through iterations 
                    khgasx(ispg,:) = knh3*(1d0+prox/k1nh3)
                    
                    dkhgas_dpro(ispg,:) = knh3*(1d0/k1nh3)
                    
                    ! obtain previous data 
                    call get_maqgasx_all( &
                        & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
                        & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
                        & ,maq,mgas,maqc,mgasc &
                        & ,maqf_loc_prev,mgas_loc  &! output
                        & )
                        
                    call get_maqgasx_all( &
                        & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
                        & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
                        & ,maqx,mgasx,maqc,mgasc &
                        & ,maqf_loc,mgasx_loc  &! output
                        & )
                    
                    pnh3 = mgas_loc(findloc(chrgas_all,'pnh3',dim=1),:)
                    pnh3x= mgasx_loc(findloc(chrgas_all,'pnh3',dim=1),:)
                    
                    ! complex with NH4
                    do ispa = 1, nsp_aq_all
                        do ispa_nh3 = 1,2
                            rspa_nh3 = real(ispa_nh3,kind=8)
                            if ( keqaq_nh3(ispa,ispa_nh3) > 0d0) then 
                                ic1 = nint(abs(base_charge(ispa)))
                                ic2 = nint(abs(base_charge(ispa)+rspa_nh3))
                                if ( ic1>0 .and. ic2 > 0) then  
                                    fkeq = gamma(ic1,:)*gamma(1,:)**rspa_nh3/gamma(ic2,:)
                                    fkeq_prev = gamma_prev(ic1,:)*gamma_prev(1,:)**rspa_nh3/gamma_prev(ic2,:)
                                    dfkeq_dios = ( &
                                        & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_nh3/gamma(ic2,:) &
                                        & + gamma(ic1,:)*rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                                        &   /gamma(ic2,:) &
                                        & + gamma(ic1,:)*gamma(1,:)**rspa_nh3*(-1d0) &
                                        &   /gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                        & )
                                elseif ( ic1==0 .and. ic2 > 0) then  
                                    fkeq = gamma(1,:)**rspa_nh3/gamma(ic2,:)
                                    fkeq_prev = gamma_prev(1,:)**rspa_nh3/gamma_prev(ic2,:)
                                    dfkeq_dios = ( &
                                        & + rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                                        &   /gamma(ic2,:) &
                                        & + gamma(1,:)**rspa_nh3*(-1d0) &
                                        &   /gamma(ic2,:)**2d0*dgamma_dios(ic2,:) &
                                        & )
                                elseif ( ic1>0 .and. ic2 == 0) then  
                                    fkeq = gamma(ic1,:)*gamma(1,:)**rspa_nh3
                                    fkeq_prev = gamma_prev(ic1,:)*gamma_prev(1,:)**rspa_nh3
                                    dfkeq_dios = ( &
                                        & + dgamma_dios(ic1,:)*gamma(1,:)**rspa_nh3 &
                                        & + gamma(ic1,:)*rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                                        & )
                                elseif ( ic1==0 .and. ic2 == 0) then  
                                    fkeq = gamma(1,:)**rspa_nh3
                                    fkeq_prev = gamma_prev(1,:)**rspa_nh3
                                    dfkeq_dios = ( &
                                        & + rspa_nh3*gamma(1,:)**(rspa_nh3-1d0)*dgamma_dios(1,:) &
                                        & )
                                endif 
                                
                                khgas(ispg,:) = khgas(ispg,:) + ( &
                                    & + fkeq_prev*keqaq_nh3(ispa,ispa_nh3) &
                                    &       *maqf_loc_prev(ispa,:)*(knh3/k1nh3*pro)**rspa_nh3*pnh3**(rspa_nh3-1d0) &
                                    & )
                                khgasx(ispg,:) = khgasx(ispg,:) + ( &
                                    & + fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(knh3/k1nh3*prox)**rspa_nh3*pnh3x**(rspa_nh3-1d0) &
                                    & )
                                dkhgas_dpro(ispg,:) = dkhgas_dpro(ispg,:) + ( &
                                    & + fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(knh3/k1nh3)**rspa_nh3*pnh3x**(rspa_nh3-1d0) &
                                    & *rspa_nh3*rspa_nh3**(rspa_nh3-1d0) &
                                    & )
                                dkhgas_dmaq(ispg,ispa,:) = dkhgas_dmaq(ispg,ispa,:) + ( &
                                    & + fkeq*keqaq_nh3(ispa,ispa_nh3)*1d0*(knh3/k1nh3*prox)**rspa_nh3*pnh3x**(rspa_nh3-1d0) &
                                    & )
                                dkhgas_dmgas(ispg,ipnh3,:) = dkhgas_dmgas(ispg,ipnh3,:) + ( &
                                    & + fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(knh3/k1nh3*prox)**rspa_nh3 &
                                    & *(rspa_nh3-1d0)*pnh3x**(rspa_nh3-2d0) &
                                    & )
                                dkhgas_dios(ispg,:) = dkhgas_dios(ispg,:) + ( &
                                    & + dfkeq_dios*keqaq_nh3(ispa,ispa_nh3) &
                                    &       *maqf_loc(ispa,:)*(knh3/k1nh3*prox)**rspa_nh3*pnh3x**(rspa_nh3-1d0) &
                                    & )
                            endif 
                        enddo 
                    enddo 

                case('pn2o')
                    khgas(ispg,:) = kn2o ! previous value; should not change through iterations 
                    khgasx(ispg,:) = kn2o
            endselect 

        enddo 


    endsubroutine calc_khgas_all_v2

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    function rough_f(ref_dummy,n_dummy,r_dummy)
        implicit none
        integer n_dummy
        real(kind=8),dimension(n_dummy):: rough_f,r_dummy
        character(10) ref_dummy
        selectcase(trim(adjustl(ref_dummy)))
            case('NSB07')       ! Navarre-Sitchler and Brantley (2007)
                rough_f = 10d0**3.3d0*r_dummy**0.33d0
            case('BM00')        ! Brantley and Mellott (2000)
                rough_f = 10d0**0.7d0*r_dummy**(-0.1d0)
            case('Letal21')     ! Lewis et al. (2021) 
                ! rough_f = 10d0**(   154.25d0 * exp( 1.0219d0 * log10( r_dummy ) ) ) ! (assuming sphere)
                ! rough_f = 10d0**(   113.41d0 * exp( 1.0219d0 * log10( r_dummy ) ) ) ! (assuming cube)
                ! rough_f = 10d0**(   max( 2.02d0*log10( r_dummy ) + 10.734d0, 1d0 ) )  ! (assuming sphere)
                rough_f = 10d0**(   max( 2.02d0*log10( r_dummy ) + 10.126d0, 1d0 ) )  ! (assuming cube)
            case('smooth')
                rough_f = 1d0
            case default 
                print*, '*** error in rough_f --> stop'
                stop
        endselect
    endfunction rough_f
    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

    #ifdef no_intr_findloc
        !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
        function findloc(chrlist_in,chrspecific,dim)
            implicit none
            character(*),intent(in)::chrlist_in(:),chrspecific
            integer,intent(in)::dim
            integer findloc,i

            findloc = 0
            do i=1, size(chrlist_in,dim=dim)
                if (trim(adjustl(chrspecific)) == trim(adjustl(chrlist_in(i)))) then
                    findloc = i
                    return
                endif 

            enddo 

        endfunction findloc
        !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    #endif 

endprogram weathering