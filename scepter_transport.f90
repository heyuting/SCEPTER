!**************************************************************************************
! Module: scepter_transport
! Purpose: Calculate transport of aqueous and gaseous species in a 1D soil column
!**************************************************************************************

module scepter_transport
    use scepter_constants    ! Physical and chemical constants
    use scepter_eq_ph        ! pH equilibrium calculations
    use scepter_sld_kin ! Solid dissolution/precipitation kinetics
    use scepter_calc_rxn_ext ! Calculate reaction rates for external reactions
    use scepter_calc_khgas   ! Gas-aqueous phase equilibrium calculations
    use scepter_concentration ! Concentration calculations  
    use scepter_thermodynamics ! Thermodynamic calculations
    use scepter_kinetics ! Kinetic reaction calculations
    use scepter_findloc ! Find location of a value in an array
    implicit none
    private
    public :: alsilicate_aq_gas_1D_v3_2

    contains
    !--------------------------------------------------------------------------------------
    ! Subroutine: alsilicate_aq_gas_1D_v3_2
    ! Purpose: Calculate coupled transport and reactions of aqueous, gaseous, and solid species
    !          in a 1D soil column
    !-------------------------------------------------------------------------------------- 
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
        real(kind=8),dimension(nsp_sld,nz)::domega_dpro,dmsld,dksld_dpro,drxnsld_dmsld &
            & ,dksld_dso4f,domega_dso4f,dksld_dios,domega_dios
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
        real(kind=8),dimension(nsp_gas,nz)::khgasx,khgas,dgas,agasx,agas,rxngas,dkhgas_dpro,dprodmgas & 
            & ,dmgas,dso4fdmgas,dkhgas_dso4f,mgasx_save,dmgasx,dkhgas_dios
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

        character(5), dimension(nflx), intent(in) :: chrflx

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

        integer ieqaq_co3,ieqaq_hco3
        data ieqaq_co3,ieqaq_hco3/1,2/

        integer ieqaq_so4,ieqaq_so42
        data ieqaq_so4,ieqaq_so42/1,2/

        integer, intent(in) :: nrxn_ext_all

        character(5), dimension(nrxn_ext_all), intent(in) :: chrrxn_ext_all

        real(kind=8), dimension(nsp_gas_all), intent(in) :: mgasth_all
        real(kind=8), dimension(nsp_aq_all), intent(in) :: maqth_all
        real(kind=8), dimension(nrxn_ext_all,nz), intent(in) :: krxn1_ext_all, krxn2_ext_all

        real(kind=8), dimension(4,nflx,nz), intent(out) :: flx_co2sp

        integer iz,row,ie,ie2,iflx,isps,ispa,ispg,ispa2,ispg2,col,irxn,isps2,iiz,isps_kinspc,row_w,col_w
        integer izp,izn
        integer::itflx,iadv,idif,irain,ires
        integer::ph_iter,ph_iter2
        data itflx,iadv,idif,irain/1,2,3,4/

        integer, dimension(nsp_sld)::irxn_sld 
        integer, dimension(nrxn_ext)::irxn_ext 

        real(kind=8), dimension(nsp_sld,nz), intent(in) :: hr

        real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp &
            & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,omega_tmp_th,rxn_ext_tmp &
            & ,edif_tmp,edif_tmp_n,edif_tmp_p,khco2n_tmp,pco2n_tmp,edifn_tmp,caqsupp_tmp,kco2,k1,k2,kho,sw_red &
            & ,flx_max,flx_max_max,proi_tmp,knh3,k1nh3,kn2o,wp_tmp,w_tmp,sporo_tmp,sporop_tmp,sporoprev_tmp  &
            & ,mn_tmp,wn_tmp,sporon_tmp,caqdif_tmp_n

        real(kind=8), parameter :: infinity = huge(0d0)
        real(kind=8), parameter :: fact = 1d-3
        real(kind=8), parameter :: dconc = 1d-14
        real(kind=8), parameter :: maxfact = 1d200
        ! real(kind=8), parameter :: threshold = log(maxfact)
        real(kind=8), parameter :: threshold = 10d0
        ! real(kind=8), parameter :: threshold = 3d0
        ! real(kind=8), parameter :: corr = 1.5d0
        real(kind=8), parameter :: corr = exp(threshold)

        real(kind=8), dimension(nz)::dummy,dummy2,dummy3,kin,dkin_dmsp,dumtest,sporo,prox_save,iosx_save

        logical print_cb,ph_error,omega_error,rxnext_error,ads_error
        character(500) print_loc
        character(20) chrfmt

        integer, parameter :: iter_max = 50
        ! integer, parameter :: iter_max = 300

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
        logical, intent(in) :: ads_ON != .true.
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
        logical, intent(in) :: sld_enforce != .true.

        ! logical::aq_close = .false.
        logical, intent(in) :: aq_close != .true.
        logical, intent(in) :: act_ON 

        character(10), dimension(nsp_sld), intent(in) :: precstyle 
        real(kind=8), dimension(nsp_sld,nz), intent(in) :: solmod, fkin ! factor to modify solubility used only to implement rxn rate law as defined by Emmanuel and Ague, 2011
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
        !-----------------------------------------------

        if (aq_close) chkflx = .false.

        !! added to enable closed gas system 
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

        ! ==============================================
        ! Main iteration loop
        ! ==============================================    
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
            
            !!!  for adsorption 
            if (ads_ON) then 
                call get_msldx_all( &
                    & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
                    & ,chrsld,chrsld_all,chrsld_cnst &
                    & ,msldx,msldc &
                    & ,msldx_loc  &! output
                    & )

                call get_maqads_all_v4( &
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
                            & = dmaqfads_sld_dmsld_loc(findloc(chraq_all,chraq(ispa),dim=1), &
                            & findloc(chrsld_all,chrsld(isps),dim=1),:) 
                        
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
                    
                        domega_dmgas(isps,ispg,:) = domega_dmgas_all(findloc(chrgas_all,chrgas(ispg),dim=1),:) &
                            & + (+ domega_dpro(isps,:)*dprodmgas(ispg,:) &
                            & + domega_dios(isps,:)*diosdmgas(ispg,:) )
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
                            & + dkhgas_dmaq_all(findloc(chrgas_all,chrgas(ispg),dim=1) & 
                            & ,findloc(chraq_all,chraq(ispa),dim=1),:) &
                            & + dkhgas_dpro(ispg,:)*dprodmaq(ispa,:) &
                            & + dkhgas_dios(ispg,:)*diosdmaq(ispa,:) &
                            & )
                    enddo 
                    do ispg2=1,nsp_gas
                        dkhgas_dmgas(ispg,ispg2,:)= ( &
                            & + dkhgas_dmgas_all(findloc(chrgas_all,chrgas(ispg),dim=1), &
                            & findloc(chrgas_all,chrgas(ispg2),dim=1),:) &
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
                    ! ==============================================
                    ! Loop over solid species to construct matrix
                    ! ==============================================
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
                        
                        ! Handle boundary conditions
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

                        ! Construct diagonal elements
                        amx3(row,row) = ( &
                            & 1d0 *  sporo_tmp /merge(1d0,dt,dt_norm)     &
                            ! & + adf(iz)*up(iz)*sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    &
                            ! & - adf(iz)*dwn(iz)*sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    &
                            & + sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    &
                            & + drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) &
                            & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & * merge(1.0d0,m_tmp,m_tmp<mth_tmp*sw_red)

                        ! Construct RHS vector
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
                            
                        ! Construct off-diagonal elements for next layer
                        if (iz/=nz) then
                            amx3(row,row+nsp3) = ( &
                            & (- sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            & ) &
                            & *merge(1.0d0,mp_tmp,m_tmp<mth_tmp*sw_red)
                        endif
                        ! if (iz/=1) amx3(row,row-nsp3) = ( &
                            ! & (+ adf(iz)*dwn(iz)* sporon_tmp*wn_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & +(+ adf(iz)*cnr(iz)* sporon_tmp*wn_tmp/dz(iz))*merge(dt,1d0,dt_norm) &
                            ! & ) &
                            ! & *merge(1.0d0,mn_tmp,m_tmp<mth_tmp*sw_red)
                        
                        ! Aqueous species coupling
                        do ispa = 1, nsp_aq
                            col = nsp3*(iz-1) + nsp_sld + ispa
                            
                            amx3(row,col ) = ( &
                                & + drxnsld_dmaq(isps,ispa,iz)*merge(dt,1d0,dt_norm) &
                                & - sum(stsld_ext(:,isps)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *maqx(ispa,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        enddo 
                        ! Gas species coupling
                        do ispg = 1, nsp_gas 
                            col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                            amx3(row,col) = ( &
                                & + drxnsld_dmgas(isps,ispg,iz)*merge(dt,1d0,dt_norm) &
                                & - sum(stsld_ext(:,isps)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                                & ) &
                                & *mgasx(ispg,iz) &
                                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        enddo 
                        ! Other solid species coupling
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
                        ! Calculate flux components
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
                        ! Total residual flux
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
            

            do iz = 1, nz   ! ==============================
                        
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
                            & + poro(iz)*sat(iz)*1d3*v(iz)*(maqx(ispa,iz) &
                            & *dmaqft_dmaqf(ispa,ispa2,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
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
                            & + poro(iz)*sat(iz)*1d3*v(iz)*(maqx(ispa,iz) &
                            & *dmaqft_dmgas(ispa,ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
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
                                & + maqx(ispa,iz)*dmaqfads_dmaqf(ispa,ispa2,iz) /merge(1d0,dt,dt_norm) &
                                & + w_tmp *maqx(ispa,iz)*dmaqfads_dmaqf(ispa,ispa2,iz) /dz(iz)*merge(dt,1d0,dt_norm) &
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
                                        & - trans(iiz,iz,isps)*maqx(ispa,iiz)& 
                                        & *dmaqfads_sld_dmsld(ispa,isps,iiz)* merge(dt,1d0,dt_norm) &
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
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp) &
                        & /dz(iz)*merge(dt,1d0,dt_norm) &
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
                            & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg,izn) &
                            & *mgasx(ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
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
                                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,izn) &
                                & *mgasx(ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
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
                            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg2,iz) &
                            & *mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
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
                                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg2,izn) &
                                & *mgasx(ispg,izn))/dz(iz)*merge(dt,1d0,dt_norm) &
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
        ! ==============================================
        ! End of main iteration loop
        ! ==============================================

        ! just adding flx calculation at the end of the iteration loop
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
            
            call get_maqads_all_v4( &
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

        do iz = 1, nz   ! ==============================
                        
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
                    & (poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*mgasx(ispg,iz) &
                    & -poroprev(iz)*sat(iz)*kco2*k1/pro(iz)*1d3*mgas(ispg,iz))/dt &
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
            
        ! ==============================================
        ! Display iteration results for debugging/monitoring
        ! This section prints out saturation, pH, and flux information
        !==============================================
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

end module scepter_transport