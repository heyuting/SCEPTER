!************************************************************************
! Module: scepter_calc_rxn_ext
! Purpose: Calculate reaction rates for external reactions
!************************************************************************
module scepter_calc_rxn_ext
    use scepter_constants
    use scepter_variables
    use scepter_findloc
    use scepter_concentration
    implicit none
    private
    public :: calc_rxn_ext_dev_3

    contains

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
                            & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                            &       *k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                            & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                            &       *k1fe2hco3*k1*k2*kco2*pco2x/prox &
                            & )*po2x &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg) &
                                    &       *k1fe2*(-1d0)/prox**2d0 &
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
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2hco3*k1*k2*kco2*pco2x/prox &
                                    & )*1d0 &
                                    & )
                            case('pco2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2co3*k1*k2*kco2*1d0/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2hco3*k1*k2*kco2*1d0/prox &
                                    & )*po2x &
                                    & )
                            case('fe2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*1d0*( &
                                    & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg  ) &
                                    &       *k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2hco3*k1*k2*kco2*pco2x/prox &
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
endmodule scepter_calc_rxn_ext