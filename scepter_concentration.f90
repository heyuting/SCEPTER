module scepter_concentration
    use scepter_constants
    use scepter_thermodynamics
    use scepter_findloc

    implicit none

    public :: get_base_charge, get_mgasx_all, get_msldx_all, get_maqgasx_all, get_maqt_all, calc_omega_v5

    contains
    subroutine get_base_charge( &
        & nsp_aq_all & 
        & ,chraq_all & 
        & ,base_charge &! output 
        & )
        implicit none
        integer,intent(in)::nsp_aq_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        real(kind=8),dimension(nsp_aq_all),intent(out)::base_charge

        integer ispa

        do ispa = 1, nsp_aq_all
            selectcase(trim(adjustl(chraq_all(ispa))))
                ! case('so4','oxa')
                case('so4')
                    base_charge(ispa) = -2d0
                case('no3','oxa','cl','ac','mes','glp')
                    base_charge(ispa) = -1d0
                case('si','im','tea')
                    base_charge(ispa) = 0d0
                case('na','k')
                    base_charge(ispa) = 1d0
                case('fe2','mg','ca')
                    base_charge(ispa) = 2d0
                case('fe3','al')
                    base_charge(ispa) = 3d0
                case default 
                    print*,'error in charge assignment'
                    stop
            endselect 
        enddo
    endsubroutine get_base_charge

    subroutine get_mgasx_all( &
        & nz,nsp_gas_all,nsp_gas,nsp_gas_cnst &
        & ,chrgas,chrgas_all,chrgas_cnst &
        & ,mgasx,mgasc &
        & ,mgasx_loc  &! output
        & )
        implicit none

        integer,intent(in)::nz,nsp_gas_all,nsp_gas,nsp_gas_cnst
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc

        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::mgasx_loc

        integer ispg

        mgasx_loc = 0d0

        do ispg = 1, nsp_gas_all
            if (any(chrgas==chrgas_all(ispg))) then 
                mgasx_loc(ispg,:) =  mgasx(findloc(chrgas,chrgas_all(ispg),dim=1),:)
            elseif (any(chrgas_cnst==chrgas_all(ispg))) then 
                mgasx_loc(ispg,:) =  mgasc(findloc(chrgas_cnst,chrgas_all(ispg),dim=1),:)
            endif 
        enddo 

    endsubroutine get_mgasx_all

    subroutine get_msldx_all( &
        & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
        & ,chrsld,chrsld_all,chrsld_cnst &
        & ,msldx,msldc &
        & ,msldx_loc  &! output
        & )
        implicit none

        integer,intent(in)::nz,nsp_sld_all,nsp_sld,nsp_sld_cnst
        character(5),dimension(nsp_sld),intent(in)::chrsld
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst
        real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx
        real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc

        real(kind=8),dimension(nsp_sld_all,nz),intent(out)::msldx_loc

        integer isps

        msldx_loc = 0d0

        do isps = 1, nsp_sld_all
            if (any(chrsld==chrsld_all(isps))) then 
                msldx_loc(isps,:) =  msldx(findloc(chrsld,chrsld_all(isps),dim=1),:)
            elseif (any(chrsld_cnst==chrsld_all(isps))) then 
                msldx_loc(isps,:) =  msldc(findloc(chrsld_cnst,chrsld_all(isps),dim=1),:)
            endif 
        enddo 


    endsubroutine get_msldx_all

    subroutine get_maqgasx_all( &
        & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
        & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
        & ,maqx,mgasx,maqc,mgasc &
        & ,maqx_loc,mgasx_loc  &! output
        & )
        implicit none

        integer,intent(in)::nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc

        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::maqx_loc
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::mgasx_loc

        integer ispa,ispg

        maqx_loc = 0d0
        mgasx_loc = 0d0

        do ispa = 1, nsp_aq_all
            if (any(chraq==chraq_all(ispa))) then 
                maqx_loc(ispa,:) =  maqx(findloc(chraq,chraq_all(ispa),dim=1),:)
            elseif (any(chraq_cnst==chraq_all(ispa))) then 
                maqx_loc(ispa,:) =  maqc(findloc(chraq_cnst,chraq_all(ispa),dim=1),:)
            endif 
        enddo 

        do ispg = 1, nsp_gas_all
            if (any(chrgas==chrgas_all(ispg))) then 
                mgasx_loc(ispg,:) =  mgasx(findloc(chrgas,chrgas_all(ispg),dim=1),:)
            elseif (any(chrgas_cnst==chrgas_all(ispg))) then 
                mgasx_loc(ispg,:) =  mgasc(findloc(chrgas_cnst,chrgas_all(ispg),dim=1),:)
            endif 
        enddo 


    endsubroutine get_maqgasx_all

    subroutine get_maqt_all( &
        & nz,nsp_aq_all,nsp_gas_all &
        & ,chraq_all,chrgas_all &
        & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
        & ,mgasx_loc,maqf_loc,prox,iosx,tc &
        & ,dmaqft_dpro,dmaqft_dmaqf,dmaqft_dmgas,dmaqft_dios &! output
        & ,maqft_loc  &! output
        & )
        ! calculating ratio of total dissolved species relative to maqf_loc
        implicit none
        integer,intent(in)::nz,nsp_aq_all,nsp_gas_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        real(kind=8),intent(in)::tc
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl
        real(kind=8),dimension(nsp_aq_all,nz),intent(in)::maqf_loc
        real(kind=8),dimension(nsp_gas_all,nz),intent(in)::mgasx_loc
        real(kind=8),dimension(nz),intent(in)::prox,iosx

        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::maqft_loc
        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::dmaqft_dpro
        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::dmaqft_dios
        real(kind=8),dimension(nsp_aq_all,nsp_aq_all,nz),intent(out)::dmaqft_dmaqf
        real(kind=8),dimension(nsp_aq_all,nsp_gas_all,nz),intent(out)::dmaqft_dmgas

        integer ispa,ispa_h,ispa_c,ispa_s,ispa_no3,ispa_nh3,ispg,iso4,ipco2,ino3,ipnh3,ispa2 &
            & ,ioxa,ispa_oxa,icl,ispa_cl,icharge,ic1,ic2

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

        real(kind=8) kco2,k1,k2,k1no3,rspa_h,rspa_s,rspa_no3,rspa_nh3,knh3,k1nh3,rspa_oxa,rspa_oxa_2,rspa_oxa_3 &
            & ,rspa_cl,rcharge
        real(kind=8),dimension(nz)::pco2x,so4f,no3f,pnh3x,oxaf,clf,fkeq,dfkeq_dios,gamma_tmp,dgamma_dios_tmp
        real(kind=8),dimension(4,nz)::gamma,dgamma_dios
        real(kind=8),dimension(nsp_aq_all)::base_charge

        iso4    = findloc(chraq_all,'so4',dim=1)
        ino3    = findloc(chraq_all,'no3',dim=1)
        ioxa    = findloc(chraq_all,'oxa',dim=1)
        icl     = findloc(chraq_all,'cl',dim=1)
        ipco2   = findloc(chrgas_all,'pco2',dim=1)
        ipnh3   = findloc(chrgas_all,'pnh3',dim=1)

        kco2    = keqgas_h(ipco2,ieqgas_h0)
        k1      = keqgas_h(ipco2,ieqgas_h1)
        k2      = keqgas_h(ipco2,ieqgas_h2)
        knh3    = keqgas_h(ipnh3,ieqgas_h0)
        k1nh3   = keqgas_h(ipnh3,ieqgas_h1)

        pnh3x   = mgasx_loc(ipnh3,:)
        pco2x   = mgasx_loc(ipco2,:)
        so4f    = maqf_loc(iso4,:)
        no3f    = maqf_loc(ino3,:)
        oxaf    = maqf_loc(ioxa,:)
        clf     = maqf_loc(icl,:)

        maqft_loc    = 0d0

        dmaqft_dpro  = 0d0
        dmaqft_dios  = 0d0
        dmaqft_dmaqf = 0d0
        dmaqft_dmgas = 0d0

        do icharge=1,4
            rcharge = 1d0*icharge
            call calc_gamma_davies(  &
                & nz,iosx,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma(icharge,:)=gamma_tmp(:)
            dgamma_dios(icharge,:)=dgamma_dios_tmp(:)
        enddo
            
        call get_base_charge( &
            & nsp_aq_all & 
            & ,chraq_all & 
            & ,base_charge &! output 
            & )

        do ispa = 1, nsp_aq_all
            
            maqft_loc(ispa,:) = maqft_loc(ispa,:) + maqf_loc(ispa,:)
            dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + 1d0
            
            ! complex with NH4
            do ispa_nh3 = 1,2
                rspa_nh3 = real(ispa_nh3,kind=8)
                if ( keqaq_nh3(ispa,ispa_nh3) > 0d0) then 
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
                    maqft_loc(ispa,:) = maqft_loc(ispa,:) + ( &
                        & + fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(pnh3x*knh3/k1nh3*prox)**rspa_nh3 & 
                        & )
                    dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + ( &
                        & + fkeq*keqaq_nh3(ispa,ispa_nh3)*1d0*(pnh3x*knh3/k1nh3*prox)**rspa_nh3 &
                        & )
                    dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) &
                        & + fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(pnh3x*knh3/k1nh3)&
                        & **rspa_nh3*rspa_nh3*prox **(rspa_nh3-1d0)
                    dmaqft_dmgas(ispa,ipnh3,:) = dmaqft_dmgas(ispa,ipnh3,:) &
                        & + fkeq*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(knh3/k1nh3*prox)&
                        & **rspa_nh3*rspa_nh3*pnh3x**(rspa_nh3-1d0)
                    dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + ( &
                        & + dfkeq_dios*keqaq_nh3(ispa,ispa_nh3)*maqf_loc(ispa,:)*(pnh3x*knh3/k1nh3*prox)**rspa_nh3 & 
                        & )
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
                ! maqft_loc(ispa,:) = 1d0
                ! account for hydrolysis speces
                do ispa_h = 1,2
                    rspa_h = real(ispa_h,kind=8)
                    if ( keqaq_h(ispa,ispa_h) > 0d0) then 
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
                        maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**rspa_h
                        dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**rspa_h
                        dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) + &
                                & fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*rspa_h*prox**(rspa_h-1d0)
                        dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**rspa_h
                    endif 
                enddo 
            ! oxalic acid 
            elseif ( &
                & trim(adjustl(chraq_all(ispa)))=='oxa' &
                & .or. trim(adjustl(chraq_all(ispa)))=='glp' &
                & ) then 
                do ispa_h = 1,2
                    if (ispa_h==1)then
                        rspa_h = real(ispa_h,kind=8)
                        if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                            fkeq = 1d0/gamma(2,:)
                            dfkeq_dios = -1d0/gamma(2,:)**2d0*dgamma_dios(2,:)
                            maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)/prox**rspa_h
                            dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*1d0/prox**rspa_h
                            dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) + ( &
                                & + fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(-rspa_h)/prox**(1d0+rspa_h) &
                                & )
                            dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + &
                                &  dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)/prox**rspa_h
                        endif 
                    elseif(ispa_h==2)then
                        rspa_h = real(ispa_h-1,kind=8)
                        if ( keqaq_h(ispa,ispa_h) > 0d0) then 
                            fkeq = gamma(1,:)**2d0
                            dfkeq_dios = 2d0*gamma(1,:)*dgamma_dios(1,:)
                            maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**rspa_h
                            dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*1d0*prox**rspa_h
                            dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) + &
                            & fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*rspa_h*prox**(rspa_h-1d0)
                            dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + &
                            & dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*prox**rspa_h
                        endif 
                    endif 
                enddo 
            ! cations
            else 
                ! maqft_loc(ispa,:) = 1d0
                ! account for hydrolysis speces
                do ispa_h = 1,4
                    rspa_h = real(ispa_h,kind=8)
                    if ( keqaq_h(ispa,ispa_h) > 0d0) then 
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
                        maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)/prox**rspa_h
                        dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_h(ispa,ispa_h)*1d0/prox**rspa_h
                        dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) + &
                        &  fkeq*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)*(-rspa_h)/prox**(1d0+rspa_h)
                        dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + dfkeq_dios*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)/prox**rspa_h
                    endif 
                enddo 
                ! account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                do ispa_c = 1,2
                    if ( keqaq_c(ispa,ispa_c) > 0d0) then 
                        if (ispa_c == 1) then ! with CO3--
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
                            maqft_loc(ispa,:) = maqft_loc(ispa,:) +&
                                & fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x/prox**2d0
                            dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + &
                                & fkeq*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x/prox**2d0
                            dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) &
                                & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(-2d0)/prox**3d0
                            dmaqft_dmgas(ispa,ipco2,:) = dmaqft_dmgas(ispa,ipco2,:) &
                                & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0/prox**2d0
                            dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + ( & 
                                & + dfkeq_dios*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x/prox**2d0 & 
                                & )
                        elseif (ispa_c == 2) then ! with HCO3- ( CO32- + H+)
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
                            maqft_loc(ispa,:) = maqft_loc(ispa,:) + &
                                & fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x/prox
                            dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + &
                                & fkeq*keqaq_c(ispa,ispa_c)*1d0*k1*k2*kco2*pco2x/prox
                            dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) &
                                & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x*(-1d0)/prox**2d0
                            dmaqft_dmgas(ispa,ipco2,:) = dmaqft_dmgas(ispa,ipco2,:) &
                                & + fkeq*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0/prox
                            dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + ( & 
                                & + dfkeq_dios*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*pco2x/prox &
                                & )
                        endif 
                    endif 
                enddo 
                ! account for complexation with free SO4
                do ispa_s = 1,2
                    rspa_s = real(ispa_s,kind=8)
                    if ( keqaq_s(ispa,ispa_s) > 0d0) then 
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
                        maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s
                        dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_s(ispa,ispa_s)*1d0*so4f**rspa_s
                        dmaqft_dmaqf(ispa,iso4,:) = dmaqft_dmaqf(ispa,iso4,:) &
                            & + fkeq*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*rspa_s*so4f**(rspa_s-1d0)
                        dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + dfkeq_dios*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s
                        
                        maqft_loc(iso4,:) = maqft_loc(iso4,:) + rspa_s*fkeq*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s
                        dmaqft_dmaqf(iso4,iso4,:) = dmaqft_dmaqf(iso4,iso4,:) + ( &
                            & + rspa_s*fkeq*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*rspa_s*so4f**(rspa_s-1d0) &
                            & )
                        dmaqft_dmaqf(iso4,ispa,:) = dmaqft_dmaqf(iso4,ispa,:) + ( &
                            & + rspa_s*fkeq*keqaq_s(ispa,ispa_s)*1d0*so4f**rspa_s &
                            & )
                        dmaqft_dios(iso4,:) = dmaqft_dios(iso4,:) + &
                                &rspa_s*dfkeq_dios*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)*so4f**rspa_s
                    endif 
                enddo 
                ! account for complexation with free NO3
                do ispa_no3 = 1,2
                    rspa_no3 = real(ispa_no3,kind=8)
                    if ( keqaq_no3(ispa,ispa_no3) > 0d0) then 
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
                        maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3
                        dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_no3(ispa,ispa_no3)*1d0*no3f**rspa_no3
                        dmaqft_dmaqf(ispa,ino3,:) = dmaqft_dmaqf(ispa,ino3,:) &
                            & + fkeq*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*rspa_no3*no3f**(rspa_no3-1d0)
                        dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + &
                            & dfkeq_dios*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3
                        
                        maqft_loc(ino3,:) = maqft_loc(ino3,:) + &
                            & rspa_no3*fkeq*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3
                        dmaqft_dmaqf(ino3,ino3,:) = dmaqft_dmaqf(ino3,ino3,:) + ( &
                            & + rspa_no3*fkeq*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*rspa_no3*no3f**(rspa_no3-1d0) &
                            & )
                        dmaqft_dmaqf(ino3,ispa,:) = dmaqft_dmaqf(ino3,ispa,:) + ( &
                            & + rspa_no3*fkeq*keqaq_no3(ispa,ispa_no3)*1d0*no3f**rspa_no3 &
                            & )
                        dmaqft_dios(ino3,:) = dmaqft_dios(ino3,:) + ( &
                            & + rspa_no3*dfkeq_dios*keqaq_no3(ispa,ispa_no3)*maqf_loc(ispa,:)*no3f**rspa_no3 &
                            & ) 
                    endif 
                enddo 
                ! account for complexation with free Cl
                do ispa_cl = 1,2
                    rspa_cl = real(ispa_cl,kind=8)
                    if ( keqaq_cl(ispa,ispa_cl) > 0d0) then 
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
                        maqft_loc(ispa,:) = maqft_loc(ispa,:) + fkeq*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl
                        dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + fkeq*keqaq_cl(ispa,ispa_cl)*1d0*clf**rspa_cl
                        dmaqft_dmaqf(ispa,icl,:) = dmaqft_dmaqf(ispa,icl,:) &
                            & + fkeq*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*rspa_cl*clf**(rspa_cl-1d0)
                        dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + &
                            & dfkeq_dios*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl
                        
                        maqft_loc(icl,:) = maqft_loc(icl,:) + rspa_cl*fkeq*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl
                        dmaqft_dmaqf(icl,icl,:) = dmaqft_dmaqf(icl,icl,:) + ( &
                            & + rspa_cl*fkeq*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*rspa_cl*clf**(rspa_cl-1d0) &
                            & )
                        dmaqft_dmaqf(icl,ispa,:) = dmaqft_dmaqf(icl,ispa,:) + ( &
                            & + rspa_cl*fkeq*keqaq_cl(ispa,ispa_cl)*1d0*clf**rspa_cl &
                            & )
                        dmaqft_dios(icl,:) = dmaqft_dios(icl,:) + &
                            &rspa_cl*dfkeq_dios*keqaq_cl(ispa,ispa_cl)*maqf_loc(ispa,:)*clf**rspa_cl
                    endif 
                enddo 
                ! account for complexation with HOxa-
                do ispa_oxa = 1,2
                    rspa_oxa   = real(ispa_oxa,kind=8)
                    rspa_oxa_2 = real(ispa_oxa,kind=8)
                    rspa_oxa_3 = real(ispa_oxa,kind=8)
                    if (trim(adjustl(chraq_all(ispa)))=='al') then
                        rspa_oxa   = real(ispa_oxa,kind=8) + 1d0
                        rspa_oxa_2 = 1d0
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
                        maqft_loc(ispa,:) = maqft_loc(ispa,:) + ( &
                            & + fkeq*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*oxaf**rspa_oxa_3/prox**rspa_oxa & 
                            & )
                        dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:) + ( & 
                            & + fkeq*keqaq_oxa(ispa,ispa_oxa)*1d0*oxaf**rspa_oxa_3/prox**rspa_oxa &
                            & )
                        dmaqft_dmaqf(ispa,ioxa,:) = dmaqft_dmaqf(ispa,ioxa,:) &
                            & + fkeq*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*rspa_oxa_3*oxaf**(rspa_oxa_3-1d0)/prox**rspa_oxa
                        dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:) &
                            & + fkeq*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*oxaf**rspa_oxa_3*(-rspa_oxa)/prox**(rspa_oxa+1d0)
                        dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:) + ( &
                            & + dfkeq_dios*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*oxaf**rspa_oxa_3/prox**rspa_oxa & 
                            & )
                        
                        maqft_loc(ioxa,:) = maqft_loc(ioxa,:) &
                            & + rspa_oxa_2*fkeq*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*oxaf**rspa_oxa_3/prox**rspa_oxa
                        dmaqft_dmaqf(ioxa,ioxa,:) = dmaqft_dmaqf(ioxa,ioxa,:) + ( &
                            & + rspa_oxa_2*fkeq*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)&
                            &*rspa_oxa_3*oxaf**(rspa_oxa_3-1d0)/prox**rspa_oxa &
                            & )
                        dmaqft_dmaqf(ioxa,ispa,:) = dmaqft_dmaqf(ioxa,ispa,:) + ( &
                            & + rspa_oxa_2*fkeq*keqaq_oxa(ispa,ispa_oxa)*1d0*oxaf**rspa_oxa_3/prox**rspa_oxa &
                            & )
                        dmaqft_dpro(ioxa,:) = dmaqft_dpro(ioxa,:) &
                            & + rspa_oxa_2*fkeq*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)&
                            & *oxaf**rspa_oxa_3*(-rspa_oxa)/prox**(rspa_oxa+1d0)
                        dmaqft_dios(ioxa,:) = dmaqft_dios(ioxa,:) &
                            & + rspa_oxa_2*dfkeq_dios*keqaq_oxa(ispa,ispa_oxa)*maqf_loc(ispa,:)*oxaf**rspa_oxa_3/prox**rspa_oxa
                    endif 
                enddo 
            endif 
            
            ! needs to devide total conc with primary conc. 
            ! maqft_loc(ispa,:) = maqft_loc(ispa,:)/maqf_loc(ispa,:)
            
            dmaqft_dpro(ispa,:) = dmaqft_dpro(ispa,:)/maqf_loc(ispa,:)
            dmaqft_dios(ispa,:) = dmaqft_dios(ispa,:)/maqf_loc(ispa,:)
            
            dmaqft_dmaqf(ispa,ispa,:) = dmaqft_dmaqf(ispa,ispa,:)/maqf_loc(ispa,:) + &
            & maqft_loc(ispa,:)*(-1d0)/maqf_loc(ispa,:)**2d0
            do ispa2 = 1,nsp_aq_all
                if (ispa2==ispa) cycle 
                dmaqft_dmaqf(ispa,ispa2,:) = dmaqft_dmaqf(ispa,ispa2,:)/maqf_loc(ispa,:)
            enddo
            
            do ispg=1,nsp_gas_all
                dmaqft_dmgas(ispa,ispg,:) = dmaqft_dmgas(ispa,ispg,:)/maqf_loc(ispa,:)
            enddo 
            
            
            maqft_loc(ispa,:) = maqft_loc(ispa,:)/maqf_loc(ispa,:)
        enddo     

    endsubroutine get_maqt_all

    subroutine get_maqads_all_v4( &
        & nz,nsp_aq_all,nsp_sld_all &
        & ,chraq_all,chrsld_all &
        & ,keqcec_all,keqiex_all,cec_pH_depend,beta_all &
        & ,msldx_loc,maqf_loc,prox &
        & ,dmaqfads_sld_dpro,dmaqfads_sld_dmaqf,dmaqfads_sld_dmsld &! output
        & ,msldf_loc,maqfads_sld_loc,beta_loc,ads_error  &! output
        & )
        ! calculating ratio of adsorbed species relative to maqf_loc
        ! (1) First calculate exposed negatively-charged sites (S-O-) (mol/m3)
        ! (2) Summing up occupied sites e.g. [S-O-Na] = K*[S-O-]*[Na] using K for S-O- + Na+ = S-O-Na
        ! *** Make sure sum([S-O-X]) = K_CEC*msld where K_CEC is in units of charged mol per unit mol of mineral
        ! *** updated version tring to implement adsorption of mono-, di- and tri-charged cations and involvement of no anions (complexes) 
        ! 
        implicit none
        integer,intent(in)::nz,nsp_aq_all,nsp_sld_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_sld_all),intent(in)::keqcec_all,beta_all
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::keqiex_all
        real(kind=8),dimension(nsp_aq_all,nz),intent(in)::maqf_loc
        real(kind=8),dimension(nsp_sld_all,nz),intent(in)::msldx_loc
        real(kind=8),dimension(nz),intent(in)::prox
        logical,dimension(nsp_sld_all),intent(in)::cec_pH_depend

        real(kind=8),dimension(nsp_sld_all,nz),intent(out)::msldf_loc,beta_loc

        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz),intent(out)::maqfads_sld_loc
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nsp_aq_all,nz),intent(out)::dmaqfads_sld_dmaqf
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz),intent(out)::dmaqfads_sld_dmsld
        real(kind=8),dimension(nsp_aq_all,nsp_sld_all,nz),intent(out)::dmaqfads_sld_dpro
        logical,intent(out)::ads_error

        ! local
        integer isps,ispa,ispa2,iter

        real(kind=8),dimension(nsp_sld_all,nz)::dmsldf_dmsld,dmsldf_dpro
        real(kind=8),dimension(nsp_sld_all,nz)::gamma_loc,dgamma_dmsld,dgamma_dmsldf,dgamma_dpro  
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all,nz)::dmsldf_dmaqf
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all,nz)::dgamma_dmaqf
        real(kind=8),dimension(nz)::f,f_chk,x,dx
        real(kind=8),dimension(nz)::a,da_dpro,da_dmsld,da
        real(kind=8),dimension(nz)::gamma,dgamma,beta,dbeta
        real(kind=8),dimension(nsp_aq_all,nz)::da_dmaqf
        real(kind=8),dimension(nsp_aq_all)::base_charge
        real(kind=8) c1_gamma,c0_gamma 
        ! real(kind=8) :: tol_dum = 1d-9
        real(kind=8) :: tol_dum = 1d-12 ! desparate for convergence 6/8/2023
        real(kind=8) :: tol_dum_2 = 1d-8
        ! real(kind=8) :: tol_dum_2 = 1d-6 ! desparate for convergence 6/8/2023
        ! real(kind=8) :: tol_dum = 1d-12   ! when beta_ON = .true.
        ! real(kind=8) :: tol_dum_2 = 1d-8  ! when beta_ON = .true.
        real(kind=8) :: low_lim = 1d-20 
        real(kind=8) :: fact = 1d5
        real(kind=8) error
        ! logical :: low_lim_ON = .true.
        logical :: low_lim_ON = .false. 
        ! logical :: beta_ON = .true. 
        logical :: beta_ON = .false.  
        logical :: gamma_ON = .true. 
        ! logical :: gamma_ON = .false. 

        ! (1) First getting fraction of negatively charged sites occupied with H+ (f[X-H]) (defined as msldf_loc)
        ! 1 = f[X-H]*beta + f[X-Na] + f[X-K] + f[X2-Ca] + f[X2-Mg] + f[X3-Mg]
        ! where: 
        !       f[X-Na]  = [X-Na]/CEC  = K * [Na+] * f[X-H] / [H+] * gamma
        !       f[X-K]   = [X-K] /CEC  = K * [K+]  * f[X-H] / [H+] * gamma
        !       f[X2-Mg] = 2[X2-Mg]/CEC = K * [Mg++] * f[X-H]^2 / [H+]^2  * gamma^2
        !       f[X2-Ca] = 2[X2-Ca]/CEC = K * [Ca++] * f[X-H]^2 / [H+]^2  * gamma^2
        !       f[X3-Al] = 3[X3-Al]/CEC = K * [Al+++] * f[X-H]^3 / [H+]^3 * gamma^3
        ! and 
        !       gamma = 10^(  3.4 * f[X-H] ) 
        !       beta  = 10^( -3.4 * ( 1 - f[X-H] ) ) = 10^-3.4 * gamma
        !       (from Appelo 1994)
        !       *** note that 10^-3.4 for gamma is accounted for in K, i.e., 
        !           log KHX = 2.5 + 3.4( 1 - f(H) ) = 5.9 - 3.4f(H)
        !           and 
        !           log KIH = log KINa - log KHNa 
        !           where  log KHNa = 5.9 in default
        ! f[X-H] ( or msldf_loc) is solved numerically considering Na+, K+, Mg++, Ca++ and Al+++

        msldf_loc = 0d0
        dmsldf_dpro = 0d0
        dmsldf_dmsld = 0d0
        dmsldf_dmaqf = 0d0

        gamma_loc = 0d0
        dgamma_dmsldf = 0d0
        dgamma_dpro = 0d0
        dgamma_dmsld = 0d0
        dgamma_dmaqf = 0d0

        c1_gamma = 3.4d0 ! between 3.1 to 3.7, average 3.4 from Appelo 1994
        ! c1_gamma = 3.1d0
        ! c1_gamma = 3.7d0
        ! c1_gamma = 5.0d0
        ! c1_gamma = 1.0d0
        ! c1_gamma = 0.0d0
        ! c1_gamma = 6.0d0
        ! c1_gamma = 8.0d0
        ! c1_gamma = 2.0d0
        c0_gamma = 0.005d0

        beta_loc = 0d0

        ads_error = .false.

        call get_base_charge( &
            & nsp_aq_all & 
            & ,chraq_all & 
            & ,base_charge &! output 
            & )

        do isps = 1, nsp_sld_all
            
            if (keqcec_all(isps) == 0d0) cycle
            
            select case(trim(adjustl(chrsld_all(isps))))
                case('ka','cabd','mgbd','kbd','nabd','g1','g2','g3','inrt')
                    ! do nothing 
                case default 
                    ! cycle
                    ! do nothing 
            endselect 
            
            ! equation to be solved:
            !       f = 1 - msldf_loc - sum ( keqiex_all(isps,ispa)* maqf_loc(ispa,:)* (msldf_loc /prox(:) ) **base_charge(ispa) ) = 0
            ! seek a solution of msldf_loc (between 0 and 1).   
            
            x = 1d0
            error = 1d4
            iter = 0
            
            c1_gamma = beta_all(isps)
            
            do while (error > tol_dum)
            
                a = 0d0
                da = 0d0
                da_dpro = 0d0
                da_dmaqf = 0d0
                da_dmsld = 0d0
                
                gamma = 10d0**(c1_gamma*x)
                dgamma = 10d0**(c1_gamma*x)*c1_gamma*log(10d0)
                
                if (.not. gamma_ON) then
                    gamma = 10d0**(c1_gamma*c0_gamma)
                    dgamma = 0d0
                endif 
                
                beta = 10d0**( -c1_gamma*( 1d0 - x )  ) 
                dbeta = 10d0**( -c1_gamma*( 1d0 - x )  ) *(c1_gamma)*log(10d0)
                
                if (.not. beta_ON) then
                    beta = 1d0
                    dbeta = 0d0
                endif 
                
                a = a + 1d0 - x * beta
                da = da     - 1d0  * beta - x * dbeta
                
                if (cec_pH_depend(isps)) then 
                    do ispa=1,nsp_aq_all
                        selectcase(trim(adjustl(chraq_all(ispa))))
                            ! case('na','k','mg','ca')
                            case('na','k','mg','ca','al')
                                a = a - keqiex_all(isps,ispa)* maqf_loc(ispa,:)*&
                                    & (x/prox)**base_charge(ispa)*gamma**base_charge(ispa)
                                da = da - keqiex_all(isps,ispa)* maqf_loc(ispa,:)*(1d0/prox)**base_charge(ispa) &
                                    &   *base_charge(ispa)*x**(base_charge(ispa)-1d0)*gamma**base_charge(ispa) &
                                    & - keqiex_all(isps,ispa)* maqf_loc(ispa,:)*(x/prox)**base_charge(ispa) &
                                    &   *base_charge(ispa)*gamma**(base_charge(ispa)-1d0)*dgamma
                                da_dmaqf(ispa,:) = da_dmaqf(ispa,:) &
                                    & - keqiex_all(isps,ispa)*1d0*(x/prox)**base_charge(ispa)*gamma**base_charge(ispa)
                                da_dpro = da_dpro -&
                                    & keqiex_all(isps,ispa)*maqf_loc(ispa,:)*x**base_charge(ispa)*gamma**base_charge(ispa) &
                                    & *(-base_charge(ispa))*(1d0/prox)**(base_charge(ispa)+1d0)
                            case default 
                                ! do nothing
                        endselect
                    enddo
                else 
                    do ispa=1,nsp_aq_all
                        selectcase(trim(adjustl(chraq_all(ispa))))
                            case('na','k','mg','ca','al')
                                a = a - fact*keqiex_all(isps,ispa)* maqf_loc(ispa,:)*x**base_charge(ispa)
                                da = da - fact*keqiex_all(isps,ispa)* maqf_loc(ispa,:)*base_charge(ispa)*x**(base_charge(ispa)-1d0)
                                da_dmaqf(ispa,:) = da_dmaqf(ispa,:) - fact*keqiex_all(isps,ispa)*1d0*x**base_charge(ispa)
                            case default 
                                ! do nothing
                        endselect
                    enddo
                endif 
                
                if (all(abs(a)<tol_dum)) exit 
                
                where (x -a/da>0d0)
                    x = x -a/da
                elsewhere 
                    x = x*exp( -a/da/x )
                endwhere
                error = maxval(abs(exp( -a/da/x )-1d0))
                iter = iter + 1
                
                ! print *, iter,error,maxval(abs(a))
            
            enddo 
            
            if (any(x > 1d0) ) then
                print *, 'solution exceeds 1: get_maqads_all_v4 ',chrsld_all(isps)
                print *,x
                ads_error = .true.
                exit
                stop
            endif 
            
            msldf_loc(isps,:) = x
            f_chk = a
            if (any(abs(f_chk)>tol_dum_2)) then 
                print *, 'mass basalnce not satisfied: get_maqads_all_v4 ',chrsld_all(isps)
                print *,f_chk
                ads_error = .true.
                exit
                stop
            endif 
            
            ! solving deviations analytically:
            ! da/dx + da/dph * dph/dx  = 0          <==> dx/dph = - (da/dph) / (da/dx)
            ! da/dx + da/dmaqf * dmaqf/dx  = 0      <==> dx/dmaqf = - (da/dmaqf) / (da/dx)
            
            if (cec_pH_depend(isps)) then 
                dmsldf_dpro(isps,:) = - da_dpro /da
            endif 
            
            do ispa=1,nsp_aq_all
                dmsldf_dmaqf(isps,ispa,:) =  - da_dmaqf(ispa,:) /da
            enddo
            
            
            if (cec_pH_depend(isps)) then 
                gamma_loc(isps,:) = 10d0**(c1_gamma*x)
                dgamma_dmsldf(isps,:) = 10d0**(c1_gamma*x)*c1_gamma*log(10d0)
                dgamma_dpro(isps,:) = dgamma_dmsldf(isps,:) * dmsldf_dpro(isps,:)
                dgamma_dmsld(isps,:) = dgamma_dmsldf(isps,:) * dmsldf_dmsld(isps,:)
                
                do ispa=1,nsp_aq_all
                    dgamma_dmaqf(isps,ispa,:) = dgamma_dmsldf(isps,:) * dmsldf_dmaqf(isps,ispa,:)
                enddo
                
                if (.not. gamma_ON) then 
                    gamma_loc(isps,:) = 10d0**(c1_gamma*c0_gamma)
                    dgamma_dmsldf(isps,:) = 0d0
                    dgamma_dpro(isps,:) = 0d0
                    dgamma_dmsld(isps,:) = 0d0
                    dgamma_dmaqf(isps,:,:) = 0d0
                endif 
                
                beta_loc(isps,:) = 10d0**(-c1_gamma* (1d0 -  x ) )  
                
                if (.not. beta_ON) beta_loc(isps,:) = 1d0
            else
                gamma_loc(isps,:) = 1d0
                beta_loc(isps,:) = 1d0
            endif
            
        enddo

        if ( ads_error ) return

        ! (2) Then getting concs of adsorbed ion concs. relative to magf (defined here as maqfads_loc)
        ! adsorbed species concs are: 
        !       [X-Na]  (mol/m3) = CEC*f[X-Na]       =       CEC * K * [Na+] * f[X-H] / [H+] * gamma
        !       [X-K]   (mol/m3) = CEC*f[X-K]        =       CEC * K * [K+]  * f[X-H] / [H+] * gamma
        !       [X2-Mg] (mol/m3) = (1/2)*CEC*f[X-Mg] = (1/2)*CEC * K * [Mg++] * f[X-H]^2 / [H+]^2 * gamma^2
        !       [X2-Ca] (mol/m3) = (1/2)*CEC*f[X-Ca] = (1/2)*CEC * K * [Ca++] * f[X-H]^2 / [H+]^2 * gamma^2
        !       [X3-Al] (mol/m3) = (1/3)*CEC*f[X-Al] = (1/3)*CEC * K * [Al+++] * f[X-H]^3 / [H+]^3 * gamma^3
        ! where 
        !       CEC (eq/m3) = msld(isps,:)*keqcec_all(isps)
        !       f[X-H] = msldf_loc
        !       [Na+] = magf_loc(isp == 'na')
        !       [H+]  = prox    
        !       gamma = 10**(3.4*msldf_loc)
        !       etc...    

        maqfads_sld_loc = 0d0
        dmaqfads_sld_dpro = 0d0
        dmaqfads_sld_dmaqf = 0d0
        dmaqfads_sld_dmsld = 0d0

        do ispa=1,nsp_aq_all
            selectcase(trim(adjustl(chraq_all(ispa))))
                case('na','k','mg','ca','al')
                    do isps=1,nsp_sld_all
                        
                        if (keqcec_all(isps) == 0d0) cycle

                        select case(trim(adjustl(chrsld_all(isps))))
                            case('ka','cabd','mgbd','kbd','nabd','g1','g2','g3','inrt')
                                ! do nothing 
                            case default 
                                ! cycle
                                ! do nothing 
                        endselect 
                        
                        if (cec_pH_depend(isps)) then
                                
                            maqfads_sld_loc(ispa,isps,:) = maqfads_sld_loc(ispa,isps,:)  &
                                & + (1d0/base_charge(ispa)) * ( &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*&
                                & keqiex_all(isps,ispa)*(msldf_loc(isps,:)/prox)**base_charge(ispa) &
                                &   *gamma_loc(isps,:)**base_charge(ispa) &
                                & )
                            dmaqfads_sld_dpro(ispa,isps,:) = dmaqfads_sld_dpro(ispa,isps,:) &
                                & + (1d0/base_charge(ispa)) * ( &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*&
                                & msldf_loc(isps,:)**base_charge(ispa) &
                                &   *(-1d0*base_charge(ispa))/(prox**(base_charge(ispa)+1d0))  &
                                &   *gamma_loc(isps,:)**base_charge(ispa) &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*(1d0/prox)**base_charge(ispa) &
                                &   *base_charge(ispa)*msldf_loc(isps,:)**(base_charge(ispa)-1d0) * dmsldf_dpro(isps,:) &
                                &   *gamma_loc(isps,:)**base_charge(ispa) &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*&
                                & (msldf_loc(isps,:)/prox)**base_charge(ispa) &
                                &   *base_charge(ispa)*gamma_loc(isps,:)**(base_charge(ispa)-1d0)*dgamma_dpro(isps,:) &
                                & )
                            do ispa2=1,nsp_aq_all
                                dmaqfads_sld_dmaqf(ispa,isps,ispa2,:) = dmaqfads_sld_dmaqf(ispa,isps,ispa2,:) &
                                & + (1d0/base_charge(ispa)) * ( &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*(1d0/prox)**base_charge(ispa) &
                                &   *base_charge(ispa)*msldf_loc(isps,:)**(base_charge(ispa)-1d0) * dmsldf_dmaqf(isps,ispa2,:)  &
                                &   *gamma_loc(isps,:)**base_charge(ispa) &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*&
                                & (msldf_loc(isps,:)/prox)**base_charge(ispa) &
                                &   *base_charge(ispa)*gamma_loc(isps,:)**(base_charge(ispa)-1d0)*dgamma_dmaqf(isps,ispa2,:) &
                                & )
                            enddo 
                            
                            dmaqfads_sld_dmsld(ispa,isps,:) = dmaqfads_sld_dmsld(ispa,isps,:) &
                                & + (1d0/base_charge(ispa)) * ( &
                                & + keqcec_all(isps)*1d0*keqiex_all(isps,ispa)*(msldf_loc(isps,:)/prox) **base_charge(ispa) &
                                &   *gamma_loc(isps,:)**base_charge(ispa) &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa) *(1d0/prox)**base_charge(ispa) &
                                &   *base_charge(ispa)*msldf_loc(isps,:)**(base_charge(ispa)-1d0) * dmsldf_dmsld(isps,:)   &
                                &   *gamma_loc(isps,:)**base_charge(ispa) &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*&
                                & (msldf_loc(isps,:)/prox)**base_charge(ispa) &
                                &   *base_charge(ispa)*gamma_loc(isps,:)**(base_charge(ispa)-1d0)*dgamma_dmsld(isps,:) &
                                & )
                        else 
                                
                            maqfads_sld_loc(ispa,isps,:) = maqfads_sld_loc(ispa,isps,:) &
                                & + (1d0/base_charge(ispa)) *fact * ( &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa)*&
                                & msldf_loc(isps,:)**base_charge(ispa) &
                                & )
                            dmaqfads_sld_dpro(ispa,isps,:) = dmaqfads_sld_dpro(ispa,isps,:) &
                                & + (1d0/base_charge(ispa)) *fact * ( &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa) &
                                &   *base_charge(ispa)*msldf_loc(isps,:)**(base_charge(ispa)-1d0)*dmsldf_dpro(isps,:)   &
                                & )
                            do ispa2=1,nsp_aq_all
                                dmaqfads_sld_dmaqf(ispa,isps,ispa2,:) = dmaqfads_sld_dmaqf(ispa,isps,ispa2,:) &
                                & + (1d0/base_charge(ispa)) *fact * ( &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa) &
                                &   *base_charge(ispa)*msldf_loc(isps,:)**(base_charge(ispa)-1d0)*dmsldf_dmaqf(isps,ispa2,:)  &
                                & )
                            enddo 
                                
                            
                            dmaqfads_sld_dmsld(ispa,isps,:) = dmaqfads_sld_dmsld(ispa,isps,:)  &
                                & + (1d0/base_charge(ispa)) *fact * ( &
                                & + keqcec_all(isps)*1d0*keqiex_all(isps,ispa)*msldf_loc(isps,:)**base_charge(ispa)   &
                                & + keqcec_all(isps)*msldx_loc(isps,:)*keqiex_all(isps,ispa) &
                                &   *base_charge(ispa)*msldf_loc(isps,:)**(base_charge(ispa)-1d0)*dmsldf_dmsld(isps,:)   &
                                & )
                        
                        endif 
                        
                    enddo
                    
                case default
                    ! do nothing
            endselect 
                            
            if (low_lim_ON) then 
                where(keqcec_all(isps)*msldx_loc(isps,:) < low_lim) 
                    maqfads_sld_loc(ispa,isps,:) = 0d0
                    dmaqfads_sld_dpro(ispa,isps,:) = 0d0
                    dmaqfads_sld_dmsld(ispa,isps,:) = 0d0
                endwhere 

                do ispa2=1,nsp_aq_all
                    where(keqcec_all(isps)*msldx_loc(isps,:) < low_lim) 
                        dmaqfads_sld_dmaqf(ispa,isps,ispa2,:) = 0d0
                    endwhere 
                enddo 
            endif 
        enddo

    endsubroutine get_maqads_all_v4

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


end module scepter_concentration 
