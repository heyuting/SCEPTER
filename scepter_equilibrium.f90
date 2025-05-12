!**************************************************************************************
! Module: scepter_equilibrium
! Purpose: Calculate equilibrium concentrations of aqueous, gaseous, and solid species
!**************************************************************************************        

module scepter_equilibrium
    use scepter_eq_coefs
    use scepter_eq_pH
    use scepter_eq_charge
    implicit none
    private
    public :: calc_charge_balance_point
    
contains
    !-----------------------------------------------------------------------
    ! Subroutine: calc_charge_balance_point
    ! Purpose: Calculate charge balance at a point
    !-----------------------------------------------------------------------
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