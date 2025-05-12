!**************************************************************************************
! Module: scepter_eq_charge
! Purpose: Calculate charge balance
!**************************************************************************************

module scepter_eq_charge
    implicit none
    private
    public :: calc_charge_balance

    contains

    !-----------------------------------------------------------------------
    ! Subroutine: calc_charge_balance
    ! Purpose: Calculate charge balance
    !-----------------------------------------------------------------------
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

endmodule scepter_eq_charge