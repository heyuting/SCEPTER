!**************************************************************************************
! Module: scepter_calc_khgas
! Purpose: Calculate equilibrium constants for gaseous species
!**************************************************************************************

module scepter_calc_khgas
    use scepter_constants
    use scepter_thermodynamics
    use scepter_concentration
    use scepter_findloc
    
    implicit none
    private
    public :: calc_khgas_all_v2

    contains

    !--------------------------------------------------------------------------------------
    ! Subroutine: calc_khgas_all_v2
    ! Purpose: Calculate equilibrium constants for gaseous species
    !--------------------------------------------------------------------------------------
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

    
endmodule scepter_calc_khgas