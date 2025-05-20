!***********************************************************************
! Module: scepter_eq_pH
! Purpose: Calculate pH and related quantities
!***********************************************************************

module scepter_eq_ph
    use scepter_constants ! Constants
    use scepter_concentration ! Concentration calculations
    use scepter_equilibrium ! Equilibrium calculations
    use scepter_eq_charge ! Charge balance calculations
    use scepter_findloc ! Find location of a value in an array
    implicit none
    private
    public :: calc_pH_v7_4

    contains
    !-----------------------------------------------------------------------
    !Subroutine: calc_pH_v7_4
    !Purpose: Calculate pH of the system
    !-----------------------------------------------------------------------
    subroutine calc_pH_v7_4( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input
        & ,poro,sat,tc &! input  
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all, keqaq_no3,keqaq_nh3 &! input
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

        external DGESV 

        integer,intent(in)::nz
        real(kind=8),intent(in)::kw,tc
        real(kind=8),dimension(nz)::so4x,prox_save,error_save,prox_save_newton,prox_init
        real(kind=8),dimension(nz)::iosx_save,ios_new
        real(kind=8),dimension(nz),intent(in)::z,poro,sat
        real(kind=8),dimension(nz),intent(inout)::prox
        logical,intent(out)::ph_error

        real(kind=8),dimension(nz)::ph_add_order,prox_tmp1,prox_tmp2
        real(kind=8),dimension(nz)::df1,f1,f2,df2,df12,d2f1
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
        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc
        real(kind=8),dimension(nsp_aq_all,nz)::df1dmaq,df1dmaqf,d2f1dmaqf
        real(kind=8),dimension(nsp_gas_all,nz)::df1dmgas,df2dmgas,d2f1dmgas

        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::dprodmaq_all
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::dprodmgas_all

        real(kind=8),dimension(nsp_aq_all,nz),intent(out)::diosdmaq_all
        real(kind=8),dimension(nsp_gas_all,nz),intent(out)::diosdmgas_all

        real(kind=8),dimension(nsp_aq_all,nz)::maqtmp_loc
        real(kind=8),dimension(nsp_gas_all,nz)::dmgas,mgastmp_loc
        real(kind=8),dimension(nz)::df1_dum,f1_dum,d2f1_dum,f1_tmp
        real(kind=8),dimension(nz)::f1_tmp1,df1_tmp1,d2f1_tmp1
        real(kind=8),dimension(nz)::f1_tmp2,df1_tmp2,d2f1_tmp2
        real(kind=8),dimension(nsp_aq_all,nz)::df1dmaqf_dum,d2f1dmaqf_dum,df1dmaqf_tmp,d2f1dmaqf_tmp
        real(kind=8),dimension(nsp_gas_all,nz)::df1dmgas_dum,d2f1dmgas_dum,df1dmgas_tmp,d2f1dmgas_tmp
        real(kind=8),dimension(nz)::df1df2,df2df1
        real(kind=8),dimension(nsp_aq_all,nz)::df2dmaqf

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
        logical print_res
        logical bisec_chk,bisec_chk_ON,bisec_only,mod_ph_order,calc_simple,halley,first_chk_done

        real(kind=8),allocatable::amx(:,:),ymx(:)
        integer,allocatable::ipiv(:)
        integer info,nmx

        real(kind=8),parameter :: threshold = 10d0

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
                        ph_tmp_min = 10d0**(-ph_max) 
                        ph_tmp_max = 10d0**(-ph_min)

                        ! initially give slp a random negative value to be saved to slp_save
                        slp = -100d0
                        ! print *,'start from alkaline pH'
                        do iph2=1,nph2 ! start from alkaline pH
                            ph_tmp = ph_max + (ph_min - ph_max) &
                                & * (real(iph2,kind=8)-1d0)/(real(nph2,kind=8)-1d0) 
                            ph_tmp = 10d0**(-ph_tmp) 
                            
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
                            ph_tmp = 10d0**(-ph_tmp) 
                            
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
                        error = abs( 10d0**(-ph_min) -  10d0**(-ph_max))/10d0**(-ph_max)
                        
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
    
endmodule scepter_eq_ph
