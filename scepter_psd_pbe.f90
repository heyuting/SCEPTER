!************************************************************************
!Name: scepter_psd_pbe
!Purpose: solve for the evolution of the particle size distribution (PSD)
! in a multi-layer (vertical) system, accounting for advection, diffusion, 
! reaction (dissolution/precipitation), and possibly rain input
!************************************************************************
module scepter_psd_pbe
    use scepter_constants
    use scepter_variables
    use scepter_psd
    use lapack95
    implicit none

    private
    public :: psd_diss_pbe

    contains
    !--------------------------------------------------------------------------------------------------
    ! Solve the population balance equation reflecting imposed dissolution rate
    !-------------------------------------------------------------------------------------------------- 
    subroutine psd_diss_pbe( &
        & nz,nps &! in
        & ,z,DV,dt,pi,tol,poro &! in 
        & ,incld_rough,roughref &! in
        & ,psd,ps,dps,ps_min,ps_max &! in 
        & ,chrsp &! in 
        & ,dpsd,psd_error_flg &! inout
        & )
        implicit none 

        integer,intent(in)::nz,nps
        real(kind=8),intent(in)::dt,ps_min,ps_max,pi,tol 
        real(kind=8),dimension(nz),intent(in)::z,poro
        real(kind=8),dimension(nps),intent(in)::ps,dps
        real(kind=8),dimension(nz),intent(in)::DV
        real(kind=8),dimension(nps,nz),intent(in)::psd
        character(5),intent(in)::chrsp
        character(10),intent(in)::roughref
        logical,intent(in)::incld_rough
        real(kind=8),dimension(nps,nz),intent(inout)::dpsd
        logical,intent(inout)::psd_error_flg
        ! local 
        real(kind=8),dimension(nps,nz)::dVd,psd_old,psd_new,dpsd_tmp,psdx,psdxx
        real(kind=8),dimension(nps)::psd_tmp,dvd_tmp,dpsx,lambda
        real(kind=8),dimension(nz)::kpsd,kpsdx,DV_chk,DV_exist
        real(kind=8) ps_new,ps_newp,dvd_res,error,vol,fact,surf
        real(kind=8),parameter::infinity = huge(0d0)
        real(kind=8),parameter::threshold = 20d0
        real(kind=8),parameter::corr = exp(threshold)
        real(kind=8),parameter::threshold_k = 2d0
        real(kind=8),parameter::corr_k = exp(threshold_k)
        integer,parameter :: iter_max = 50
        integer ips,iips,ips_new,iz,isps,row,col,ie,ie2,iter,iiz

        logical :: logcalc = .true.
        ! logical :: logcalc = .false.

        ! logical :: safe_mode = .true.
        logical :: safe_mode = .false.

        logical,dimension(nz) :: ms_not_ok! = .false.

        real(kind=8) amx3(nps+1,nps+1),ymx3(nps+1),emx3(nps+1),ymx3_pre(nps+1)
        integer ipiv3(nps+1) 
        integer info 

        external DGESV

        ! attempt to calculate psd change by dissolution ( defined with particle number / bulk m3 / log (r) )
        ! assumptions/formulations: 
        ! 1. Solve population balance equation with shrinking particle: 
        !       [ ( f(r,t+dt) - f(r,t))/dt ] * dr 
        !               = mv * k * ( lambda(r+dr,t) * f(r+dr,t) - lambda(r,t) * f(r,t) ) 
        !               = k' * ( lambda(r+dr,t) * f(r+dr,t) - lambda(r,t) * f(r,t) ) 
        !               where k' = mv * k
        ! 2. Constraint from main reaction-transport scheme: dV
        !       sigma [ 4 * pi /3 * r^3 * [f(r,t+dt) - f(r,t)]/dt *dr * dt ]  =  dV
        !       sigma [ 4 * pi * r^2 * f(r,t) *dr * dt ]  =  dV
        ! 3. Solve for f(r,t+dt) (r=1, .., nps) and k from above (nps + 1) equations
        !
        ! *** Caution must be paid to dr where in the rest of calculation it is defined as d log(r)

        ! check whether there is enough material to dissolve DV 
        ms_not_ok = .false.
        DV_exist = 0d0
        do iz = 1,nz
            DV_exist(iz) =  sum( 4d0/3d0*pi*(10d0**ps(:))**3d0 * psd(:,iz) * dps(:) )
            if (DV(iz)>0d0 .and.  DV(iz) > DV_exist(iz) ) then 
                ms_not_ok(iz) = .true.
                print *, '*** not enough stuff to dissolve material ',chrsp,iz
                ! exit
            endif 
        enddo 

        if (any(isnan(DV)) .or. any(isnan(psd)) .or. dt==0d0) then  
            print *,any(isnan(DV)),any(isnan(psd)),dt
            stop
        endif 

        ! if ( trim(adjustl(chrsp)) == 'blk') then 
            ! safe_mode = .true.
        ! else 
            ! safe_mode = .false.
        ! endif 
        select case(trim(adjustl(chrsp)))
            case('cc','dlm','arg','gps')
                safe_mode = .false.
            case default
                safe_mode = .true.
        endselect 

        if (safe_mode .and. any(ms_not_ok)) then
            psd_error_flg = .true.
            return
        endif 

        ! if time step is large, return to the main loop with reduced time step
        if (trim(adjustl(chrsp)) == 'blk' .and. dt >= 10d0 .and. any(ms_not_ok)) then
            psd_error_flg = .true.
            return
        endif 

        ! roughness factor as functin of radius
        lambda = 1d0
        if (incld_rough)  lambda = rough_f( roughref, nps, (10d0**ps(:)) )

        ! R = log10 r 
        ! dR/dr =  1/(r * log10)
        ! dr = dR * r * log10
        ! note that dps is dR; we need dr denoted here as dpsx (?)
        ! dpsx = dps * 10d0**(ps) * log(10d0)
        do ips = 1,nps
            dpsx(ips) = 10d0**(ps(ips)+0.5d0*dps(ips)) - 10d0**(ps(ips)-0.5d0*dps(ips))
        enddo

        ! psd is given as number of particles/m3/log(m)
        ! converting to number of particles/m3/m
        ! note that psd*dps = psdx*dpsx
        do iz=1,nz
            psdx(:,iz) = psd(:,iz) *dps(:) / dpsx(:)
        enddo 

        kpsd = 0d0
        do iz = 1, nz
            if (dV(iz)>=0d0) then
                do ips=1,nps
                    if (ips==nps) then
                        kpsd(iz)  = kpsd(iz) - ( - lambda(ips)*psdx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    else
                        kpsd(iz)  = kpsd(iz) - ( lambda(ips+1)*psdx(ips+1,iz) - lambda(ips)*psdx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    endif
                enddo 
            else 
                do ips=1,nps
                    if (ips==1) then
                        kpsd(iz)  = kpsd(iz) - ( lambda(ips)*psdx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    elseif (ips==nps) then
                        kpsd(iz)  = kpsd(iz) - ( - lambda(ips-1)*psdx(ips-1,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    else
                        kpsd(iz)  = kpsd(iz) - ( lambda(ips)*psdx(ips,iz) - lambda(ips-1)*psdx(ips-1,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    endif
                enddo 
            endif
            ! if (kpsd(iz) * dV(iz) < 0d0) then 
                ! print *,' *** inconsistent between kpsd and DV ' , chrsp
                ! psd_error_flg = .true.
                ! return
                ! kpsd(iz) = sum( 4d0*pi*(10d0**ps(:))**2d0 * lambda(:) * psdx(:,iz) * dpsx(:) )
            ! endif 
            ! kpsd(iz) = dV(iz)/dt/ kpsd(iz) 
            ! kpsd(iz) = dV(iz)/dt/ abs( kpsd(iz) )
            
            ! kpsd(iz) = sum( 4d0*pi*(10d0**ps(:))**2d0 * lambda(:) * psdx(:,iz) * dpsx(:) )

            kpsd(iz) = dV(iz)/dt/ kpsd(iz) 
        enddo 

        kpsdx = kpsd
        psdxx = psdx

        do iz = 1, nz
            
            if (.not.safe_mode .and. ms_not_ok(iz)) then 
                print *,'****** DV is larger than the amount existing at ',iz, ' for ' ,chrsp
                ! simply limiting the dissolution by existing particles
                ! dpsd(:,iz) = - psd(:,iz)
                ! cycle
                !!
                ! calling previous subroutine but only at a depth 
                dpsd_tmp = 0d0
                call psd_diss_iz( &
                    & nz,nps,iz &! in
                    & ,z,DV,dt,pi,tol,poro &! in 
                    & ,incld_rough,roughref &! in
                    & ,psd,ps,dps,ps_min,ps_max &! in 
                    & ,chrsp &! in 
                    & ,dpsd_tmp,psd_error_flg &! inout
                    & )
                if (psd_error_flg) return
                dpsd(:,iz) = dpsd_tmp(:,iz)
                cycle
                !!
                ! scaling so that total dissolution volume change is conserved while assuming immediate dissolution 
                ! dpsd(:,iz) = - psd(:,iz) * ( DV(iz) ) &
                    ! & / sum(4d0/3d0*pi*(10d0**ps(:))**3d0 * psd(:,iz) * dps(:))
                ! cycle
                !! 
                ! dvd = 0d0
                ! dVd(:,iz) = ( psd (:,iz) * (10d0**ps(:))**2d0 *lambda(:) )
                ! dVd(:,iz) = dVd(:,iz)*( DV(iz) )/sum(dVd(:,iz) * dps(:))
                ! dpsd(:,iz) = - dVd(:,iz)/(4d0/3d0*pi*(10d0**ps(:))**3d0) 
                ! cycle
                ! !!! 
                ! dpsd(:,iz) = -k *psd(:,iz)
                ! int{dpsd(:,iz)*dps) = -k * int{ psd(:,iz)*dps} = DV
                ! k = DV/int{ psd(:,iz)*dps}
                ! dpsd(:,iz) = - psd(:,iz)*DV(iz)/sum(psd(:,iz)*dps(:))
                ! cycle
                ! or scaling dissolution with reflecting the surface area and trying keeping the mass balance
                ! DV(iz) = sum( 4d0/3*pi*(10d0**ps(:))**3d0 * lambda(:) * dpsd(:,iz) * dps(:) )
                ! dvd = 0d0
                ! dVd(:,iz) = ( psd (:,iz) * (10d0**ps(:))**2d0 *lambda(:) )
                ! dVd(:,iz) = dVd(:,iz)*( DV(iz) -DV_exist(iz) )/sum(dVd(:,iz) * dps(:))
                ! dpsd(:,iz) = dpsd(:,iz) - dVd(:,iz)/(4d0/3d0*pi*(10d0**ps(:))**3d0) 
                ! cycle
                ! or do explicit calculation based on the previous PSD
                do iips = 1, nps
                    vol  = 4d0/3d0*pi*(10d0**ps(iips))**3d0
                    surf  = 4d0*pi*(10d0**ps(iips))**2d0
                    if ( dV(iz) >= 0d0) then 
                        if (iips==nps) then 
                            ! & + ( psdxx(iips,iz) -  psdx(iips,iz) ) /dt * dpsx(iips) &
                            ! & - kpsdx(iz) * ( - lambda(iips)*psdxx(iips,iz) )  &
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) * ( - lambda(iips)*psdx(iips,iz) ) *dt/dpsx(iips)
                        else
                            ! & + ( psdxx(iips,iz) -  psdx(iips,iz) ) /dt * dpsx(iips) &
                            ! & - kpsdx(iz) * ( lambda(iips+1)*psdxx(iips+1,iz) - lambda(iips)*psdxx(iips,iz) )  &
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) & 
                            & * ( lambda(iips+1)*psdx(iips+1,iz) - lambda(iips)*psdx(iips,iz) ) *dt/dpsx(iips)
                        endif
                    
                    else
                    
                        if (iips==1) then
                            ! & + ( psdxx(iips,iz) -  psdx(iips,iz) ) /dt * dpsx(iips) &
                            ! & - kpsdx(iz) * ( lambda(iips)*psdxx(iips,iz) )  &
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) * ( lambda(iips)*psdx(iips,iz) )*dt/dpsx(iips)
                        elseif (iips==nps) then 
                            ! & + ( psdxx(iips,iz) -  psdx(iips,iz) ) /dt * dpsx(iips) &
                            ! & - kpsdx(iz) * ( lambda(iips)*psdxx(iips,iz) - lambda(iips-1)*psdxx(iips-1,iz) )  &
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) * ( - lambda(iips-1)*psdx(iips-1,iz) ) &
                                & *dt/dpsx(iips)
                        else
                            ! & + ( psdxx(iips,iz) -  psdx(iips,iz) ) /dt * dpsx(iips) &
                            ! & - kpsdx(iz) * ( lambda(iips)*psdxx(iips,iz) - lambda(iips-1)*psdxx(iips-1,iz) )  &
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) & 
                            & * (lambda(iips)*psdx(iips,iz) - lambda(iips-1)*psdx(iips-1,iz))*dt/dpsx(iips)
                        endif
                        
                    endif 
                    
                enddo ! end of do-loop for iips
                
                dpsd(:,iz) = psdxx(:,iz)*dpsx(:)/dps(:) - psd(:,iz)
                
                cycle
            endif 

            error = 1d4 
            iter = 1
            ! print*, DV(iz)
            if (DV(Iz) == 0d0) then 
                ! print *, 'DV is zero so dpsd = 0; return' 
                dpsd(:,iz) =0d0
                cycle
            endif 
                    
            do while (error > tol) 
            
                amx3 = 0d0
                ymx3 = 0d0
            
                do ips = 1, nps
                
                    row =  ips 
                    
                    vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    surf  = 4d0*pi*(10d0**ps(ips))**2d0
                    
                    ! if ( dV(iz) >= 0d0) then 
                    if ( kpsdx(iz) >= 0d0) then 
                    
                        if (ips==nps) then
                            ymx3(row) = ( & 
                                & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( - lambda(ips)*psdxx(ips,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( - lambda(ips)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                            
                            col = nps + 1
                            amx3(row,col) = ( & 
                                & - 1d0 * ( - lambda(ips)*psdxx(ips,iz) )  &
                                & ) &
                                & * 1d0
                                
                            ymx3(col) = ymx3(col) + ( &
                                & - kpsdx(iz) * ( - lambda(ips)*psdxx(ips,iz) )*vol  &
                                & )
                                
                            amx3(col,col) = amx3(col,col) + ( &
                                & - 1d0 * ( - lambda(ips)*psdxx(ips,iz) )*vol  &
                                & )&
                                & * 1d0
                                
                            amx3(col,row) = amx3(col,row) + ( &
                                & - kpsdx(iz) * ( - lambda(ips)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                        else
                            ymx3(row) = ( & 
                                & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( - lambda(ips)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                                
                            col = ips + 1
                            amx3(row,col) = ( & 
                                & - kpsdx(iz) * ( lambda(ips+1)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips+1,iz),1d0,logcalc)
                            
                            col = nps + 1
                            amx3(row,col) = ( & 
                                & - 1d0 * (lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz))  &
                                & ) &
                                & * 1d0
                                
                            ! volume conserv.
                            ymx3(col) = ymx3(col) + ( &
                                & - kpsdx(iz) * (lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) )*vol  &
                                & )
                                
                            amx3(col,col) = amx3(col,col) + ( &
                                & - 1d0 * ( lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) )*vol  &
                                & )&
                                & * 1d0
                                
                            amx3(col,row) = amx3(col,row) + ( &
                                & - kpsdx(iz) * ( - lambda(ips)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                                
                            amx3(col,row+1) = amx3(col,row+1) + ( &
                                & - kpsdx(iz) * ( lambda(ips+1)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips+1,iz),1d0,logcalc)
                        endif
                    
                    else
                    
                        if (ips==1) then
                            ymx3(row) = ( & 
                                & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                            
                            col = nps + 1
                            amx3(row,col) = ( & 
                                & - 1d0 * ( lambda(ips)*psdxx(ips,iz) )  &
                                & ) &
                                & * 1d0
                                
                            ! volume conserv.
                            ymx3(col) = ymx3(col) + ( &
                                & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) )*vol  &
                                & )
                                
                            amx3(col,col) = amx3(col,col) + ( &
                                & - 1d0 * ( lambda(ips)*psdxx(ips,iz) )*vol  &
                                & )&
                                & * 1d0
                                
                            amx3(col,row) = amx3(col,row) + ( &
                                & - kpsdx(iz) * ( lambda(ips)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                                
                        elseif (ips==nps) then
                            ymx3(row) = ( & 
                                & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( - lambda(ips-1)*psdxx(ips-1,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & ) &
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                                
                            col = ips - 1
                            amx3(row,col) = ( & 
                                & - kpsdx(iz) * ( - lambda(ips-1)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips-1,iz),1d0,logcalc)
                            
                            col = nps + 1
                            amx3(row,col) = ( & 
                                & - 1d0 * ( - lambda(ips-1)*psdxx(ips-1,iz) )  &
                                & ) &
                                & * 1d0
                                
                            ! volume conserv.
                            ymx3(col) = ymx3(col) + ( &
                                & - kpsdx(iz) * ( - lambda(ips-1)*psdxx(ips-1,iz) )*vol  &
                                & )
                                
                            amx3(col,col) = amx3(col,col) + ( &
                                & - 1d0 * (  - lambda(ips-1)*psdxx(ips-1,iz) )*vol  &
                                & )&
                                & * 1d0
                                
                            amx3(col,row-1) = amx3(col,row-1) + ( &
                                & - kpsdx(iz) * ( - lambda(ips-1)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips-1,iz),1d0,logcalc)
                        else
                            ymx3(row) = ( & 
                                & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                                
                            col = ips - 1
                            amx3(row,col) = ( & 
                                & - kpsdx(iz) * ( - lambda(ips-1)*1d0 )  &
                                & ) &
                                & * merge( psdxx(ips-1,iz),1d0,logcalc)
                            
                            col = nps + 1
                            amx3(row,col) = ( & 
                                & - 1d0 * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )  &
                                & ) &
                                & * 1d0
                            
                            ! volume conserv.
                            ymx3(col) = ymx3(col) + ( &
                                & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )*vol  &
                                & )
                                
                            amx3(col,col) = amx3(col,col) + ( &
                                & - 1d0 * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )*vol  &
                                & )&
                                & * 1d0
                                
                            amx3(col,row) = amx3(col,row) + ( &
                                & - kpsdx(iz) * ( lambda(ips)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips,iz),1d0,logcalc)
                                
                            amx3(col,row-1) = amx3(col,row-1) + ( &
                                & - kpsdx(iz) * ( - lambda(ips-1)*1d0 )*vol  &
                                & )&
                                & * merge( psdxx(ips-1,iz),1d0,logcalc)
                        endif
                        
                    endif 
                    
                    ! fact = max( abs( ymx3(row) ), maxval( abs( amx3(row,:) ) ) )
                    
                    ! ymx3(row) = ymx3(row) / fact
                    ! amx3(row,:) = amx3(row,:) / fact
                    
                enddo ! end of do-loop for ips
                
                row = nps+1
                ymx3(row) = ymx3(row) + ( &
                    & - dV(iz)/dt  &
                    & )
                
                ! do ips=1,nps
                    ! vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    ! surf  = 4d0*pi*(10d0**ps(ips))**2d0
                                
                    ! ymx3(row) = ymx3(row) + ( &
                        ! & - kpsdx(iz)* surf * lambda(ips)*psdxx(ips,iz) *dpsx(ips)   &
                        ! & )
                        
                    ! amx3(row,row) = amx3(row,row) + ( &
                        ! & - 1d0 * surf * lambda(ips)*psdxx(ips,iz) *dpsx(ips)   &
                        ! & )&
                        ! & * 1d0
                    
                    ! col = ips
                    ! amx3(row,col) = amx3(row,col) + ( &
                        ! & - kpsdx(iz) * surf * lambda(ips)*1d0 *dpsx(ips)   &
                        ! & )&
                        ! & * psdxx(ips,iz)
                        
                                
                    ! ymx3(row) = ymx3(row) + ( &
                        ! & - vol * ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & )
                    
                    ! col = ips
                    ! amx3(row,col) = amx3(row,col) + ( &
                        ! & - vol * ( 1d0 -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & )&
                        ! & * psdxx(ips,iz)
                        
                        
                ! enddo 
                    
                fact = 1d1
                
                ymx3(row) = ymx3(row) * fact
                amx3(row,:) = amx3(row,:) * fact
            
                ymx3=-1.0d0*ymx3
                ymx3_pre = ymx3

                if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
                ! if (.true.) then 
                    print*,'PBE--'//chrsp//': error in mtx'
                    print*,'PBE--'//chrsp//': any(isnan(amx3)),any(isnan(ymx3))'
                    print*,any(isnan(amx3)),any(isnan(ymx3))
                    psd_error_flg = .true.
                    ! pause
                    exit

                    if (any(isnan(ymx3))) then 
                        do ie = 1, nps+1
                            if (isnan(ymx3(ie))) then 
                                print*,'PBE: NAN is here...',iz,ie
                            endif
                        enddo 
                    endif


                    if (any(isnan(amx3))) then 
                        do ie = 1,nps+1
                            do ie2 = 1,nps+1
                                if (isnan(amx3(ie,ie2))) then 
                                    print*,'PBE: NAN is here...',iz,ie,ie2
                                endif
                            enddo
                        enddo
                    endif
                    

                    open(unit=11,file='amx.txt',status = 'replace')
                    open(unit=12,file='ymx.txt',status = 'replace')
                    do ie = 1,nps+1
                        write(11,*) (amx3(ie,ie2),ie2 = 1,nps+1)
                        write(12,*) ymx3(ie)
                    enddo 
                    close(11)
                    close(12)         
                    stop
                    
                endif

                call DGESV(nps+1,int(1),amx3,nps+1,IPIV3,ymx3,nps+1,INFO) 

                if (any(isnan(ymx3)) .or. info/=0 ) then
                    print*,'*** PBE--'//chrsp//': error in soultion',any(isnan(ymx3)),info, kpsdx(iz),ms_not_ok(iz)
                    psd_error_flg = .true.
                    ! pause
                    exit
                    

                    ! open(unit=11,file='amx.txt',status = 'replace')
                    ! open(unit=12,file='ymx.txt',status = 'replace')
                    ! open(unit=10,file='ymx_pre.txt',status = 'replace')
                    ! do ie = 1,nps+1
                        ! write(11,*) (amx3(ie,ie2),ie2 = 1,nps+1)
                        ! write(12,*) ymx3(ie)
                        ! write(10,*) ymx3_pre(ie)
                    ! enddo 
                    ! close(11)
                    ! close(12)         
                    ! close(10)         
                    ! stop
                    
                endif
                
                do ips = 1,nps
                    row =  ips 

                    if (isnan(ymx3(row))) then 
                        print *,'PBE--'//chrsp//': nan at', iz,ips
                        stop
                    endif
                    
                    
                    if (logcalc) then
                        ! emx3(row) = dpsx(ips)*psdxx(ips,iz)*exp(ymx3(row)) - dpsx(ips)*psdxx(ips,iz)
                        emx3(row) = exp(abs(ymx3(row))) - 1d0
                        
                        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                            psdxx(ips,iz) = psdxx(ips,iz)*corr
                        else if (ymx3(row) < -threshold) then 
                            psdxx(ips,iz) = psdxx(ips,iz)/corr
                        else   
                            psdxx(ips,iz) = psdxx(ips,iz)*exp(ymx3(row))
                        endif
                    else
                        emx3(row) = abs(ymx3(row)/psdxx(ips,iz))
                        psdxx(ips,iz) = psdxx(ips,iz) + ymx3(row)
                    endif 
                enddo 
                
                row = nps+1
                
                if (isnan(ymx3(row))) then 
                    print *,'PBE--'//chrsp//': nan at', iz
                    stop
                endif
                
                ! linear case
                
                ! if ( ( kpsdx(iz) + ymx3(row) ) *kpsdx(iz) >= 0d0) then ! signs are the same 
                    ! emx3(row) = log(1d0 + ymx3(row)/kpsdx(iz))
                    
                    ! if ((.not.isnan(emx3(row))).and.emx3(row) >threshold_k) then 
                        ! kpsdx(iz) = kpsdx(iz) * corr_k
                    ! else if (emx3(row) < -threshold_k) then 
                        ! kpsdx(iz) = kpsdx(iz) / corr_k
                    ! else   
                        ! kpsdx(iz) = kpsdx(iz)*exp(emx3(row))
                    ! endif
                ! else ! signs are different 
                    ! emx3(row) = abs(ymx3(row)/kpsdx(iz))
                    ! kpsdx(iz) = kpsdx(iz) + ymx3(row)
                ! endif 
                emx3(row) = abs(ymx3(row)/kpsdx(iz))
                kpsdx(iz) = kpsdx(iz) + ymx3(row)
                
                error = maxval(emx3)

                if ( isnan(error) .or. any(isnan(psdxx)) .or. any(isnan(kpsdx)) ) then 
                    error = 1d3
                    print*, 'PBE--'//chrsp//': !! error is NaN; values are returned to those before iteration with reducing dt'
                    print*, 'PBE--'//chrsp//': isnan(error), info/=0,any(isnan(pdsx))'
                    print*, isnan(error), any(isnan(psdx)), any(isnan(kpsdx)) 
                    
                    psd_error_flg = .true.
                    stop
                    exit
                endif
#ifdef show_PSDiter
                print '(a,E11.3,a,i0,a,E11.3,a,E11.3,a,E11.3,a,E11.3)' &
                    & , 'PBE--'//chrsp//': iteration error = ',error, ', iteration = ',iter &
                    & ,', time step [yr] = ',dt &
                    & , ', max psd = ',maxval(psdxx(:,iz)) &
                    & , ', min psd = ',minval(psdxx(:,iz))   &
                    & , ', diss-rate [m/yr] = ',kpsdx(iz) 
                ! print *, error > tol
#endif 
                iter = iter + 1 
                
                if (iter > iter_Max ) then
                    if (dt==0d0) then 
                        print *, chrsp,'dt==0d0; stop'
                        stop
                    endif 
                    psd_error_flg = .true.
                    exit 
                end if
                
                if (psd_error_flg) exit
                
            enddo 
            
            dpsd(:,iz) = psdxx(:,iz)*dpsx(:)/dps(:) - psd(:,iz)
            ! dpsd(:,iz) = -dpsd(:,iz)

            DV_chk(iz) = 0d0
            if (dV(iz)>=0d0) then
                do ips=1,nps
                    if (ips==nps) then
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( - lambda(ips)*psdxx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    else
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    endif
                enddo 
            else 
                do ips=1,nps
                    if (ips==1) then
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    elseif (ips==nps) then
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * (  - lambda(ips-1)*psdxx(ips-1,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    else
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    endif
                enddo 
            endif
            ! DV_chk(iz) = - sum( 4d0/3d0*pi*(10d0**ps(:))**3d0 * (psdxx(:,iz) - psdx(:,iz))/dt * dpsx(:))
            
            ! DV_chk(iz) = sum(kpsdx(iz) * 4d0*pi*(10d0**ps(:))**2d0 *lambda(:) * psdxx(:,iz) * dpsx(:))
            
            
            DV_chk(iz) = DV_chk(iz) *dt
            
            if (abs(DV(iz))>tol .and. abs( (DV(iz) - DV_chk(iz))/DV(iz) ) > tol) then 
            ! if (abs(DV(iz))>0d0 .and. abs( (DV(iz) - DV_chk(iz))/DV(iz) ) > tol) then 
                print *, 'PBE--volume conservation not satisfied',iz,DV(iz),DV_chk(iz)
                psd_error_flg = .true.
            endif 
                
            if (psd_error_flg) exit

        enddo 

    endsubroutine psd_diss_pbe

endmodule scepter_psd_pbe