module scepter_psd
    use scepter_constants
    use scepter_kinetics
    implicit none
    public :: calc_psd_pr
    public :: calc_p80
    public :: psd_diss_iz


contains
!--------------------------------------------------------------------------------------------------
! Calculate parent rock particle size distribution
!--------------------------------------------------------------------------------------------------
    subroutine calc_psd_pr( &
        & nps &! input
        & ,pi,p80,ps_sigma_std,poroi,volsld,tol &! input
        & ,ps,dps &! input
        & ,msldunit &! input
        & ,psd_pr &! output 
        & )

        implicit none
        ! input 
        integer,intent(in)::nps
        real(kind=8),intent(in)::pi,p80,ps_sigma_std,poroi,volsld,tol
        real(kind=8),dimension(nps),intent(in)::ps,dps
        character(3),intent(in)::msldunit
        ! output 
        real(kind=8),dimension(nps),intent(out)::psd_pr
        ! local
        real(kind=8) psu_pr,pssigma_pr


        psu_pr = log10(p80)
        pssigma_pr = 1d0
        pssigma_pr = ps_sigma_std

        ! calculate parent rock particle size distribution 
        psd_pr = 1d0/pssigma_pr/sqrt(2d0*pi)*exp( -0.5d0*( (ps - psu_pr)/pssigma_pr )**2d0 )

        ! to ensure sum is 1
        ! print *, sum(psd_pr*dps)
        psd_pr = psd_pr/sum(psd_pr*dps)  
        ! print *, sum(psd_pr*dps)
        ! stop

        ! balance for volumes
        ! sum(msldi*mv*1d-6) (m3/m3) must be equal to sum( 4/3(pi)r3 * psd_pr * dps) 
        ! where psd is number / bulk m3 / log r 
        ! (if msld is defined as mol/sld m3 then msldi needs to be multiplied by (1 - poro)
        ! volsld = sum(msldi*mv*1d-6) + mblki*mvblk*1d-6 ! bulk case
        psd_pr = psd_pr*( volsld )  &
            & /sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_pr(:)*dps(:))

        if ( msldunit == 'sld') then 
            psd_pr = psd_pr*(1d0-poroi)

            if ( abs( ( ( volsld ) &
                & * (1d0-poroi) & 
                & -sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_pr(:)*dps(:))) &
                & / ( ( volsld ) &
                & * (1d0-poroi) & 
                & )  ) > tol) then 
                print *,( volsld ) &
                    & * (1d0-poroi) & 
                    & ,sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_pr(:)*dps(:))
                stop
            endif 
        elseif ( msldunit == 'blk') then 

            if ( abs( ( ( volsld ) &
                & -sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_pr(:)*dps(:))) &
                & / ( ( volsld ) &
                & )  ) > tol) then 
                print *,( volsld ) &
                    & ,sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_pr(:)*dps(:))
                stop
            endif 
        endif 

    endsubroutine calc_psd_pr

!--------------------------------------------------------------------------------------------------
!  calculate the P80 value from a given PSD
!--------------------------------------------------------------------------------------------------
    subroutine calc_p80( &
        & nps,ps,intpsd &! input 
        & ,p80_tmp &! output
        & )

        implicit none
        ! input 
        integer,intent(in)::nps
        real(kind=8),dimension(nps),intent(in)::ps,intpsd
        ! output 
        real(kind=8),intent(out)::p80_tmp
        ! local
        real(kind=8),dimension(nps)::dm
        real(kind=8) slp
        integer ips
        logical found_p80

        dm = 10d0**ps*1d6*2d0 ! diameter in um

        p80_tmp = 0d0
        found_p80 = .false.
        do ips=1,nps-1
            if ( (intpsd(ips)<=0.80d0) .and. (intpsd(ips+1)>=0.80d0) ) then
                slp = ( intpsd(ips+1) - intpsd(ips) ) / ( dm(ips+1) - dm(ips) )
                p80_tmp = dm(ips) + ( 0.80d0 - intpsd(ips) ) / slp
                found_p80 = .true.
            endif 
        enddo

        if (.not. found_p80) then 
            print*,'p80 could not be found'
            print*,dm
            print*,intpsd
            stop
        endif 

    endsubroutine calc_p80

!--------------------------------------------------------------------------------------------------
! Dissolution/precipitation of particles for a specific layer (i.e., iz)    
!-------------------------------------------------------------------------------------------------- 
    subroutine psd_diss_iz( &
        & nz,nps,iz &! in
        & ,DV,pi,tol &! in 
        & ,incld_rough,roughref &! in
        & ,psd,ps,dps,ps_min,ps_max &! in 
        & ,chrsp &! in 
        & ,dpsd,psd_error_flg &! inout
        & )
        implicit none 

        integer,intent(in)::nz,nps,iz
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
        real(kind=8),dimension(nps,nz)::dVd,psd_old,psd_new,dpsd_tmp
        real(kind=8),dimension(nps)::dvd_tmp,rough_tmp
        real(kind=8) ps_new,ps_newp,dvd_res
        integer ips,iips,ips_new
        logical :: safe_mode = .false.
        ! logical :: safe_mode = .true.

        ! attempt to do psd ( defined with particle number / bulk m3 / log (r) )
        ! assumptions: 
        ! 1. particle numbers are only affected by transport (including raining/dusting) 
        ! 2. dissolution does not change particle numbers: it only affect particle distribution 
        ! unless particle is the minimum radius. in this case particle can be lost via dissolution 
        ! 3. when a mineral precipitates, it is assumed to increase particle radius?
        ! e.g., when a 1 um of particle is dissolved by X m3, its radius is changed and this particle is put into a different bin of (smaller) radius 

        ! sum of volume change of minerals at iz is DV = sum(flx_sld(5:5+nsp_sld,iz)*mv(:)*1d-6)*dt (m3 / m3) 
        ! this must be distributed to different particle size bins (dV(r)) in proportion to psd * (4*pi*r^2)
        ! dV(r) = DV/(psd*4*pi*r^2) where dV is m3 / bulk m3 / log(r) 
        ! has to modify so that sum( dV * dps ) = DV 
        ! new psd is obtained by dV(r) = psd(r)*( 4/3 * pi * r^3 - 4/3 * pi * r'^3 ) where r' is the new radius as a result of dissolution
        ! if r' is exactly one of ps value (ps(ips) == r'), then  psd(r') = psd(r)
        ! else: 
        ! first find the closest r* value which is one of ps values.
        ! then DV(r) = psd(r)* 4/3 * pi * r^3 - psd(r*)* 4/3 * pi * r*^3 
        ! i.e., psd(r*) = [ psd(r)* 4/3 * pi * r^3 - DV(r)]  /( 4/3 * pi * r*^3)
        !               = [ psd(r)* 4/3 * pi * r^3 - psd(r)*( 4/3 * pi * r^3 - 4/3 * pi * r'^3 ) ] /( 4/3 * pi * r*^3)
        !               = psd(r) * 4/3 * pi * r'^3 /( 4/3 * pi * r*^3)
        !               = psd(r) * (r'/r*)^3
        ! in this way volume is conservative? 
        ! check: sum( psd(r) * 4/3 * pi * r^3 * dps) - sum( psd(r') * 4/3 * pi * r'^3 * dps) = DV 

        dpsd_tmp = 0d0 
        psd_old = psd
        psd_new = psd
        ! do iz=1,nz

        ! the following should be removed
        do ips = 1, nps
            if ( psd (ips,iz) /= 0d0 ) then 
                dVd(ips,iz) = DV(iz)/ ( psd (ips,iz) * (10d0**ps(ips))**2d0 )
            else 
                dVd(ips,iz) = 0d0
            endif 
        enddo 

        ! correct one?
        dVd = 0d0
        rough_tmp = 1d0
        if (.not.incld_rough) then 
            rough_tmp(:) = rough_f( 'smooth', nps, (10d0**ps(:)) )
        else
            rough_tmp(:) = rough_f( roughref, nps, (10d0**ps(:)) )
        endif 
        dVd(:,iz) = ( psd (:,iz) * (10d0**ps(:))**2d0 *rough_tmp(:) )

        if (all(dVd == 0d0)) then 
            print *,chrsp,'all dissolved loc1?'
            psd_error_flg = .true.
            return
            stop
        endif 

        ! scale with DV
        dVd(:,iz) = dVd(:,iz)*DV(iz)/sum(dVd(:,iz) * dps(:))

        if ( abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz)) > tol ) then
            print *, chrsp,' vol. balance failed somehow ',abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz))
            print *, iz, sum(dVd(:,iz) * dps(:)), DV(iz)
            stop
        endif 

        if (any(isnan(dVd(:,iz))) ) then 
            print *,chrsp, 'nan in dVd loc 1'
            stop
        endif 

        do ips = 1, nps
            
            if ( psd(ips,iz) == 0d0) cycle
            
            if ( ips == 1 .and. dVd(ips,iz) > 0d0 ) then 
                ! this is the minimum size dealt within the model 
                ! so if dissolved (dVd > 0), particle number must reduce 
                ! (revised particle volumes) = (initial particle volumes) - (volume change) 
                ! psd'(ips,iz) * 4d0/3d0*pi*(10d0**ps(ips))**3d0 =  psd(ips,iz) * 4d0/3d0*pi*(10d0**ps(ips))**3d0 - dVd(ips,iz) 
                ! [ psd'(ips,iz) - psd(ips,iz) ] * 4d0/3d0*pi*(10d0**ps(ips))**3d0 = - dVd(ips,iz) 
                if ( .not. safe_mode ) then  ! do not care about producing negative particle number
                    dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0) 
                elseif ( safe_mode ) then   ! never producing negative particle number
                    if ( dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0) < psd(ips,iz) ) then ! when dissolution does not consume existing particles 
                        dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0) 
                    else ! when dissolution exceeds potential consumption of existing particles 
                        ! dvd_res is defined as residual volume to be dissolved 
                        dvd_res = dVd(ips,iz) - psd(ips,iz)*(4d0/3d0*pi*(10d0**ps(ips))**3d0)  ! residual 
                        
                        ! distributing the volume to whole radius 
                        ! the wrong one?
                        dvd_tmp = 0d0
                        do iips = ips+1,nps
                            if ( psd (iips,iz) /= 0d0 ) then 
                                dVd_tmp(iips) = dvd_res *dps(ips)/ ( psd (iips,iz) * (10d0**ps(iips))**2d0 )
                            else 
                                dVd_tmp(iips) = 0d0
                            endif 
                        enddo 
                        
                        ! correct one?
                        dVd_tmp = 0d0
                        ! if (.not.incld_rough) then 
                            ! dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 )
                        ! else
                            ! dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 *rough_c0*(10d0**ps(ips+1:))**rough_c1 )
                        ! endif 
                        dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 *rough_tmp(ips+1:) )
                        
                        if (all(dVd_tmp == 0d0)) then 
                            print *,chrsp,'all dissolved loc2?',ips, psd(ips+1:,iz)
                            psd_error_flg = .true.
                            return
                            stop
                        endif 
                        
                        ! scale with dvd_res*dps
                        dVd_tmp(ips+1:) = dVd_tmp(ips+1:)*dvd_res*dps(ips)/sum(dVd_tmp(ips+1:) * dps(ips+1:))
                        
                        ! dVd(ips+1:,iz) = dVd(ips+1:,iz) + dvd_res/(nps - ips)
                        dVd(ips+1:,iz) = dVd(ips+1:,iz) + dVd_tmp(ips+1:)
                        dVd(ips,iz) = dVd(ips,iz) - dvd_res
                        
                        if ( abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz)) > tol ) then
                            print *,chrsp, ' vol. balance failed somehow loc2 ',abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz))
                            print *, iz, sum(dVd(:,iz) * dps(:)), DV(iz)
                            stop
                        endif 
                        
                        if (any(isnan(dVd(:,iz))) ) then 
                            print *, chrsp,'nan in dVd loc 2'
                            stop
                        endif 
                        
                        ! if (dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0) > psd(ips,iz)) then 
                            ! print *, 'error: stop',psd(ips,iz),dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0)
                            ! stop
                        ! endif 
                        
                        ! dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0) 
                        dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - psd(ips,iz) 
                    endif 
                endif 
            
            elseif ( ips == nps .and. dVd(ips,iz) < 0d0 ) then 
                ! this is the max size dealt within the model 
                ! so if precipirated (dVd < 0), particle number must increase  
                ! (revised particle volumes) = (initial particle volumes) - (volume change) 
                ! psd'(ips,iz) * 4d0/3d0*pi*(10d0**ps(ips))**3d0 =  psd(ips,iz) * 4d0/3d0*pi*(10d0**ps(ips))**3d0 - dVd(ips,iz) 
                ! [ psd'(ips,iz) - psd(ips,iz) ] * 4d0/3d0*pi*(10d0**ps(ips))**3d0 = - dVd(ips,iz) 
                dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0**ps(ips))**3d0) 
            
            else 
                ! new r*^3 after dissolution/precipitation 
                ps_new =  ( 4d0/3d0*pi*(10d0**ps(ips))**3d0 - dVd(ips,iz) /psd(ips,iz) )/(4d0/3d0*pi) 
                
                if (ps_new <= 0d0) then  ! too much dissolution so removing all particles cannot explain dVd at a given bin ips
                    ps_new = ps_min
                    ps_newp = 10d0**ps(1)
                    dpsd_tmp(1,iz) =  dpsd_tmp(1,iz) + psd(ips,iz)
                    dpsd_tmp(ips,iz) =  dpsd_tmp(ips,iz) - psd(ips,iz)
                    
                    dvd_res = dVd(ips,iz) - psd(ips,iz)*(4d0/3d0*pi*(10d0**ps(ips))**3d0)  ! residual
                    
                    ! distributing the volume to whole radius 
                    ! wrong one ?
                    dvd_tmp = 0d0
                    do iips = ips+1,nps
                        if ( psd (iips,iz) /= 0d0 ) then 
                            dVd_tmp(iips) = dvd_res *dps(ips)/ ( psd (iips,iz) * (10d0**ps(iips))**2d0 )
                        else 
                            dVd_tmp(iips) = 0d0
                        endif 
                    enddo 
                    
                    ! correct one?
                    dVd_tmp = 0d0
                    ! if (.not.incld_rough) then 
                        ! do iips=1,nps
                            ! if (iips == ips) cycle
                            ! dVd_tmp(iips) = ( psd (iips,iz) * (10d0**ps(iips))**2d0 )
                        ! enddo 
                    ! else
                        ! do iips = 1,nps
                            ! if (iips == ips) cycle
                            ! dVd_tmp(iips) = ( psd (iips,iz) * (10d0**ps(iips))**2d0 * rough_c0*(10d0**ps(iips))**rough_c1 )
                        ! enddo 
                    ! endif 
                    do iips = 1,nps
                        if (iips == ips) cycle
                        dVd_tmp(iips) = ( psd (iips,iz) * (10d0**ps(iips))**2d0 * rough_tmp(iips) )
                    enddo 
                    
                    if (all(dVd_tmp == 0d0)) then 
                        print *,chrsp,'all dissolved loc3?',ips, psd(:,iz)
                        psd_error_flg = .true.
                        return
                        stop
                    endif 
                    
                    ! dVd_tmp(ips+1:) = dVd_tmp(ips+1:)*dvd_res*dps(ips)/sum(dVd_tmp(ips+1:) * dps(ips+1:))
                    dVd_tmp(:) = dVd_tmp(:)*dvd_res*dps(ips)/sum(dVd_tmp(:) * dps(:))
                    dVd_tmp(ips) = - dvd_res
                    
                    ! dVd(ips+1:,iz) = dVd(ips+1:,iz) + dVd_tmp(ips+1:)
                    ! dVd(ips,iz) = dVd(ips,iz) - dvd_res
                    dVd(:,iz) = dVd(:,iz) + dVd_tmp(:)
                    
                    if ( abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz)) > tol ) then
                        print *,chrsp, ' vol. balance failed somehow loc4 ',abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz))
                        print *, 'going to stop'
                        print *, iz, sum(dVd(:,iz) * dps(:)), DV(iz)
                        stop
                    endif 
                    
                    if (any(isnan(dVd(:,iz))) ) then 
                        print *,chrsp, 'nan in dVd loc 4'
                        stop
                    endif 
                    
                else 
                    ips_new = -1
                    ! new r*
                    ps_new =  ps_new**(1d0/3d0) 
                    if (ps_new <= ps_min) then 
                        ips_new = 1
                    elseif (ps_new >= ps_max) then 
                        ips_new = nps
                    else 
                        do iips = 1, nps -1
                            if ( ( ps_new - 10d0**ps(iips) ) *  ( ps_new - 10d0**ps(iips+1) ) <= 0d0 ) then 
                                if ( log10(ps_new) <= 0.5d0*( ps(iips) + ps(iips+1) ) ) then 
                                    ips_new = iips
                                else 
                                    ips_new = iips + 1 
                                endif 
                                exit 
                            endif 
                        enddo 
                    endif 
                    if (ips_new == -1) then 
                        print *,chrsp,'error: ips_new is not initialized'
                        stop
                    endif 
                    ps_newp = 10d0**ps(ips_new)  ! closest binned particle radius to r*
                    dpsd_tmp(ips_new,iz) = dpsd_tmp(ips_new,iz) + psd(ips,iz)*(ps_new/ps_newp)**3d0
                    dpsd_tmp(ips,iz) =  dpsd_tmp(ips,iz) - psd(ips,iz)
                    ! print *,iz,ips,dVd(ips,iz), psd(ips,iz)* 4d0/3d0 * pi * (10d0**ps(ips) - 10d0**ps(ips_new))**3d0 
                endif 
            
            endif 
            ! print *, iz, ips,ps_new,ps_newp
        enddo 
        ! enddo 

        if (any(isnan(psd))) then 
            print *,chrsp, 'nan in psd'
            stop
        endif 
        if (any(psd<0d0)) then 
            print *,chrsp, 'negative psd'
            stop
        endif 
        if (any(isnan(dvd))) then 
            print *,chrsp, 'nan in dvd' 
            ! do iz = 1, nz
            do ips=1,nps
                if (isnan(dvd(ips,iz))) then 
                    print *,chrsp, 'ips,iz,dvd,psd',ips,iz,dvd(ips,iz),psd(ips,iz)
                endif 
            enddo 
            ! enddo 
            stop
        endif 

        psd_new = psd_new + dpsd_tmp
        ! do iz = 1, nz
        ! if ( abs(DV(iz)) > tol  &
            ! & .and. abs ( ( sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) &
            ! & - sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) - DV(iz) ) / DV(iz) ) > tol ) then  
        if ( abs(DV(iz)) > tol  &
            & .and. abs ( ( - sum( dpsd_tmp(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) &
            &  - DV(iz) ) / DV(iz) ) > tol ) then  
            print *, chrsp,'checking the vol. balance and failed ... ' &
                ! & , abs ( ( sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) &
                ! & - sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) - DV(iz) ) / DV(iz) )
                & , abs ( ( - sum( dpsd_tmp(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) &
                &  - DV(iz) ) / DV(iz) )
            print *, iz, sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) &
                & ,sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) & 
                & ,sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:)) &
                &   -sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0**ps(:))**3d0 * dps(:))  &
                & ,DV(iz)
            ! stop 
            ! pause
            psd_error_flg = .true.
        endif 
        ! enddo 

        if (any(isnan(dpsd_tmp))) then 
            print *,chrsp, 'nan in dpsd _rxn' 
            ! do iz = 1, nz
            do ips=1,nps
                if (isnan(dpsd_tmp(ips,iz))) then 
                    print *, 'ips,iz,dpsd_tmp,psd',ips,iz,dpsd_tmp(ips,iz),psd(ips,iz)
                endif 
            enddo 
            ! enddo 
            stop
        endif 

        dpsd = dpsd + dpsd_tmp
    
    endsubroutine psd_diss_iz

end module scepter_psd 