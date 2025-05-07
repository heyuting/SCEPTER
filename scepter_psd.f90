module scepter_psd
    use scepter_constants
    use scepter_variables
    implicit none

contains

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

    subroutine psd_diss( &
        & nz,nps &! in
        & ,z,DV,dt,pi,tol,poro &! in 
        & ,incld_rough,rough_c0,rough_c1 &! in
        & ,psd,ps,dps,ps_min,ps_max &! in 
        & ,chrsp &! in 
        & ,dpsd,psd_error_flg &! inout
        & )
        implicit none 

        integer,intent(in)::nz,nps
        real(kind=8),intent(in)::dt,ps_min,ps_max,pi,tol,rough_c0,rough_c1 
        real(kind=8),dimension(nz),intent(in)::z,poro
        real(kind=8),dimension(nps),intent(in)::ps,dps
        real(kind=8),dimension(nz),intent(in)::DV
        real(kind=8),dimension(nps,nz),intent(in)::psd
        character(5),intent(in)::chrsp
        logical,intent(in)::incld_rough
        real(kind=8),dimension(nps,nz),intent(inout)::dpsd
        logical,intent(inout)::psd_error_flg
        ! local 
        real(kind=8),dimension(nps,nz)::dVd,psd_old,psd_new,dpsd_tmp
        real(kind=8),dimension(nps)::psd_tmp,dvd_tmp
        real(kind=8) ps_new,ps_newp,dvd_res
        integer ips,iips,ips_new,iz,isps
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
        do iz=1,nz

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
            if (.not.incld_rough) then 
                dVd(:,iz) = ( psd (:,iz) * (10d0**ps(:))**2d0 )
            else
                dVd(:,iz) = ( psd (:,iz) * (10d0**ps(:))**2d0 *rough_c0*(10d0**ps(:))**rough_c1)
            endif 
            
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
                            if (.not.incld_rough) then 
                                dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 )
                            else
                                dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 *rough_c0*(10d0**ps(ips+1:))**rough_c1 )
                            endif 
                            
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
                        if (.not.incld_rough) then 
                            ! dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 )
                            do iips=1,nps
                                if (iips == ips) cycle
                                dVd_tmp(iips) = ( psd (iips,iz) * (10d0**ps(iips))**2d0 )
                            enddo 
                        else
                            ! dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0**ps(ips+1:))**2d0 * rough_c0*(10d0**ps(ips+1:))**rough_c1 )
                            do iips = 1,nps
                                if (iips == ips) cycle
                                dVd_tmp(iips) = ( psd (iips,iz) * (10d0**ps(iips))**2d0 * rough_c0*(10d0**ps(iips))**rough_c1 )
                            enddo 
                        endif 
                        
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
                        ps_newp = 10d0**ps(ips_new)  ! closest binned particle radius to r*
                        dpsd_tmp(ips_new,iz) = dpsd_tmp(ips_new,iz) + psd(ips,iz)*(ps_new/ps_newp)**3d0
                        dpsd_tmp(ips,iz) =  dpsd_tmp(ips,iz) - psd(ips,iz)
                        ! print *,iz,ips,dVd(ips,iz), psd(ips,iz)* 4d0/3d0 * pi * (10d0**ps(ips) - 10d0**ps(ips_new))**3d0 
                    endif 
                
                endif 
                ! print *, iz, ips,ps_new,ps_newp
            enddo 
        enddo 

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
            do iz = 1, nz
                do ips=1,nps
                    if (isnan(dvd(ips,iz))) then 
                        print *,chrsp, 'ips,iz,dvd,psd',ips,iz,dvd(ips,iz),psd(ips,iz)
                    endif 
                enddo 
            enddo 
            stop
        endif 

        psd_new = psd_new + dpsd_tmp
        do iz = 1, nz
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
        enddo 

        if (any(isnan(dpsd_tmp))) then 
            print *,chrsp, 'nan in dpsd _rxn' 
            do iz = 1, nz
                do ips=1,nps
                    if (isnan(dpsd_tmp(ips,iz))) then 
                        print *, 'ips,iz,dpsd_tmp,psd',ips,iz,dpsd_tmp(ips,iz),psd(ips,iz)
                    endif 
                enddo 
            enddo 
            stop
        endif 

        dpsd = dpsd + dpsd_tmp
    
    endsubroutine psd_diss

    subroutine psd_diss_iz( &
        & nz,nps,iz &! in
        & ,z,DV,dt,pi,tol,poro &! in 
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
        real(kind=8),dimension(nps)::psd_tmp,dvd_tmp,rough_tmp
        real(kind=8) ps_new,ps_newp,dvd_res
        integer ips,iips,ips_new,isps
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

    subroutine psd_diss_pbe( &
        & nz,nps &! in
        & ,z,DV,dt,pi,tol,poro &! in 
        & ,incld_rough,roughref &! in
        & ,psd,ps,dps,ps_min,ps_max &! in 
        & ,chrsp &! in 
        & ,dpsd,psd_error_flg &! inout
        & )
        ! an attempt to solve population balance equation reflecting imposed dissolution rate
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
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) * ( lambda(iips+1)*psdx(iips+1,iz) - lambda(iips)*psdx(iips,iz) ) &
                                & *dt/dpsx(iips)
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
                            psdxx(iips,iz) = psdx(iips,iz) + kpsdx(iz) * ( lambda(iips)*psdx(iips,iz) - lambda(iips-1)*psdx(iips-1,iz) ) &
                                & *dt/dpsx(iips)
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

    subroutine psd_diss_pbe_exp( &
        & nz,nps &! in
        & ,z,DV,dt,pi,tol,poro &! in 
        & ,incld_rough,rough_c0,rough_c1 &! in
        & ,psd,ps,dps,ps_min,ps_max &! in 
        & ,chrsp &! in 
        & ,dpsd,psd_error_flg &! inout
        & )
        ! an attempt to solve population balance equation reflecting imposed dissolution rate
        implicit none 

        integer,intent(in)::nz,nps
        real(kind=8),intent(in)::dt,ps_min,ps_max,pi,tol,rough_c0,rough_c1 
        real(kind=8),dimension(nz),intent(in)::z,poro
        real(kind=8),dimension(nps),intent(in)::ps,dps
        real(kind=8),dimension(nz),intent(in)::DV
        real(kind=8),dimension(nps,nz),intent(in)::psd
        character(5),intent(in)::chrsp
        logical,intent(in)::incld_rough
        real(kind=8),dimension(nps,nz),intent(inout)::dpsd
        logical,intent(inout)::psd_error_flg
        ! local 
        real(kind=8),dimension(nps,nz)::dVd,psd_old,psd_new,dpsd_tmp,psdx,psdxx
        real(kind=8),dimension(nps)::psd_tmp,dvd_tmp,dpsx,lambda
        real(kind=8),dimension(nz)::kpsd,kpsdx,DV_chk,kpsdxx
        real(kind=8) ps_new,ps_newp,dvd_res,error,vol,fact,surf
        real(kind=8),parameter::infinity = huge(0d0)
        real(kind=8),parameter::threshold = 20d0
        real(kind=8),parameter::corr = exp(threshold)
        integer,parameter :: iter_max = 50
        integer ips,iips,ips_new,iz,isps,row,col,ie,ie2,iter
        logical :: safe_mode = .false.

        real(kind=8) amx3(nps,nps),ymx3(nps),emx3(nps)
        integer ipiv3(nps) 
        integer info 
                

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


        ! roughness factor as functin of radius
        lambda = 1d0
        if (incld_rough)  lambda = rough_c0*(10d0**ps(:))**rough_c1 

        ! R = log10 r 
        ! dR/dr =  1/(r * log10)
        ! dr = dR * r * log10
        ! note that dps is dR; we need dr denoted here as dpsx (?)
        dpsx = dps * 10d0**(ps) * log(10d0)

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
                ! print *,' *** inconsistent between kpsd and DV ' , chrsld
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
        kpsdxx = kpsd
        psdxx = psdx

        do iz = 1, nz

            error = 1d4 
            iter = 1
            
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
                                & + (  -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                ! & - kpsdx(iz) * ( - lambda(ips)*psdxx(ips,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( - lambda(ips)*1d0 )  &
                                & ) 
                        else
                            ymx3(row) = ( & 
                                & + (  -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                ! & - kpsdx(iz) * ( lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( - lambda(ips)*1d0 )  &
                                & ) 
                                
                            col = ips + 1
                            amx3(row,col) = ( & 
                                & - kpsdx(iz) * ( lambda(ips+1)*1d0 )  &
                                & ) 
                        endif
                    
                    else
                    
                        if (ips==1) then
                            ymx3(row) = ( & 
                                & + (  -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                ! & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips)*1d0 )  &
                                & ) 
                            
                        elseif (ips==nps) then
                            ymx3(row) = ( & 
                                & + (  -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                ! & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                ! & - kpsdx(iz) * ( lambda(ips)*1d0 )  &
                                & ) 
                                
                            col = ips - 1
                            amx3(row,col) = ( & 
                                & - kpsdx(iz) * ( - lambda(ips-1)*1d0 )  &
                                & ) 
                            
                        else
                            ymx3(row) = ( & 
                                & + (  -  psdx(ips,iz) ) /dt * dpsx(ips) &
                                ! & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )  &
                                & )
                            amx3(row,row) = ( & 
                                & + ( 1d0 ) /dt * dpsx(ips) &
                                & - kpsdx(iz) * ( lambda(ips)*1d0 )  &
                                & ) 
                                
                            col = ips - 1
                            amx3(row,col) = ( & 
                                & - kpsdx(iz) * ( - lambda(ips-1)*1d0 )  &
                                & ) 
                        endif
                        
                    endif 
                    
                    ! fact = max( abs( ymx3(row) ), maxval( abs( amx3(row,:) ) ) )
                    
                    ! ymx3(row) = ymx3(row) / fact
                    ! amx3(row,:) = amx3(row,:) / fact
                    
                enddo ! end of do-loop for ips
            
                ymx3=-1.0d0*ymx3

                if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
                ! if (.true.) then 
                    print*,'PBE--'//chrsp//': error in mtx'
                    print*,'PBE--'//chrsp//': any(isnan(amx3)),any(isnan(ymx3))'
                    print*,any(isnan(amx3)),any(isnan(ymx3))

                    if (any(isnan(ymx3))) then 
                        do ie = 1, nps+1
                            if (isnan(ymx3(ie))) then 
                                print*,'NAN is here...',iz,ie
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
                    stop
                    
                endif

                call DGESV(nps,int(1),amx3,nps,IPIV3,ymx3,nps,INFO) 

                if (any(isnan(ymx3)) .or. info/=0 ) then
                    print*,'PBE--'//chrsp//': error in soultion',any(isnan(ymx3)),info
                    psd_error_flg = .true.
                    ! pause
                    exit
                endif
                
                psdxx(:,iz) = ymx3
                
                kpsdxx(iz) = 0d0
                ! if (dV(iz)>=0d0) then
                if (kpsdx(iz)>=0d0) then
                    do ips=1,nps
                        if (ips==nps) then
                            kpsdxx(iz)  = kpsdxx(iz) - ( - lambda(ips)*psdxx(ips,iz) ) &
                                & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        else
                            kpsdxx(iz)  = kpsdxx(iz) - ( lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) ) &
                                & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        endif
                    enddo 
                else 
                    do ips=1,nps
                        if (ips==1) then
                            kpsdxx(iz)  = kpsdxx(iz) - ( lambda(ips)*psdxx(ips,iz) ) &
                                & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        elseif (ips==nps) then
                            kpsdxx(iz)  = kpsdxx(iz) - ( - lambda(ips-1)*psdxx(ips-1,iz) ) &
                                & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        else
                            kpsdxx(iz)  = kpsdxx(iz) - ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) ) &
                                & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        endif
                    enddo 
                endif

                kpsdxx(iz) = dV(iz)/dt/ kpsdxx(iz) 
                
                error = 0
                if (abs( kpsdx(iz) ) > 1d-12) error = abs( (kpsdxx(iz)-kpsdx(iz))/kpsdx(iz) )
                
                kpsdx(iz) = kpsdxx(iz)

                if ( isnan(error) .or. any(isnan(psdxx)) ) then 
                    error = 1d3
                    print*, 'PBE--'//chrsp//': !! error is NaN; values are returned to those before iteration with reducing dt'
                    print*, 'PBE--'//chrsp//': isnan(error), info/=0,any(isnan(pdsx))'
                    print*, isnan(error), any(isnan(psdx)) 
                    
                    psd_error_flg = .true.
                    ! stop
                    exit
                endif

                print '(a,E11.3,a,i0,a,E11.3)', 'PBE--'//chrsp//': iteration error = ',error, ', iteration = ',iter,', time step [yr] = ',dt
                ! print *, error > tol
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
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( - lambda(ips-1)*psdxx(ips-1,iz) ) &
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

    endsubroutine psd_diss_pbe_exp

    subroutine psd_diss_pbe_expall( &
        & nz,nps &! in
        & ,z,DV,dt,pi,tol,poro &! in 
        & ,incld_rough,rough_c0,rough_c1 &! in
        & ,psd,ps,dps,ps_min,ps_max &! in 
        & ,chrsp &! in 
        & ,dpsd,psd_error_flg &! inout
        & )
        ! an attempt to solve population balance equation reflecting imposed dissolution rate
        implicit none 

        integer,intent(in)::nz,nps
        real(kind=8),intent(in)::dt,ps_min,ps_max,pi,tol,rough_c0,rough_c1 
        real(kind=8),dimension(nz),intent(in)::z,poro
        real(kind=8),dimension(nps),intent(in)::ps,dps
        real(kind=8),dimension(nz),intent(in)::DV
        real(kind=8),dimension(nps,nz),intent(in)::psd
        character(5),intent(in)::chrsp
        logical,intent(in)::incld_rough
        real(kind=8),dimension(nps,nz),intent(inout)::dpsd
        logical,intent(inout)::psd_error_flg
        ! local 
        real(kind=8),dimension(nps,nz)::dVd,psd_old,psd_new,dpsd_tmp,psdx,psdxx
        real(kind=8),dimension(nps)::psd_tmp,dvd_tmp,dpsx,lambda
        real(kind=8),dimension(nz)::kpsd,kpsdx,DV_chk
        real(kind=8) ps_new,ps_newp,dvd_res,error,vol,fact,surf
        real(kind=8),parameter::infinity = huge(0d0)
        real(kind=8),parameter::threshold = 20d0
        real(kind=8),parameter::corr = exp(threshold)
        integer,parameter :: iter_max = 50
        integer ips,iips,ips_new,iz,isps,row,col,ie,ie2,iter
        logical :: safe_mode = .false.

        real(kind=8) amx3(nps+1,nps+1),ymx3(nps+1),emx3(nps+1)
        integer ipiv3(nps+1) 
        integer info 
                

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


        ! roughness factor as functin of radius
        lambda = 1d0
        if (incld_rough)  lambda = rough_c0*(10d0**ps(:))**rough_c1 

        ! R = log10 r 
        ! dR/dr =  1/(r * log10)
        ! dr = dR * r * log10
        ! note that dps is dR; we need dr denoted here as dpsx (?)
        dpsx = dps * 10d0**(ps) * log(10d0)

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
            ! kpsdx(iz) = kpsd(iz)
            ! do while (kpsdx(iz) * dV(iz) < 0d0) 
                ! kpsd = 0d0
                ! if (kpsdx(iz)>=0d0) then
                    ! do ips=1,nps
                        ! if (ips==nps) then
                            ! kpsd(iz)  = kpsd(iz) - ( - lambda(ips)*psdx(ips,iz) ) &
                                ! & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        ! else
                            ! kpsd(iz)  = kpsd(iz) - ( lambda(ips+1)*psdx(ips+1,iz) - lambda(ips)*psdx(ips,iz) ) &
                                ! & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        ! endif
                    ! enddo 
                ! else 
                    ! do ips=1,nps
                        ! if (ips==1) then
                            ! kpsd(iz)  = kpsd(iz) - ( lambda(ips)*psdx(ips,iz) ) &
                                ! & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        ! elseif (ips==nps) then 
                            ! kpsd(iz)  = kpsd(iz) - ( - lambda(ips-1)*psdx(ips-1,iz) ) &
                                ! & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        ! else
                            ! kpsd(iz)  = kpsd(iz) - ( lambda(ips)*psdx(ips,iz) - lambda(ips-1)*psdx(ips-1,iz) ) &
                                ! & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                        ! endif
                    ! enddo 
                ! endif
                ! kpsdx(iz) = kpsd(iz)
            ! enddo 
            ! if (kpsd(iz) * dV(iz) < 0d0) then 
                ! print *,' *** inconsistent between kpsd and DV ' , chrsld
                ! print *,' *** maybe because previously dissolving phase starts to precipitate or vice versa ...?'
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
            
            do ips = 1, nps
            
                row =  ips 
                
                vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                surf  = 4d0*pi*(10d0**ps(ips))**2d0
                
                if ( dV(iz) >= 0d0) then 
                ! if ( kpsdx(iz) >= 0d0) then 
                ! if ( .true.) then 
                
                    if (ips==nps) then 
                        ! & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & - kpsdx(iz) * ( - lambda(ips)*psdxx(ips,iz) )  &
                        psdxx(ips,iz) = psdx(ips,iz) + kpsdx(iz) * ( - lambda(ips)*psdx(ips,iz) ) *dt/dpsx(ips)
                    else
                        ! & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & - kpsdx(iz) * ( lambda(ips+1)*psdxx(ips+1,iz) - lambda(ips)*psdxx(ips,iz) )  &
                        psdxx(ips,iz) = psdx(ips,iz) + kpsdx(iz) * ( lambda(ips+1)*psdx(ips+1,iz) - lambda(ips)*psdx(ips,iz) ) &
                            & *dt/dpsx(ips)
                    endif
                
                else
                
                    if (ips==1) then
                        ! & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) )  &
                        psdxx(ips,iz) = psdx(ips,iz) + kpsdx(iz) * ( lambda(ips)*psdx(ips,iz) )*dt/dpsx(ips)
                    elseif (ips==nps) then 
                        ! & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )  &
                        psdxx(ips,iz) = psdx(ips,iz) + kpsdx(iz) * ( - lambda(ips-1)*psdx(ips-1,iz) ) &
                            & *dt/dpsx(ips)
                    else
                        ! & + ( psdxx(ips,iz) -  psdx(ips,iz) ) /dt * dpsx(ips) &
                        ! & - kpsdx(iz) * ( lambda(ips)*psdxx(ips,iz) - lambda(ips-1)*psdxx(ips-1,iz) )  &
                        psdxx(ips,iz) = psdx(ips,iz) + kpsdx(iz) * ( lambda(ips)*psdx(ips,iz) - lambda(ips-1)*psdx(ips-1,iz) ) &
                            & *dt/dpsx(ips)
                    endif
                    
                endif 
                
            enddo ! end of do-loop for ips
            
            dpsd(:,iz) = psdxx(:,iz)*dpsx(:)/dps(:) - psd(:,iz)
            ! dpsd(:,iz) = -dpsd(:,iz)

            DV_chk(iz) = 0d0
            if (dV(iz)>=0d0) then
            ! if (kpsdx(iz)>=0d0) then
                do ips=1,nps
                    if (ips==nps) then
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( - lambda(ips)*psdx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    else
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( lambda(ips+1)*psdx(ips+1,iz) - lambda(ips)*psdx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    endif
                enddo 
            else 
                do ips=1,nps
                    if (ips==1) then
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( lambda(ips)*psdx(ips,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    elseif (ips==nps) then
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( - lambda(ips-1)*psdx(ips-1,iz) ) &
                            & * 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    else
                        DV_chk(iz) = DV_chk(iz) - kpsdx(iz) * ( lambda(ips)*psdx(ips,iz) - lambda(ips-1)*psdx(ips-1,iz) ) &
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

    endsubroutine psd_diss_pbe_expall

    subroutine psd_implicit_all_v2( &
        & nz,nsp_sld,nps,nflx_psd &! in
        & ,z,dz,dt,pi,tol,w0,w,poro,poroi,poroprev &! in
        & ,trans &! in
        & ,psd,psd_pr,ps,dps,dpsd,psd_rain &! in   
        & ,chrsp &! in 
        & ,flgback,flx_max_max &! inout
        & ,psdx,flx_psd &! out
        & )
        implicit none 

        integer,intent(in)::nz,nsp_sld,nps,nflx_psd
        real(kind=8),intent(in)::dt,pi,tol,w0,poroi
        real(kind=8),dimension(nz),intent(in)::z,dz,w,poro,poroprev
        real(kind=8),dimension(nz,nz,nsp_sld),intent(in)::trans
        real(kind=8),dimension(nps),intent(in)::ps,dps,psd_pr
        logical,intent(inout)::flgback
        real(kind=8),intent(inout)::flx_max_max
        character(5),intent(in)::chrsp
        real(kind=8),dimension(nps,nz),intent(in)::psd,dpsd,psd_rain
        real(kind=8),dimension(nps,nz),intent(out)::psdx
        ! local 
        real(kind=8),dimension(nps,nz)::psd_old,dpsd_tmp
        real(kind=8),dimension(nz)::DV,kpsd,sporo
        integer iz,isps,ips,iiz,row,col,ie,ie2,iips
        real(kind=8) vol,surf,m_tmp,mp_tmp,mi_tmp,mprev_tmp,rxn_tmp,drxn_tmp,w_tmp,wp_tmp,trans_tmp,msupp_tmp  &
            & ,sporo_tmp, sporop_tmp,sporoprev_tmp,dtinv,dzinv

        logical::dt_norm = .true.
        real(kind=8),parameter::infinity = huge(0d0)
        real(kind=8),parameter::threshold = 20d0
        ! real(kind=8),parameter::threshold = 3d0
        ! real(kind=8),parameter::corr = 1.5d0
        real(kind=8),parameter::corr = exp(threshold)
        integer,parameter :: iter_max = 50
        ! integer,parameter :: nflx_psd = 6
        real(kind=8) error,fact,flx_max! ,flx_max_max
        real(kind=8) :: flx_tol = 1d-3
        ! real(kind=8) :: flx_tol = 1d-4
        real(kind=8) :: flx_max_tol = 1d-6
        ! real(kind=8) :: flx_max_tol = 1d-5
        real(kind=8) :: dt_th = 1d-6
        real(kind=8) :: fact_tol(nps) 
        integer iter,iflx
        real(kind=8),dimension(nps,nflx_psd,nz),intent(out) :: flx_psd ! itflx,iadv,idif,irain,irxn,ires
        integer  itflx_psd,iadv_psd,idif_psd,irain_psd,irxn_psd,ires_psd
        data itflx_psd,iadv_psd,idif_psd,irain_psd,irxn_psd,ires_psd/1,2,3,4,5,6/
        character(5),dimension(nflx_psd)::chrflx_psd
        character(20) chrfmt
        ! logical :: chkflx = .false.
        logical :: chkflx = .true.

        real(kind=8) amx3(nz,nz),ymx3(nz),emx3(nps),emx3_loc(nz)
        integer ipiv3(nz) 
        integer info 

        external DGESV
                
        ! do iz=1,nz
            ! hr(iz) = sum( 4d0*pi*(10d0**ps(:))**2d0*psd(:,iz)*dps(:) )
        ! enddo 
        ! do iz=1,nz
            ! hr(iz) = sum( 4d0*pi*(10d0**ps(:))**2d0*rough_c0*(10d0**ps(:))**rough_c1*psd(:,iz)*dps(:) )
        ! enddo 
        ! 
        ! attempt to do psd ( defined with particle number / bulk m3 / log (r) )
        ! solve as for msld 
        ! one of particle size equations are used to give massbalance constraint  


        chrflx_psd = (/'tflx ','adv  ','dif  ','rain ','rxn  ','res  '/)

        dtinv = 1d0/dt

        sporo = 1d0 - poro
        sporo = 1d0 

        psdx = psd

        kpsd = 0d0
        ! do iz = 1, nz
            ! do isps = 1,nsp_sld 
                ! DV(iz) = DV(iz) + flx_sld(isps, 4 + isps,iz)*mv(isps)*1d-6
            ! enddo 
            ! kpsd(iz) = DV(iz) / hr(iz) !/ sum( msldx(:,iz) * mv(:) * 1d-6 )
        ! enddo 


        error = 1d4 
        iter = 1
        emx3 = error

        do ips =1,nps
            fact_tol(ips) = maxval(psd(ips,:)) * 1d-12
            ! fact_tol(ips) = maxval(psd(ips,:)) * 1d-15
        enddo 
        ! fact_tol = 1d0

        ! do while (error > tol*fact_tol) 
        do while (error > 1d0) 
        ! do while ( any (emx3 > fact_tol )  ) 
            
            
            do ips = 1, nps
            
                if (emx3(ips) <= fact_tol(ips)) cycle

                amx3 = 0d0
                ymx3 = 0d0
                
                do iz = 1, nz
                
                    row =  iz   
                    
                    vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    surf = 4d0*pi*(10d0**ps(ips))**2d0
                            
                    m_tmp = vol * psdx(ips,iz) * dps(ips)
                    mprev_tmp = vol * psd(ips,iz) * dps(ips)        
                    ! rxn_tmp = vol * psdx(ips,iz)*dps(ips) &
                        ! & * surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                    ! rxn_tmp = surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                    ! rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) / dt 
                    rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) * dtinv
                    
                    ! msupp_tmp = vol * psd_rain(ips,iz) * dps(ips)  / dt
                    msupp_tmp = vol * psd_rain(ips,iz) * dps(ips) * dtinv
                    
                    mi_tmp = vol * psd_pr(ips) * dps(ips)
                    mp_tmp = vol * psdx(ips,min(iz+1,nz)) * dps(ips)
                    
                    dzinv = 1d0/dz(iz)
                    
                    ! drxn_tmp = & 
                        ! & vol * 1d0 * dps(ips) &
                        ! & * surf * psdx(ips,iz) * dps(ips) * kpsd(iz) &
                        ! & + vol * psdx(ips,iz) * dps(ips) &
                        ! & * su
                    ! drxn_tmp = surf * 1d0 * dps(ips) * kpsd(iz) 
                    drxn_tmp = 0d0 
                    
                    w_tmp = w(iz)
                    wp_tmp = w(min(nz,iz+1))

                    sporo_tmp = 1d0-poro(iz)
                    sporop_tmp = 1d0-poro(min(nz,iz+1))
                    sporoprev_tmp = 1d0-poroprev(iz)
                    
                    if (iz==nz) then 
                        mp_tmp = mi_tmp
                        wp_tmp = w0
                        sporop_tmp = 1d0- poroi
                    endif 
                    
                    sporo_tmp = 1d0
                    sporop_tmp = 1d0
                    sporoprev_tmp = 1d0

                    amx3(row,row) = ( &
                        ! & sporo_tmp * vol * dps(ips) * 1d0  /  merge(1d0,dt,dt_norm)     &
                        & sporo_tmp * vol * dps(ips) * 1d0  *  merge(1d0,dtinv,dt_norm)     &
                        ! & + sporo_tmp * vol * dps(ips) * w_tmp / dz(iz)  *merge(dt,1d0,dt_norm)    &
                        & + sporo_tmp * vol * dps(ips) * w_tmp * dzinv  *merge(dt,1d0,dt_norm)    &
                        & + sporo_tmp * drxn_tmp * merge(dt,1d0,dt_norm) &
                        & ) &
                        & * psdx(ips,iz)

                    ymx3(row) = ( &
                        ! & ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) / merge(1d0,dt,dt_norm) &
                        & ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) * merge(1d0,dtinv,dt_norm) &
                        ! & -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) / dz(iz) * merge(dt,1d0,dt_norm)  &
                        & -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) * dzinv * merge(dt,1d0,dt_norm)  &
                        & + sporo_tmp* rxn_tmp * merge(dt,1d0,dt_norm) &
                        & - sporo_tmp* msupp_tmp * merge(dt,1d0,dt_norm) &
                        & ) &
                        & * 1d0
                                
                    if (iz/=nz) then 
                        col = iz + 1  
                        amx3(row,col) = &
                            ! & (- sporop_tmp * vol * dps(ips) * wp_tmp / dz(iz)) * merge(dt,1d0,dt_norm) * psdx(ips,min(iz+1,nz))
                            & (- sporop_tmp * vol * dps(ips) * wp_tmp * dzinv) * merge(dt,1d0,dt_norm) * psdx(ips,min(iz+1,nz))
                    endif 
                    
                    do iiz = 1, nz
                        col = iiz   
                        trans_tmp = sum(trans(iiz,iz,:))/nsp_sld
                        if (trans_tmp == 0d0) cycle
                            
                        amx3(row,col) = amx3(row,col) - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips)  * merge(dt,1d0,dt_norm) 
                        ymx3(row) = ymx3(row) - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips)  * merge(dt,1d0,dt_norm) 
                    enddo
                    
                    
                    ! if (ips == nps) then 
                    
                        ! amx3(row,:) = 0d0
                        ! ymx3(row  ) = 0d0
                        
                        
                        ! ymx3(row) = ( &
                            ! &  sum( msldx(:,iz) * mv(:) * 1d-6 ) - sum ( 4d0/3d0*pi*(10d0**ps(:))**3d0 * psdx(:,iz) * dps(:) ) &
                            ! & ) &
                            ! & * 1d0
                            
                        ! do iips = 1, nps
                            ! col = iz  + ( iips - 1 ) * nz 
                            ! amx3(row,col) = ( &
                                ! &   -  ( 4d0/3d0*pi*(10d0**ps(iips))**3d0 * psdx(iips,iz) * dps(iips) )  &
                                ! & ) &
                                ! & * 1d0
                        
                        ! enddo 
                    
                    ! endif 
                    
                    ! fact = max( abs( ymx3(row) ), maxval( abs( amx3(row,:) ) ) )
                    
                    ! ymx3(row) = ymx3(row) / fact
                    ! amx3(row,:) = amx3(row,:) / fact
                    
                enddo
            
                ymx3=-1.0d0*ymx3

                if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
                ! if (.true.) then 
                    print*,'PSD--'//chrsp//': error in mtx'
                    print*,'PSD--'//chrsp//': any(isnan(amx3)),any(isnan(ymx3))'
                    print*,any(isnan(amx3)),any(isnan(ymx3))

                    if (any(isnan(ymx3))) then 
                        do iz = 1, nz
                            if (isnan(ymx3(iz))) then 
                                print*,'NAN is here...',ips,iz
                            endif
                        enddo 
                    endif


                    if (any(isnan(amx3))) then 
                        do ie = 1,(nz)
                            do ie2 = 1,(nz)
                                if (isnan(amx3(ie,ie2))) then 
                                    print*,'PSD: NAN is here...',ips,ie,ie2
                                endif
                            enddo
                        enddo
                    endif

                    open(unit=11,file='amx.txt',status = 'replace')
                    open(unit=12,file='ymx.txt',status = 'replace')
                    do ie = 1,nz
                        write(11,*) (amx3(ie,ie2),ie2 = 1,nz)
                        write(12,*) ymx3(ie)
                    enddo 
                    close(11)
                    close(12)        
                    
                    stop
                    
                endif
            
                call DGESV(Nz,int(1),amx3,Nz,IPIV3,ymx3,Nz,INFO) 
            
                if (any(isnan(ymx3)) .or. info/=0 ) then
                    print*,'PSD--'//chrsp//': error in soultion',any(isnan(ymx3)),info
                    flgback = .true.
                    ! pause
                    exit
                endif
            
                do iz = 1, nz
                    
                    row =  iz   

                    if (isnan(ymx3(row))) then 
                        print *,'PSD--'//chrsp//': nan at', iz,ips
                        stop
                    endif
                    
                    emx3_loc(row) = dps(ips)*psdx(ips,iz)*exp(ymx3(row)) - dps(ips)*psdx(ips,iz)
                    
                    if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                        psdx(ips,iz) = psdx(ips,iz)*corr
                    else if (ymx3(row) < -threshold) then 
                        psdx(ips,iz) = psdx(ips,iz)/corr
                    else   
                        psdx(ips,iz) = psdx(ips,iz)*exp(ymx3(row))
                    endif
                enddo 

                if (all(fact_tol == 1d0)) then 
                    emx3(ips) = maxval(exp(abs(ymx3))) - 1.0d0
                else 
                    emx3(ips) = maxval(abs(emx3_loc))
                endif 
            
            enddo 
            
            ! error = maxval(emx3)
            error = maxval(emx3/fact_tol)
            
            ! if (isnan(error)) error = 1d4

            if ( isnan(error) .or. any(isnan(psdx)) ) then 
                error = 1d3
                print*, 'PSD--'//chrsp//': !! error is NaN; values are returned to those before iteration with reducing dt'
                print*, 'PSD--'//chrsp//': isnan(error), info/=0,any(isnan(pdsx))'
                print*, isnan(error), any(isnan(psdx)) 
                
                flgback = .true.
                stop
                exit
            endif
        #ifdef show_PSDiter
            print '(a,E11.3,a,i0,a,E11.3)', 'PSD--'//chrsp//': iteration error = ',error, ', iteration = ',iter,', time step [yr] = ',dt
        #endif 
            iter = iter + 1 
            
            if (iter > iter_Max ) then
                ! dt = dt/1.01d0
                ! dt = dt/10d0
                if (dt==0d0) then 
                    print *, chrsp,'dt==0d0; stop'
                
        ! #ifdef errmtx_printout
                    ! open(unit=11,file='amx.txt',status = 'replace')
                    ! open(unit=12,file='ymx.txt',status = 'replace')
                    ! do ie = 1,nsp3*(nz)
                        ! write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                        ! write(12,*) ymx3(ie)
                    ! enddo 
                    ! close(11)
                    ! close(12)      
        ! #endif 
                    stop
                endif 
                flgback = .true.
                
                exit 
            end if
            
            if (flgback) exit

        enddo 

        ! calculating flux 
        flx_psd = 0d0
        do ips = 1, nps
            
            do iz = 1, nz
            
                row =  iz   
                
                vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                surf = 4d0*pi*(10d0**ps(ips))**2d0
                        
                m_tmp = vol * psdx(ips,iz) * dps(ips)
                mprev_tmp = vol * psd(ips,iz) * dps(ips)        
                ! rxn_tmp = vol * psdx(ips,iz)*dps(ips) &
                    ! & * surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                ! rxn_tmp = surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) / dt 
                
                msupp_tmp = vol * psd_rain(ips,iz) * dps(ips)  / dt
                
                mi_tmp = vol * psd_pr(ips) * dps(ips)
                mp_tmp = vol * psdx(ips,min(iz+1,nz)) * dps(ips)
                
                dzinv = 1d0/dz(iz)
                
                ! drxn_tmp = & 
                    ! & vol * 1d0 * dps(ips) &
                    ! & * surf * psdx(ips,iz) * dps(ips) * kpsd(iz) &
                    ! & + vol * psdx(ips,iz) * dps(ips) &
                    ! & * su
                ! drxn_tmp = surf * 1d0 * dps(ips) * kpsd(iz) 
                drxn_tmp = 0d0 
                
                w_tmp = w(iz)
                wp_tmp = w(min(nz,iz+1))

                sporo_tmp = 1d0-poro(iz)
                sporop_tmp = 1d0-poro(min(nz,iz+1))
                sporoprev_tmp = 1d0-poroprev(iz)
                
                if (iz==nz) then 
                    mp_tmp = mi_tmp
                    wp_tmp = w0
                    sporop_tmp = 1d0- poroi
                endif 
                
                sporo_tmp = 1d0
                sporop_tmp = 1d0
                sporoprev_tmp = 1d0
                
                flx_psd(ips,itflx_psd,iz) = ( &
                    & ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) * dtinv  &
                    & )
                flx_psd(ips,iadv_psd,iz) = ( &
                    & -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) * dzinv &
                    & )
                flx_psd(ips,irxn_psd,iz) = ( &
                    & + sporo_tmp* rxn_tmp  &
                    & )
                flx_psd(ips,irain_psd,iz) = ( &
                    & - sporo_tmp* msupp_tmp  &
                    & )
                
                do iiz = 1, nz  
                    trans_tmp = sum(trans(iiz,iz,:))/nsp_sld
                    if (trans_tmp == 0d0) cycle
                    
                    flx_psd(ips,idif_psd,iz) = flx_psd(ips,idif_psd,iz) + ( &
                        & - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips) &
                        & )
                enddo
                
                
                flx_psd(ips,ires_psd,iz) = sum(flx_psd(ips,:,iz))
                
            enddo
        enddo
            
        #ifdef dispPSDiter

        write(chrfmt,'(i0)') nflx_psd
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'

        print *
        print *,' [fluxes -- PSD] '
        print *, '<'//chrsp//'>'
        print trim(adjustl(chrfmt)),'rad',(chrflx_psd(iflx),iflx=1,nflx_psd)

        write(chrfmt,'(i0)') nflx_psd
        chrfmt = '(f5.2,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        do ips = 1, nps
            print trim(adjustl(chrfmt)), ps(ips), (sum(flx_psd(ips,iflx,:)*dz(:)),iflx=1,nflx_psd)
        enddo 
        print *

        #endif     
        

                
        if ( chkflx .and. dt > dt_th) then 
            ! flx_max_max = 0d0
            do ips = 1, nps
                flx_max = 0d0
                do iflx=1,nflx_psd 
                    flx_max= max( flx_max, abs( sum(flx_psd(ips,iflx,:)*dz(:)) ) )
                enddo 
                flx_max_max = max( flx_max_max, flx_max)
            enddo 
            do ips = 1, nps
            
                if ( flx_max > flx_max_max*flx_max_tol .and. abs( sum(flx_psd(ips,ires_psd,:)*dz(:)) ) > flx_max * flx_tol ) then 
                    
                    print *, chrsp,' too large error in PSD flx?'
                    print *, 'flx_max, flx_max_max,tol = ', flx_max,flx_max_max,flx_max_tol
                    print *, 'res, max = ', abs( sum(flx_psd(ips,ires_psd,:)*dz(:)) ), flx_max
                    print *, 'res/max, target = ', abs( sum(flx_psd(ips,ires_psd,:)*dz(:)) )/flx_max, flx_tol
                    
                    flgback = .true.
                
                endif 
                
            enddo 
        endif  
    
    endsubroutine psd_implicit_all_v2

    subroutine psd_implicit_all_v4( &
        & nz,nsp_sld,nps,nflx_psd &! in
        & ,z,dz,dt,pi,tol,w0,w,poro,poroi,poroprev &! in 
        & ,incld_rough,roughref &! in
        & ,trans &! in
        & ,psd,psd_pr,ps,dps,dpsd,psd_rain,DV,psd_norm_fact &! in   
        & ,chrsp &! in 
        & ,flgback,flx_max_max &! inout
        & ,psdx,flx_psd &! out
        & )
        implicit none 

        integer,intent(in)::nz,nsp_sld,nps,nflx_psd
        real(kind=8),intent(in)::dt,pi,tol,w0,poroi
        real(kind=8),dimension(nz),intent(in)::z,dz,w,poro,poroprev,DV
        real(kind=8),dimension(nz,nz,nsp_sld),intent(in)::trans
        real(kind=8),dimension(nps),intent(in)::ps,dps,psd_pr,psd_norm_fact
        logical,intent(in)::incld_rough
        logical,intent(inout)::flgback
        real(kind=8),intent(inout)::flx_max_max
        character(5),intent(in)::chrsp
        character(10),intent(in)::roughref
        real(kind=8),dimension(nps,nz),intent(in)::psd,dpsd,psd_rain
        real(kind=8),dimension(nps,nz),intent(out)::psdx
        ! local 
        real(kind=8),dimension(nps,nz)::psd_old,dpsd_tmp
        real(kind=8),dimension(nz)::kpsd,sporo,kpsdx
        real(kind=8),dimension(nps)::rough_tmp
        integer iz,isps,ips,iiz,row,col,ie,ie2,iips
        real(kind=8) vol,surf,m_tmp,mp_tmp,mi_tmp,mprev_tmp,rxn_tmp,drxn_tmp,w_tmp,wp_tmp,trans_tmp,msupp_tmp  &
            & ,sporo_tmp, sporop_tmp,sporoprev_tmp,dtinv,dzinv,drxndk_tmp

        logical::dt_norm = .true.
        real(kind=8),parameter::infinity = huge(0d0)
        real(kind=8),parameter::threshold = 20d0
        ! real(kind=8),parameter::threshold = 3d0
        ! real(kind=8),parameter::corr = 1.5d0
        real(kind=8),parameter::corr = exp(threshold)
        integer,parameter :: iter_max = 50
        ! integer,parameter :: nflx_psd = 6
        real(kind=8) error,fact,flx_max! ,flx_max_max
        real(kind=8) error_k
        real(kind=8) :: flx_tol = 1d-3
        ! real(kind=8) :: flx_tol = 1d-4
        real(kind=8) :: flx_max_tol = 1d-6
        ! real(kind=8) :: flx_max_tol = 1d-5
        real(kind=8) :: dt_th = 1d-6
        real(kind=8) :: k_tol = 1d-4
        ! real(kind=8) :: k_tol = 1d-6
        real(kind=8) :: fact_tol(nps) 
        integer iter,iflx
        integer iter_k
        real(kind=8),dimension(nps,nflx_psd,nz),intent(out) :: flx_psd ! itflx,iadv,idif,irain,irxn,ires
        integer  itflx_psd,iadv_psd,idif_psd,irain_psd,irxn_psd,ires_psd
        data itflx_psd,iadv_psd,idif_psd,irain_psd,irxn_psd,ires_psd/1,2,3,4,5,6/
        character(5),dimension(nflx_psd)::chrflx_psd
        character(20) chrfmt
        ! logical :: chkflx = .false.
        logical :: chkflx = .true.
        logical :: explicit = .false.
        ! logical :: explicit = .true.

        real(kind=8) amx3(nz,nz),ymx3(nz),emx3(nps),emx3_loc(nz)
        integer ipiv3(nz) 
        integer info 
                

        ! attempt to do psd ( defined with particle number / bulk m3 / log (r) )
        ! assumptions/formulations: 
        ! 1. d/dt|diss [4 * pi /3 * r^3 * f(r,t)] = k(r) * 4 * pi * r^2 * f(r,t)  
        !      ( < == >  d/dt [f(r,t)] = k(r) * 3 * r^(-1) * f(r,t) )
        ! 2. integral for r (or summation numerically) :
        !    int [ d/dt|diss [4 * pi /3 * r^3 * f(r,t)] ] dr =  int [ k(r) * 4 * pi * r^2 * f(r,t) ] dr
        ! 3. Both sides must equal dV calculated from `bulk' RTM 
        !    int [ k(r) * 4 * pi * r^2 * f(r,t) ] dr = dV
        !    Attempt here is to obtain k(r) using obtained dV  
        !    This can be easily done if one assumes k(r) is constant and independent r; 
        !    Even if k(r) is dependent on r (e.g., Emmanuel and Ague, 2011) 
        !    one would be able to obtain k(r) as long as r-dependent part is mechanistically known: k(r) = k' * gamma(r) 


        chrflx_psd = (/'tflx ','adv  ','dif  ','rain ','rxn  ','res  '/)

        dtinv = 1d0/dt

        sporo = 1d0 - poro
        sporo = 1d0 

        psdx = psd

        kpsd = 0d0
        ! dV = 0d0
        do iz = 1, nz
            if (.not. incld_rough) then 
                rough_tmp(:) = rough_f( 'smooth', nps, (10d0**ps(:)) )
            else
                rough_tmp(:) = rough_f( roughref, nps, (10d0**ps(:)) )
            endif 
            kpsd(iz) = sum( psd_norm_fact(:) * psd(:,iz) * dps(:) * 4d0*pi*(10d0**ps(:))**2d0 *rough_tmp(:) ) ! available surface area
            kpsd(iz) = dV(iz)/kpsd(iz) 
            ! kpsd(iz) = abs(dV(iz))/kpsd(iz) 
        enddo 

        do ips =1,nps
            fact_tol(ips) = maxval(psd(ips,:)) * 1d-12
            ! fact_tol(ips) = maxval(psd(ips,:)) * 1d-15
        enddo 
        ! fact_tol = 1d0


        ! do loop for kpsd  
        error_k = 1d4
        iter_k = 1
        kpsdx = kpsd
                                do while (error_k > k_tol)


        error = 1d4 
        iter = 1
        emx3 = error

        ! do while (error > tol*1d-3) 
        ! do while (error > tol) 
        do while (error > 1d0) 
        ! do while ( any (emx3 > fact_tol )  ) 
            
            do ips = 1, nps
            
                fact = maxval( psd(ips,:) ) 
            
                if (emx3(ips) <= fact_tol(ips)) cycle

                amx3 = 0d0
                ymx3 = 0d0
                
                do iz = 1, nz
                
                    row =  iz !+ (ips - 1)*nz
                    
                    vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                    surf = 4d0*pi*(10d0**ps(ips))**2d0*rough_tmp(ips)
                            
                    m_tmp = vol * psdx(ips,iz) * dps(ips)
                    mprev_tmp = vol * psd(ips,iz) * dps(ips)        
                    ! rxn_tmp = vol * psdx(ips,iz)*dps(ips) &
                        ! & * surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                    ! rxn_tmp = surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                    ! rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) / dt 
                    ! rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) * dtinv
                    ! rxn_tmp =  - dsign(1d0,dV(iz))*  kpsdx(iz)*surf * psdx(ips,iz) * dps(ips) * dtinv
                    rxn_tmp =  kpsdx(iz)*surf * psdx(ips,iz) * dps(ips) * dtinv
                    
                    ! msupp_tmp = vol * psd_rain(ips,iz) * dps(ips)  / dt
                    msupp_tmp = vol * psd_rain(ips,iz) * dps(ips) * dtinv
                    
                    mi_tmp = vol * psd_pr(ips) * dps(ips)
                    mp_tmp = vol * psdx(ips,min(iz+1,nz)) * dps(ips)
                    
                    dzinv = 1d0/dz(iz)
                    
                    ! drxn_tmp = & 
                        ! & vol * 1d0 * dps(ips) &
                        ! & * surf * psdx(ips,iz) * dps(ips) * kpsd(iz) &
                        ! & + vol * psdx(ips,iz) * dps(ips) &
                        ! & * su
                    ! drxn_tmp = - dsign(1d0,dV(iz))* kpsdx(iz)* surf * 1d0 * dps(ips) *  dtinv
                    drxn_tmp =  kpsdx(iz)* surf * 1d0 * dps(ips) *  dtinv
                    ! drxn_tmp = 0d0 
                    drxndk_tmp =  1d0* surf * psdx(ips,iz) * dps(ips) *  dtinv
                    
                    if (explicit) then 
                        rxn_tmp =  kpsd(iz)*surf * psd(ips,iz) * dps(ips) * dtinv
                        drxn_tmp = 0d0
                        drxndk_tmp = 0d0
                    endif 
                    
                    w_tmp = w(iz)
                    wp_tmp = w(min(nz,iz+1))

                    sporo_tmp = 1d0-poro(iz)
                    sporop_tmp = 1d0-poro(min(nz,iz+1))
                    sporoprev_tmp = 1d0-poroprev(iz)
                    
                    if (iz==nz) then 
                        mp_tmp = mi_tmp
                        wp_tmp = w0
                        sporop_tmp = 1d0- poroi
                    endif 
                    
                    sporo_tmp = 1d0
                    sporop_tmp = 1d0
                    sporoprev_tmp = 1d0

                    amx3(row,row) = ( &
                        ! & sporo_tmp * vol * dps(ips) * 1d0  /  merge(1d0,dt,dt_norm)     &
                        & sporo_tmp * vol * dps(ips) * 1d0  *  merge(1d0,dtinv,dt_norm)     &
                        ! & + sporo_tmp * vol * dps(ips) * w_tmp / dz(iz)  *merge(dt,1d0,dt_norm)    &
                        & + sporo_tmp * vol * dps(ips) * w_tmp * dzinv  *merge(dt,1d0,dt_norm)    &
                        & + sporo_tmp * drxn_tmp * merge(dt,1d0,dt_norm) &
                        & ) &
                        & * psdx(ips,iz)

                    ymx3(row) = ( &
                        ! & ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) / merge(1d0,dt,dt_norm) &
                        & ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) * merge(1d0,dtinv,dt_norm) &
                        ! & -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) / dz(iz) * merge(dt,1d0,dt_norm)  &
                        & -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) * dzinv * merge(dt,1d0,dt_norm)  &
                        & + sporo_tmp* rxn_tmp * merge(dt,1d0,dt_norm) &
                        & - sporo_tmp* msupp_tmp * merge(dt,1d0,dt_norm) &
                        & ) &
                        & * 1d0
                                
                    if (iz/=nz) then 
                        col = iz + 1  !+ (ips - 1)*nz
                        amx3(row,col) = &
                            ! & (- sporop_tmp * vol * dps(ips) * wp_tmp / dz(iz)) * merge(dt,1d0,dt_norm) * psdx(ips,min(iz+1,nz))
                            & (- sporop_tmp * vol * dps(ips) * wp_tmp * dzinv) * merge(dt,1d0,dt_norm) * psdx(ips,min(iz+1,nz))
                    endif 
                    
                    do iiz = 1, nz
                        col = iiz   ! + (ips - 1)*nz
                        trans_tmp = sum(trans(iiz,iz,:))/nsp_sld
                        if (trans_tmp == 0d0) cycle
                            
                        amx3(row,col) = amx3(row,col) - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips)  * merge(dt,1d0,dt_norm) 
                        ymx3(row) = ymx3(row) - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips)  * merge(dt,1d0,dt_norm) 
                    enddo
                    
                    ! col = iz + nz + (nps - 1)*nz
                    ! amx3(row,col) = amx3(row,col) + ( &
                        ! & + sporo_tmp * drxndk_tmp * merge(dt,1d0,dt_norm) &
                        ! & )&
                        ! & * 1d0
                        ! & * kpsdx(iz)
                    
                    
                    ! if (ips == nps) then 
                    
                        ! amx3(row,:) = 0d0
                        ! ymx3(row  ) = 0d0
                        
                        
                        ! ymx3(row) = ( &
                            ! &  sum( msldx(:,iz) * mv(:) * 1d-6 ) - sum ( 4d0/3d0*pi*(10d0**ps(:))**3d0 * psdx(:,iz) * dps(:) ) &
                            ! & ) &
                            ! & * 1d0
                            
                        ! do iips = 1, nps
                            ! col = iz  + ( iips - 1 ) * nz 
                            ! amx3(row,col) = ( &
                                ! &   -  ( 4d0/3d0*pi*(10d0**ps(iips))**3d0 * psdx(iips,iz) * dps(iips) )  &
                                ! & ) &
                                ! & * 1d0
                        
                        ! enddo 
                    
                    ! endif 
                    
                    ! fact = max( abs( ymx3(row) ), maxval( abs( amx3(row,:) ) ) )
                    
                    ! ymx3(row) = ymx3(row) / fact
                    ! amx3(row,:) = amx3(row,:) / fact
                    
                enddo
                    
                ymx3 = ymx3 / fact
                amx3 = amx3 / fact

                ymx3=-1.0d0*ymx3

                if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
                ! if (.true.) then 
                    print*,'PSD--'//chrsp//': error in mtx'
                    print*,'PSD--'//chrsp//': any(isnan(amx3)),any(isnan(ymx3))'
                    print*,any(isnan(amx3)),any(isnan(ymx3))

                    if (any(isnan(ymx3))) then 
                        do iz = 1, nz!*nps+nz
                            if (isnan(ymx3(iz))) then 
                                print*,'NAN is here...',ips,iz
                            endif
                        enddo 
                    endif


                    if (any(isnan(amx3))) then 
                        do ie = 1,nz!*nps+nz)
                            do ie2 = 1,nz!*nps+nz)
                                if (isnan(amx3(ie,ie2))) then 
                                    print*,'PSD: NAN is here...',ips,ie,ie2
                                endif
                            enddo
                        enddo
                    endif
                    stop
                    
                endif

                call DGESV(Nz,int(1),amx3,Nz,IPIV3,ymx3,Nz,INFO) 

                if (any(isnan(ymx3)) .or. info/=0 ) then
                    print*,'PSD--'//chrsp//': error in soultion',any(isnan(ymx3)),info
                    flgback = .true.
                    ! pause
                    exit
                endif
            
                do iz = 1, nz
                    
                    row =  iz   !+ (ips - 1)* nz

                    if (isnan(ymx3(row))) then 
                        print *,'PSD--'//chrsp//': nan at', iz,ips
                        stop
                    endif
                    
                    emx3_loc(row) = dps(ips)*psdx(ips,iz)*exp(ymx3(row)) - dps(ips)*psdx(ips,iz)
                    
                    if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                        psdx(ips,iz) = psdx(ips,iz)*corr
                    else if (ymx3(row) < -threshold) then 
                        psdx(ips,iz) = psdx(ips,iz)/corr
                    else   
                        psdx(ips,iz) = psdx(ips,iz)*exp(ymx3(row))
                    endif
                enddo
                
                if (all(fact_tol == 1d0)) then 
                    emx3(ips) = maxval(exp(abs(ymx3))) - 1.0d0
                else 
                    emx3(ips) = maxval(abs(emx3_loc))
                endif 
            
            enddo 
            
            ! error = maxval(emx3)
            error = maxval(emx3/fact_tol)
            
            ! if (isnan(error)) error = 1d4

            if ( isnan(error) .or. any(isnan(psdx)) ) then 
                error = 1d3
                print*, 'PSD--'//chrsp//': !! error is NaN; values are returned to those before iteration with reducing dt'
                print*, 'PSD--'//chrsp//': isnan(error), info/=0,any(isnan(pdsx))'
                print*, isnan(error), any(isnan(psdx)) 
                
                flgback = .true.
                stop
                exit
            endif

            print '(a,E11.3,a,i0,a,E11.3)', 'PSD--'//chrsp//': iteration error = ',error, ', iteration = ',iter,', time step [yr] = ',dt
            iter = iter + 1 
            
            if (iter > iter_Max ) then
                ! dt = dt/1.01d0
                ! dt = dt/10d0
                if (dt==0d0) then 
                    print *, chrsp,'dt==0d0; stop'
                
        ! #ifdef errmtx_printout
                    ! open(unit=11,file='amx.txt',status = 'replace')
                    ! open(unit=12,file='ymx.txt',status = 'replace')
                    ! do ie = 1,nsp3*(nz)
                        ! write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                        ! write(12,*) ymx3(ie)
                    ! enddo 
                    ! close(11)
                    ! close(12)      
        ! #endif 
                    stop
                endif 
                flgback = .true.
                
                exit 
            end if
            
            if (flgback) exit

        enddo 

        if (explicit) exit

        kpsd = kpsdx

        do iz = 1, nz
            kpsdx(iz) = dV(iz)/sum( psd_norm_fact(:) * psdx(:,iz) * dps(:)* 4d0*pi*(10d0**ps(:))**2d0 *rough_tmp(:) ) ! available surface area
            ! kpsdx(iz) = dV(iz)/kpsdx(iz) 
            ! kpsd(iz) = abs(dV(iz))/kpsd(iz) 
        enddo 

        do iz = 1, nz
            if ( abs(DV(iz)) > tol  &
                & .and. abs ( ( kpsdx(iz) * sum( psd_norm_fact(:) * psdx(:,iz) * dps(:) &
                & * 4d0*pi*(10d0**ps(:))**2d0 *rough_tmp(:) ) &
                &  - DV(iz) ) / DV(iz) ) > tol ) then  
                print *, chrsp,'checking the vol. balance within do-loop and failed ... ' &
                    & , abs ( ( kpsdx(iz) * sum( psd_norm_fact(:) * psdx(:,iz) * dps(:) &
                    & * 4d0*pi*(10d0**ps(:))**2d0 *rough_tmp(:) ) &
                    &  - DV(iz) ) / DV(iz) )
            endif 
        enddo 

        error_k = 0d0
        do iz=1,nz
            if (kpsdx(iz)/=0d0) then 
                error_k = max(error_k, abs((kpsdx(iz) - kpsd(iz))/kpsdx(iz)) )
                ! error_k = max(error_k, abs((kpsdx(iz) - kpsd(iz)))*dtinv )
            endif 
        enddo

        print *, '----- kpsd iteration', iter_k, error_k
        iter_k = iter_k + 1 

        if (iter_k > iter_Max ) then
            ! dt = dt/1.01d0
            ! dt = dt/10d0
            if (dt==0d0) then 
                print *, chrsp,'dt==0d0; stop'
                stop
            endif 
            flgback = .true.
            
            exit 
        end if

        if (flgback) exit

        ! end of do-loop for k
                                enddo

        if (.not.explicit) then 
            do iz = 1, nz
                if ( abs(DV(iz)) > tol  &
                    & .and. abs ( ( kpsdx(iz) * sum( psd_norm_fact(:) * psdx(:,iz) * dps(:) &
                    & * 4d0*pi*(10d0**ps(:))**2d0 *rough_tmp(:) ) &
                    &  - DV(iz) ) / DV(iz) ) > tol ) then  
                    print *, chrsp,'checking the vol. balance and failed ... ' &
                        & , abs ( ( kpsdx(iz) * sum( psd_norm_fact(:) * psdx(:,iz) * dps(:) &
                        & * 4d0*pi*(10d0**ps(:))**2d0 *rough_tmp(:) ) &
                        &  - DV(iz) ) / DV(iz) )
                    flgback = .true.
                endif 
            enddo 
        endif 

        ! calculating flux 
        flx_psd = 0d0
        do ips = 1, nps
            
            do iz = 1, nz
            
                row =  iz   
                
                vol  = 4d0/3d0*pi*(10d0**ps(ips))**3d0
                surf = 4d0*pi*(10d0**ps(ips))**2d0 * rough_tmp(ips)
                        
                m_tmp = vol * psdx(ips,iz) * dps(ips)
                mprev_tmp = vol * psd(ips,iz) * dps(ips)        
                ! rxn_tmp = vol * psdx(ips,iz)*dps(ips) &
                    ! & * surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                ! rxn_tmp = surf * psdx(ips,iz)*dps(ips) * kpsd(iz) 
                ! rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) / dt 
                rxn_tmp =   kpsdx(iz)*surf * psdx(ips,iz) * dps(ips) * dtinv
                    
                if (explicit) then 
                    rxn_tmp =  kpsd(iz)*surf * psd(ips,iz) * dps(ips) * dtinv
                endif 
                
                msupp_tmp = vol * psd_rain(ips,iz) * dps(ips)  / dt
                
                mi_tmp = vol * psd_pr(ips) * dps(ips)
                mp_tmp = vol * psdx(ips,min(iz+1,nz)) * dps(ips)
                
                dzinv = 1d0/dz(iz)
                
                ! drxn_tmp = & 
                    ! & vol * 1d0 * dps(ips) &
                    ! & * surf * psdx(ips,iz) * dps(ips) * kpsd(iz) &
                    ! & + vol * psdx(ips,iz) * dps(ips) &
                    ! & * su
                ! drxn_tmp = surf * 1d0 * dps(ips) * kpsd(iz) 
                drxn_tmp = 0d0 
                
                w_tmp = w(iz)
                wp_tmp = w(min(nz,iz+1))

                sporo_tmp = 1d0-poro(iz)
                sporop_tmp = 1d0-poro(min(nz,iz+1))
                sporoprev_tmp = 1d0-poroprev(iz)
                
                if (iz==nz) then 
                    mp_tmp = mi_tmp
                    wp_tmp = w0
                    sporop_tmp = 1d0- poroi
                endif 
                
                sporo_tmp = 1d0
                sporop_tmp = 1d0
                sporoprev_tmp = 1d0
                
                flx_psd(ips,itflx_psd,iz) = ( &
                    & ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) * dtinv  &
                    & )
                flx_psd(ips,iadv_psd,iz) = ( &
                    & -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) * dzinv &
                    & )
                flx_psd(ips,irxn_psd,iz) = ( &
                    & + sporo_tmp* rxn_tmp  &
                    & )
                flx_psd(ips,irain_psd,iz) = ( &
                    & - sporo_tmp* msupp_tmp  &
                    & )
                
                do iiz = 1, nz  
                    trans_tmp = sum(trans(iiz,iz,:))/nsp_sld
                    if (trans_tmp == 0d0) cycle
                    
                    flx_psd(ips,idif_psd,iz) = flx_psd(ips,idif_psd,iz) + ( &
                        & - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips) &
                        & )
                enddo
                
                
                flx_psd(ips,ires_psd,iz) = sum(flx_psd(ips,:,iz))
                
            enddo
        enddo
            
        #ifdef dispPSDiter

        write(chrfmt,'(i0)') nflx_psd
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'

        print *
        print *,' [fluxes -- PSD] '
        print *, '<'//chrsp//'>'
        print trim(adjustl(chrfmt)),'rad',(chrflx_psd(iflx),iflx=1,nflx_psd)

        write(chrfmt,'(i0)') nflx_psd
        chrfmt = '(f5.2,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        do ips = 1, nps
            print trim(adjustl(chrfmt)), ps(ips), (sum(flx_psd(ips,iflx,:)*dz(:)),iflx=1,nflx_psd)
        enddo 
        print *

        #endif     
        

                
        if ( chkflx .and. dt > dt_th) then 
            ! flx_max_max = 0d0
            do ips = 1, nps
                flx_max = 0d0
                do iflx=1,nflx_psd 
                    flx_max= max( flx_max, abs( sum(flx_psd(ips,iflx,:)*dz(:)) ) )
                enddo 
                flx_max_max = max( flx_max_max, flx_max)
            enddo 
            do ips = 1, nps
            
                if ( flx_max > flx_max_max*flx_max_tol .and. abs( sum(flx_psd(ips,ires_psd,:)*dz(:)) ) > flx_max * flx_tol ) then 
                    
                    print *, chrsp,' too large error in PSD flx?'
                    print *, 'flx_max, flx_max_max,tol = ', flx_max,flx_max_max,flx_max_tol
                    print *, 'res, max = ', abs( sum(flx_psd(ips,ires_psd,:)*dz(:)) ), flx_max
                    print *, 'res/max, target = ', abs( sum(flx_psd(ips,ires_psd,:)*dz(:)) )/flx_max, flx_tol
                    
                    flgback = .true.
                
                endif 
                
            enddo 
        endif  
    
    endsubroutine psd_implicit_all_v4


end module scepter_psd 