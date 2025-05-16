!***********************************************************************
!Name: scepter_psd_implicit
!Purpose: solve for the evolution of the particle size distribution (PSD)
! in a multi-layer (vertical) system, accounting for advection, diffusion, 
! reaction (dissolution/precipitation), and possibly rain input
!***********************************************************************

module scepter_psd_implicit
    use scepter_constants
    use scepter_variables
    implicit none
    public :: psd_implicit_all_v2
    public :: psd_implicit_all_v4


    contains

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

    !--------------------------------------------------------------------------------------------------
    ! Subroutine to do psd ( defined with particle number / bulk m3 / log (r) )
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
    !--------------------------------------------------------------------------------------------------
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

        external DGESV

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

endmodule scepter_psd_implicit