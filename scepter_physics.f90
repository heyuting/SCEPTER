module scepter_physics
    use scepter_constants
    use scepter_variables
    implicit none
    private
    public :: make_transmx, calc_poro   

    contains

    subroutine make_transmx(  &
        & nsp_sld,imix,dz,poro,nz,z,zml,dbl_ref,tol,save_trans  &! input
        & ,trans  &! output 
        & )
        implicit none
        integer,intent(in)::nsp_sld,nz
        integer,dimension(nsp_sld),intent(in)::imix
        real(kind=8),intent(in)::dz(nz),poro(nz),z(nz),zml(nsp_sld),dbl_ref,tol
        real(kind=8)::sporo(nz)
        logical,intent(in)::save_trans
        real(kind=8),intent(out)::trans(nz,nz,nsp_sld)
        integer izml
        integer iz,isp,iiz,izdbl
        real(kind=8) :: translabs(nz,nz),dbio(nz),transdbio(nz,nz),transturbo2(nz,nz),transtill(nz,nz)
        real(kind=8) :: probh,dbl
        character(10) chr
        ! following must synch with main 
        integer,parameter :: imixtype_nobio     = 0
        integer,parameter :: imixtype_fick      = 1
        integer,parameter :: imixtype_turbo2    = 2
        integer,parameter :: imixtype_till      = 3
        integer,parameter :: imixtype_labs      = 4

        sporo = 1d0 - poro
        trans = 0d0
        !~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
        if (any(imix==imixtype_labs)) then
            translabs = 0d0

            open(unit=88,file='../input/labs-mtx.txt',action='read',status='unknown')
            do iz=1,nz
                read(88,*) translabs(iz,:)  ! writing 
            enddo
            close(88)

        endif

        if (.true.) then  ! devided by the time duration when transition matrices are created in LABS and weakening by a factor
        ! if (.false.) then 
            ! translabs = translabs *365.25d0/10d0*1d0/3d0  
            ! translabs = translabs *365.25d0/10d0*1d0/15d0  
            translabs = translabs *365.25d0/10d0*1d0/10d0
        endif
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! zml=zml_ref ! mixed layer depth assumed to be a reference value at first 

        dbl = dbl_ref

        do isp=1,nsp_sld
            
            dbio=0d0
            izdbl=0
            do iz = 1, nz
                if (z(iz) <= dbl) then 
                    dbio(iz) = 0d0
                    izdbl = iz
                elseif (dbl < z(iz) .and. z(iz) <=zml(isp)) then
                    ! dbio(iz) =  0.15d-4   !  within mixed layer 150 cm2/kyr (Emerson, 1985) 
                    dbio(iz) =  2d-4   !  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
                    ! dbio(iz) =  3d-5   !  within mixed layer (Jarvis et al., 2010) 
                    ! dbio(iz) =  2d-4*exp(z(iz)/0.1d0)   !  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
                    ! dbio(iz) =  2d-7*exp(z(iz)/1d0)   !  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
                    ! dbio(iz) =  2d-10   !  just a small value 
                    ! dbio(iz) =  2d-3   !  just a value changed 
                    izml = iz   ! determine grid of bottom of mixed layer 
                else
                    dbio(iz) =  0d0 ! no biodiffusion in deeper depths 
                endif
            enddo

            transdbio = 0d0   ! transition matrix to realize Fickian mixing with biodiffusion coefficient dbio which is defined just above 
            do iz = max(1,izdbl), izml
                if (iz==max(1,izdbl)) then 
                    transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                    transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                elseif (iz==izml) then 
                    transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
                    transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
                else 
                    transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                        + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                    transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
                    transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                endif
            enddo
            ! do iz = max(1,izdbl), izml
                ! if (iz==max(1,izdbl)) then 
                    ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                    ! transdbio(iz+1,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                ! elseif (iz==izml) then 
                    ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
                    ! transdbio(iz-1,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
                ! else 
                    ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                        ! + 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                    ! transdbio(iz-1,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
                    ! transdbio(iz+1,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
                ! endif
            ! enddo
            
            ! Added; changes have been made here rather than in solving governing eqs.
            do iz=1,nz
                transdbio(:,iz) = transdbio(:,iz)/dz(iz)
            enddo 

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            ! transition matrix for random mixing 
            transturbo2 = 0d0
            ! ending up in upward mixing 
            probh = 0.0010d0
            transturbo2(max(1,izdbl):izml,max(1,izdbl):izml) = probh  ! arbitrary assumed probability 
            do iz=1,izml  ! when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
            transturbo2(iz,iz)=-probh*(izml-max(1,izdbl))  
            enddo
            ! trying real homogeneous 
            transturbo2 = 0d0
            probh = 0.001d0 ! def used in IMP
            ! probh = 0.01d0 ! strong mixing
            ! probh = 0.05d0 ! strong mixing
            probh = 0.1d0 ! strong mixing
            probh = 0.5d0 ! strong mixing
            probh = 1.0d0 ! strong mixing
            ! probh = 1.5d0 ! strong mixing
            ! probh = 2d0 ! strong mixing
            ! probh = 5d0 ! strong mixing
            ! probh = 10d0 ! strong mixing
            probh = 20d0 ! strong mixing
            ! probh = 0.0005d0 ! just testing smaller mixing (used for tuning)
            ! probh = 0.0001d0 ! just testing smaller mixing for PSDs
            do iz=1,izml 
                do iiz=1,izml
                    if (iiz/=iz) then 
                        transturbo2(iiz,iz) = probh!*dz(iz)/dz(iiz)
                        transturbo2(iiz,iiz) = transturbo2(iiz,iiz) - transturbo2(iiz,iz)
                    endif 
                enddo
            enddo
            
            ! trying inverse mixing 
            transtill = 0d0
            probh = 0.010d0
            probh = 0.10d0
            ! probh = 1.00d0
            ! probh = 10.0d0
            do iz=1,izml  ! when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
                ! transtill(iz,iz)=-probh*dz(iz)/dz(izml+1-iz) !*(iz - izml*0.5d0)**2d0/(izml**2d0*0.25d0)
                ! transtill(izml+1-iz,iz)=probh*dz(iz)/dz(izml+1-iz) ! *(iz - izml*0.5d0)**2d0/(izml**2d0*0.25d0)
                ! do iiz = izml+1-iz,iz+1,-1
                    ! transtill(iiz,iz)= probh*dz(iz)/dz(iiz)*(iiz/real(izml+1-iz,kind=8))
                ! enddo 
                ! transtill(iz,iz) = -sum(transtill(:,iz))
                do iiz=1,izml
                    if (iiz/=iz) then 
                        if (iiz==iz-1 .or. iiz == iz+1 .or. iiz == izml + 1 - iz ) then 
                            transtill(iiz,iz)= probh !*dz(iz)/dz(iiz) 
                            ! transtill(iiz,iiz) = transtill(iiz,iiz) - transtill(iiz,iz)
                        endif 
                    endif 
                enddo 
                transtill(iz,iz) = -sum(transtill(:,iz))
            enddo
            

            ! if (turbo2(isp)) translabs = transturbo2   ! translabs temporarily used to represents nonlocal mixing 
            
            ! added 
            do iz =1,nz
                do iiz= 1,nz
                    translabs(iiz,iz) = translabs(iiz,iz)/dz(iz)*dz(iiz)
                    transturbo2(iiz,iz) = transturbo2(iiz,iz)/dz(iz)*dz(iiz)
                    transtill(iiz,iz) = transtill(iiz,iz)/dz(iz)*dz(iiz)
                enddo 
            enddo 
            
            trans(:,:,isp) = 0d0 
            
            if ( imix(isp)==imixtype_nobio ) cycle
            
            if ( imix(isp)==imixtype_fick ) then 
                trans(:,:,isp) = trans(:,:,isp) + transdbio(:,:)
            endif 
            
            if ( imix(isp)==imixtype_turbo2 ) then 
                trans(:,:,isp) = trans(:,:,isp) + transturbo2(:,:)
            endif 
            
            if ( imix(isp)==imixtype_labs ) then 
                trans(:,:,isp) = trans(:,:,isp) + translabs(:,:)
            endif 
            
            if ( imix(isp)==imixtype_till ) then 
                trans(:,:,isp) = trans(:,:,isp) + transtill(:,:)
            endif 
            
            ! if (any(abs(sum(trans(:,:,isp),dim=1))>tol)) then 
                ! print *, 'transition matrix can be non-conservative'
                ! write(chr,'(i3.3)') isp
                ! open(unit=88,file='./mtx-'//trim(adjustl(chr))//'.txt',action='write',status='replace')
                ! do iz=1,nz
                    ! write(88,*) (trans(iz,iiz,isp),iiz=1,nz)  ! writing 
                ! enddo
                ! close(88)
                ! stop
            ! endif 
            
            if (save_trans) then 
                write(chr,'(i3.3)') isp
                open(unit=88,file='./mtx-'//trim(adjustl(chr))//'.txt',action='write',status='replace')
                do iz=1,nz
                    write(88,*) (trans(iz,iiz,isp),iiz=1,nz)  ! writing 
                enddo
                close(88)
            endif 
            
            ! trans(:,:,isp) = transdbio(:,:)  !  firstly assume local mixing implemented by dbio 

            ! if (nonlocal(isp)) trans(:,:,isp) = translabs(:,:)  ! if nonlocal, replaced by either turbo2 mixing or labs mixing 
            ! if (nobio(isp)) trans(:,:,isp) = 0d0  ! if assuming no bioturbation, transition matrix is set at zero  
        enddo

    endsubroutine make_transmx

    subroutine calc_poro( &
        & nz,nsp_sld,nflx,idif,irain &! in
        & ,flx_sld,mv,poroprev,w,poroi,w_btm,dz,tol,dt &! in
        & ,poro &! inout
        & )
        implicit none 

        integer,intent(in)::nz,nsp_sld,nflx,idif,irain
        real(kind=8),intent(in)::poroi,w_btm,tol,dt
        real(kind=8),dimension(nz),intent(in)::dz,w,poroprev
        real(kind=8),dimension(nz),intent(inout)::poro
        real(kind=8),dimension(nsp_sld,nflx,nz),intent(in)::flx_sld
        real(kind=8),dimension(nsp_sld),intent(in)::mv
        ! local 
        real(kind=8),dimension(nz)::DV,resi_poro
        integer iz,isps,row,ie,ie2
        real(kind=8) w_tmp,wp_tmp,sporo_tmp,sporop_tmp,sporoprev_tmp

        real(kind=8),parameter::infinity = huge(0d0)

        real(kind=8) amx3(nz,nz),ymx3(nz)
        integer ipiv3(nz) 
        integer info 

        external DGESV

        ! 
        ! attempt to solve porosity under any kind of porosity - uplift rate relationship 
        ! based on equation:
        ! d(1-poro)/dt = d(1-poro)*w/dz - mv*1d-6*sum( flx_sld(mixing, dust, rxns) ) 

            
        ymx3 = 0d0
        amx3 = 0d0
        DV = 0d0

        do iz=1,nz
            do isps = 1,nsp_sld 
                DV(iz) = DV(iz) + ( flx_sld(isps, 4 + isps,iz) + flx_sld(isps, idif ,iz) + flx_sld(isps, irain ,iz) ) &
                    & *mv(isps)*1d-6 
            enddo 
            
            row = iz
            
            w_tmp = w(iz)
            wp_tmp = w(min(nz,iz+1))
            sporo_tmp = 1d0-poro(iz)
            sporop_tmp = 1d0-poro(min(nz,iz+1)) 
            sporoprev_tmp = 1d0-poroprev(iz)
                    
            if (iz==nz) then 
                wp_tmp = w_btm
                sporop_tmp = 1d0 - poroi
            endif 
            
            if (iz/=nz) then 
            
                ymx3(row) = ( &
                    & + (1d0 - sporoprev_tmp)/dt    &
                    & - ( 1d0*wp_tmp - 1d0*w_tmp)/dz(iz)  &
                    & + DV(iz) &
                    & )
                    
                amx3(row,row) = ( &
                    & + (-1d0 )/dt    &
                    & - ( - (-1d0)*w_tmp)/dz(iz)  &
                    & )
                    
                amx3(row,row+1) = ( &
                    & - ( -1d0*wp_tmp )/dz(iz)  &
                    & )
                
            else 
            
                ymx3(row) = ( &
                    & + (1d0 - sporoprev_tmp)/dt    &
                    & - ( sporop_tmp*wp_tmp - 1d0*w_tmp)/dz(iz)  &
                    & + DV(iz) &
                    & )
                    
                amx3(row,row) = ( &
                    & + (-1d0 )/dt    &
                    & - (- (-1d0)*w_tmp)/dz(iz)  &
                    & )
            
            
            endif 
            
        enddo 
            
        ymx3=-1.0d0*ymx3

        if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
        ! if (.true.) then 
            print*,'porocalc: error in mtx'
            print*,'porocalc: any(isnan(amx3)),any(isnan(ymx3))'
            print*,any(isnan(amx3)),any(isnan(ymx3))

            if (any(isnan(ymx3))) then 
                do iz = 1, nz
                    if (isnan(ymx3(iz))) then 
                        print*,'porocalc: NAN is here...',iz
                    endif
                enddo 
            endif


            if (any(isnan(amx3))) then 
                do ie = 1,(nz)
                    do ie2 = 1,(nz)
                        if (isnan(amx3(ie,ie2))) then 
                            print*,'porocalc: NAN is here...',ie,ie2
                        endif
                    enddo
                enddo
            endif
            stop
            
        endif

        call DGESV(Nz,int(1),amx3,Nz,IPIV3,ymx3,Nz,INFO) 

        poro = ymx3

        ! resi_poro = 0d0

        ! do iz=1,nz
            
            ! w_tmp = w(iz)
            ! wp_tmp = w(min(nz,iz+1))
            ! sporo_tmp = 1d0-poro(iz)
            ! sporop_tmp = 1d0-poro(min(nz,iz+1)) 
            ! sporoprev_tmp = 1d0-poroprev(iz)
                    
            ! if (iz==nz) then 
                ! wp_tmp = w_btm
                ! sporop_tmp = 1d0 - poroi
            ! endif 
            
            ! resi_poro(iz) = ( &
                    ! & + (sporo_tmp - sporoprev_tmp)/dt    &
                    ! & - ( sporop_tmp*wp_tmp - sporo_tmp*w_tmp)/dz(iz)  &
                    ! & + DV(iz)  &
                    ! & )
        ! enddo 

        ! if ( maxval(resi_poro) > tol) then 
            ! print *, 'porosity calculation is wierd!?' 
            ! print *, info
            ! print *, poro
            ! print *, resi_poro
            ! stop
        ! endif 
    
    endsubroutine calc_poro

    subroutine calc_uplift( &
        & nz,nsp_sld,nflx,idif,irain &! IN
        & ,iwtype &! in
        & ,flx_sld,mv,poroi,w_btm,dz,poro,poroprev,dt &! in
        & ,w &! inout
        & )
        implicit none 

        integer,intent(in)::nz,nsp_sld,nflx,idif,irain
        integer,intent(in)::iwtype
        real(kind=8),intent(in)::poroi,w_btm,dt
        real(kind=8),dimension(nz),intent(in)::dz,poro,poroprev
        real(kind=8),dimension(nz),intent(inout)::w
        real(kind=8),dimension(nsp_sld,nflx,nz),intent(in)::flx_sld
        real(kind=8),dimension(nsp_sld),intent(in)::mv
        ! local 
        real(kind=8),dimension(nz)::DV,wsporo
        integer iz,isps
        integer,parameter :: iwtype_cnst = 0
        integer,parameter :: iwtype_pwcnst = 1
        integer,parameter :: iwtype_spwcnst = 2
        integer,parameter :: iwtype_flex = 3



        select case(iwtype)
            case(iwtype_cnst) ! default case with constant uplift rate
                
                w = w_btm
                
            case(iwtype_pwcnst) ! poro * w = const
                
                wsporo = w_btm*poroi
                w = wsporo/poro
                
            case(iwtype_spwcnst) ! (1 - poro) * w = const
            
                wsporo = w_btm*(1d0 - poroi)
                w = wsporo/(1d0-poro)
                
            case(iwtype_flex) ! flexible w (including const porosity)
                ! in this case, porosity is not calculated but given so equation for porosity is used to solve w instead  
                ! based on equation:
                ! d(1-poro)/dt = d(1-poro)*w/dz - mv*1d-6*sum( flx_sld(mixing, dust, rxns) ) 
                ! now porosity is cont.
                ! 0 = (1-poro)*dw/dz - mv*1d-6*sum( flx_sld(mixing, dust, rxns) ) 
                
                DV = 0d0

                do iz=1,nz
                    do isps = 1,nsp_sld 
                        DV(iz) = DV(iz) + ( flx_sld(isps, 4 + isps,iz) + flx_sld(isps, idif ,iz) + flx_sld(isps, irain ,iz) ) &
                            & *mv(isps)*1d-6 
                    enddo 
                enddo 
                
                do iz=nz,1,-1
                    ! 0d0 = (1d0 - poro(iz)) * (w(iz+1) - w(iz))/dz(iz) + DV(iz)
                    ! 0d0 = (1-poro(iz)) * (w(iz+1) - w(iz)) + DV(iz)*dz(iz)
                    ! 0d0 = w(iz+1) - w(iz) + DV(iz)*dz(iz)/(1d0 - poro(iz)) 
                    ! w(iz) = w(iz+1)  + DV(iz)*dz(iz)/(1d0 - poro(iz)) 
                    ! ... more generally ... 
                    ! ( (1d0 - poro(iz)) - (1d0 - poroprev(iz)) )/dt = ( (1d0 - poro(iz+1)) * w(iz+1) - (1d0 - poro(iz)) * w(iz) )/dz(iz) + DV(iz)
                    ! 0d0 = ( (1d0 - poro(iz+1)) * w(iz+1) - (1d0 - poro(iz)) * w(iz) )/dz(iz) + DV(iz) - ( (1d0 - poro(iz)) - (1d0 - poroprev(iz)) )/dt
                    ! 0d0 = ( (1d0 - poro(iz+1)) * w(iz+1) - (1d0 - poro(iz)) * w(iz) ) + DV(iz)*dz(iz) - ( (1d0 - poro(iz)) - (1d0 - poroprev(iz)) )/dt*dz(iz)
                    ! (1d0 - poro(iz)) * w(iz) =  (1d0 - poro(iz+1)) * w(iz+1)  + DV(iz)*dz(iz) - ( (1d0 - poro(iz)) - (1d0 - poroprev(iz)) )/dt*dz(iz)
                    ! w(iz) =  (1d0 - poro(iz+1))/(1d0 - poro(iz))  * w(iz+1)  + DV(iz)*dz(iz)/(1d0 - poro(iz))  - ( 1d0 - (1d0 - poroprev(iz))/(1d0 - poro(iz)) )/dt*dz(iz)
                    if (iz==nz) then 
                        w(iz) = w_btm  + DV(iz)*dz(iz)/(1d0 - poro(iz)) 
                        ! general version
                        w(iz) =  (1d0 - poroi)/(1d0 - poro(iz))  * w_btm  &
                            & + DV(iz)*dz(iz)/(1d0 - poro(iz))  - ( 1d0 - (1d0 - poroprev(iz))/(1d0 - poro(iz)) )/dt*dz(iz)
                    else
                        w(iz) = w(iz+1)  + DV(iz)*dz(iz)/(1d0 - poro(iz)) 
                        ! general version
                        w(iz) =  (1d0 - poro(iz+1))/(1d0 - poro(iz))  * w(iz+1)  &
                            & + DV(iz)*dz(iz)/(1d0 - poro(iz))  - ( 1d0 - (1d0 - poroprev(iz))/(1d0 - poro(iz)) )/dt*dz(iz)
                    endif 
                enddo 
                
            case default 
                
                w = w_btm
                
        endselect 
    
    endsubroutine calc_uplift

    subroutine calcupwindscheme(  &
        up,dwn,cnr,adf & ! output 
        ,w,nz   & ! input &
        )
        implicit none
        integer,intent(in)::nz
        real(kind=8),intent(in)::w(nz)
        real(kind=8),dimension(nz),intent(out)::up,dwn,cnr,adf
        real(kind=8) corrf
        real(kind=8) :: cnr_save(nz)
        integer iz
        ! copied and pasted from iMP code and modified for weathering 

        ! ------------ determine variables to realize advection 
        !  upwind scheme 
        !  up  ---- burial advection at grid i = sporo(i)*w(i)*(some conc. at i) - sporo(i-1)*w(i-1)*(some conc. at i - 1) 
        !  dwn ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i)*w(i)*(some conc. at i) 
        !  cnr ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i-1)*w(i-1)*(some conc. at i - 1) 
        !  when burial rate is positive, scheme need to choose up, i.e., up = 1.  
        !  when burial rate is negative, scheme need to choose dwn, i.e., dwn = 1.  
        !  where burial change from positive to negative or vice versa, scheme chooses cnr, i.e., cnr = 1. for the mass balance sake 

        up = 0
        dwn=0
        cnr =0
        adf=1d0
        do iz=1,nz 
            if (iz==1) then 
                if (w(iz)>=0d0 .and. w(iz+1)>=0d0) then  ! positive burial 
                    up(iz) = 1
                elseif (w(iz)<=0d0 .and. w(iz+1)<=0d0) then  ! negative burial 
                    dwn(iz) = 1
                else   !  where burial sign changes  
                    if (.not.(w(iz)*w(iz+1) <=0d0)) then 
                        print*,'error'
                        stop
                    endif
                    cnr(iz) = 1
                endif
            elseif (iz==nz) then 
                if (w(iz)>=0d0 .and. w(iz-1)>=0d0) then
                    up(iz) = 1
                elseif (w(iz)<=0d0 .and. w(iz-1)<=0d0) then
                    dwn(iz) = 1
                else 
                    if (.not.(w(iz)*w(iz-1) <=0d0)) then 
                        print*,'error'
                        stop
                    endif
                    cnr(iz) = 1
                endif
            else 
                ! if iz-1 and iz+1 have the same sign, then it can be assigned either as up or dwn
                ! else cnr whose neighbor has a different sign 
                if (w(iz) >=0d0) then 
                    if (w(iz+1)>=0d0 .and. w(iz-1)>=0d0) then
                        up(iz) = 1
                    else
                        cnr(iz) = 1
                    endif
                else  
                    if (w(iz+1)<=0d0 .and. w(iz-1)<=0d0) then
                        dwn(iz) = 1
                    else
                        cnr(iz) = 1
                    endif
                endif
            endif
        enddo        

        if (sum(up(:)+dwn(:)+cnr(:))/=nz) then
            print*,'error',sum(up),sum(dwn),sum(cnr)
            stop
        endif

        ! try to make sure mass balance where advection direction changes 
        ! 
        ! case (i)
        !       :           w         direction    
        !     iz - 2        +             ^              w(iz-1) - w(iz-2)
        !     iz - 1        +             ^              w(iz  ) - w(iz-1)
        !     iz            +             ^            a[w(iz+1) - w(iz  )] + b[w(iz+1) - w(iz-1)]  
        !     iz + 1        -             v            c[w(iz+1) - w(iz  )] + d[w(iz+2) - w(iz  )]
        !     iz + 2        -             v              w(iz+2) - w(iz+1)
        !     iz + 3        -             v              w(iz+3) - w(iz+2) 
        ! layers [iz] & [iz+1] must yield [w(iz+1) - w(iz  )]
        ! and calculated as (a+b+c)w(iz+1) - (a+c+d)w(iz  ) - b w(iz-1) + d w(iz+1) 
        ! thus b = d = 0 and   a + b + c = 1 and a + c + d = 1
        ! a and c can be arbitrary as long as satisfying a + c = 1 
        ! --------------------------------------------------------------------------------------------
        ! case (ii)
        !       :           w         direction    
        !     iz - 2        -             v              w(iz-2) - w(iz-3)
        !     iz - 1        -             v              w(iz-1) - w(iz-2)
        !     iz            -             v            a[w(iz  ) - w(iz-1)] + b[w(iz+1) - w(iz-1)]  
        !     iz + 1        +             ^            c[w(iz+2) - w(iz+1)] + d[w(iz+2) - w(iz  )]
        !     iz + 2        +             ^              w(iz+3) - w(iz+2)
        !     iz + 3        +             ^              w(iz+4) - w(iz+3) 
        ! layers [iz] & [iz+1] must yield [w(iz+2) - w(iz-1)]
        ! and calculated as (c+d)w(iz+2) - (a+b)w(iz-1) + (a-d)w(iz  ) + (b-c)w(iz+1) 
        ! thus c + d = 1, a + b = 1, a - d = 0, and b - c = 0
        ! these can be satisfied by b = c = 1 - a and d = a and a can be arbitrary as long as 0 <= a <= 1
        cnr_save = cnr
        do iz=1,nz-1
            if (cnr_save(iz)==1 .and. cnr_save(iz+1)==1) then 
            ! if (cnr(iz)==1 .and. cnr(iz+1)==1) then 
                if (w(iz) < 0d0 .and. w(iz+1) >= 0d0) then
                    corrf = 5d0  !  This assignment of central advection term helps conversion especially when assuming turbo2 mixing 
                    cnr(iz+1)=abs(w(iz)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
                    cnr(iz)=abs(w(iz+1)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
                    dwn(iz+1)=1d0-cnr(iz+1)
                    up(iz)=1d0-cnr(iz)
                endif 
            endif 
            if (cnr_save(iz)==1 .and. cnr_save(iz+1)==1) then 
            ! if (cnr(iz)==1 .and. cnr(iz+1)==1) then 
                if (w(iz)>= 0d0 .and. w(iz+1) < 0d0) then
                    cnr(iz+1)=0
                    cnr(iz)=0
                    up(iz+1)=1
                    dwn(iz)=1
                    adf(iz)=abs(w(iz+1))/(abs(w(iz+1))+abs(w(iz)))
                    adf(iz+1)=abs(w(iz))/(abs(w(iz+1))+abs(w(iz)))
                endif 
            endif 
        enddo       

    endsubroutine calcupwindscheme

end module scepter_physics 