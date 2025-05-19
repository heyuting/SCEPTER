module scepter_kinetics
    use scepter_constants
    use scepter_variables
    
    implicit none
    private
    public :: sld_rxn, rough_f
    ! Constants
    real(kind=8), parameter :: cal2j = 4.184d0

    contains

    !-----------------------------------------------------------------------
    ! Subroutine: sld_rxn
    ! Purpose: Calculate reaction rates for solid dissolution/precipitation
    !-----------------------------------------------------------------------
    subroutine sld_rxn( &
        & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx &! input 
        & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle,solmod &! input
        & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas,chrsld &! input
        & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
        & ) 
        implicit none 

        integer,intent(in)::nz,nsp_sld,nsp_aq,nsp_gas
        real(kind=8),intent(in)::msld_seed,dt
        real(kind=8),dimension(nz),intent(in)::poro,sat,dz
        real(kind=8),dimension(nsp_sld,nz),intent(in)::hr
        real(kind=8),dimension(nsp_sld),intent(in)::mv,msldth
        real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
        real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas
        real(kind=8),dimension(nsp_aq),intent(in)::maqth
        real(kind=8),dimension(nsp_aq,nz),intent(in)::maq
        real(kind=8),dimension(nsp_gas),intent(in)::mgasth
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,agas
        real(kind=8),dimension(nsp_sld,nz),intent(in)::ksld,omega,nonprec,msldx,msld,solmod
        real(kind=8),dimension(nsp_sld,nsp_aq,nz),intent(in)::dksld_dmaq,domega_dmaq
        real(kind=8),dimension(nsp_sld,nsp_gas,nz),intent(in)::dksld_dmgas,domega_dmgas
        character(10),dimension(nsp_sld),intent(in)::precstyle
        character(5),dimension(nsp_sld),intent(in)::chrsld
        real(kind=8),dimension(nsp_sld,nz),intent(out)::rxnsld,drxnsld_dmsld
        real(kind=8),dimension(nsp_sld,nsp_aq,nz),intent(out)::drxnsld_dmaq
        real(kind=8),dimension(nsp_sld,nsp_gas,nz),intent(out)::drxnsld_dmgas

        integer ispa,isps,ispg,iz
        real(kind=8),dimension(nsp_sld,nz)::maxdis,maxprec
        real(kind=8),parameter::infinity = huge(0d0)
            
            

        ! *** sanity check
        if (any(isnan(ksld)) .or. any(ksld>infinity)) then 
            print *,' *** found insanity in ksld (in sld_rxn): listing below -- '
            do isps=1,nsp_sld
                do iz=1,nz
                    if (isnan(ksld(isps,iz)) .or. ksld(isps,iz)>infinity) print*,chrsld(isps),iz,ksld(isps,iz)
                enddo
            enddo 
            stop
        endif 
            
        rxnsld = 0d0
        drxnsld_dmsld = 0d0
        drxnsld_dmaq = 0d0
        drxnsld_dmgas = 0d0

        maxdis = 1d200
        maxprec = -1d200
        do isps=1,nsp_sld
            maxdis(isps,:) = min(maxdis(isps,:),(msld(isps,:)-msldth(isps))/dt)
            do ispa = 1,nsp_aq
                if (staq(isps,ispa)<0d0) then 
                    maxdis(isps,:) = min(maxdis(isps,:), -1d0/staq(isps,ispa)*poro*sat*1d3*(maq(ispa,:)-maqth(ispa))/dt )
                endif 
            enddo 
            do ispg = 1,nsp_gas
                if (stgas(isps,ispg)<0d0) then 
                    maxdis(isps,:) = min(maxdis(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:)*(mgas(ispg,:)-mgasth(ispg))/dt )
                endif 
            enddo 
            do ispa = 1,nsp_aq
                if (staq(isps,ispa)>0d0) then 
                    maxprec(isps,:) = max(maxprec(isps,:), -1d0/staq(isps,ispa)*poro*sat*1d3*(maq(ispa,:)-maqth(ispa))/dt )
                endif 
            enddo 
            do ispg = 1,nsp_gas
                if (stgas(isps,ispg)>0d0) then 
                    maxprec(isps,:) = max(maxprec(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:)*(mgas(ispg,:)-mgasth(ispg))/dt )
                endif 
            enddo 
        enddo 

        do isps = 1,nsp_sld
            select case(trim(adjustl(precstyle(isps))))
            
                case ('full_lim') 
                    
                    do iz = 1,nz
                        if (1d0-omega(isps,iz) > 0d0) then 
                            rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) 
                            if (rxnsld(isps,iz)> maxdis(isps,iz)) then 
                                rxnsld(isps,iz) = maxdis(isps,iz)
                                drxnsld_dmsld(isps,iz) = 0d0
                                drxnsld_dmaq(isps,:,iz) = 0d0
                                drxnsld_dmgas(isps,:,iz) = 0d0
                            else 
                                drxnsld_dmsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz)) 
                                drxnsld_dmaq(isps,:,iz) = ( +ksld(isps,iz)*poro(iz)*hr(isps,iz)* &
                                    & mv(isps)*1d-6*msldx(isps,iz)* (-domega_dmaq(isps,:,iz)) + & 
                                    & dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps) &
                                    & *1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)))
                                drxnsld_dmgas(isps,:,iz) = ( +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps) & 
                                    &*1d-6*msldx(isps,iz)*(-domega_dmgas(isps,:,iz)) + dksld_dmgas(isps,:,iz)* &
                                    & poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz) *(1d0-omega(isps,iz)))
                            endif 
                        elseif (1d0-omega(isps,iz) < 0d0) then 
                            if (nonprec(isps,iz)==1d0) then 
                                rxnsld(isps,iz) = 0d0
                                drxnsld_dmsld(isps,iz) = 0d0
                                drxnsld_dmaq(isps,:,iz) = 0d0
                                drxnsld_dmgas(isps,:,iz) = 0d0
                            elseif (nonprec(isps,iz)==0d0) then 
                                rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz))
                                if (rxnsld(isps,iz) < maxprec(isps,iz)) then 
                                    rxnsld(isps,iz) = maxprec(isps,iz)
                                    drxnsld_dmsld(isps,iz) = 0d0
                                    drxnsld_dmaq(isps,:,iz) = 0d0
                                    drxnsld_dmgas(isps,:,iz) = 0d0
                                else
                                    rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz))
                                    drxnsld_dmsld(isps,iz) = 0d0
                                    drxnsld_dmaq(isps,:,iz) = ( &
                                        & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*(-domega_dmaq(isps,:,iz)) &
                                        & +dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz)) &
                                        & )
                                    drxnsld_dmgas(isps,:,iz) = ( &
                                        & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*(-domega_dmgas(isps,:,iz)) &
                                        & +dksld_dmgas(isps,:,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz)) &
                                        & )
                                endif 
                            endif 
                        endif 
                    enddo 
                    
            
                case ('full') 
                    rxnsld(isps,:) = ( &
                        & + min(ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:) &
                        & *(1d0-omega(isps,:))/(1d0-poro),maxdis(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                        &  + max(ksld(isps,:)*poro*hr(isps,:)*(1d0-omega(isps,:)),maxprec(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                        & )
                    
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                        & )
                        
                    do ispa = 1,nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            &  + ksld(isps,:)*poro*hr(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            &  + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            & )
                    enddo 
                        
                    do ispg = 1,nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            & /(1d0-poro) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                            & /(1d0-poro) &
                            &  + ksld(isps,:)*poro*hr(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            &  + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                            & )
                    enddo 
                    
                    
                case('seed')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6 &
                            & *(msldx(isps,:)+msld_seed)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6 &
                            & *(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6 &
                            & *(msldx(isps,:)+msld_seed)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6 &
                            & *(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('decay')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*msldx(isps,:)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('2/3noporo')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6)**(2d0/3d0) &
                        &       *(2d0/3d0)*msldx(isps,:)**(-1d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('2/3')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:)*(mv(isps)*1d-6 &
                        & *msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:)*(mv(isps)*1d-6)**(2d0/3d0) &
                        &       *(2d0/3d0)*msldx(isps,:)**(-1d0/3d0)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro**(2d0/3d0)*hr(isps,:) &
                            & *(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('psd_full')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*hr(isps,:)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*hr(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*hr(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case('emmanuel')
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)*solmod(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6 &
                            & *msldx(isps,:)*(-domega_dmaq(isps,ispa,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps) &
                            & *1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:) &
                            & *(-domega_dmgas(isps,ispg,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps) &
                            & *1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                
                
                case default
                    rxnsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    if (any(isnan(rxnsld(isps,:)))) then 
                        print *, 'NAN in rxnsld: ',chrsld(isps)
                        print *, 'ksld(isps,:)',ksld(isps,:)
                        print *, 'hr(isps,:)',hr(isps,:)
                        print *, 'msldx(isps,:)',msldx(isps,:)
                        print *, 'omega(isps,:)',omega(isps,:)
                        stop
                    endif 
                
                    drxnsld_dmsld(isps,:) = ( &
                        & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                        & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                        & )
                    
                    do ispa = 1, nsp_aq
                        drxnsld_dmaq(isps,ispa,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    ! correcting for solid fraction available to porewater 
                    
                    ! rxnsld(isps,:) = rxnsld(isps,:)/(1d0-poro)
                    ! drxnsld_dmsld(isps,:) = drxnsld_dmsld(isps,:)/(1d0-poro)
                    ! drxnsld_dmaq(isps,:,:) = drxnsld_dmaq(isps,:,:)/(1d0-poro)
                    ! drxnsld_dmgas(isps,:,:) = drxnsld_dmgas(isps,:,:)/(1d0-poro)
                    
                    
                    ! attempt to add authigenesis above some threshould for omega
                    ! do iz=1,nz 
                        ! if (nonprec(isps,iz)==0d0 .and. omega(isps,iz) > auth_th) then 
                            ! rxnsld(isps,iz) = rxnsld(isps,iz) + ( &
                                ! & + ksld(isps,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                ! & )
                            
                            ! do ispa = 1, nsp_aq
                                ! drxnsld_dmaq(isps,ispa,iz) = drxnsld_dmaq(isps,ispa,iz) + ( &
                                    ! & + ksld(isps,iz)*poro(iz)*hr(iz)*(-domega_dmaq(isps,ispa,iz)) &
                                    ! & + dksld_dmaq(isps,ispa,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                    ! & )
                            ! enddo 
                            
                            ! do ispg = 1, nsp_gas
                                ! drxnsld_dmgas(isps,ispg,iz) = drxnsld_dmgas(isps,ispg,iz) + ( &
                                    ! & + ksld(isps,iz)*poro(iz)*hr(iz)*(-domega_dmgas(isps,ispg,iz)) &
                                    ! & + dksld_dmgas(isps,ispg,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                    ! & )
                            ! enddo 
                        ! endif 
                    ! enddo 
                    
                    ! print *, 'max-rxnflx', isps, sum (ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*dz)
            
            endselect
        enddo   
    
    endsubroutine sld_rxn

    !-----------------------------------------------------------------------
    ! Function: rough_f
    ! Purpose: Calculate roughness factor for solid dissolution/precipitation
    !-----------------------------------------------------------------------
    function rough_f(ref_dummy,n_dummy,r_dummy)
        implicit none
        integer n_dummy
        real(kind=8),dimension(n_dummy):: rough_f,r_dummy
        character(10) ref_dummy
        selectcase(trim(adjustl(ref_dummy)))
            case('NSB07')       ! Navarre-Sitchler and Brantley (2007)
                rough_f = 10d0**3.3d0*r_dummy**0.33d0
            case('BM00')        ! Brantley and Mellott (2000)
                rough_f = 10d0**0.7d0*r_dummy**(-0.1d0)
            case('Letal21')     ! Lewis et al. (2021) 
                ! rough_f = 10d0**(   154.25d0 * exp( 1.0219d0 * log10( r_dummy ) ) ) ! (assuming sphere)
                ! rough_f = 10d0**(   113.41d0 * exp( 1.0219d0 * log10( r_dummy ) ) ) ! (assuming cube)
                ! rough_f = 10d0**(   max( 2.02d0*log10( r_dummy ) + 10.734d0, 1d0 ) )  ! (assuming sphere)
                rough_f = 10d0**(   max( 2.02d0*log10( r_dummy ) + 10.126d0, 1d0 ) )  ! (assuming cube)
            case('smooth')
                rough_f = 1d0
            case default 
                print*, '*** error in rough_f --> stop'
                stop
        endselect
    endfunction rough_f
    
end module scepter_kinetics