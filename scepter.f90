program weathering
    use scepter_constants
    use scepter_variables
    use scepter_weathering_main
    use scepter_IO
    use scepter_input
    use scepter_physics
    use scepter_psd
    use scepter_concentration
    use scepter_reactions
    implicit none 

    CALL getcwd(cwd)
    WRITE(*,*) TRIM(cwd)

    call getarg(0,path)
    WRITE(*,*) TRIM(path)
    ! call get_command_argument(0, cmd)
    ! WRITE(*,*) TRIM(cmd)
    ! path2 = path(:index(path,'weathering')-2)
    path2 = path(:len(trim(path))-len('scepter')-1)
    WRITE(*,*) TRIM(path2)

    CALL chdir(TRIM(path2))
    CALL getcwd(path)
    WRITE(*,*) TRIM(path)

    call get_variables_num( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! output
        & )

    print *,nsp_sld,nsp_aq,nsp_gas,nrxn_ext

    allocate(chraq(nsp_aq),chrsld(nsp_sld),chrgas(nsp_gas),chrrxn_ext(nrxn_ext))
    allocate(chrsld_kinspc(nsld_kinspc),kin_sld_spc(nsld_kinspc))
        
    call get_variables( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! input
        & ,chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc &! output
        & ) 
        
    print *,chraq
    print *,chrsld 
    print *,chrgas 
    print *,chrrxn_ext 
    print *,chrsld_kinspc 
    print *,kin_sld_spc 
    
    ! pause
    sim_name = 'chkchk'

    call get_bsdvalues( &
        & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
        & ,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau &
        & )
        
    call weathering_main( &
        & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,ttot,plant_rain,rainpowder_2nd  &! input
        & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save &! input
        & ,count_dtunchanged_Max,tc,step_tau &! input 
        & ,nsld_kinspc,chrsld_kinspc,kin_sld_spc &! input
        & )

    contains 

    subroutine makegrid(beta,nz,ztot,dz,z,regular_grid)  !  making grid, after Hoffmann & Chiang, 2000
        implicit none
        integer(kind=4),intent(in) :: nz
        logical,intent(in)::regular_grid
        real(kind=8),intent(in)::beta,ztot
        real(kind=8),intent(out)::dz(nz),z(nz)
        integer(kind=4) iz

        do iz = 1, nz 
            z(iz) = iz*ztot/nz  ! regular grid 
            if (iz==1) then
                dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
            endif
            if (iz/=1) then 
                dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(:iz-1))
            endif
        enddo

        if (regular_grid) then 
            dz = ztot/nz  ! when implementing regular grid
        endif 

        do iz=1,nz  ! depth is defined at the middle of individual layers 
            if (iz==1) z(iz)=dz(iz)*0.5d0  
            if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
        enddo

    endsubroutine makegrid

    #ifdef no_intr_findloc
        function findloc(chrlist_in,chrspecific,dim)
            implicit none
            character(*),intent(in)::chrlist_in(:),chrspecific
            integer,intent(in)::dim
            integer findloc,i

            findloc = 0
            do i=1, size(chrlist_in,dim=dim)
                if (trim(adjustl(chrspecific)) == trim(adjustl(chrlist_in(i)))) then
                    findloc = i
                    return
                endif 

            enddo 

        endfunction findloc
    #endif 

endprogram weathering