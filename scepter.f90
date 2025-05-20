!***********************************************************************
! File: scepter.f90
! Purpose: Main program for SCEPTER (Soil Chemical Evolution and Physical Transport of Elements in Regolith)
! This program integrates chemical reactions, transport processes, and physical evolution
! of soil systems over time.
!***********************************************************************

program weathering
    ! Module imports for various components of the weathering simulation
    use scepter_constants    ! Physical and chemical constants
    use scepter_IO          ! Input/Output operations
    use scepter_input       ! Input parameter handling
    use scepter_weathering_main  ! Main weathering simulation routines
    
    implicit none 
    integer nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nz,nsld_kinspc
    character(5),dimension(:),allocatable::chraq,chrsld,chrgas,chrrxn_ext,chrsld_kinspc 
    real(kind=8),dimension(:),allocatable::kin_sld_spc
    character(500) sim_name,runname_save,cwd,path,path2,cmd
    real(kind=8) ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,plant_rain,zml_ref,tc,rainpowder_2nd &
        & ,step_tau
    integer count_dtunchanged_Max
    ! Get and set up working directory
    CALL getcwd(cwd)
    WRITE(*,*) TRIM(cwd)

    ! Get program path and set up simulation directory
    call getarg(0,path)
    WRITE(*,*) TRIM(path)
    ! call get_command_argument(0, cmd)
    ! WRITE(*,*) TRIM(cmd)
    ! path2 = path(:index(path,'weathering')-2)
    path2 = path(:len(trim(path))-len('scepter')-1)
    WRITE(*,*) TRIM(path2)

    ! Change to simulation directory
    CALL chdir(TRIM(path2))
    CALL getcwd(path)
    WRITE(*,*) TRIM(path)

    ! Initialize simulation parameters
    ! Get number of species and reactions
    call get_variables_num( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! output
        & )

    print *,nsp_sld,nsp_aq,nsp_gas,nrxn_ext

    ! Allocate arrays for chemical species and reactions
    allocate(chraq(nsp_aq),chrsld(nsp_sld),chrgas(nsp_gas),chrrxn_ext(nrxn_ext))
    allocate(chrsld_kinspc(nsld_kinspc),kin_sld_spc(nsld_kinspc))
    
    ! Get chemical species names and reaction information
    call get_variables( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! input
        & ,chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc &! output
        & ) 
        
    ! Print species information for verification
    print *,chraq
    print *,chrsld 
    print *,chrgas 
    print *,chrrxn_ext 
    print *,chrsld_kinspc 
    print *,kin_sld_spc 
    
    ! pause
    sim_name = 'chkchk'

    ! Get soil and simulation parameters
    call get_bsdvalues( &
        & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
        & ,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau &
        & )
        
    ! Run main weathering simulation
    call weathering_main( &
        & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,ttot,plant_rain,rainpowder_2nd  &! input
        & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save &! input
        & ,count_dtunchanged_Max,tc,step_tau &! input 
        & ,nsld_kinspc,chrsld_kinspc,kin_sld_spc &! input
        & )

endprogram weathering