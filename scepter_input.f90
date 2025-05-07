module scepter_input
    use scepter_constants
    use scepter_variables
    implicit none

    contains

    subroutine get_clim_num( &
        & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
        & ,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau &
        & )
        implicit none

        integer,intent(out):: nz,count_dtunchanged_Max
        real(kind=8),intent(out):: ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,tc,rainpowder_2nd,step_tau
        character(256),intent(out):: sim_name,plant_rain,runname_save

        character(500) file_name
        integer i

        file_name = './clim.in'
        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        read(50,*) nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau
        close(50)

    endsubroutine get_clim_num

    subroutine get_bsdvalues( &
        & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
        & ,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau &
        & )
        implicit none

        integer,intent(out):: nz,count_dtunchanged_Max
        real(kind=8),intent(out):: ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,tc,rainpowder_2nd,step_tau
        character(256),intent(out):: sim_name,plant_rain,runname_save

        call get_clim_num( &
            & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
            & ,count_dtunchanged_Max,tc,rainpowder_2nd,step_tau &
            & )

    endsubroutine get_bsdvalues

    subroutine get_rainwater( &
        & nsp_aq_all,chraq_all,def_rain &! input
        & ,rain_all &! output
        & )
        implicit none

        integer,intent(in):: nsp_aq_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        real(kind=8),dimension(nsp_aq_all),intent(out)::rain_all
        real(kind=8),intent(in)::def_rain 
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i,n_tmp

        file_name = './rain.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1

        ! in default 
        rain_all = def_rain

        if (n_tmp <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,n_tmp
            read(50,*) chr_tmp,val_tmp
            if (any(chraq_all == chr_tmp)) then 
                rain_all(findloc(chraq_all,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_rainwater

    subroutine get_dust( &
        & nsp_sld_all,chrsld_all,def_dust &! input
        & ,dust_frct_all &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_sld_all),intent(out)::dust_frct_all
        real(kind=8),intent(in)::def_dust 
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i,n_tmp

        file_name = './dust.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1

        ! in default 
        dust_frct_all = def_dust

        if (n_tmp <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,n_tmp
            read(50,*) chr_tmp,val_tmp
            if (any(chrsld_all == chr_tmp)) then 
                dust_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_dust

    subroutine get_dust_2nd( &
        & nsp_sld_all,chrsld_all,def_dust &! input
        & ,dust_frct_all &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_sld_all),intent(out)::dust_frct_all
        real(kind=8),intent(in)::def_dust 
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i,n_tmp

        file_name = './dust_2nd.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1

        ! in default 
        dust_frct_all = def_dust

        if (n_tmp <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,n_tmp
            read(50,*) chr_tmp,val_tmp
            if (any(chrsld_all == chr_tmp)) then 
                dust_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_dust_2nd

    subroutine get_OM_rain( &
        & nsp_sld_all,chrsld_all,def_OM_frc &! input
        & ,OM_frct_all &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_sld_all),intent(out)::OM_frct_all
        real(kind=8),intent(in)::def_OM_frc 
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i,n_tmp

        file_name = './OM_rain.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1

        ! in default 
        OM_frct_all = def_OM_frc

        if (n_tmp <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,n_tmp
            read(50,*) chr_tmp,val_tmp
            if (any(chrsld_all == chr_tmp)) then 
                OM_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_OM_rain

    subroutine get_parentrock( &
        & nsp_sld_all,chrsld_all,def_pr &! input
        & ,parentrock_frct_all &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        real(kind=8),dimension(nsp_sld_all),intent(out)::parentrock_frct_all
        real(kind=8),intent(in)::def_pr 
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i,n_tmp

        file_name = './parentrock.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1

        ! in default 
        parentrock_frct_all = def_pr

        if (n_tmp <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,n_tmp
            read(50,*) chr_tmp,val_tmp
            if (any(chrsld_all == chr_tmp)) then 
                parentrock_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_parentrock

    subroutine get_atm( &
        & nsp_gas_all,chrgas_all &! input
        & ,atm_all &! output
        & )
        implicit none

        integer,intent(in):: nsp_gas_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        real(kind=8),dimension(nsp_gas_all),intent(out)::atm_all
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i,n_tmp

        file_name = './atm.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1

        ! in default 
        atm_all(findloc(chrgas_all,'po2',dim=1)) = 0.21d0
        atm_all(findloc(chrgas_all,'pco2',dim=1)) = 10d0**(-3.5d0)
        atm_all(findloc(chrgas_all,'pnh3',dim=1)) = 1d-9
        atm_all(findloc(chrgas_all,'pn2o',dim=1)) = 270d-9

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,n_tmp
            read(50,*) chr_tmp,val_tmp
            if (any(chrgas_all == chr_tmp)) then 
                atm_all(findloc(chrgas_all,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_atm

    subroutine get_switches( &
        & iwtype,imixtype,poroiter_in,display,display_lim_in,read_data,incld_rough &
        & ,act_ON,timestep_fixed,ads_ON,regular_grid,aq_close &! inout
        & ,poroevol,surfevol1,surfevol2,do_psd,lim_minsld_in,do_psd_full,season &! inout
        & )
        implicit none

        character(100) chr_tmp
        logical,intent(inout):: poroiter_in,display,display_lim_in,read_data,incld_rough &
            & ,act_ON,timestep_fixed,ads_ON,regular_grid,aq_close &
            & ,poroevol,surfevol1,surfevol2,do_psd,lim_minsld_in,do_psd_full,season
        integer,intent(out) :: imixtype,iwtype

        character(500) file_name
        integer i,n_tmp

        file_name = './switches.in'

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')

        read(50,*) iwtype,chr_tmp
        read(50,*) imixtype,chr_tmp
        read(50,*) poroiter_in,chr_tmp
        read(50,*) lim_minsld_in,chr_tmp
        read(50,*) display,chr_tmp
        read(50,*) display_lim_in,chr_tmp
        read(50,*) read_data,chr_tmp
        read(50,*) incld_rough,chr_tmp
        read(50,*) act_ON,chr_tmp
        read(50,*) timestep_fixed,chr_tmp
        read(50,*) ads_ON,chr_tmp
        read(50,*) regular_grid,chr_tmp
        read(50,*) aq_close,chr_tmp
        read(50,*) poroevol,chr_tmp
        read(50,*) surfevol1,chr_tmp
        read(50,*) surfevol2,chr_tmp
        read(50,*) do_psd,chr_tmp
        read(50,*) do_psd_full,chr_tmp
        read(50,*) season,chr_tmp

        close(50)


    endsubroutine get_switches

    subroutine get_sa_num(nsld_sa_dum)
        implicit none

        integer,intent(out):: nsld_sa_dum

        character(500) file_name
        integer n_tmp

        file_name = './sa.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1
        nsld_sa_dum = n_tmp


    endsubroutine get_sa_num

    subroutine get_sa( &
        & nsp_sld,chrsld,def_hr,nsld_sa &! input
        & ,hrii,chrsld_sa_dum &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld,nsld_sa
        character(5),dimension(nsp_sld),intent(in)::chrsld
        character(5),dimension(nsld_sa),intent(out)::chrsld_sa_dum
        real(kind=8),dimension(nsp_sld),intent(out)::hrii
        real(kind=8),intent(in)::def_hr 
        character(5) chr_tmp
        real(kind=8) val_tmp

        character(500) file_name
        integer i

        file_name = './sa.in'
        ! in default 
        hrii = def_hr

        if (nsld_sa <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,nsld_sa
            read(50,*) chr_tmp,val_tmp
            chrsld_sa_dum(i) = chr_tmp
            if (any(chrsld == chr_tmp)) then 
                hrii(findloc(chrsld,chr_tmp,dim=1)) = val_tmp
            endif 
        enddo 
        close(50)


    endsubroutine get_sa

    subroutine get_psdrain_num(nps_rain_char)
        implicit none

        integer,intent(out):: nps_rain_char

        character(500) file_name
        integer n_tmp

        file_name = './psdrain.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1
        nps_rain_char = n_tmp


    endsubroutine get_psdrain_num

    subroutine get_psdrain( &
        & nps_rain_char &! input
        & ,psu_rain_list,pssigma_rain_list,psw_rain_list &! output
        & )
        implicit none

        integer,intent(in):: nps_rain_char
        real(kind=8),dimension(nps_rain_char),intent(out)::psu_rain_list,pssigma_rain_list,psw_rain_list
        real(kind=8),dimension(3)::val_tmp

        character(500) file_name
        integer i

        file_name = './psdrain.in'

        if (nps_rain_char <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,nps_rain_char
            read(50,*) val_tmp(:)
            psu_rain_list(i)        = val_tmp(1)
            pssigma_rain_list(i)    = val_tmp(2)
            psw_rain_list(i)        = val_tmp(3)
        enddo 
        close(50)


    endsubroutine get_psdrain

    subroutine get_cec_num(nsld_cec_dum)
        implicit none

        integer,intent(out):: nsld_cec_dum

        character(500) file_name
        integer n_tmp

        file_name = './cec.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1
        nsld_cec_dum = n_tmp


    endsubroutine get_cec_num

    subroutine get_cec( &
        & nsp_sld_all,chrsld_all,mcec_def,nsld_cec,nsp_aq_all,chraq_all,logkhaq_def,beta_def &! input
        & ,mcec,chrsld_cec_dum,logkhaq,beta  &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld_all,nsld_cec,nsp_aq_all
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsld_cec),intent(out)::chrsld_cec_dum
        real(kind=8),dimension(nsp_sld_all),intent(in)::mcec_def,beta_def
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::logkhaq_def
        real(kind=8),dimension(nsp_sld_all),intent(out)::mcec,beta
        real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(out)::logkhaq
        character(5) chr_tmp
        real(kind=8) val_tmp,val_tmp3
        real(kind=8),dimension(5):: val_tmp2

        character(500) file_name
        integer i

        file_name = './cec.in'
        ! in default 
        mcec = mcec_def
        logkhaq = logkhaq_def
        beta = beta_def

        if (nsld_cec <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,nsld_cec
            read(50,*) chr_tmp,val_tmp,val_tmp2(:),val_tmp3
            chrsld_cec_dum(i) = chr_tmp
            if (any(chrsld_all == chr_tmp)) then 
                mcec(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
                logkhaq(findloc(chrsld_all,chr_tmp,dim=1),findloc(chraq_all,'na',dim=1)) = val_tmp2(1) ! Na
                logkhaq(findloc(chrsld_all,chr_tmp,dim=1),findloc(chraq_all,'k',dim=1))  = val_tmp2(2) ! K
                logkhaq(findloc(chrsld_all,chr_tmp,dim=1),findloc(chraq_all,'ca',dim=1)) = val_tmp2(3) ! Ca
                logkhaq(findloc(chrsld_all,chr_tmp,dim=1),findloc(chraq_all,'mg',dim=1)) = val_tmp2(4) ! Mg
                logkhaq(findloc(chrsld_all,chr_tmp,dim=1),findloc(chraq_all,'al',dim=1)) = val_tmp2(5) ! Al
                beta(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp3
            endif 
        enddo 
        close(50)


    endsubroutine get_cec

    subroutine get_nopsd_num(nsld_nopsd_dum)
        implicit none

        integer,intent(out):: nsld_nopsd_dum

        character(500) file_name
        integer n_tmp

        file_name = './nopsd.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1
        nsld_nopsd_dum = n_tmp


    endsubroutine get_nopsd_num

    subroutine get_nopsd( &
        & nsp_sld,chrsld,nsld_nopsd &! input
        & ,chrsld_nopsd_dum &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld,nsld_nopsd
        character(5),dimension(nsp_sld),intent(in)::chrsld
        character(5),dimension(nsld_nopsd),intent(out)::chrsld_nopsd_dum
        character(5) chr_tmp

        character(500) file_name
        integer i

        file_name = './nopsd.in'

        if (nsld_nopsd <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,nsld_nopsd
            read(50,*) chr_tmp
            chrsld_nopsd_dum(i) = chr_tmp
        enddo 
        close(50)


    endsubroutine get_nopsd

    subroutine get_2ndsld_num(nsp_sld_2)
        implicit none

        integer,intent(out):: nsp_sld_2

        character(500) file_name
        integer n_tmp

        file_name = './2ndslds.in'
        call Console4(file_name,n_tmp)

        n_tmp = n_tmp - 1
        nsp_sld_2 = n_tmp


        endsubroutine get_2ndsld_num

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        subroutine get_2ndsld( &
            & nsp_sld_2 &! input
            & ,chrsld_2 &! output
            & )
        implicit none

        integer,intent(in):: nsp_sld_2
        character(5),dimension(nsp_sld_2),intent(out)::chrsld_2
        character(5) chr_tmp

        character(500) file_name
        integer i

        file_name = './2ndslds.in'

        if (nsp_sld_2 <= 0) return

        open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
        read(50,'()')
        do i =1,nsp_sld_2
            read(50,*) chr_tmp
            chrsld_2(i) = chr_tmp
        enddo 
        close(50)


    endsubroutine get_2ndsld

end module scepter_input 