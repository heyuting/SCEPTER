module scepter_IO
    use scepter_constants
    use scepter_variables
    implicit none

    subroutine get_variables_num( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! output
        & )
        implicit none

        integer,intent(out):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nsld_kinspc
        character(500) file_name

        file_name = './slds.in'
        call Console4(file_name,nsp_sld)
        file_name = './solutes.in'
        call Console4(file_name,nsp_aq)
        file_name = './gases.in'
        call Console4(file_name,nsp_gas)
        file_name = './extrxns.in'
        call Console4(file_name,nrxn_ext)
        file_name = './kinspc.in'
        call Console4(file_name,nsld_kinspc)

        nsp_sld = nsp_sld - 1
        nsp_aq = nsp_aq - 1
        nsp_gas = nsp_gas - 1
        nrxn_ext = nrxn_ext - 1
        nsld_kinspc = nsld_kinspc - 1

    endsubroutine get_variables_num

    subroutine get_variables( &
        & nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc &! input
        & ,chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc &! output
        & )
        implicit none

        integer,intent(in):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nsld_kinspc
        character(5),dimension(nsp_sld),intent(out)::chrsld 
        character(5),dimension(nsp_aq),intent(out)::chraq 
        character(5),dimension(nsp_gas),intent(out)::chrgas 
        character(5),dimension(nrxn_ext),intent(out)::chrrxn_ext 
        character(5),dimension(nsld_kinspc),intent(out)::chrsld_kinspc
        real(kind=8),dimension(nsld_kinspc),intent(out)::kin_sld_spc

        character(500) file_name
        integer ispa,ispg,isps,irxn,isldspc

        if (nsp_aq>=1) then 
            file_name = './solutes.in'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do ispa =1,nsp_aq
                read(50,*) chraq(ispa) 
            enddo 
            close(50)
        endif 

        if (nsp_sld>=1) then 
            file_name = './slds.in'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do isps =1,nsp_sld
                read(50,*) chrsld(isps) 
            enddo  
            close(50)
        endif 

        if (nsp_gas>=1) then 
            file_name = './gases.in'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do ispg =1,nsp_gas
                read(50,*) chrgas(ispg) 
            enddo 
            close(50)
        endif 

        if (nrxn_ext>=1) then 
            file_name = './extrxns.in'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do irxn =1,nrxn_ext
                read(50,*) chrrxn_ext(irxn) 
            enddo 
            close(50)
        endif 

        if (nsld_kinspc>=1) then 
            file_name = './kinspc.in'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do isldspc =1,nsld_kinspc
                read(50,*) chrsld_kinspc(isldspc), kin_sld_spc(isldspc) 
            enddo 
            close(50)
        endif 

    endsubroutine get_variables


    subroutine get_saved_variables_num( &
        & workdir,runname_save &! input 
        & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc,nsld_sa_save &! output
        & )
        implicit none

        integer,intent(out):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nsld_kinspc,nsld_sa_save
        character(256),intent(in):: workdir,runname_save
        character(500) file_name

        file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/slds.save'
        call Console4(file_name,nsp_sld)
        file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/solutes.save'
        call Console4(file_name,nsp_aq)
        file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/gases.save'
        call Console4(file_name,nsp_gas)
        file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/extrxns.save'
        call Console4(file_name,nrxn_ext)
        file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/kinspc.save'
        call Console4(file_name,nsld_kinspc)
        file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/sa.save'
        call Console4(file_name,nsld_sa_save)

        nsp_sld = nsp_sld - 1
        nsp_aq = nsp_aq - 1
        nsp_gas = nsp_gas - 1
        nrxn_ext = nrxn_ext - 1
        nsld_kinspc = nsld_kinspc - 1
        nsld_sa_save = nsld_sa_save - 1

    endsubroutine get_saved_variables_num

    subroutine get_saved_variables( &
        & workdir,runname_save &! input 
        & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc,nsld_sa_dum &! input
        & ,chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc &! output
        & ,chrsld_sa_dum,hrii_dum &! output 
        & )
        implicit none

        integer,intent(in):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nsld_kinspc,nsld_sa_dum
        character(5),dimension(nsp_sld),intent(out)::chrsld 
        character(5),dimension(nsp_aq),intent(out)::chraq 
        character(5),dimension(nsp_gas),intent(out)::chrgas 
        character(5),dimension(nrxn_ext),intent(out)::chrrxn_ext 
        character(5),dimension(nsld_kinspc),intent(out)::chrsld_kinspc 
        character(5),dimension(nsld_sa_dum),intent(out)::chrsld_sa_dum 
        character(256),intent(in):: workdir,runname_save
        real(kind=8),dimension(nsld_kinspc),intent(out)::kin_sld_spc
        real(kind=8),dimension(nsld_sa_dum),intent(out)::hrii_dum

        character(500) file_name
        integer ispa,ispg,isps,irxn,isldspc,isldsa

        if (nsp_aq>=1) then 
            file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/solutes.save'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do ispa =1,nsp_aq
                read(50,*) chraq(ispa) 
            enddo 
            close(50)
        endif 

        if (nsp_sld>=1) then 
            file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/slds.save'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do isps =1,nsp_sld
                read(50,*) chrsld(isps) 
            enddo  
            close(50)
        endif 

        if (nsp_gas>=1) then 
            file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/gases.save'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do ispg =1,nsp_gas
                read(50,*) chrgas(ispg) 
            enddo 
            close(50)
        endif 

        if (nrxn_ext>=1) then 
            file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/extrxns.save'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do irxn =1,nrxn_ext
                read(50,*) chrrxn_ext(irxn) 
            enddo 
            close(50)
        endif 

        if (nsld_kinspc>=1) then 
            file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/kinspc.save'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do isldspc =1,nsld_kinspc
                read(50,*) chrsld_kinspc(isldspc), kin_sld_spc(isldspc) 
            enddo 
            close(50)
        endif 

        if (nsld_sa_dum>=1) then 
            file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/sa.save'
            open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
            read(50,'()')
            do isldsa =1,nsld_sa_dum
                read(50,*) chrsld_sa_dum(isldsa), hrii_dum(isldsa) 
            enddo 
            close(50)
        endif 

    endsubroutine get_saved_variables

    subroutine Console4(file_name,i)

        implicit none

        integer,intent(out) :: i
        character(500),intent(in)::file_name

        open(9, file =trim(adjustl(file_name)))

        i = 0
        do 
            read(9, *, end = 99)
            i = i + 1
        end do 

        ! 99 print *, i
        99 continue
        close(9)

    end subroutine Console4

end module scepter_IO
