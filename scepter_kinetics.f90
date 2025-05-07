module scepter_kinetics
    use scepter_constants
    use scepter_variables
    use scepter_equilibrium
    use scepter_transport
    implicit none
    private
    public :: sld_kin, sld_rxn, rough_f, calc_rxn_ext_dev_3
    ! Constants
    real(kind=8), parameter :: cal2j = 4.184d0

    contains
    subroutine sld_kin( &
        & nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv_tmp &! input
        & ,nsp_gas_all,chrgas_all,mgas_loc &! input
        & ,nsp_aq_all,chraq_all,maqf_loc &! input
        & ,mineral,dev_sp &! input 
        & ,kin,dkin_dmsp &! output
        & ) 
        implicit none

        integer,intent(in)::nz
        real(kind=8),intent(in)::rg,tc,sec2yr,tempk_0,mv_tmp,kw,kho
        real(kind=8),dimension(nz),intent(in)::prox

        character(5),intent(in)::mineral,dev_sp
        real(kind=8),dimension(nz),intent(out)::kin,dkin_dmsp
        real(kind=8) mh,moh,kinn_ref,kinh_ref,kinoh_ref,ean,eah,eaoh,tc_ref

        integer,intent(in)::nsp_gas_all,nsp_aq_all
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        real(kind=8),dimension(nsp_gas_all,nz),intent(in)::mgas_loc
        real(kind=8),dimension(nsp_aq_all,nz),intent(in)::maqf_loc

        real(kind=8),dimension(nz) :: pco2,alf
        real(kind=8) mco2,kinco2_ref,eaco2,q10,kref

        ! real(kind=8) k_arrhenius


        kin = 0d0
        dkin_dmsp = 0d0

        select case(trim(adjustl(mineral)))
            ! case('ka','al2o3') ! corundum dissolution rate is assumed to be the same as kaolinite (cf., Carroll-Webb and Walther, 1988)
            case('ka') ! corundum dissolution rate is assumed to be the same as kaolinite (cf., Carroll-Webb and Walther, 1988)
                mh = 0.777d0
                moh = -0.472d0
                kinn_ref = 10d0**(-13.18d0)*sec2yr
                kinh_ref = 10d0**(-11.31d0)*sec2yr
                kinoh_ref = 10d0**(-17.05d0)*sec2yr
                ean = 22.2d0
                eah = 65.9d0
                eaoh = 17.9d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('ab')
                mh = 0.457d0
                moh = -0.572d0
                kinn_ref = 10d0**(-12.56d0)*sec2yr
                kinh_ref = 10d0**(-10.16d0)*sec2yr
                kinoh_ref = 10d0**(-15.6d0)*sec2yr
                ean = 69.8d0
                eah = 65d0
                eaoh = 71d0
                ! above is from table 13 ( this could be a confused mixture of linear and non-linear regression parameters in Table 1)
                ! following is linear regression result of Table 1
                mh = 0.457d0
                moh = -0.572d0
                kinn_ref = 10d0**(-12.04d0)*sec2yr
                kinh_ref = 10d0**(-9.87d0)*sec2yr
                kinoh_ref = 10d0**(-16.98d0)*sec2yr
                ean = 69.8d0
                eah = 65d0
                eaoh = 71d0
                ! then non-linear regression
                ! mh = 0.317d0
                ! moh = -0.471d0
                ! kinn_ref = 10d0**(-12.56d0)*sec2yr
                ! kinh_ref = 10d0**(-10.16d0)*sec2yr
                ! kinoh_ref = 10d0**(-15.6d0)*sec2yr
                ! ean = 65d0
                ! eah = 65d0
                ! eaoh = 66.5d0
                !!! 
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('kfs','sdn') ! sanidine assumed to follow rate law of K-feldspar
                mh = 0.5d0
                moh = -0.823d0
                kinn_ref = 10d0**(-12.41d0)*sec2yr
                kinh_ref = 10d0**(-10.06d0)*sec2yr
                ! kinoh_ref = 10d0**(-9.68d0)*sec2yr*kw**(-moh)
                kinoh_ref = 10d0**(-21.2d0)*sec2yr
                ! ean = 9.08*cal2j
                ! eah = 12.4d0*cal2j
                ! eaoh = 22.5d0*cal2j
                ean = 38d0
                eah = 51.7d0
                eaoh = 94.1d0
                tc_ref = 25d0
                ! Brantley et al 2008 used data from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('leu') 
                mh = 0.7d0
                moh = -0.2d0
                kinn_ref = 10d0**(-9.20d0)*sec2yr
                kinh_ref = 10d0**(-6.00d0)*sec2yr
                kinoh_ref = 10d0**(-10.66d0)*sec2yr
                ean = 75.5d0
                eah = 132.2d0
                eaoh = 56.6d0
                tc_ref = 25d0
                ! Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('fo')
                mh = 0.47d0
                moh = 0d0
                kinn_ref = 10d0**(-10.64d0)*sec2yr
                kinh_ref = 10d0**(-6.85d0)*sec2yr
                kinoh_ref = 0d0
                ean = 79d0
                eah = 67.2d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('fa')
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-12.80d0)*sec2yr
                kinh_ref = 10d0**(-4.80d0)*sec2yr
                kinoh_ref = 0d0
                ean = 94.4d0
                eah = 94.4d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('an')
                mh = 1.411d0
                moh = 0d0
                kinn_ref = 10d0**(-9.12d0)*sec2yr
                kinh_ref = 10d0**(-3.5d0)*sec2yr
                kinoh_ref = 0d0
                ean = 17.8d0
                eah = 16.6d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('la')
                mh = 0.626d0
                moh = 0d0
                moh = -0.57d0
                kinn_ref = 10d0**(-10.91d0)*sec2yr
                kinh_ref = 10d0**(-7.87d0)*sec2yr
                kinoh_ref = 0d0 ! original data 
                ! kinoh_ref = 10d0**(-15.57d0)*sec2yr ! added by Beerling et al. 2020 from albite data of Palandri and Kharaka 2004
                ean = 45.2d0
                eah = 42.1d0
                eaoh = 0d0 ! original data
                ! eaoh = 71d0 ! Beerling
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('and')
                mh = 0.541d0
                moh = 0d0
                kinn_ref = 10d0**(-11.47d0)*sec2yr
                kinh_ref = 10d0**(-8.88d0)*sec2yr
                kinoh_ref = 0d0
                ean = 57.4d0
                eah = 53.5d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('olg')
                mh = 0.457d0
                moh = 0d0
                kinn_ref = 10d0**(-11.84d0)*sec2yr
                kinh_ref = 10d0**(-9.67d0)*sec2yr
                kinoh_ref = 0d0
                ean = 69.8d0
                eah = 65.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('by')
                mh = 1.018d0
                moh = 0d0
                kinn_ref = 10d0**(-9.82d0)*sec2yr
                kinh_ref = 10d0**(-5.85d0)*sec2yr
                kinoh_ref = 0d0
                ean = 31.5d0
                eah = 29.3d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('cc')
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-5.81d0)*sec2yr
                kinh_ref = 10d0**(-0.3d0)*sec2yr
                kinoh_ref = 0d0
                ean = 23.5d0
                eah = 14.4d0
                eaoh = 0d0
                tc_ref = 25d0
                ! adding co2 mechanism
                mco2 = 1d0
                kinco2_ref = 10d0**(-3.48d0)*sec2yr
                eaco2 = 35.4d0
                pco2 = mgas_loc(findloc(chrgas_all,'pco2',dim=1),:)
                ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & + pco2**mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case('pco2')
                        dkin_dmsp = ( & 
                            & + mco2*pco2**(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('arg')
                ! assumed to be the same as those for cc
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-5.81d0)*sec2yr
                kinh_ref = 10d0**(-0.3d0)*sec2yr
                kinoh_ref = 0d0
                ean = 23.5d0
                eah = 14.4d0
                eaoh = 0d0
                tc_ref = 25d0
                ! adding co2 mechanism
                mco2 = 1d0
                kinco2_ref = 10d0**(-3.48d0)*sec2yr
                eaco2 = 35.4d0
                pco2 = mgas_loc(findloc(chrgas_all,'pco2',dim=1),:)
                ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & + pco2**mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case('pco2')
                        dkin_dmsp = ( & 
                            & + mco2*pco2**(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('dlm') ! for disordered dolomite
                mh = 0.500d0
                moh = 0d0
                kinn_ref = 10d0**(-7.53d0)*sec2yr
                kinh_ref = 10d0**(-3.19d0)*sec2yr
                kinoh_ref = 0d0
                ean = 52.2d0
                eah = 36.1d0
                eaoh = 0d0
                tc_ref = 25d0
                ! adding co2 mechanism
                mco2 = 0.5d0
                kinco2_ref = 10d0**(-5.11d0)*sec2yr
                eaco2 = 34.8d0
                pco2 = mgas_loc(findloc(chrgas_all,'pco2',dim=1),:)
                ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & + pco2**mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case('pco2')
                        dkin_dmsp = ( & 
                            & + mco2*pco2**(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('gb','amal')
                mh = 0.992d0
                moh = -0.784d0
                kinn_ref = 10d0**(-11.50d0)*sec2yr
                kinh_ref = 10d0**(-7.65d0)*sec2yr
                kinoh_ref = 10d0**(-16.65d0)*sec2yr
                ean = 61.2d0
                eah = 47.5d0
                eaoh = 80.1d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('amsi')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-12.23d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 74.5d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('phsi')
                mh = 1d0
                moh = -0.33d0
                kinn_ref = 5d-18*1d4*sec2yr ! given as mol/cm2/sec; converted to mol/m2/yr
                kinh_ref = 6d-16*1d4*sec2yr ! given as mol/cm2/sec; converted to mol/m2/yr
                kinoh_ref = 3.5d-13*1d4*sec2yr *kw**(-moh) ! given as mol/cm2/sec; converted to mol/m2/yr
                ean = 74.5d0 
                eah = 74.5d0 
                eaoh = 74.5d0 
                tc_ref = 25d0
                ! from Fraysse et al. 2009 except for Ea which is taken from value for Palandri and Kharaka 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('qtz')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-13.40d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 90.9d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                mh = 0.309d0 ! brantley 2008
                moh = -0.411d0 ! brantley 2008
                kinn_ref = 10d0**(-13.40d0)*sec2yr
                kinh_ref = 4.34d-12 *sec2yr ! brantley 2008
                kinoh_ref = 6.06d-10 * kw**(-moh) *sec2yr ! brantley 2008
                ean = 90.9d0
                eah = 90.9d0 ! assumed to be the same as ean
                eaoh = 90.9d0 ! assumed to be the same as ean
                tc_ref = 25d0
                ! pH neutral range from Palandri and Kharaka, 2004 pH dependence from Brantley et al 2008
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('gt','amfe3')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-7.94d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 86.5d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('hm')
                mh = 1d0
                moh = 0d0
                kinn_ref = 10d0**(-14.60d0)*sec2yr
                kinh_ref = 10d0**(-9.39d0)*sec2yr
                kinoh_ref = 0d0
                ean = 66.2d0
                eah = 66.2d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                dkin_dmsp = 0d0

            case('ct')
                mh = 0d0
                moh = -0.23d0
                kinn_ref = 10d0**(-12d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 10d0**(-13.58d0)*sec2yr
                ean = 73.5d0
                eah = 0d0
                eaoh = 73.5d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('mscv')
                mh = 0.370d0
                moh = -0.22d0
                kinn_ref = 10d0**(-13.55d0)*sec2yr
                kinh_ref = 10d0**(-11.85d0)*sec2yr
                kinoh_ref = 10d0**(-13.55d0)*sec2yr
                ean = 22d0
                eah = 22d0
                eaoh = 22d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('plgp')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-12.4d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 29d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('cabd','ill','kbd','nabd','mgbd','casp','ksp','nasp','mgsp') ! illite kinetics is assumed to be the same as smectite (Bibi et al., 2011)
                mh = 0.34d0
                moh = -0.4d0
                kinn_ref = 10d0**(-12.78d0)*sec2yr
                kinh_ref = 10d0**(-10.98d0)*sec2yr
                kinoh_ref = 10d0**(-16.52d0)*sec2yr
                ean = 35d0
                eah = 23.6d0
                eaoh = 58.9d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('nph','anl') ! analcime kinetics is assumed to be the same as nepherine (cf. Ragnarsdottir, GCA, 1993)
                mh = 1.130d0
                moh = -0.200d0
                kinn_ref = 10d0**(-8.56d0)*sec2yr
                kinh_ref = 10d0**(-2.73d0)*sec2yr
                kinoh_ref = 10d0**(-10.76d0)*sec2yr
                ean = 65.4d0
                eah = 62.9d0
                eaoh = 37.8d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 

            case('dp')
                mh = 0.71d0
                moh = 0d0
                kinn_ref = 10d0**(-11.11d0)*sec2yr
                kinh_ref = 10d0**(-6.36d0)*sec2yr
                kinoh_ref = 0d0
                ean = 50.6d0
                eah = 96.1d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('hb','cpx','agt')
                mh = 0.70d0
                moh = 0d0
                kinn_ref = 10d0**(-11.97d0)*sec2yr
                kinh_ref = 10d0**(-6.82d0)*sec2yr
                kinoh_ref = 0d0
                ean = 78.0d0
                eah = 78.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for augite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('en','opx','fer')
                mh = 0.60d0
                moh = 0d0
                kinn_ref = 10d0**(-12.72d0)*sec2yr
                kinh_ref = 10d0**(-9.02d0)*sec2yr
                kinoh_ref = 0d0
                ean = 80.0d0
                eah = 80.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for enstatite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('jd')
                mh = 0.70d0
                moh = 0d0
                kinn_ref = 10d0**(-9.50d0)*sec2yr
                kinh_ref = 10d0**(-6.00d0)*sec2yr
                kinoh_ref = 0d0
                ean = 94.4d0
                eah = 132.2d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for jadeite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('wls')
                mh = 0.40d0
                moh = 0d0
                kinn_ref = 10d0**(-8.88d0)*sec2yr
                kinh_ref = 10d0**(-5.37d0)*sec2yr
                kinoh_ref = 0d0
                ean = 54.7d0
                eah = 54.7d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for wollastonite from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('tm')
                mh = 0.70d0
                moh = 0d0
                kinn_ref = 10d0**(-10.60d0)*sec2yr
                kinh_ref = 10d0**(-8.40d0)*sec2yr
                kinoh_ref = 0d0
                ean = 94.4d0
                eah = 18.9d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('ep')
                mh = 0.338d0
                moh = -0.556d0
                kinn_ref = 10d0**(-11.99d0)*sec2yr
                kinh_ref = 10d0**(-10.60d0)*sec2yr
                kinoh_ref = 10d0**(-17.33d0)*sec2yr
                ean = 70.7d0
                eah = 71.1d0
                eaoh = 79.1d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('clch')
                mh = 0.5d0
                moh = 00
                kinn_ref = 10d0**(-12.52d0)*sec2yr
                kinh_ref = 10d0**(-11.11d0)*sec2yr
                kinoh_ref = 0d0
                ean = 88.0d0
                eah = 88.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('cdr')
                mh = 1.0d0
                moh = 00
                kinn_ref = 10d0**(-11.20d0)*sec2yr
                kinh_ref = 10d0**(-3.8d0)*sec2yr
                kinoh_ref = 0d0
                ean = 28.3d0
                eah = 113.3d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('antp','splt')
                mh = 0.440d0
                moh = 0d0
                kinn_ref = 10d0**(-14.24d0)*sec2yr
                kinh_ref = 10d0**(-11.94d0)*sec2yr
                kinoh_ref = 0d0
                ean = 51.0d0
                eah = 51.0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! for anthophyllite from Palandri and Kharaka, 2004
                ! sepiolite dissolution is assumed to be close to anthophyllite dissolution (Mulders and Oelkers 2020 GCA)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('gps')
                mh = 0d0
                moh = 0d0
                kinn_ref = 10d0**(-2.79d0)*sec2yr
                kinh_ref = 0d0
                kinoh_ref = 0d0
                ean = 0d0
                eah = 0d0
                eaoh = 0d0
                tc_ref = 25d0
                ! from Palandri and Kharaka, 2004
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            ! case('caso4') ! from Palandri and Kharaka, 2004
                ! mh = 0d0
                ! moh = 0d0
                ! kinn_ref = 10d0**(-3.19d0)*sec2yr
                ! kinh_ref = 0d0
                ! kinoh_ref = 0d0
                ! ean = 14.3d0
                ! eah = 0d0
                ! eaoh = 0d0
                ! tc_ref = 25d0
                ! kin = ( & 
                    ! & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    ! & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    ! & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    ! & ) 
                ! select case(trim(adjustl(dev_sp)))
                    ! case('pro')
                        ! dkin_dmsp = ( & 
                            ! & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            ! & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            ! & ) 
                    ! case default 
                        ! dkin_dmsp = 0d0
                ! endselect 
                
            case('fe2o','mgo','k2o','cao','na2o','al2o3','sio2','caso4')
                kin = ( &
                    & 1d0/0.01d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 1 year as in Chen et al. (2010, AFM) 
                    & )
                ! assuming randomly fast dissolution rate for oxides
                dkin_dmsp = 0d0
            
            case('gbas')
                mh = 1.013d0
                moh = -0.258d0
                kinn_ref = 0d0
                kinh_ref = 10d0**(-4.27d0)*sec2yr
                kinoh_ref = 10d0**(-11.00d0)*sec2yr
                ean = 0d0
                eah = 39.7d0
                eaoh = 38.4d0
                tc_ref = 25d0
                ! from Pollyea and Rimstidt, 2017 (this is normalized to geometric surface area)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                ! from Brantley et al. (2008)
                mh = 1.16d0
                moh = -0.16d0
                eah = 47.5d0
                eaoh = 47.5d0
                kinh_ref = 588d0*sec2yr
                kinoh_ref = 0.0822d0*kw**(-moh)*sec2yr
                kin = ( &
                    & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**mh &
                    & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**moh &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                        ! from Brantley et al. (2008)
                        dkin_dmsp = ( &
                            & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**(mh-1d0)*mh &
                            & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**(moh-1d0)*moh &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
            
            case('cbas')
                mh = 0.68d0
                moh = -0.286d0
                kinn_ref = 0d0
                kinh_ref = 10d0**(-6.15d0)*sec2yr
                kinoh_ref = 10d0**(-11.83d0)*sec2yr
                ean = 0d0
                eah = 40.1d0
                eaoh = 32.9d0
                tc_ref = 25d0
                ! from Pollyea and Rimstidt, 2017 (this is normalized to geometric surface area)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                ! ----- from Brantley et al. (2008)
                mh = 1.16d0
                moh = -0.16d0
                eah = 47.5d0
                eaoh = 47.5d0
                kinh_ref = 588d0*sec2yr
                kinoh_ref = 0.0822d0*kw**(-moh)*sec2yr
                kin = ( &
                    & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**mh &
                    & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**moh &
                    & ) 
                ! ----- from Oelkers and Gislason (2001) & Gislason & Oelkers (2003)
                alf = maqf_loc(findloc(chraq_all,'al',dim=1),:)
                eah = 25.5d0 ! from Gislason and Oelkers (2003) (but for geometric surface based equation)
                tc_ref = 25d0
                mh = 0.35d0 ! OG01
                ! mh = 1d0/3d0 ! GO03
                kinh_ref = 10d0**(-11.64d0)*1d4*sec2yr ! units converted from mol/cm2/sec to mol/m2/yr (OG01)
                ! kinh_ref = 10d0**(-5.6d0)*1d4/92d0*sec2yr ! units converted from mol/cm2/sec to mol/m2/yr (GO03); 92 is surface roughness
                kin = (&
                    & k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg)*(prox**3d0/alf)**mh &
                    ! & kinh_ref*exp(-eah/(rg*(tc+tempk_0)))*(prox**3d0/alf)**mh &
                    & )
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                        ! ----- from Brantley et al. (2008)
                        dkin_dmsp = ( &
                            & kinh_ref *exp(-eah/(rg*(tc+tempk_0)))*prox**(mh-1d0)*mh &
                            & + kinoh_ref *exp(-eaoh/(rg*(tc+tempk_0)))*prox**(moh-1d0)*moh &
                            & ) 
                        ! ----- from Oelkers and Gislason (2001)
                        dkin_dmsp = ( &
                            & k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg)*prox**(3d0*mh-1d0)*3d0*mh *(1d0/alf)**mh &
                            ! & kinh_ref*exp(-eah/(rg*(tc+tempk_0)))*prox**(3d0*mh-1d0)*3d0*mh *(1d0/alf)**mh &
                            & ) 
                    case('al')
                        ! ----- from Oelkers and Gislason (2001)
                        dkin_dmsp = ( &
                            & k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg)*prox**(3d0*mh)*(-mh) *(1d0/alf)**(mh+1d0) &
                            ! & kinh_ref*exp(-eah/(rg*(tc+tempk_0)))*prox**(3d0*mh)*(-mh) *(1d0/alf)**(mh+1d0) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
                
            case('py')
                mh = 0d0
                moh = -0.11d0
                kinn_ref = 0d0
                kinh_ref = 0d0
                kinoh_ref = 10.0d0**(-8.19d0)*sec2yr*kho**0.5d0
                ean = 0d0
                eah = 0d0
                eaoh = 57d0
                tc_ref = 15d0
                ! from Williamson and Rimstidt (1994)
                kin = ( & 
                    & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
                    & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
                select case(trim(adjustl(dev_sp)))
                    case('pro')
                        dkin_dmsp = ( & 
                            & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                            & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                            & ) 
                    case default 
                        dkin_dmsp = 0d0
                endselect 
                
            case('g1')
                kin = ( &
                    & 1d0/1d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 1 year as in Chen et al. (2010, AFM) 
                    & )
                ! adding temperature dependence in the form of Q10
                kref = 1d0/1d0
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                dkin_dmsp = 0d0
                
            case('g2')
                kin = ( &
                    & 1d0/8d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 8 year as in Chen et al. (2010, AFM) 
                    & )
                ! adding temperature dependence in the form of Q10
                kref = 1d0/8d0
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                dkin_dmsp = 0d0
                
            case('g3')
                kin = ( &
                    & 1d0/1d3 &! mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 1k year  
                    & )
                ! adding temperature dependence in the form of Q10
                kref = 1d0/1d3
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                dkin_dmsp = 0d0
                
            case('amnt')
                kin = ( &
                    & 1d0/0.01d0 &! just a value assumed; turnover time of 0.1 year for NH4NO3 
                    & )
                ! adding temperature dependence in the form of Q10
                ! kref = 0.01d0/1d0
                kref = 1d0/0.1d0
                tc_ref = 15d0
                q10 = 3d0
                kin = k_q10(kref,tc,tc_ref,q10)
                ! kin = kref
                dkin_dmsp = 0d0
                
            case('kcl','gac','mesmh','ims','teas','naoh','naglp','cacl2','nacl')
                kin = ( &
                    & 1d0/0.01d0 &! just a value assumed; turnover time of 0.1 year for NH4NO3 
                    & )
                ! adding temperature dependence in the form of Q10
                ! kref = 0.01d0/1d0
                kref = 1d0/0.0001d0
                tc_ref = 15d0
                q10 = 1d0
                kin = k_q10(kref,tc,tc_ref,q10)
                ! kin = kref
                dkin_dmsp = 0d0
                
            case default 
                kin =0d0
                dkin_dmsp = 0d0

        endselect  


    endsubroutine sld_kin

    subroutine sld_rxn( &
        & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
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

        real(kind=8)::auth_th = 1d2
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
        ! print *, 'in sld_rxn'
        ! print *, ksld(findloc(chrsld,'kfs',dim=1),:)

            
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
                                drxnsld_dmaq(isps,:,iz) = ( &
                                    & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmaq(isps,:,iz)) &
                                    & +dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                                    & )
                                drxnsld_dmgas(isps,:,iz) = ( &
                                    & +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmgas(isps,:,iz)) &
                                    & +dksld_dmgas(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                                    & )
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
                        & + min(ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:))/(1d0-poro),maxdis(isps,:)) &
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
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(-domega_dmaq(isps,ispa,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(-domega_dmgas(isps,ispg,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
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
                        & + ksld(isps,:)*poro**(2d0/3d0)*hr(isps,:)*(mv(isps)*1d-6*msldx(isps,:))**(2d0/3d0)*(1d0-omega(isps,:)) &
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
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmaq(isps,ispa,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & )
                    enddo 
                    
                    do ispg = 1, nsp_gas
                        drxnsld_dmgas(isps,ispg,:) = ( &
                            & + ksld(isps,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)*solmod(isps,:)) &
                            & *merge(0d0,1d0,1d0-omega(isps,:)*solmod(isps,:)*nonprec(isps,:) < 0d0) &
                            & + dksld_dmgas(isps,ispg,:)*poro*hr(isps,:)*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)*solmod(isps,:)) &
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

    subroutine calc_rxn_ext_dev_3( &
        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
        & ,rg,tempk_0,tc,iosx &!input
        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
        & ,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl  &! input 
        & ,rxn_name,sp_name &! input 
        & ,rxn_ext,drxnext_dmsp,rxnext_error &! output
        & )
        implicit none
        integer,intent(in)::nz
        real(kind=8),dimension(nz),intent(in):: poro,sat,prox,iosx
        real(kind=8),dimension(nz),intent(out):: drxnext_dmsp
        real(kind=8),dimension(nz),intent(out):: rxn_ext
        character(5),intent(in)::rxn_name,sp_name
        logical,intent(out)::rxnext_error

        integer,intent(in)::nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst 

        character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all
        character(5),dimension(nsp_gas),intent(in)::chrgas
        character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
        character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
        character(5),dimension(nsp_aq),intent(in)::chraq
        character(5),dimension(nsp_aq_all),intent(in)::chraq_all
        character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst

        real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
        real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
        real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
        real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
        real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
        real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
        real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all
        real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
        real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl
        real(kind=8),dimension(nrxn_ext_all,nz),intent(in)::krxn1_ext_all,krxn2_ext_all

        integer,intent(in)::nsp_sld,nsp_sld_cnst,nsp_sld_all

        character(5),dimension(nsp_sld),intent(in)::chrsld
        character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst
        character(5),dimension(nsp_sld_all),intent(in)::chrsld_all

        real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx
        real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc
        real(kind=8),dimension(nsp_sld_all),intent(in)::msldth_all,mv_all

        real(kind=8),intent(in)::rho_grain,kw,rg,tempk_0,tc

        real(kind=8),dimension(nsp_sld,nz),intent(in):: hr

        ! local variables
        real(kind=8):: po2th,fe2th,mwtom,g1th,g2th,g3th,mvpy,fe3th,knh3,k1nh3,ko2,v_tmp,km_tmp1,km_tmp2,km_tmp3  &
            & ,kn2o,k1fe2,k1fe2co3,k1fe2hco3,k1fe2so4,kco2,k1,k2
        real(kind=8),dimension(nz):: po2x,vmax,mo2,fe2x,koxa,vmax2,mom2,komb,beta,omx,ombx &
            & ,mo2g1,mo2g2,mo2g3,kg1,kg2,kg3,g1x,g2x,g3x,pyx,fe3x,koxpy,pnh3x,nh4x,dnh4_dpro,dnh4_dpnh3 &
            & ,no3x,pn2ox,dv_dph_tmp,fe2f,dfe2f_dfe2,dfe2f_dpco2,dfe2f_dpro,dfe2f_dso4f,pco2x,hrpy,oxax &
            & ,vmax_tmp
        real(kind=8),dimension(nsp_aq_all,nz)::maqx_loc,maqft_loc,dmaqft_dpro_loc,dmaqft_dios_loc
        real(kind=8),dimension(nsp_aq_all,nsp_aq_all,nz)::dmaqft_dmaqf_loc
        real(kind=8),dimension(nsp_aq_all,nsp_gas_all,nz)::dmaqft_dmgas_loc
        real(kind=8),dimension(nsp_gas_all,nz)::mgasx_loc
        real(kind=8),dimension(nsp_sld_all,nz)::msldx_loc

        ! real(kind=8):: thon = 1d0
        real(kind=8):: thon = -1d100

        integer ieqgas_h0,ieqgas_h1,ieqgas_h2
        data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

        integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
        data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

        integer ieqaq_co3,ieqaq_hco3
        data ieqaq_co3,ieqaq_hco3/1,2/

        integer ieqaq_so4,ieqaq_so42
        data ieqaq_so4,ieqaq_so42/1,2/

        character(5) sp_tmp
        character(25) scheme

        ! ... need to clean up the following mess at some day ...

        vmax = krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)
        mo2 = krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)

        po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

        koxa = krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) 

        fe2th = maqth_all(findloc(chraq_all,'fe2',dim=1))

        fe3th = maqth_all(findloc(chraq_all,'fe3',dim=1))

        g1th = msldth_all(findloc(chrsld_all,'g1',dim=1))
        g2th = msldth_all(findloc(chrsld_all,'g2',dim=1))
        g3th = msldth_all(findloc(chrsld_all,'g3',dim=1))

        ko2 = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)

        kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
        k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
        k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)

        knh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h0)
        k1nh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h1)

        kn2o = keqgas_h(findloc(chrgas_all,'pn2o',dim=1),ieqgas_h0)

        k1fe2 = keqaq_h(findloc(chraq_all,'fe2',dim=1),ieqaq_h1)
        k1fe2co3 = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_co3)
        k1fe2hco3  = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_hco3)
        k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)

        vmax2 = krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)
        mom2 = krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)

        komb = krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)
        beta = krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)

        koxpy = krxn1_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1),:)

        kg2 = krxn1_ext_all(findloc(chrrxn_ext_all,'g2k',dim=1),:)

        mvpy = mv_all(findloc(chrsld_all,'py',dim=1))

        hrpy = 0d0
        if (any(chrsld=='py')) then 
            hrpy = hr(findloc(chrsld,'py',dim=1),:)
        endif 


        ! getting mgasx_loc & maqx_loc
        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maqx,mgasx,maqc,mgasc &
            & ,maqx_loc,mgasx_loc  &! output
            & )

        ! getting maqft_loc and its derivatives
        call get_maqt_all( &
        ! call get_maqt_all_v2( &
            & nz,nsp_aq_all,nsp_gas_all &
            & ,chraq_all,chrgas_all &
            & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
            & ,mgasx_loc,maqx_loc,prox,iosx,tc &
            & ,dmaqft_dpro_loc,dmaqft_dmaqf_loc,dmaqft_dmgas_loc,dmaqft_dios_loc &! output
            & ,maqft_loc  &! output
            & )

        call get_msldx_all( &
            & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
            & ,chrsld,chrsld_all,chrsld_cnst &
            & ,msldx,msldc &
            & ,msldx_loc  &! output
            & )
            
        po2x    = mgasx_loc(findloc(chrgas_all,'po2',dim=1),:)
        pco2x   = mgasx_loc(findloc(chrgas_all,'pco2',dim=1),:)
        pnh3x   = mgasx_loc(findloc(chrgas_all,'pnh3',dim=1),:)
        pn2ox   = mgasx_loc(findloc(chrgas_all,'pn2o',dim=1),:)

        oxax    = maqx_loc(findloc(chraq_all,'oxa',dim=1),:)    *maqft_loc(findloc(chraq_all,'oxa',dim=1),:)
        fe2x    = maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    *maqft_loc(findloc(chraq_all,'fe2',dim=1),:)
        fe3x    = maqx_loc(findloc(chraq_all,'fe3',dim=1),:)    *maqft_loc(findloc(chraq_all,'fe3',dim=1),:)
        no3x    = maqx_loc(findloc(chraq_all,'no3',dim=1),:)    *maqft_loc(findloc(chraq_all,'no3',dim=1),:)
        fe2f    = maqx_loc(findloc(chraq_all,'fe2',dim=1),:)

        omx     = msldx_loc(findloc(chrsld_all,'om',dim=1),:)
        ombx    = msldx_loc(findloc(chrsld_all,'omb',dim=1),:)
        g1x     = msldx_loc(findloc(chrsld_all,'g1',dim=1),:)
        g2x     = msldx_loc(findloc(chrsld_all,'g2',dim=1),:)
        g3x     = msldx_loc(findloc(chrsld_all,'g3',dim=1),:)
        pyx     = msldx_loc(findloc(chrsld_all,'py',dim=1),:)


        nh4x = pnh3x*knh3*prox/k1nh3
        dnh4_dpro = pnh3x*knh3*1d0/k1nh3
        dnh4_dpnh3 = 1d0*knh3*prox/k1nh3

        select case(trim(adjustl(rxn_name)))

            case('resp')
                rxn_ext = vmax*po2x/(po2x+mo2)
                ! rxn_ext = vmax*merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
                
                select case(trim(adjustl(sp_name)))
                    case('po2')
                        drxnext_dmsp = (&
                            & vmax*1d0/(po2x+mo2) &
                            & +vmax*po2x*(-1d0)/(po2x+mo2)**2d0 &
                            & )
                    case default
                        drxnext_dmsp = 0d0
                endselect 
            
            case('oxao2')
                ! parameterization adopted by Lawrence et al. 2014; Perez-Fodich and Derry, 2019
                ! decay const of 0.5 yr-1
                rxn_ext = ( &
                    & + poro*sat*1d3*0.5d0 &
                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                    & *maqft_loc(findloc(chraq_all,'oxa',dim=1),:) &
                    & )
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & *dmaqft_dpro_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & )
                            case('ios')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & *dmaqft_dios_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & ) 
                            case('oxa')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & *dmaqft_dmaqf_loc( &
                                    &       findloc(chraq_all,'oxa',dim=1) &
                                    &       ,findloc(chraq_all,'oxa',dim=1),:) &
                                    & + poro*sat*1d3*0.5d0 &
                                    & *1d0 &
                                    & *maqft_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                    & ) 
                            case default
                                if (any( trim(adjustl(sp_name)) == chraq ) ) then
                                    drxnext_dmsp = ( &
                                        & + poro*sat*1d3*0.5d0 &
                                        & * maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                        & * dmaqft_dmaqf_loc( &
                                        &       findloc(chraq_all,'oxa',dim=1)  &
                                        &       ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                        & ) 
                                elseif (any( trim(adjustl(sp_name)) == chrgas ) ) then
                                    drxnext_dmsp = ( &
                                        & + poro*sat*1d3*0.5d0 &
                                        & * maqx_loc(findloc(chraq_all,'oxa',dim=1),:) &
                                        & * dmaqft_dmgas_loc( &
                                        &       findloc(chraq_all,'oxa',dim=1)  &
                                        &       ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:) &
                                        & ) 
                                else
                                    drxnext_dmsp = 0d0
                                endif 
                        endselect 
                
            case('fe2o2')
                ! scheme = 'full' ! reflecting individual rate consts for different Fe2+ species (after Kanzaki and Murakami 2016)
                scheme = 'default' ! as a function of pH and pO2
                
                selectcase(trim(adjustl(scheme)))
                    case('full')
                        rxn_ext = ( &
                            & + poro*sat*1d3*fe2f*( &
                            & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                            & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                            & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                            & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x/prox &
                            & )*po2x &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2*(-1d0)/prox**2d0 &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) &
                                    &       *k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0 &
                                    & )*po2x &
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x/prox &
                                    & )*1d0 &
                                    & )
                            case('pco2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*fe2f*( &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*1d0/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*1d0/prox &
                                    & )*po2x &
                                    & )
                            case('fe2')
                                drxnext_dmsp = ( &
                                    & + poro*sat*1d3*1d0*( &
                                    & + k_arrhenius(10d0**(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) &
                                    & + k_arrhenius(10d0**(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2/prox &
                                    & + k_arrhenius(10d0**(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x/prox**2d0 &
                                    & + k_arrhenius(10d0**(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x/prox &
                                    & )*po2x &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect
                        
                    case default
                        rxn_ext = ( &
                            & poro*sat*1d3*fe2x*po2x &
                            & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                            & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & poro*sat*1d3*fe2x*po2x &
                                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)*2.0d0*(kw*(-1d0)/prox**2d0)) &
                                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & poro*sat*1d3*fe2x*1d0 &
                                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                    & )
                            case('fe2')
                                drxnext_dmsp = ( &
                                    & poro*sat*1d3*po2x &
                                    & * ( & 
                                    & + 1d0    & 
                                    & * maqft_loc(findloc(chraq_all,'fe2',dim=1),:)   &
                                    & + maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    & 
                                    & * dmaqft_dmaqf_loc( &
                                    &       findloc(chraq_all,'fe2',dim=1) &
                                    &       ,findloc(chraq_all,'fe2',dim=1),:)   &
                                    & ) & 
                                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                    & )
                            case default
                                ! drxnext_dmsp = 0d0
                                if (any( trim(adjustl(sp_name)) == chraq ) ) then
                                    drxnext_dmsp = ( &
                                        & poro*sat*1d3*po2x &
                                        & * ( & 
                                        & + maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    & 
                                        & * dmaqft_dmaqf_loc( &
                                        &       findloc(chraq_all,'fe2',dim=1) &
                                        &       ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:)   &
                                        & ) & 
                                        & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                        & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                        & )
                                elseif (any( trim(adjustl(sp_name)) == chrgas ) ) then
                                    drxnext_dmsp = ( &
                                        & poro*sat*1d3*po2x &
                                        & * ( & 
                                        & + maqx_loc(findloc(chraq_all,'fe2',dim=1),:)    & 
                                        & * dmaqft_dmgas_loc( &
                                        &       findloc(chraq_all,'fe2',dim=1) &
                                        &       ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:)   &
                                        & ) & 
                                        & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                                        & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                                        & )
                                else
                                    drxnext_dmsp = 0d0
                                endif 
                        endselect
                endselect 
            
            case('omomb')
                rxn_ext = vmax2 & ! mg C / soil g /yr
                    & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                    & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                    & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
                
                select case(trim(adjustl(sp_name)))
                    case('om')
                        drxnext_dmsp = ( &
                            & vmax2 & ! mg C / soil g /yr
                            & *1d0*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                            & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & + vmax2 & ! mg C / soil g /yr
                            & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                            & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & *(-1d0)/(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3)**2d0 &
                            & *(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & ) 
                    case('omb')
                        drxnext_dmsp = ( &
                            & vmax2 & ! mg C / soil g /yr
                            & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                            & *1d0*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                            & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & ) 
                    case default
                        drxnext_dmsp = 0d0
                endselect
            
            case('ombto')
                rxn_ext = komb*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**beta &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
                
                select case(trim(adjustl(sp_name)))
                    case('omb')
                        drxnext_dmsp = ( &
                            & komb*beta*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**(beta-1d0) &
                            & *(1d0-poro)*rho_grain*rho_grain*1d6 &
                            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                            & ) 
            
                    case default
                        drxnext_dmsp = 0d0
                        
                endselect 
                        
            
            case('pyfe3') 
                rxn_ext = ( &
                    & koxpy*poro*hrpy*mvpy*pyx*fe3x**0.93d0*fe2x**(-0.40d0) &
                    & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                    & /(1d0 - poro) &
                    & )
                
                select case(trim(adjustl(sp_name)))
                    case('py')
                        drxnext_dmsp = (&
                            & koxpy*poro*hrpy*mvpy*1d0*fe3x**0.93d0*fe2x**(-0.40d0) &
                            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                            & )
                    case('fe3')
                        drxnext_dmsp = (&
                            & koxpy*poro*hrpy*mvpy*pyx*(0.93d0)*fe3x**(0.93d0-1d0)*fe2x**(-0.40d0) &
                            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                            & )
                    case('fe2')
                        drxnext_dmsp = (&
                            & koxpy*poro*hrpy*mvpy*pyx*fe3x**0.93d0*(-0.4d0)*fe2x**(-0.40d0-1d0) &
                            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                            & )
                    case default
                        drxnext_dmsp = 0d0
                endselect 
                
            case('amo2o')
                scheme = 'maggi08' ! Maggi et al. (2008) wihtout baterial, pH and water saturation functions
                ! scheme = 'Fennel' ! from biogem_box_geochem.f90 in GENIE model referring to Fennel et al. 2005 with a correction 
                ! scheme = 'FennelOLD' ! from biogem_box_geochem.f90 in GENIE model referring to Fennel et al. 2005 without a correction 
                ! scheme = 'Ozaki' ! from biogem_box_geochem.f90 in GENIE model referring to Ozaki et al. [EPSL ... ?]
                
                select case(trim(adjustl(scheme)))
                    case('maggi08')
                        v_tmp = 9.53d-6*60d0*60d0*24d0*365d0 ! (~300 /yr)
                        ! v_tmp = v_tmp/100d0 ! (~3 /yr; default value produces too much nitrate (pH goes down to ~1)
                        km_tmp1 = 14d-5
                        km_tmp2 = 2.41d-5
                        rxn_ext = ( &
                            & v_tmp &
                            & *nh4x/(nh4x + km_tmp1 ) &
                            & *po2x*ko2/(po2x*ko2 + km_tmp2 ) &
                            & *min(2d0*sat,1d0) &
                            ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                            & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                            & )
                        
                        dv_dph_tmp = 0d0
                        where (3d0 < -log10(prox) .and. -log10(prox) < 7d0)
                            dv_dph_tmp = 0.25d0
                        elsewhere (7d0 < -log10(prox) .and. -log10(prox) < 11d0)
                            dv_dph_tmp = -0.25d0
                        elsewhere (-log10(prox) == 7d0)
                            dv_dph_tmp = 0d0
                        elsewhere (-log10(prox) == 3d0)
                            dv_dph_tmp = 0.125d0
                        elsewhere (-log10(prox) == 11d0)
                            dv_dph_tmp = -0.125d0
                        elsewhere 
                            dv_dph_tmp = 0d0
                        endwhere 
                        
                        ! when using modified version using normal distribution with sigma = 1 
                        dv_dph_tmp = exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &
                            & *-0.5d0*2d0*((-log10(prox)-7d0)/1d0)  &
                            & *(-1d0) &
                            & *1d0/log(10d0)/prox
                        
                        select case(trim(adjustl(sp_name)))
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *nh4x/(nh4x + km_tmp1 ) &
                                    & *( &
                                    & 1d0*ko2/(po2x*ko2 + km_tmp2) &
                                    & + po2x*ko2*(-1d0)/(po2x*ko2 + km_tmp2)**2d0 * ko2 &
                                    & ) &
                                    & *min(2d0*sat,1d0) &
                                    ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pnh3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & dnh4_dpnh3/(nh4x + km_tmp1 ) &
                                    & + nh4x*(-1d0)/(nh4x + km_tmp1 )**2d0 * dnh4_dpnh3 &
                                    & ) &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & *min(2d0*sat,1d0) &
                                    ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & dnh4_dpro/(nh4x + km_tmp1 ) &
                                    & + nh4x*(-1d0)/(nh4x + km_tmp1 )**2d0 * dnh4_dpro &
                                    & ) &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & *min(2d0*sat,1d0) &
                                    ! & *max( min( 0.25d0*(-log10(prox))-0.75d0, -0.25d0*(-log10(prox))+2.75d0 ), 0d0 ) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & + &
                                    & v_tmp &
                                    & *nh4x/(nh4x + km_tmp1 ) &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *dv_dph_tmp &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                    case('Fennel','FennelOLD')
                        if (trim(adjustl(scheme))== 'Fennel') v_tmp = 6.0d0 ! /yr
                        if (trim(adjustl(scheme))== 'FennelOLD') v_tmp = 0.16667d0 ! /yr
                        km_tmp2 = 2.0D-05
                        rxn_ext = ( &
                            & v_tmp &
                            & *nh4x &
                            & *po2x*ko2/(po2x*ko2 + km_tmp2 ) &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *nh4x &
                                    & *( &
                                    & 1d0*ko2/(po2x*ko2 + km_tmp2) &
                                    & + po2x*ko2*(-1d0)/(po2x*ko2 + km_tmp2)**2d0 * ko2 &
                                    & ) &
                                    & )
                            case('pnh3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *dnh4_dpnh3 &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * dnh4_dpro &
                                    & *po2x*ko2/(po2x*ko2 + km_tmp2) &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                    case('Ozaki')
                        v_tmp = 18250.0d0/1027.649d0 ! /yr
                        rxn_ext = ( &
                            & v_tmp &
                            & *nh4x &
                            & *po2x*ko2 &
                            & )
                        
                        select case(trim(adjustl(sp_name)))
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *nh4x &
                                    & * 1d0*ko2 &
                                    & )
                            case('pnh3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *dnh4_dpnh3 &
                                    & *po2x*ko2 &
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * dnh4_dpro &
                                    & *po2x*ko2 &
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                    endselect
                        
            case('g2n0','g2n21') 
                ! overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
                ! first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O)   
                ! (assuming that oxidation by N2O governs overall denitrification)
                scheme = 'maggi08' ! Maggi et al. (2008) wihtout baterial, pH and water saturation functions; vmax from oxidation by N2O (rate-limiting)
                
                select case(trim(adjustl(scheme)))
                    case('maggi08')
                        v_tmp = 1.23d-7*60d0*60d0*24d0*365d0
                        km_tmp1 = 10d-5 * 1d6 ! mol L-1 converted to mol m-3
                        km_tmp2 = 11.3d-5
                        km_tmp3 = 2.52d-5
                        rxn_ext = ( &
                            & v_tmp &
                            & *g2x/(g2x + km_tmp1 ) &
                            & *no3x/(no3x + km_tmp2 ) &
                            & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                            & *min(2d0*sat,1d0) &
                            & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                            & )
                        
                        ! when using modified version using normal distribution with sigma = 1 
                        dv_dph_tmp = exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &
                            & *-0.5d0*2d0*((-log10(prox)-7d0)/1d0)  &
                            & *(-1d0) &
                            & *1d0/log(10d0)/prox
                            
                        select case(trim(adjustl(sp_name)))
                            case('g2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & 1d0/(g2x + km_tmp1 ) &
                                    & + g2x*(-1d0)/(g2x + km_tmp1 )**2d0 * 1d0 &
                                    & ) &
                                    & *no3x/(no3x + km_tmp2 ) &
                                    & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('no3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & * ( & 
                                    & 1d0/(no3x + km_tmp2 ) &
                                    & + no3x*(-1d0)/(no3x + km_tmp2 )**2d0 * 1d0 &
                                    & ) &
                                    & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *no3x/(no3x + km_tmp2 ) &
                                    & *km_tmp3*(-1d0)/(po2x*ko2 + km_tmp3 )**2d0 * ko2 &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *no3x/(no3x + km_tmp2 ) &
                                    & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *dv_dph_tmp &! modified version using normal distribution with sigma = 1 
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                endselect 
                        
            case('g2n22') ! 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O)  
                scheme = 'maggi08' ! Maggi et al. (2008) wihtout baterial, pH and water saturation functions
                
                select case(trim(adjustl(scheme)))
                    case('maggi08')
                        v_tmp = 1.23d-7*60d0*60d0*24d0*365d0
                        km_tmp1 = 10d-5 * 1d6 ! mol L-1 converted to mol m-3
                        km_tmp2 = 11.3d-5
                        km_tmp3 = 2.52d-5
                        rxn_ext = ( &
                            & v_tmp &
                            & *g2x/(g2x + km_tmp1 ) &
                            & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                            ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                            & *km_tmp3/(no3x + km_tmp3 ) &
                            & *min(2d0*sat,1d0) &
                            & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                            & )
                        
                        ! when using modified version using normal distribution with sigma = 1 
                        dv_dph_tmp = exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &
                            & *-0.5d0*2d0*((-log10(prox)-7d0)/1d0)  &
                            & *(-1d0) &
                            & *1d0/log(10d0)/prox
                            
                        select case(trim(adjustl(sp_name)))
                            case('g2')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & * ( & 
                                    & 1d0/(g2x + km_tmp1 ) &
                                    & + g2x*(-1d0)/(g2x + km_tmp1 )**2d0 * 1d0 &
                                    & ) &
                                    & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *km_tmp3/(no3x + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pn2o')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & * ( & 
                                    & kn2o/(kn2o*pn2ox + km_tmp2 ) &
                                    & + kn2o*pn2ox*(-1d0)/(kn2o*pn2ox + km_tmp2 )**2d0 * kn2o &
                                    & ) &
                                    ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *km_tmp3/(no3x + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            ! case('po2')
                                ! drxnext_dmsp = ( &
                                    ! & v_tmp &
                                    ! & *g2x/(g2x + km_tmp1 ) &
                                    ! & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3*(-1d0)/(po2x*ko2 + km_tmp3 )**2d0 * ko2 &
                                    ! & *min(2d0*sat,1d0) &
                                    ! & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    ! & )
                            case('no3')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3*(-1d0)/(po2x*ko2 + km_tmp3 )**2d0 * ko2 &
                                    & *km_tmp3*(-1d0)/(no3x + km_tmp3 )**2d0 * 1d0 &
                                    & *min(2d0*sat,1d0) &
                                    & *exp(-0.5d0*((-log10(prox)-7d0)/1d0)**2d0) &! modified version using normal distribution with sigma = 1 
                                    & )
                            case('pro')
                                drxnext_dmsp = ( &
                                    & v_tmp &
                                    & *g2x/(g2x + km_tmp1 ) &
                                    & *kn2o*pn2ox/(kn2o*pn2ox + km_tmp2 ) &
                                    ! & *km_tmp3/(po2x*ko2 + km_tmp3 ) &
                                    & *km_tmp3/(no3x + km_tmp3 ) &
                                    & *min(2d0*sat,1d0) &
                                    & *dv_dph_tmp &! modified version using normal distribution with sigma = 1 
                                    & )
                            case default
                                drxnext_dmsp = 0d0
                        endselect 
                        
                endselect 
            
            case('g2k','g2ca','g2mg')
                ! cation uptake by g2 
                ! 0.0395 Eq mol-1 C harvested from Kantzas et al. 2022 
                ! for now OM is parameterized as NPP x 1.5 
                ! harvest = NPP * HI * RS /DF / C (Monfreda et al. 2008)
                ! where 
                !   HI = harvest index (0.85 as assumed by Kantzas et al. 2022)
                !   RS = root:shoot ratio (0.25 assumed) 
                !   DF = dry proportion of economic yield (0.9 assumed, cf. Kroodsma and Field, 2006)  
                !   C  = carbon content (0.45 g C / g dry matter) 
                !   
                ! So, vmax = 0.0395 * kresp * g2 / 1.5 * 0.85 * 0.25 / 0.9 / 0.45 
                ! And assume ion uptake kinetics by Fageria et al. 2010 for wheat (data from Barber 1995)
                ! Michaelis-Menten formulation 
                !   maq/(Kmaq + maq) 
                !       where 
                !           maq is ion conc. (M)
                !           Kmaq is Michaeris constants (M)
                ! for now, only consider K, Ca and Mg uptake
                ! Maximum uptake is distributed between these ions by K:Ca:Mg = 7:1.6:0.4
                ! and Kmaq for K, Ca and Mg are 7e-6, 5e-6 and 1e-6 M, respectively.
                
                select case(trim(adjustl(rxn_name)))
                    case('g2k')
                        sp_tmp = 'k    '
                        km_tmp1 = 7d0/11d0
                        km_tmp2 = 7d-6
                    case('g2ca')
                        sp_tmp = 'ca   '
                        km_tmp1 = 1.6d0/11d0
                        km_tmp2 = 5d-6
                    case('g2mg')
                        sp_tmp = 'mg   '
                        km_tmp1 = 0.4d0/11d0
                        km_tmp2 = 1d-6
                    case default
                        print *, '*** FUNDAMENTAL ERROR in extra reaction: cation uptake'
                        stop
                endselect  
                
                vmax_tmp = 0.0395d0 * kg2  / 1.5d0 * 0.85d0 * 0.25d0 / 0.9d0 / 0.45d0
                
                rxn_ext = ( &
                    & + vmax_tmp * g2x  &
                    & *po2x/(mo2 + po2x)  &
                    & *km_tmp1  &
                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & /( km_tmp2 &
                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                    & ) &
                    & )
                    
                        select case(trim(adjustl(sp_name)))
                            case('pro')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dpro_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & +  ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *(-1d0)/( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & )**2d0 &
                                    & ) &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dpro_loc(findloc(chraq_all,sp_tmp,dim=1),:) 
                            case('ios')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dios_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & +  ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *(-1d0)/( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & )**2d0 &
                                    & ) &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dios_loc(findloc(chraq_all,sp_tmp,dim=1),:) 
                            case('g2')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * 1d0   &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & )
                            case('po2')
                                drxnext_dmsp = ( &
                                    & + vmax_tmp * g2x  &
                                    & *1d0/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & + ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x*(-1d0)/(mo2 + po2x)**2d0  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) &
                                    & + ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dmgas_loc( &
                                        & findloc(chraq_all,sp_tmp,dim=1) &
                                        & ,findloc(chrgas_all,'po2',dim=1),:) &
                                    & /( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & ) &
                                    & ) & 
                                    & + ( &
                                    & + vmax_tmp * g2x  &
                                    & *po2x/(mo2 + po2x)  &
                                    & *km_tmp1  &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *(-1d0)/( km_tmp2 &
                                    & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & )**2d0 &
                                    & ) &
                                    & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                    & *dmaqft_dmgas_loc( &
                                    & findloc(chraq_all,sp_tmp,dim=1) &
                                    & ,findloc(chrgas_all,'po2',dim=1),:) 
                            case default
                                if (any( trim(adjustl(sp_name)) == chraq ) ) then
                                    if ( sp_name == sp_tmp ) then 
                                        drxnext_dmsp = ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                                & findloc(chraq_all,sp_tmp,dim=1) &
                                                & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                            & /( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & ) &
                                            & ) & 
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *(-1d0)/( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )**2d0 &
                                            & ) &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                            & findloc(chraq_all,sp_tmp,dim=1) &
                                            & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *1d0 &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & /( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & ) &
                                            & ) &
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *(-1d0)/( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )**2d0 &
                                            & *1d0 &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )
                                    else
                                        drxnext_dmsp = ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                                & findloc(chraq_all,sp_tmp,dim=1) &
                                                & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:) &
                                            & /( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & ) &
                                            & ) & 
                                            & + ( &
                                            & + vmax_tmp * g2x  &
                                            & *po2x/(mo2 + po2x)  &
                                            & *km_tmp1  &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *(-1d0)/( km_tmp2 &
                                            & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & )**2d0 &
                                            & ) &
                                            & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                            & *dmaqft_dmaqf_loc( &
                                            & findloc(chraq_all,sp_tmp,dim=1) &
                                            & ,findloc(chraq_all,trim(adjustl(sp_name)),dim=1),:)
                                    endif 
                                elseif (any( trim(adjustl(sp_name)) == chrgas ) ) then
                                    drxnext_dmsp = ( &
                                        & + vmax_tmp * g2x  &
                                        & *po2x/(mo2 + po2x)  &
                                        & *km_tmp1  &
                                        & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *dmaqft_dmgas_loc( &
                                            & findloc(chraq_all,sp_tmp,dim=1) &
                                            & ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:) &
                                        & /( km_tmp2 &
                                        & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & ) &
                                        & ) & 
                                        & + ( &
                                        & + vmax_tmp * g2x  &
                                        & *po2x/(mo2 + po2x)  &
                                        & *km_tmp1  &
                                        & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *(-1d0)/( km_tmp2 &
                                        & +maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *maqft_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & )**2d0 &
                                        & ) &
                                        & *maqx_loc(findloc(chraq_all,sp_tmp,dim=1),:) &
                                        & *dmaqft_dmgas_loc( &
                                        & findloc(chraq_all,sp_tmp,dim=1) &
                                        & ,findloc(chrgas_all,trim(adjustl(sp_name)),dim=1),:) 
                                else
                                    drxnext_dmsp = 0d0
                                endif 
                        endselect 
                
                
            case default 
                rxn_ext = 0d0
                drxnext_dmsp = 0d0
                
        endselect

        rxnext_error = .false.
        if (any(isnan(rxn_ext)) .or. any(isnan(drxnext_dmsp))) then 
            print *,'nan in calc_rxn_ext_dev_3'
            print *,'any(isnan(rxn_ext)) | any(isnan(drxnext_dmsp))'
            print *,rxn_name,sp_name,any(isnan(rxn_ext)),any(isnan(drxnext_dmsp)) &
                & ,any(isnan(dmaqft_dios_loc(findloc(chraq_all,sp_tmp,dim=1),:) ))
            rxnext_error = .true.
            ! stop
        endif 

    endsubroutine calc_rxn_ext_dev_3

    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
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
    !ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


end module scepter_kinetics