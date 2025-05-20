!************************************************************************
! Module for solid dissolution/precipitation kinetics
!************************************************************************   
module scepter_sld_kin
    use scepter_constants
    use scepter_thermodynamics
    implicit none
    private
    public :: sld_kin

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
        
endmodule scepter_sld_kin