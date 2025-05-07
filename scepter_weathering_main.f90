module scepter_weathering_main
use scepter_constants
use scepter_variables
implicit none

subroutine weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup0,zsat,zml_ref,w0,q0,p80,ttot,plant_rain,rainpowder_2nd  &! input
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save &! input
    & ,count_dtunchanged_Max,tcin,step_tau &! input 
    & ,nsld_kinspc_in,chrsld_kinspc_in,kin_sld_spc_in &! input 
    & )

    implicit none

    !-----------------------------
    #ifdef mod_basalt_cmp
    #include <basalt_defines.h>
    #endif                                    
    !-------------------------

    tc = tcin
    qin = q0
    satup = satup0

    nsp_sld_cnst = nsp_sld_all - nsp_sld
    nsp_aq_cnst = nsp_aq_all - nsp_aq
    nsp_gas_cnst = nsp_gas_all - nsp_gas
    nsp3 = nsp_sld + nsp_aq + nsp_gas

    #ifdef calcw_full
    nsp3 = nsp3 + 1
    #endif 

    isldprof    = idust + nsp_sld + nsp_gas + nsp_aq + 1
    isldprof2   = idust + nsp_sld + nsp_gas + nsp_aq + 2
    isldprof3   = idust + nsp_sld + nsp_gas + nsp_aq + 3
    iaqprof     = idust + nsp_sld + nsp_gas + nsp_aq + 4
    iaqprof2    = idust + nsp_sld + nsp_gas + nsp_aq + 5
    iaqprof3    = idust + nsp_sld + nsp_gas + nsp_aq + 6
    iaqprof4    = idust + nsp_sld + nsp_gas + nsp_aq + 7
    iaqprof5    = idust + nsp_sld + nsp_gas + nsp_aq + 8
    iaqprof6    = idust + nsp_sld + nsp_gas + nsp_aq + 9
    igasprof    = idust + nsp_sld + nsp_gas + nsp_aq + 10
    isldsat     = idust + nsp_sld + nsp_gas + nsp_aq + 11
    ibsd        = idust + nsp_sld + nsp_gas + nsp_aq + 12
    irate       = idust + nsp_sld + nsp_gas + nsp_aq + 13
    ipsd        = idust + nsp_sld + nsp_gas + nsp_aq + 14
    ipsdv       = idust + nsp_sld + nsp_gas + nsp_aq + 15
    ipsds       = idust + nsp_sld + nsp_gas + nsp_aq + 16
    ipsdflx     = idust + nsp_sld + nsp_gas + nsp_aq + 17
    isa         = idust + nsp_sld + nsp_gas + nsp_aq + 18
    isa2        = idust + nsp_sld + nsp_gas + nsp_aq + 19

    ! species whose flux is saved all time
    ! chrsp_saveall = (/'pco2 '/)
    chrsp_saveall = (/'xxxxx'/)

    nflx = 5 + nrxn_ext + nsp_sld

    do isps=1,nsp_sld
        irxn_sld(isps) = 4+isps
    enddo 

    do irxn=1,nrxn_ext
        irxn_ext(irxn) = 4+nsp_sld+irxn
    enddo 

    ires = nflx

    chrflx(1:4) = (/'tflx ','adv  ','dif  ','rain '/)
    if (nsp_sld > 0) chrflx(irxn_sld(:)) = chrsld
    if (nrxn_ext > 0) chrflx(irxn_ext(:)) = chrrxn_ext
    chrflx(nflx) = 'res  '

    ! print *,chrflx,irxn_sld,chrsld

    ! pause

    ! define all species and rxns definable in the model 
    ! note that rxns here exclude diss(/prec) of mineral 
    ! which are automatically included when associated mineral is chosen

    chrsld_all = (/'fo   ','ab   ','an   ','cc   ','ka   ','gb   ','py   ','ct   ','fa   ','gt   ','cabd ' &
        & ,'dp   ','hb   ','kfs  ','om   ','omb  ','amsi ','arg  ','dlm  ','hm   ','ill  ','anl  ','nph  ' &
        & ,'qtz  ','gps  ','tm   ','la   ','by   ','olg  ','and  ','cpx  ','en   ','fer  ','opx  ','kbd  ' &
        & ,'mgbd ','nabd ','mscv ','plgp ','antp ','agt  ','jd   ','wls  ','phsi ','splt ','casp ','ksp  ' &
        & ,'nasp ','mgsp ','fe2o ','mgo  ','k2o  ','cao  ','na2o ','al2o3','gbas ','cbas ','ep   ','clch ' &
        & ,'sdn  ','cdr  ','leu  ','amal ','amfe3' &
        & ,'g1   ','g2   ','g3   ','amnt ','kcl  ','gac  ','mesmh','ims  ','teas ','naoh ','naglp','cacl2' &
        & ,'nacl ','sio2 ','caso4' &
        & ,'inrt '/)
    chraq_all  = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    ','no3  ','oxa  ' &
        & ,'cl   ','ac   ','mes  ','im   ','tea  ','glp  '/)
    chrgas_all = (/'pco2 ','po2  ','pnh3 ','pn2o '/)
    chrrxn_ext_all = (/'resp ','fe2o2','omomb','ombto','pyfe3','amo2o','g2n0 ','g2n21','g2n22','oxao2' &
        & ,'g2k  ','g2ca ','g2mg '/)

    ! define the species and rxns explicitly simulated in the model in a fully coupled way
    ! should be chosen from definable species & rxn lists above 

    ! chrsld = (/'fo   ','ab   ','an   ','cc   ','ka   '/)
    ! chraq = (/'mg   ','si   ','na   ','ca   ','al   '/)
    ! chrgas = (/'pco2 ','po2  '/)
    ! chrrxn_ext = (/'resp '/)


    ! define solid species which can precipitate
    ! in default, all minerals only dissolve 
    ! should be chosen from the chrsld list
    ! #ifdef diss_only
    ! chrsld_2(:) = '     '
    ! #else
    ! chrsld_2 = (/'cc   ','ka   ','gb   ','ct   ','gt   ','cabd ','amsi ','hm   ','ill  ','anl  ','gps  '  &
        ! & ,'arg  ','dlm  ','qtz  ','mgbd ','nabd ','kbd  ','phsi ','casp ','ksp  ','nasp ','mgsp ','al2o3'  &
        ! & ,'amal ','amfe3' /) 
        
    call get_2ndsld_num(nsp_sld_2)

    if (allocated(chrsld_2)) deallocate(chrsld_2)
    allocate(chrsld_2(nsp_sld_2))

    call get_2ndsld( &
        & nsp_sld_2 &! input
        & ,chrsld_2 &! output
        & )
        
    ! version that removes dolomite from 2ndary minerals
    ! chrsld_2 = (/'cc   ','ka   ','gb   ','ct   ','gt   ','cabd ','amsi ','hm   ','ill  ','anl  ','gps  '  &
        ! & ,'arg  ','qtz  ','mgbd ','nabd ','kbd  ','phsi ','casp ','ksp  ','nasp ','mgsp ','al2o3','amal '  &
        ! & ,'amfe3' /) 
        
    ! version that removes all carbonates from 2ndary minerals
    ! chrsld_2 = (/'ka   ','gb   ','ct   ','gt   ','cabd ','amsi ','hm   ','ill  ','anl  ','gps  '  &
        ! ,'qtz  ','mgbd ','nabd ','kbd  ','phsi ','casp ','ksp  ','nasp ','mgsp ','al2o3'/) 
        
    ! version that removes all base-cation bearer from 2ndary minerals
    ! chrsld_2 = (/'ka   ','gb   ','ct   ','gt   ','amsi ','hm   ','gps  '  &
        ! & ,'qtz  ','al2o3','amal '  &
        ! & ,'amfe3' /) 
    ! #endif 
    ! below are species which are sensitive to pH 
    chraq_ph   = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    ','no3  ','oxa  ' &
        & ,'cl   ','ac   ','mes  ','im   ','tea  ','glp  '/)
    chrgas_ph = (/'pco2 ','pnh3 '/)

    chrco2sp = (/'co2g ','co2aq','hco3 ','co3  ','DIC  ','ALK  '/)

    if (nsp_gas_cnst .ne. 0) then 
        do ispg = 1, nsp_gas_cnst
            do ispg2=1,nsp_gas_all
                if (.not.any(chrgas==chrgas_all(ispg2)) .and. .not.any(chrgas_cnst==chrgas_all(ispg2))) then 
                    chrgas_cnst(ispg) = chrgas_all(ispg2)
                    exit 
                endif 
            enddo
        enddo 
        print *, chrgas_cnst
        ! pause 
    endif 

    if (nsp_sld_cnst .ne. 0) then 
        do isps = 1, nsp_sld_cnst
            do isps2=1,nsp_sld_all
                if (.not.any(chrsld==chrsld_all(isps2)) .and. .not.any(chrsld_cnst==chrsld_all(isps2))) then 
                    chrsld_cnst(isps) = chrsld_all(isps2)
                    exit 
                endif 
            enddo
        enddo 
        print *, chrsld_cnst
        ! pause 
    endif 

    ! molar volume 

    mv_all = (/mvfo,mvab,mvan,mvcc,mvka,mvgb,mvpy,mvct,mvfa,mvgt,mvcabd,mvdp,mvhb,mvkfs,mvom,mvomb,mvamsi &
        & ,mvarg,mvdlm,mvhm,mvill,mvanl,mvnph,mvqtz,mvgps,mvtm,mvla,mvby,mvolg,mvand,mvcpx,mven,mvfer,mvopx &
        & ,mvkbd,mvmgbd,mvnabd,mvmscv,mvplgp,mvantp,mvagt,mvjd,mvwls,mvphsi,mvsplt,mvcasp,mvksp,mvnasp,mvmgsp &
        & ,mvfe2o,mvmgo,mvk2o,mvcao,mvna2o,mval2o3,mvgbas,mvcbas,mvep,mvclch,mvsdn,mvcdr,mvleu,mvamal,mvamfe3 &
        & ,mvg1,mvg2,mvg3,mvamnt,mvkcl,mvgac,mvmesmh,mvims,mvteas,mvnaoh,mvnaglp,mvcacl2,mvnacl,mvsio2,mvcaso4  &
        & ,mvinrt/)
    mwt_all = (/mwtfo,mwtab,mwtan,mwtcc,mwtka,mwtgb,mwtpy,mwtct,mwtfa,mwtgt,mwtcabd,mwtdp,mwthb,mwtkfs,mwtom,mwtomb,mwtamsi &
        & ,mwtarg,mwtdlm,mwthm,mwtill,mwtanl,mwtnph,mwtqtz,mwtgps,mwttm,mwtla,mwtby,mwtolg,mwtand,mwtcpx,mwten,mwtfer,mwtopx &
        & ,mwtkbd,mwtmgbd,mwtnabd,mwtmscv,mwtplgp,mwtantp,mwtagt,mwtjd,mwtwls,mwtphsi,mwtsplt,mwtcasp,mwtksp,mwtnasp,mwtmgsp &
        & ,mwtfe2o,mwtmgo,mwtk2o,mwtcao,mwtna2o,mwtal2o3,mwtgbas,mwtcbas,mwtep,mwtclch,mwtsdn,mwtcdr,mwtleu,mwtamal,mvamfe3 &
        & ,mwtg1,mwtg2,mwtg3,mwtamnt,mwtkcl,mwtgac,mwtmesmh,mwtims,mwtteas,mwtnaoh,mwtnaglp,mwtcacl2,mwtnacl,mwtsio2,mwtcaso4 &
        & ,mwtinrt/)

    do isps = 1, nsp_sld 
        mv(isps) = mv_all(findloc(chrsld_all,chrsld(isps),dim=1))
        mwt(isps) = mwt_all(findloc(chrsld_all,chrsld(isps),dim=1))
    enddo 

    mwtaq_all = (/ mwtaqmg,mwtaqsi,mwtaqna,mwtaqca,mwtaqal,mwtaqfe2,mwtaqfe3,mwtaqso4,mwtaqk,mwtaqno3,mwtaqoxa  &
        & ,mwtaqcl,mwtaqac,mwtaqmes,mwtaqim,mwtaqtea,mwtaqglp /)
        
    do ispa = 1, nsp_aq 
        mwtaq(ispa) = mwtaq_all(findloc(chraq_all,chraq(ispa),dim=1))
    enddo 


    ! maqi_all = 0d0
        
    def_rain = 1d-20
    ! def_rain = 1d-50
    def_pr = 1d-20
    ! def_pr = 1d-50
    ! def_pr = 1d-0
        
    call get_rainwater( &
        & nsp_aq_all,chraq_all,def_rain &! input
        & ,maqi_all &! output
        & )
        
    call get_parentrock( &
        & nsp_sld_all,chrsld_all,def_pr &! input
        & ,msldi_all &! output
        & )

    ! bulk soil concentration     
    mblki = 0d0
    incld_blk = .false.

    ! adding the case where input wt% exceeds 100% 
    if ( sum(msldi_all) > 1d0) then 
        print *, 'parent rock comp. exceeds 100% so rescale'
        msldi_all = msldi_all/sum(msldi_all)  ! now the units are g/g
    ! endif 
    ! msldi_all = msldi_all/mwt_all*rho_grain*1d6 ! converting g/g to mol/sld m3
    ! msldi_all = (1d0 - poroi) * msldi_all       ! mol/sld m3 to mol/bulk m3 

    ! when input is less than 100wt% add bulk soil 
    elseif ( sum(msldi_all) < 1d0 ) then 
        print *, 'parent rock comp. is less than 100% so add "bulk soil"'
        ! sum(msldi_all)  + mblki = 1d0
        incld_blk = .true.
        mblki = 1d0 - sum(msldi_all)
        if ( mblki < 0d0 ) mblki = 0d0
    endif 

    rho_grain_calc = rho_grain
    msldi_allx = msldi_all/mwt_all*rho_grain_calc*1d6 !  converting g/g to mol/sld m3
    mblkix = mblki/mwtblk*rho_grain_calc*1d6 !  converting g/g to mol/sld m3
    ! msldi_allx = msldi_allx/sum(msldi_allx*mv_all*1d-6)  !  try to make sure volume total must be 1 
    msldi_allx = msldi_allx/( sum(msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6 )  !  try to make sure volume total must be 1 (including bulk soil if any)
    mblkix = mblkix/( sum(msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6 )  !  try to make sure volume total must be 1 (including bulk soil if any)
    if (msldunit=='blk') then 
        msldi_allx = msldi_allx*(1d0 - poroi)  !  try to make sure volume total must be 1 - poroi
        mblkix = mblkix*(1d0 - poroi)  !  try to make sure volume total must be 1 - poroi
    endif 
    ! then the follwoing must be satisfied
    ! (1d0 - poroi)*rho_grain_calc = msldi_all*mwt_all*1d-6
    ! (1d0 - poroi) = msldi_all*mv_all*1d-6
    rho_error = 1d4
    rho_tol = 1d-6
    ! poroi_calc = 1d0 - sum(msldi_allx*mv_all*1d-6)
    poroi_calc = 1d0 - ( sum( msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6 ) ! corrected for bulk soil if any 
    do while (rho_error > rho_tol) 
        rho_grain_calcx = rho_grain_calc
        
        ! rho_grain_calc = sum(msldi_allx(:)*mwt_all(:)*1d-6) ! /(1d0-poroi_calc)
        rho_grain_calc = sum(msldi_allx(:)*mwt_all(:)*1d-6)  + mblkix*mwtblk*1d-6  ! /(1d0-poroi_calc) | corrected for bulk soil 
        if (msldunit=='blk') rho_grain_calc = rho_grain_calc / (1d0-poroi_calc)
        
        msldi_allx = msldi_all/mwt_all*rho_grain_calc*1d6 !  converting g/g to mol/sld m3
        mblkix = mblki/mwtblk*rho_grain_calc*1d6 !  converting g/g to mol/sld m3
        ! msldi_allx = msldi_allx/sum(msldi_allx*mv_all*1d-6)  !  try to make sure volume total must be 1 
        msldi_allx = msldi_allx/( sum(msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6 )  !  try to make sure volume total must be 1 | corrected for bulk soil 
        if (msldunit=='blk')  then 
            msldi_allx = (1d0 - poroi) * msldi_allx  !  converting mol/sld m3 to mol/bulk m3
            mblkix = (1d0 - poroi) * mblkix  !  converting mol/sld m3 to mol/bulk m3
        endif 
        
        ! poroi_calc = 1d0 - sum(msldi_allx*mv_all*1d-6)
        poroi_calc = 1d0 - ( sum(msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6 ) ! corrected for bulk soil
        
        rho_error = abs ((rho_grain_calc - rho_grain_calcx)/rho_grain_calc) 
        
        print*,rho_error
        
    enddo

    if (msldunit=='sld') then 
        ! print *,1d0,sum(msldi_allx*mv_all*1d-6)
        ! print *,sum(msldi_allx*mwt_all*1d-6),rho_grain_calc
        print *,1d0,sum(msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6
        print *,sum(msldi_allx*mwt_all*1d-6) + mblkix*mwtblk*1d-6 ,rho_grain_calc
    elseif (msldunit=='blk') then 
        ! print *,1d0,sum(msldi_allx*mv_all*1d-6) /(1d0 - poroi)
        ! print *,sum(msldi_allx*mwt_all*1d-6)/(1d0 - poroi),rho_grain_calc
        print *,1d0,( sum(msldi_allx*mv_all*1d-6) + mblkix*mvblk*1d-6 )/(1d0 - poroi)
        print *,( sum(msldi_allx*mwt_all*1d-6) + mblkix*mwtblk*1d-6 )/(1d0 - poroi),rho_grain_calc
    endif 
    ! print *,mblkix
    ! stop
    msldi_all = msldi_allx
    mblki = mblkix
    rho_grain = rho_grain_calc

    call get_atm( &
        & nsp_gas_all,chrgas_all &! input
        & ,mgasi_all &! output
        & )

    ! print*,maqi_all 
    ! print*,mgasi_all 
    print*,msldi_all

    ! pause

    ! constant values are taken from the boundary values specified above 
    do ispg = 1,nsp_gas_cnst
        mgasc(ispg,:) = mgasi_all(findloc(chrgas_all,chrgas_cnst(ispg),dim=1))
        ! mgasc(ispg,:) = 0d0
    enddo 
    do ispa = 1,nsp_aq_cnst
        maqc(ispa,:) = maqi_all(findloc(chraq_all,chraq_cnst(ispa),dim=1))
        ! maqc(ispa,:) = 0d0
    enddo 
    do isps = 1,nsp_sld_cnst
        msldc(isps,:) = msldi_all(findloc(chrsld_all,chrsld_cnst(isps),dim=1))
        ! msldc(isps,:) = 0d0
    enddo 

    ! threshould values 
    mgasth_all = 1d-200
    maqth_all = 1d-200
    msldth_all = 1d-200


    ! passing initial and threshold values to explcit variables 
    do isps = 1, nsp_sld    
        print *, chrsld(isps)
        if (any(chrsld_all == chrsld(isps))) then 
            msldi(isps) = msldi_all(findloc(chrsld_all,chrsld(isps),dim=1))
            msldth(isps) = msldth_all(findloc(chrsld_all,chrsld(isps),dim=1))
            print *,msldi(isps),msldi_all(findloc(chrsld_all,chrsld(isps),dim=1))
        endif 
    enddo 
    do ispa = 1, nsp_aq    
        if (any(chraq_all == chraq(ispa))) then 
            maqi(ispa) = maqi_all(findloc(chraq_all,chraq(ispa),dim=1))
            maqth(ispa) = maqth_all(findloc(chraq_all,chraq(ispa),dim=1))
        endif 
    enddo 
    do ispg = 1, nsp_gas    
        if (any(chrgas_all == chrgas(ispg))) then 
            mgasi(ispg) = mgasi_all(findloc(chrgas_all,chrgas(ispg),dim=1))
            mgasth(ispg) = mgasth_all(findloc(chrgas_all,chrgas(ispg),dim=1))
        endif 
    enddo 

    ! print*,maqi 
    ! print*,mgasi
    print*,msldi

    ! pause

    ! stoichiometry
    ! mineral dissolution(/precipitation)
    staq_all = 0d0
    stgas_all = 0d0
    ! Forsterite; Mg2SiO4
    staq_all(findloc(chrsld_all,'fo',dim=1), findloc(chraq_all,'mg',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'fo',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Albite; NaAlSi3O8
    ! staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
    ! staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    ! staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! Analcime; NaAlSi2O6*H2O
    staq_all(findloc(chrsld_all,'anl',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'anl',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'anl',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! Nepheline; NaAlSiO4
    staq_all(findloc(chrsld_all,'nph',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'nph',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'nph',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! K-feldspar; KAlSi3O8
    staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! Sanidine; KAlSi3O8
    staq_all(findloc(chrsld_all,'sdn',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'sdn',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'sdn',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! Anothite; CaAl2Si2O8
    ! staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    ! staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    ! staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
    ! Albite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_ab
    staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_ab
    staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_ab
    staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_ab
    ! Anothite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_an
    staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_an
    staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_an
    staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_an
    ! Labradorite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_la
    staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_la
    staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_la
    staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_la
    ! Andesine; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(findloc(chrsld_all,'and',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_and
    staq_all(findloc(chrsld_all,'and',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_and
    staq_all(findloc(chrsld_all,'and',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_and
    staq_all(findloc(chrsld_all,'and',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_and
    ! Oligoclase; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(findloc(chrsld_all,'olg',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_olg
    staq_all(findloc(chrsld_all,'olg',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_olg
    staq_all(findloc(chrsld_all,'olg',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_olg
    staq_all(findloc(chrsld_all,'olg',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_olg
    ! Bytownite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(findloc(chrsld_all,'by',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_by
    staq_all(findloc(chrsld_all,'by',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_by
    staq_all(findloc(chrsld_all,'by',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_by
    staq_all(findloc(chrsld_all,'by',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_by
    ! Calcite; CaCO3
    staq_all(findloc(chrsld_all,'cc',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    stgas_all(findloc(chrsld_all,'cc',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    ! Kaolinite; Al2Si2O5(OH)4
    staq_all(findloc(chrsld_all,'ka',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'ka',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
    ! Amorphous Al(OH)3
    staq_all(findloc(chrsld_all,'amal',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! Gibbsite; Al(OH)3
    staq_all(findloc(chrsld_all,'gb',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    ! Pyrite; FeS2
    staq_all(findloc(chrsld_all,'py',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'py',dim=1), findloc(chraq_all,'so4',dim=1)) = 2d0
    stgas_all(findloc(chrsld_all,'py',dim=1), findloc(chrgas_all,'po2',dim=1)) = -7d0/2d0
    ! Chrysotile; Mg3Si2O5(OH)4
    staq_all(findloc(chrsld_all,'ct',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'ct',dim=1), findloc(chraq_all,'mg',dim=1)) = 3d0
    ! Fayalite; Fe2SiO4
    staq_all(findloc(chrsld_all,'fa',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'fa',dim=1), findloc(chraq_all,'fe2',dim=1)) = 2d0
    ! Amorphous Fe(OH)3
    staq_all(findloc(chrsld_all,'amfe3',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
    ! Goethite; FeO(OH)
    staq_all(findloc(chrsld_all,'gt',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
    ! Hematite; Fe2O3
    staq_all(findloc(chrsld_all,'hm',dim=1), findloc(chraq_all,'fe3',dim=1)) = 2d0
    ! Ca-beidellite; Ca(1/6)Al(7/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0/6d0
    staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'al',dim=1)) = 7d0/3d0
    staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! Mg-beidellite; Mg(1/6)Al(7/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'mgbd',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0/6d0
    staq_all(findloc(chrsld_all,'mgbd',dim=1), findloc(chraq_all,'al',dim=1)) = 7d0/3d0
    staq_all(findloc(chrsld_all,'mgbd',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! K-beidellite; K(1/3)Al(7/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'kbd',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'kbd',dim=1), findloc(chraq_all,'al',dim=1)) = 7d0/3d0
    staq_all(findloc(chrsld_all,'kbd',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! Na-beidellite; Na(1/3)Al(7/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'nabd',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'nabd',dim=1), findloc(chraq_all,'al',dim=1)) = 7d0/3d0
    staq_all(findloc(chrsld_all,'nabd',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! Ca-saponite; Ca(1/6)Mg3Al(1/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'casp',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0/6d0
    staq_all(findloc(chrsld_all,'casp',dim=1), findloc(chraq_all,'mg',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'casp',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'casp',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! K-saponite; K(1/3)Mg3Al(1/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'ksp',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'ksp',dim=1), findloc(chraq_all,'mg',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'ksp',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'ksp',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! Na-saponite; Na(1/3)Mg3Al(1/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'nasp',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'nasp',dim=1), findloc(chraq_all,'mg',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'nasp',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'nasp',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! Mg-saponite; Mg(1/6)Mg3Al(1/3)Si(11/3)O10(OH)2
    staq_all(findloc(chrsld_all,'mgsp',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0/6d0 + 3d0
    staq_all(findloc(chrsld_all,'mgsp',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0/3d0
    staq_all(findloc(chrsld_all,'mgsp',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
    ! Illite; K0.6Mg0.25Al2.3Si3.5O10(OH)2
    staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'k',dim=1)) = 0.6d0
    staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'mg',dim=1)) = 0.25d0
    staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'al',dim=1)) = 2.3d0
    staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'si',dim=1)) = 3.5d0
    ! Diopside (MgCaSi2O6)
    staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    ! Hedenbergite (FeCaSi2O6)
    staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    ! Clinopyroxene (FexMg(1-x)CaSi2O6)
    staq_all(findloc(chrsld_all,'cpx',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'cpx',dim=1), findloc(chraq_all,'fe2',dim=1)) = fr_hb_cpx
    staq_all(findloc(chrsld_all,'cpx',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0 - fr_hb_cpx
    staq_all(findloc(chrsld_all,'cpx',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    ! Enstatite (MgSiO3)
    staq_all(findloc(chrsld_all,'en',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'en',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Ferrosilite (FeSiO3)
    staq_all(findloc(chrsld_all,'fer',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'fer',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Orthopyroxene (FexMg(1-x)SiO3)
    staq_all(findloc(chrsld_all,'opx',dim=1), findloc(chraq_all,'fe2',dim=1)) = fr_fer_opx
    staq_all(findloc(chrsld_all,'opx',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0 - fr_fer_opx
    staq_all(findloc(chrsld_all,'opx',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Jadeite (NaAlSi2O6)
    staq_all(findloc(chrsld_all,'jd',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'jd',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 
    staq_all(findloc(chrsld_all,'jd',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    ! Wollastonite (CaSiO3)
    staq_all(findloc(chrsld_all,'wls',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0 
    staq_all(findloc(chrsld_all,'wls',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Augite (Fe(xy+x)Mg(y-xy+1-x)Ca(1-y)Si2O6); x=fr_fer_agt ; y=fr_opx_agt 
    ! Augite (Fe(xy+x)(1-z)Mg(y-xy+1-x)(1-z)Ca(1-y)(1-z)NazAlzSi2O6); x=fr_fer_agt ; y=fr_opx_agt ; z=fr_napx_agt
    staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'fe2',dim=1)) &
        & = ( fr_fer_agt* (1d0 + fr_opx_agt) ) * (1d0 - fr_napx_agt) 
    staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'mg',dim=1)) &
        & = ( (1d0 - fr_fer_agt )*(fr_opx_agt + 1d0) ) * (1d0 - fr_napx_agt)
    staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'ca',dim=1)) & 
        & = ( 1d0 - fr_opx_agt ) * (1d0 - fr_napx_agt)
    staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'na',dim=1)) = fr_napx_agt
    staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'al',dim=1)) = fr_napx_agt
    staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    if ( &
        & abs( 2d0*staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'fe2',dim=1))  &
        & + 2d0*staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'mg',dim=1)) &
        & + 2d0*staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'ca',dim=1)) & 
        & + staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'na',dim=1)) &
        & + 3d0*staq_all(findloc(chrsld_all,'agt',dim=1), findloc(chraq_all,'al',dim=1)) &
        & - 4d0) &
        & /4d0 > tol &
        & ) then 
        print *, '*** ERROR in stoichiometry of augite'
        stop
    endif 
    ! print *,( fr_fer_agt* (1d0 + fr_opx_agt) ) * (1d0 - fr_napx_agt)
    ! print *,( (1d0 - fr_fer_agt )*(fr_opx_agt + 1d0) ) * (1d0 - fr_napx_agt)
    ! print *,( 1d0 - fr_opx_agt ) * (1d0 - fr_napx_agt)
    ! stop
    ! Tremolite (Ca2Mg5(Si8O22)(OH)2)
    staq_all(findloc(chrsld_all,'tm',dim=1), findloc(chraq_all,'ca',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'tm',dim=1), findloc(chraq_all,'mg',dim=1)) = 5d0
    staq_all(findloc(chrsld_all,'tm',dim=1), findloc(chraq_all,'si',dim=1)) = 8d0
    ! Anthophyllite (Mg2Mg5(Si8O22)(OH)2)
    staq_all(findloc(chrsld_all,'antp',dim=1), findloc(chraq_all,'mg',dim=1)) = 7d0
    staq_all(findloc(chrsld_all,'antp',dim=1), findloc(chraq_all,'si',dim=1)) = 8d0
    ! Muscovite; KAl2(AlSi3O10)(OH)2
    staq_all(findloc(chrsld_all,'mscv',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'mscv',dim=1), findloc(chraq_all,'al',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'mscv',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    ! Phlogopite; KMg3(AlSi3O10)(OH)2
    staq_all(findloc(chrsld_all,'plgp',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'plgp',dim=1), findloc(chraq_all,'mg',dim=1)) = 3d0
    staq_all(findloc(chrsld_all,'plgp',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'plgp',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    ! Sepiolite (Mg4Si6O15(OH)2:6H2O)
    staq_all(findloc(chrsld_all,'splt',dim=1), findloc(chraq_all,'mg',dim=1)) = 4d0
    staq_all(findloc(chrsld_all,'splt',dim=1), findloc(chraq_all,'si',dim=1)) = 6d0
    ! Amorphous silica; SiO2
    staq_all(findloc(chrsld_all,'amsi',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Phytolith silica; SiO2
    staq_all(findloc(chrsld_all,'phsi',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Quartz; SiO2
    staq_all(findloc(chrsld_all,'qtz',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Aragonite (CaCO3)
    staq_all(findloc(chrsld_all,'arg',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    stgas_all(findloc(chrsld_all,'arg',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    ! Dolomite (CaMg(CO3)2)
    staq_all(findloc(chrsld_all,'dlm',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'dlm',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
    stgas_all(findloc(chrsld_all,'dlm',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 2d0
    ! Gypsum; CaSO4*2H2O
    staq_all(findloc(chrsld_all,'gps',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'gps',dim=1), findloc(chraq_all,'so4',dim=1)) = 1d0
    ! Anhydrite; CaSO4
    staq_all(findloc(chrsld_all,'caso4',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'caso4',dim=1), findloc(chraq_all,'so4',dim=1)) = 1d0
    ! Ferrous oxide; FeO
    staq_all(findloc(chrsld_all,'fe2o',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
    ! Periclase; MgO
    staq_all(findloc(chrsld_all,'mgo',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
    ! Dipotasium monoxide; K2O
    staq_all(findloc(chrsld_all,'k2o',dim=1), findloc(chraq_all,'k',dim=1)) = 2d0
    ! Calcium oxide; CaO
    staq_all(findloc(chrsld_all,'cao',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    ! Disodium monoxide; Na2O
    staq_all(findloc(chrsld_all,'na2o',dim=1), findloc(chraq_all,'na',dim=1)) = 2d0
    ! Corundum; Al2O3
    staq_all(findloc(chrsld_all,'al2o3',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
    ! hypothetical SiO2 
    staq_all(findloc(chrsld_all,'sio2',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
    ! Glass basalt
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'si',dim=1)) = fr_si_gbas
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'al',dim=1)) = fr_al_gbas
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'na',dim=1)) = fr_na_gbas
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'k',dim=1)) = fr_k_gbas
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'mg',dim=1)) = fr_mg_gbas
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_ca_gbas
    staq_all(findloc(chrsld_all,'gbas',dim=1), findloc(chraq_all,'fe2',dim=1)) = fr_fe2_gbas
    ! crystaline basalt
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'si',dim=1)) = fr_si_cbas
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'al',dim=1)) = fr_al_cbas
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'na',dim=1)) = fr_na_cbas
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'k',dim=1)) = fr_k_cbas
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'mg',dim=1)) = fr_mg_cbas
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_ca_cbas
    staq_all(findloc(chrsld_all,'cbas',dim=1), findloc(chraq_all,'fe2',dim=1)) = fr_fe2_cbas
    ! Epidote (Ca2FeAl2Si3O12OH)
    staq_all(findloc(chrsld_all,'ep',dim=1), findloc(chraq_all,'ca',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'ep',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'ep',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'ep',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    ! Clinochlore (Mg5Al2Si3O10(OH)8)
    staq_all(findloc(chrsld_all,'clch',dim=1), findloc(chraq_all,'mg',dim=1)) = 5d0
    staq_all(findloc(chrsld_all,'clch',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'clch',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
    ! Cordierite (Mg2Al4Si5O18)
    staq_all(findloc(chrsld_all,'cdr',dim=1), findloc(chraq_all,'mg',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'cdr',dim=1), findloc(chraq_all,'al',dim=1)) = 4d0
    staq_all(findloc(chrsld_all,'cdr',dim=1), findloc(chraq_all,'si',dim=1)) = 5d0
    ! Leucite (KAlSi206)
    staq_all(findloc(chrsld_all,'leu',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'leu',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'leu',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
    ! OMs; CH2O
    stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
    stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = n2c_g1

    stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0 - 2d0*oxa2c_g2 - 2d0*ac2c_g2
    stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
    stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = n2c_g2
    staq_all(findloc(chrsld_all,'g2',dim=1), findloc(chraq_all,'oxa',dim=1)) = oxa2c_g2
    staq_all(findloc(chrsld_all,'g2',dim=1), findloc(chraq_all,'ac',dim=1)) = ac2c_g2

    stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
    stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = n2c_g3
    ! the above need to be modified to enable anoxic degradation 

    ! Fertilizers 
    ! Ammonium nitrate (NH4NO3); assuming simplified reaction implementable without enabling full N cycles 
    ! NH4NO3 = NH4+ + NO3- coupling with (NH4+ + 2O2 -> NO3- + H2O + 2 H+) 
    ! NH4NO3 + 2O2 = 2NO3- + H2O + 2H+
    staq_all(findloc(chrsld_all,'amnt',dim=1), findloc(chraq_all,'no3',dim=1)) = 2d0
    stgas_all(findloc(chrsld_all,'amnt',dim=1), findloc(chrgas_all,'po2',dim=1)) = -2d0

    ! when fully doing N cycle
    ! staq_all(findloc(chrsld_all,'amnt',dim=1), findloc(chraq_all,'no3',dim=1)) = 1d0 
    ! stgas_all(findloc(chrsld_all,'amnt',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = 1d0

    ! KCl solid
    staq_all(findloc(chrsld_all,'kcl',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'kcl',dim=1), findloc(chraq_all,'cl',dim=1)) = 1d0

    ! CH3COOH solid
    staq_all(findloc(chrsld_all,'gac',dim=1), findloc(chraq_all,'ac',dim=1)) = 1d0

    ! MES monohydrate
    staq_all(findloc(chrsld_all,'mesmh',dim=1), findloc(chraq_all,'mes',dim=1)) = 1d0

    ! Imidazole solid 
    staq_all(findloc(chrsld_all,'ims',dim=1), findloc(chraq_all,'im',dim=1)) = 1d0

    ! TriEthanoleAmine solid 
    staq_all(findloc(chrsld_all,'teas',dim=1), findloc(chraq_all,'tea',dim=1)) = 1d0

    ! NaOH solid 
    staq_all(findloc(chrsld_all,'naoh',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0

    ! Sodium glycerophosphate
    staq_all(findloc(chrsld_all,'naglp',dim=1), findloc(chraq_all,'na',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'naglp',dim=1), findloc(chraq_all,'glp',dim=1)) = 1d0

    ! hydrophilite: CaCl2 
    staq_all(findloc(chrsld_all,'cacl2',dim=1), findloc(chraq_all,'cl',dim=1)) = 2d0
    staq_all(findloc(chrsld_all,'cacl2',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0

    ! halite: NaCl
    staq_all(findloc(chrsld_all,'nacl',dim=1), findloc(chraq_all,'cl',dim=1)) = 1d0
    staq_all(findloc(chrsld_all,'nacl',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0


    staq = 0d0
    stgas = 0d0

    do isps = 1, nsp_sld
        if (any(chrsld_all == chrsld(isps))) then 
            do ispa = 1, nsp_aq 
                if (any(chraq_all == chraq(ispa))) then 
                    staq(isps,ispa) = &
                        & staq_all(findloc(chrsld_all,chrsld(isps),dim=1), findloc(chraq_all,chraq(ispa),dim=1))
                endif 
            enddo 
            do ispg = 1, nsp_gas 
                if (any(chrgas_all == chrgas(ispg))) then 
                    stgas(isps,ispg) = &
                        & stgas_all(findloc(chrsld_all,chrsld(isps),dim=1), findloc(chrgas_all,chrgas(ispg),dim=1))
                endif 
            enddo 
        endif 
    enddo 

    ! external reactions
    staq_ext_all = 0d0
    stgas_ext_all = 0d0
    stsld_ext_all = 0d0
    ! respiration 
    stgas_ext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
    ! fe2 oxidation 
    staq_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe2',dim=1)) = -1d0
    staq_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0/4d0
    ! SOC assimilation by microbes 
    stsld_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'om',dim=1)) = -1d0
    stsld_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'omb',dim=1)) = 0.31d0
    stgas_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 0.69d0
    stgas_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrgas_all,'po2',dim=1)) = -0.69d0
    ! turnover of microbes 
    stsld_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'om',dim=1)) = 1d0
    stsld_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'omb',dim=1)) = -1d0
    ! pyrite oxidation by fe3
    stsld_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chrsld_all,'py',dim=1)) = -1d0
    staq_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe3',dim=1)) = -14d0
    staq_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe2',dim=1)) = 15d0
    staq_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'so4',dim=1)) = 2d0
    ! ammonia oxidation by O2 (NH4+ + 2O2 -> NO3- + H2O + 2 H+) 
    staq_ext_all(findloc(chrrxn_ext_all,'amo2o',dim=1), findloc(chraq_all,'no3',dim=1)) = 1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'amo2o',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = -1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'amo2o',dim=1), findloc(chrgas_all,'po2',dim=1)) = -2d0
    ! overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
    staq_ext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chraq_all,'no3',dim=1)) = -4d0/5d0 ! values relative to CH2O 
    stsld_ext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chrsld_all,'g2',dim=1)) = -1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = n2c_g2
    ! stgas_ext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chrgas_all,'pn2',dim=1)) = 2d0/5d0 ! should be added after enabling pn2 
    ! first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O) 
    staq_ext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chraq_all,'no3',dim=1)) = -1d0  ! values relative to CH2O
    stsld_ext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chrsld_all,'g2',dim=1)) = -1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chrgas_all,'pn2o',dim=1)) = 0.5d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = n2c_g2
    ! 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O) 
    stsld_ext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrsld_all,'g2',dim=1)) = -1d0 ! values relative to CH2O
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrgas_all,'pn2o',dim=1)) = -2d0
    stgas_ext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = n2c_g2
    ! stgas_ext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrgas_all,'pn2',dim=1)) = 2d0 ! should be added after enabling pn2 
    ! oxalate oxidation to CO2 (2 H2C2O4 + O2 -> 4 CO2 + 2 H2O)
    staq_ext_all(findloc(chrrxn_ext_all,'oxao2',dim=1), findloc(chraq_all,'oxa',dim=1)) = -2d0 
    stgas_ext_all(findloc(chrrxn_ext_all,'oxao2',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 4d0
    stgas_ext_all(findloc(chrrxn_ext_all,'oxao2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
    ! cation uptake by OM 
    staq_ext_all(findloc(chrrxn_ext_all,'g2k',dim=1), findloc(chraq_all,'k',dim=1)) = -1d0 
    staq_ext_all(findloc(chrrxn_ext_all,'g2ca',dim=1), findloc(chraq_all,'ca',dim=1)) = -1d0 
    staq_ext_all(findloc(chrrxn_ext_all,'g2mg',dim=1), findloc(chraq_all,'mg',dim=1)) = -1d0 

    ! define 1 when a reaction is sensitive to a speces 
    stgas_dext_all = 0d0
    staq_dext_all = 0d0
    stsld_dext_all = 0d0
    ! respiration 
    stgas_dext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    ! fe2 oxidation 
    stgas_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
    staq_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
    ! SOC assimilation by microbes 
    stsld_dext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'om',dim=1)) = 1d0
    stsld_dext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'omb',dim=1)) = 1d0
    ! turnover of microbes 
    stsld_dext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'omb',dim=1)) = 1d0
    ! pyrite oxidation by fe3
    stsld_dext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chrsld_all,'py',dim=1)) = 1d0
    staq_dext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
    staq_dext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
    ! ammonia oxidation by O2 (NH4+ + 2O2 -> NO3- + H2O + 2 H+) 
    stgas_dext_all(findloc(chrrxn_ext_all,'amo2o',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'amo2o',dim=1), findloc(chrgas_all,'pnh3',dim=1)) = 1d0
    ! overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
    staq_dext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chraq_all,'no3',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    stsld_dext_all(findloc(chrrxn_ext_all,'g2n0',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
    ! first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O) 
    staq_dext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chraq_all,'no3',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    stsld_dext_all(findloc(chrrxn_ext_all,'g2n21',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
    ! 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O) 
    staq_dext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chraq_all,'no3',dim=1)) = 1d0
    ! stgas_dext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrgas_all,'pn2o',dim=1)) = 1d0
    stsld_dext_all(findloc(chrrxn_ext_all,'g2n22',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
    ! oxalate oxidation to CO2 (2 H2C2O4 + O2 -> 4 CO2 + 2 H2O)
    staq_dext_all(findloc(chrrxn_ext_all,'oxao2',dim=1), findloc(chraq_all,'oxa',dim=1)) = 1d0
    ! cation uptake by OM
    stsld_dext_all(findloc(chrrxn_ext_all,'g2k',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'g2k',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    staq_dext_all(findloc(chrrxn_ext_all,'g2k',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
    stsld_dext_all(findloc(chrrxn_ext_all,'g2ca',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'g2ca',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    staq_dext_all(findloc(chrrxn_ext_all,'g2ca',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
    stsld_dext_all(findloc(chrrxn_ext_all,'g2mg',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
    stgas_dext_all(findloc(chrrxn_ext_all,'g2mg',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
    staq_dext_all(findloc(chrrxn_ext_all,'g2mg',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0

    staq_ext = 0d0
    stgas_ext = 0d0
    stsld_ext = 0d0

    do irxn = 1, nrxn_ext
        if (any(chrrxn_ext_all == chrrxn_ext(irxn))) then 
            do ispa = 1, nsp_aq 
                if (any(chraq_all == chraq(ispa))) then 
                    staq_ext(irxn,ispa) = &
                        & staq_ext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                        &   ,findloc(chraq_all,chraq(ispa),dim=1))
                    staq_dext(irxn,ispa) = &
                        & staq_dext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                        &   ,findloc(chraq_all,chraq(ispa),dim=1))
                endif 
            enddo 
            do ispg = 1, nsp_gas 
                if (any(chrgas_all == chrgas(ispg))) then 
                    stgas_ext(irxn,ispg) = &
                        & stgas_ext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                        &   ,findloc(chrgas_all,chrgas(ispg),dim=1))
                    stgas_dext(irxn,ispg) = &
                        & stgas_dext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                        &   ,findloc(chrgas_all,chrgas(ispg),dim=1))
                endif 
            enddo 
            do isps = 1, nsp_sld 
                if (any(chrsld_all == chrsld(isps))) then 
                    stsld_ext(irxn,isps) = &
                        & stsld_ext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                        &   ,findloc(chrsld_all,chrsld(isps),dim=1))
                    stsld_dext(irxn,isps) = &
                        & stsld_dext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                        &   ,findloc(chrsld_all,chrsld(isps),dim=1))
                endif 
            enddo 
        endif 
    enddo 


    def_dust = 0d0
        
    call get_dust( &
        & nsp_sld_all,chrsld_all,def_dust &! input
        & ,rfrc_sld_all &! output
        & )

    ! added 2nd type of dust (fertilizer/manure etc.)
    call get_dust_2nd( &
        & nsp_sld_all,chrsld_all,def_dust &! input
        & ,rfrc_sld_all_2nd &! output
        & )

    rfrc_sld_all = rfrc_sld_all/mwt_all
    rfrc_sld_all_2nd = rfrc_sld_all_2nd/mwt_all ! added 
    ! rfrc_sld_all = rfrc_sld_all/sum(rfrc_sld_all)




    ! rfrc_sld_plant_all(findloc(chrsld_all,'om',dim=1)) = 1d0

    ! rfrc_sld_plant_all(findloc(chrsld_all,'g1',dim=1)) = 0.1d0
    ! rfrc_sld_plant_all(findloc(chrsld_all,'g2',dim=1)) = 0.8d0
    ! rfrc_sld_plant_all(findloc(chrsld_all,'g3',dim=1)) = 0.1d0

    def_OM_frc = 0d0

    call get_OM_rain( &
        & nsp_sld_all,chrsld_all,def_OM_frc &! input
        & ,rfrc_sld_plant_all &! output
        & )

    ! rfrc_sld_plant_all = rfrc_sld_plant_all/mwt_all
    ! rfrc_sld_plant_all = rfrc_sld_plant_all/sum(rfrc_sld_plant_all)

    dust_step = .false.
    if (step_tau >0d0 .and. step_tau <1d0)dust_step = .true.


    do isps = 1, nsp_sld 
        rfrc_sld(isps) = rfrc_sld_all(findloc(chrsld_all,chrsld(isps),dim=1))
        rfrc_sld_2nd(isps) = rfrc_sld_all_2nd(findloc(chrsld_all,chrsld(isps),dim=1))
        rfrc_sld_plant(isps) = rfrc_sld_plant_all(findloc(chrsld_all,chrsld(isps),dim=1))
    enddo


    call get_switches( &
        & iwtype,imixtype,poroiter_in,display,display_lim_in,read_data,incld_rough &
        & ,act_ON,timestep_fixed,ads_ON,regular_grid,aq_close &! inout
        & ,poroevol,surfevol1,surfevol2,do_psd,lim_minsld_in,do_psd_full,season &!
        & )

    select case(imixtype)
        case(imixtype_nobio)
            print*, 'no bioturbation'
        case(imixtype_fick)
            print*, 'Fickian mixing'
        case(imixtype_turbo2)
            print*, 'homogeneous mixing'
        case(imixtype_till)
            print*, 'tilling'
        case(imixtype_labs)
            print*, 'LABS mixing'
        case default 
            print *, '***| chosen number is not available for mixing styles (choose between 0 to 4)'
            print *, '***| thus choose default |---- > no mixing'
            imixtype = imixtype_nobio
    endselect 

    imix                = imixtype
    imixtype_OM         = imixtype_OM_in
    imixtype_background = imixtype_background_in     ! in default it is fickian assuming driver of background mixing is the same as driver of OM mixing

    zml_background  = zsupp                 ! assume background mixing is equivalent as the driver of OM mixing
    zml_OM          = zsupp                 ! zsupp is mixing depth for OM 
    zml_dust        = zml_ref               ! zml_ref is mixing depth for dust

    select case(iwtype)
        case(iwtype_cnst)
            print *, 'const w',iwtype
        case(iwtype_flex)
            print *, 'w flex (cnst porosity profile)',iwtype
        case(iwtype_pwcnst)
            print *, 'w x porosity = cnst',iwtype
        case(iwtype_spwcnst)
            print *, 'w x (1 - porosity) = cnst',iwtype
        case default 
            print *, '***| chosen number is not available for advection styles (choose between 0 to 3)'
            print *, '***| thus choose default |---- > cnst w'
            iwtype = iwtype_cnst
    endselect 

    if (poroiter_in) then 
        print *, 'porosity iteration is ON'
    else 
        print *, 'porosity iteration is OFF'
    endif 

    if (lim_minsld_in) then 
        print *, 'limiting lowest mineral conc. is ON'
    else 
        print *, 'limiting lowest mineral conc. is OFF'
    endif 

    if (do_psd_full) do_psd = .true.

    if (display_lim_in) display_lim = .true.

    if (sld_enforce) nsp3 = nsp_aq + nsp_gas ! excluding solid phases


    ! kinetic formulation type
    precstyle = 'def'
    ! precstyle = 'full'
    ! precstyle = 'full_lim'
    ! precstyle = 'seed '
    ! precstyle = '2/3'
    ! precstyle = 'psd_full'
    ! precstyle = 'emmanuel' ! enables solubility change by a const factor or as a function of radius, T and interfacial energy (as in Emmanuel and Ague, 2011, Chem Geol)

    solmod = 1d0
    fkin   = 1d0

    do isps = 1, nsp_sld
        select case(trim(adjustl(chrsld(isps))))
            case('g1','g2','g3','amnt','inrt','kcl','gac','mesmh','ims','teas','naoh','naglp','cacl2','nacl')
                precstyle(isps) = 'decay'
            case('cc','arg','dlm') ! added to change solubility 
                precstyle(isps) = 'def'
                ! precstyle(isps) = 'emmanuel'
                ! solmod(isps,:) = 0.1d0 ! assumed factor to be multiplied with omega
            case('casp','ksp','nasp','mgsp')
                precstyle(isps) = 'def'
                ! precstyle(isps) = 'emmanuel'
                ! solmod(isps,:) = 0.05d0 ! assumed factor to be multiplied with omega
                ! fkin(isps,:) = 0.01d0
            case default 
                precstyle(isps) = 'def'
                ! precstyle(isps) = '2/3'
                ! precstyle(isps) = '2/3noporo'
                ! precstyle(isps) = 'psd_full'
                ! fkin(isps,:) = 0.01d0
        endselect
    enddo 

    cec_pH_depend=.true.
    ! cec_pH_depend=.false.

    call get_base_charge( &
        & nsp_aq_all & 
        & ,chraq_all & 
        & ,base_charge_all &! output 
        & )

    ! detault cec
    mcec_all_def = 0d0
    do isps=1,nsp_sld_all
        ! if not tracked assume zero cec
        ! for trakced species assign default values 
        if ( any(chrsld == chrsld_all(isps)) ) then
            select case(trim(adjustl(chrsld_all(isps))))
                case('ka')
                    mcec_all_def(isps) = 16.2d0 ! 16.2 cmol/kg (from Beerling et al. 2020)
                case('inrt')
                    mcec_all_def(isps) = 0d0 
                case('cabd','mgbd','kbd','nabd') 
                    mcec_all_def(isps) = 70d0 !  70 cmol/kg (from Parfitt et al. 1996)
                case('g1','g2','g3')
                    mcec_all_def(isps) = 330d0 ! 330 cmol/kg (from Parfitt et al. 1996)
                case default 
                    mcec_all_def(isps) = 0d0
            endselect
        endif 
    enddo 

    ! detault logKH\Na
    logkhaq_all_def = 0d0
    do isps=1,nsp_sld_all
        if (mcec_all_def(isps)>0d0) then 
            do ispa=1,nsp_aq_all
                select case( trim(adjustl(chraq_all(ispa))) )
                    case('na')
                        logkhaq_all_def(isps,ispa) = 5.9d0
                    case('k')
                        logkhaq_all_def(isps,ispa) = 4.8d0
                    case('ca')
                        logkhaq_all_def(isps,ispa) = 10.47d0
                    case('mg')
                        logkhaq_all_def(isps,ispa) = 10.786d0
                    case('al')
                        logkhaq_all_def(isps,ispa) = 16.47d0
                endselect
            enddo
        endif 
    enddo

    ! detault beta
    beta_all_def = 3.4d0


    do ispa = 1, nsp_aq    
        if (any(chraq_all == chraq(ispa))) base_charge(ispa) = base_charge_all(findloc(chraq_all,chraq(ispa),dim=1))
    enddo 

    rectime_flx = 0d0
    do irec_flx = 1,20
        rectime_flx(irec_flx) = irec_flx/20d0
    enddo
    do irec_flx = 21,38
        rectime_flx(irec_flx) = rectime_flx(20) + (irec_flx-20)/20d0*10d0
    enddo
    do irec_flx = 39,60
        rectime_flx(irec_flx) = rectime_flx(38) + (irec_flx-38)/20d0*100d0
    enddo

    if (linear_rectime)  rectime_prof =  (/(irec_prof*ttot/nrec_prof, irec_prof = 1,nrec_prof)/)

    if (rectime_scheme_old) then 
        do while (rectime_flx(nrec_flx)>ttot) 
            rectime_flx = rectime_flx/10d0
        enddo 
        do while (rectime_flx(nrec_flx)<ttot) 
            rectime_flx = rectime_flx*10d0
        enddo 

        do while (rectime_prof(nrec_prof)>ttot) 
            rectime_prof = rectime_prof/10d0
        enddo 
        do while (rectime_prof(nrec_prof)<ttot) 
            rectime_prof = rectime_prof*10d0
        enddo 
        
        savetime = rectime_prof(18)/100d0
        dsavetime = rectime_prof(18)/100d0
        
    else 
        rectime_flx = rectime_flx*ttot/maxval(rectime_flx)
        rectime_prof = rectime_prof*ttot/maxval(rectime_prof)
        savetime = ttot/100d0
        dsavetime = ttot/100d0
    endif 

    ! print*, rectime_flx
    ! stop

    ! write(chrq(1),'(i0)') int(qin/(10d0**(floor(log10(qin)))))
    ! write(chrq(2),'(i0)') floor(log10(qin))
    ! chrq(3) = trim(adjustl(chrq(1)))//'E'//trim(adjustl(chrq(2)))
    write(chrq(3),'(E10.2)') qin
    write(chrz(3),'(i0)') nint(zsat)
    write(chrrain,'(E10.2)') rainpowder


    ! write(workdir,*) '../pyweath_output/'     
    write(workdir,*) './'    
    write(flxdir,*) './flx'    
    write(profdir,*) './prof'     

    ! if (cplprec) then 
        ! write(base,*) 'test_cplp_test'
    ! else 
        ! write(base,*) 'test_cpl'
    ! endif 

    base = trim(adjustl(sim_name))

    if (al_inhibit) base = trim(adjustl(base))//'_alx'

    base = trim(adjustl(base))//'_rain-'//trim(adjustl(chrrain))    
    
    if (poroevol) then      
        base = trim(adjustl(base))//'_pevol'
    endif 
    if (surfevol1)then 
        base = trim(adjustl(base))//'_sevol1'
    elseif (surfevol2) then 
        base = trim(adjustl(base))//'_sevol2'
    #if defined(surfssa)
        base = trim(adjustl(base))//'_ssa'
    #endif 
    endif 

    if (.not. regular_grid) then 
        base = trim(adjustl(base))//'_irr'
    endif 

    if (dust_wave)then 
        write(chrrain,'(E10.2)') wave_tau
        base = trim(adjustl(base))//'_rwave-'//trim(adjustl(chrrain))
    endif 

    if (incld_rough)then 
        write(chrrain,'(E10.2)') p80
        base = trim(adjustl(base))//'_p80r-'//trim(adjustl(chrrain))
    else
        write(chrrain,'(E10.2)') p80
        base = trim(adjustl(base))//'_p80-'//trim(adjustl(chrrain))
    endif 

    write(runname,*) trim(adjustl(base))//'_q-'//trim(adjustl(chrq(3)))//'_zsat-'  &
        & //trim(adjustl(chrz(3)))
        
    ! directly name runname from input 
    ! write(runname,*) trim(adjustl(sim_name))
    write(runname,*) 'output'

    #ifdef full_flux_report
    do isps = 1, nsp_sld 
        do iz = 1, nz
            isldflx(isps,iz) = idust + (isps-1)*nz + iz
            ! print *,isldflx(isps,iz)
        enddo 
    enddo     
    do ispa = 1, nsp_aq 
        do iz=1,nz
            iaqflx(ispa,iz) = idust + nsp_sld*nz  + (ispa - 1)*nz + iz
            ! print *,iaqflx(ispa,iz)
        enddo
    enddo 

    do ispg = 1, nsp_gas
        do iz= 1,nz
            igasflx(ispg,iz) = idust + nsp_sld*nz + nsp_aq*nz + (ispg-1)*nz + iz
            ! print*,igasflx(ispg,iz)
        enddo 
    enddo 

    do ico2 = 1, 6
        do iz = 1, nz
            ico2flx(ico2,iz) = idust + nsp_sld*nz + nsp_aq*nz + nsp_gas*nz + (ico2 - 1)*nz + iz
            ! print*,ico2flx(ico2,iz)
        enddo 
    enddo 
    ! pause
    #else 
    do isps = 1, nsp_sld 
        isldflx(isps) = idust + isps
    enddo 
        
    do ispa = 1, nsp_aq 
        iaqflx(ispa) = idust + nsp_sld  + ispa
    enddo 

    do ispg = 1, nsp_gas
        igasflx(ispg) = idust + nsp_sld + nsp_aq + ispg
    enddo 

    do ico2 = 1, 6
        ico2flx(ico2) = idust + nsp_sld + nsp_aq + nsp_gas + ico2
    enddo 

    iphint  = idust + nsp_sld + nsp_aq + nsp_gas + 7
    iphint2 = idust + nsp_sld + nsp_aq + nsp_gas + 8

    #endif 

    ! print*,workdir
    ! print*,runname
    ! pause

    ! call system ('mkdir -p '//trim(adjustl(workdir))//trim(adjustl(runname)))
    call system ('mkdir -p '//trim(adjustl(flxdir)))
    call system ('mkdir -p '//trim(adjustl(profdir)))

    ! call system ('cp gases.in solutes.in slds.in extrxns.in '//trim(adjustl(workdir))//trim(adjustl(runname)))

    #ifdef full_flux_report

    write(chrfmt,'(i0)') nflx+2

    chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a))'

    do isps = 1,nsp_sld
        do iz = 1,nz
            write(chriz,'(i3.3)') iz
            open(isldflx(isps,iz), file=trim(adjustl(flxdir))//'/' &
                & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'-'//trim(adjustl(chriz))//'.txt' &
                & , status='replace')
            write(isldflx(isps,iz),trim(adjustl(chrfmt))) 'time','z',(chrflx(iflx),iflx=1,nflx)
            close(isldflx(isps,iz))
        enddo 
    enddo 

    do ispa = 1,nsp_aq
        do iz= 1,nz
            write(chriz,'(i3.3)') iz
            open(iaqflx(ispa,iz), file=trim(adjustl(flxdir))//'/' &
                & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'-'//trim(adjustl(chriz))//'.txt' &
                & , status='replace')
            write(iaqflx(ispa,iz),trim(adjustl(chrfmt))) 'time','z',(chrflx(iflx),iflx=1,nflx)
            close(iaqflx(ispa,iz))
        enddo 
    enddo 

    do ispg = 1,nsp_gas
        do iz=1,nz
            write(chriz,'(i3.3)') iz
            open(igasflx(ispg,iz), file=trim(adjustl(flxdir))//'/' &
                & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'-'//trim(adjustl(chriz))//'.txt' &
                & , status='replace')
            write(igasflx(ispg,iz),trim(adjustl(chrfmt))) 'time','z',(chrflx(iflx),iflx=1,nflx)
            close(igasflx(ispg,iz))
        enddo 
    enddo 

    do ico2 = 1,6
        do iz= 1,nz
            write(chriz,'(i3.3)') iz
            open(ico2flx(ico2,iz), file=trim(adjustl(flxdir))//'/' &
                & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'-'//trim(adjustl(chriz))//'.txt' &
                & , status='replace')
            write(ico2flx(ico2,iz),trim(adjustl(chrfmt))) 'time','z',(chrflx(iflx),iflx=1,nflx)
            close(ico2flx(ico2,iz))
        enddo 
    enddo 

    #else 

    write(chrfmt,'(i0)') nflx+1

    chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a))'

    do isps = 1,nsp_sld
        open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
            & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', status='replace')
        write(isldflx(isps),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(isldflx(isps))
        
        open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
            & //'int_flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', status='replace')
        write(isldflx(isps),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(isldflx(isps))
        
        if ( any ( chrsp_saveall == chrsld(isps)) ) then
            open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
                & //'flx_all_sld-'//trim(adjustl(chrsld(isps)))//'.txt', status='replace')
            write(isldflx(isps),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
            close(isldflx(isps))
        endif
    enddo 

    do ispa = 1,nsp_aq
        open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
            & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', status='replace')
        write(iaqflx(ispa),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(iaqflx(ispa))
        
        open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
            & //'int_flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', status='replace')
        write(iaqflx(ispa),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(iaqflx(ispa))
        
        if ( any ( chrsp_saveall == chraq(ispa) ) ) then
            open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
                & //'flx_all_aq-'//trim(adjustl(chraq(ispa)))//'.txt', status='replace')
            write(iaqflx(ispa),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
            close(iaqflx(ispa))
        endif 
    enddo 

    do ispg = 1,nsp_gas
        open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
            & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', status='replace')
        write(igasflx(ispg),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(igasflx(ispg))
        
        open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
            & //'int_flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', status='replace')
        write(igasflx(ispg),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(igasflx(ispg))
        
        if ( any ( chrsp_saveall == chrgas(ispg) ) ) then
            open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
                & //'flx_all_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', status='replace')
            write(igasflx(ispg),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
            close(igasflx(ispg))
        endif 
    enddo 

    do ico2 = 1,6
        open(ico2flx(ico2), file=trim(adjustl(flxdir))//'/' &
            & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', status='replace')
        write(ico2flx(ico2),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(ico2flx(ico2))
        
        open(ico2flx(ico2), file=trim(adjustl(flxdir))//'/' &
            & //'int_flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', status='replace')
        write(ico2flx(ico2),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
        close(ico2flx(ico2))
    enddo 

    #endif 

    open(idust, file=trim(adjustl(flxdir))//'/'//'dust.txt', &
        & status='replace')
    write(idust,*) ' time ', ' dust(relative_to_average) '
    close(idust)

    climate(:) = .false.
    if (season) climate(:) = .true.

    open(idust, file=trim(adjustl(flxdir))//'/'//'climate.txt', &
        & status='replace')
    write(idust,*) ' time ', ' T(oC) ', ' q(m/yr) ', ' Wet(-) '
    close(idust)

    clim_file = (/'T_temp.in  ','q_temp.in  ','Wet_temp.in'/)

    do iclim = 1,3
        if (climate(iclim)) then 
            call get_clim_num( &
                & clim_file(iclim) &! in 
                & ,nclim(iclim) &! output
                & ) 
            select case (iclim) 
                case(1)
                    if ( allocated(clim_T) ) deallocate(clim_T)
                    allocate(clim_T(2,nclim(iclim)))    
                    open(idust,file=trim(adjustl(workdir))//'/'//trim(adjustl(clim_file(iclim))),  &
                        & status ='old',action='read')
                    read (idust,'()')
                    clim_T = 0d0
                    do ict = 1, nclim(iclim)
                        read (idust,*) clim_T(1,ict),clim_T(2,ict)
                    enddo 
                    close(idust)
                    ! print *
                    ! do ict = 1, nclim(iclim)
                        ! print *, clim_T(:,ict)
                    ! enddo 
                    dct(iclim) = clim_T(1,2) - clim_T(1,1)
                    ctau(iclim) = clim_T(1,nclim(iclim)) + dct(iclim)
                case(2)
                    if ( allocated(clim_q) ) deallocate(clim_q)
                    allocate(clim_q(2,nclim(iclim)))
                    open(idust,file=trim(adjustl(workdir))//'/'//trim(adjustl(clim_file(iclim))),  &
                        & status ='old',action='read')
                    read (idust,'()')
                    clim_q = 0d0
                    do ict = 1, nclim(iclim)
                        read (idust,*) clim_q(1,ict),clim_q(2,ict)
                    enddo 
                    close(idust)
                    ! converting mm/month to m/yr 
                    clim_q(2,:) = clim_q(2,:)*12d0/1d3
                    ! print *
                    ! do ict = 1, nclim(iclim)
                        ! print *, clim_q(:,ict)
                    ! enddo 
                    dct(iclim) = clim_q(1,2) - clim_q(1,1)
                    ctau(iclim) = clim_q(1,nclim(iclim)) + dct(iclim)
                case(3)
                    if ( allocated(clim_sat) ) deallocate(clim_sat)
                    allocate(clim_sat(2,nclim(iclim)))
                    open(idust,file=trim(adjustl(workdir))//'/'//trim(adjustl(clim_file(iclim))),  &
                        & status ='old',action='read')
                    read (idust,'()')
                    clim_sat = 0d0
                    do ict = 1, nclim(iclim)
                        read (idust,*) clim_sat(1,ict),clim_sat(2,ict)
                    enddo 
                    close(idust)
                    ! converting mm/m to m/m 
                    clim_sat(2,:) = clim_sat(2,:)*1d0/1d3
                    ! print *
                    ! do ict = 1, nclim(iclim)
                        ! print *, clim_sat(:,ict)
                    ! enddo 
                    dct(iclim) = clim_sat(1,2) - clim_sat(1,1)
                    ctau(iclim) = clim_sat(1,nclim(iclim)) + dct(iclim)
                case default
                    print*, 'error in obtaining climate'
                    stop
            endselect 
        endif 
    enddo 

    ! stop

    !!!  MAKING GRID !!!!!!!!!!!!!!!!! 
    beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
    beta = 1.00005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
    call makegrid(beta,nz,ztot,dz,z,regular_grid)
        
        
    open(iphint, file=trim(adjustl(flxdir))//'/'//'int_ph.txt', status='replace')
    write(iphint,*) 'time\depth',(z(iz),iz=1,nz)
    close(iphint)

        
    open(iphint2, file=trim(adjustl(flxdir))//'/'//'ph.txt', status='replace')
    write(iphint2,*) 'time\depth',(z(iz),iz=1,nz)
    close(iphint2)


    sat = min(1.0d0,(1d0-satup)*z/zsat + satup)
    #ifdef satconvex 
    sat = min(1.0d0, satup+(1d0-satup)*(z/zsat)**2d0)
    #endif 
    #ifdef satconcave 
    sat = min(1.0d0, 1d0-(1d0-satup)*(1d0-z/zsat)**2d0)
    do iz=1,nz
        if (z(iz)>=zsat) sat(iz)=1d0
    enddo 
    #endif 

    ! getting user-defined SA

    call get_sa_num(nsld_sa)

    if (allocated(chrsld_sa)) deallocate(chrsld_sa)
    allocate(chrsld_sa(nsld_sa))

    call get_sa( &
        & nsp_sld,chrsld,p80,nsld_sa &! input
        & ,hrii,chrsld_sa &! output
        & )

    do isps = 1, nsp_sld
        hri(isps,:) = 1d0/hrii(isps)
    enddo

    ! getting user-defined PSD for dust

    call get_psdrain_num(nps_rain_char_in)

    if (nps_rain_char_in <= 0) then
        ! random default used in GMD paper
        nps_rain_char = 4
        
        if (allocated(pssigma_rain_list_in)) deallocate(pssigma_rain_list_in)
        if (allocated(psu_rain_list_in)) deallocate(psu_rain_list_in)
        if (allocated(psw_rain_list_in)) deallocate(psw_rain_list_in)
        allocate(pssigma_rain_list_in(nps_rain_char),psu_rain_list_in(nps_rain_char),psw_rain_list_in(nps_rain_char))
        
        psu_rain_list_in        = (/ log10(5d-6), log10(20d-6),  log10(50d-6), log10(70d-6) /)
        pssigma_rain_list_in(:) = 0.2d0
        psw_rain_list_in(:)     = 1d0
        
        
    else

        nps_rain_char = nps_rain_char_in
        
        if (allocated(pssigma_rain_list_in)) deallocate(pssigma_rain_list_in)
        if (allocated(psu_rain_list_in)) deallocate(psu_rain_list_in)
        if (allocated(psw_rain_list_in)) deallocate(psw_rain_list_in)
        allocate(pssigma_rain_list_in(nps_rain_char),psu_rain_list_in(nps_rain_char),psw_rain_list_in(nps_rain_char))

        call get_psdrain( &
            & nps_rain_char &! input
            & ,psu_rain_list_in,pssigma_rain_list_in,psw_rain_list_in &! output
            & )
            
        do ips = 1, nps_rain_char
            psu_rain_list_in(ips) = log10( psu_rain_list_in(ips) )
        enddo
    endif 
        
    if (allocated(pssigma_rain_list)) deallocate(pssigma_rain_list)
    if (allocated(psu_rain_list)) deallocate(psu_rain_list)
    if (allocated(psw_rain_list)) deallocate(psw_rain_list)
    allocate(pssigma_rain_list(nps_rain_char),psu_rain_list(nps_rain_char),psw_rain_list(nps_rain_char))

    ! print*,'printing psd inputs',nps_rain_char_in,nps_rain_char
    ! print*,psu_rain_list
    ! print*,pssigma_rain_list
    ! print*,psw_rain_list
    ! stop

    ! getting user-defined cec

    call get_cec_num(nsld_cec)

    if (allocated(chrsld_cec)) deallocate(chrsld_cec)
    allocate(chrsld_cec(nsld_cec))

    call get_cec( &
        & nsp_sld_all,chrsld_all,mcec_all_def,nsld_cec,nsp_aq_all,chraq_all,logkhaq_all_def,beta_all_def &! input
        & ,mcec_all,chrsld_cec,logkhaq_all,beta_all &! output
        & ) 

    call get_nopsd_num(nsld_nopsd)

    if (allocated(chrsld_nopsd)) deallocate(chrsld_nopsd)
    allocate(chrsld_nopsd(nsld_nopsd))

    call get_nopsd( &
        & nsp_sld,chrsld,nsld_nopsd &! input
        & ,chrsld_nopsd &! output
        & )


    rough = 1d0
    rough_ps = 1d0
    rough_ps_b = 1d0
    ! roughness factor parameterization references
    roughref_b = 'NSB07'      ! Navarre-Sitchler and Brantley (2007)
    ! roughref_b  = 'BM00'       ! Brantley and Mellott (2000)
    ! roughref_b  = 'Letal21'    ! Lewis et al. (2021) (assuming sphere)
    ! roughref_b  = 'smooth'     ! smooth surface 

    roughref = roughref_b 
    ! from Navarre-Sitchler and Brantley (2007)
    ! rough_c0 = 10d0**(3.3d0)
    ! rough_c1 = 0.33d0
    ! from Brantley and Mellott (2000)
    ! rough_c0 = 10d0**(0.7d0)
    ! rough_c1 = -0.1d0
    if (incld_rough) then 
        ! rough = 10d0**(3.3d0)*p80**0.33d0 ! from Navarre-Sitchler and Brantley (2007)
        do isps=1,nsp_sld
            rough(isps,:) = rough_f( roughref(isps), nz, 1d0/hri(isps,:) )
        enddo 
    endif 

    ! ssa_cmn = -4.4528d0*log10(p80*1d6) + 11.578d0 ! m2/g

    do isps=1,nsp_sld
        hr(isps,:) = hri(isps,:)*rough(isps,:)
    enddo 
    v = qin/poroi/sat
    poro = poroi
    torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
    tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)

    w_btm = w0
    w = w_btm
    ! if (noncnstw) then 
        ! w(:) = w0/(1d0- poro(:)) ! from w*(1- poro) = w0*(1 - poroi) --- isovolumetric weathering?
        ! w_btm = w0/(1d0- poroi) ! from w*(1- poro) = w0*(1 - poroi) --- isovolumetric weathering?
    ! endif 

    ! adding dispersion calculation
    if (disp_ON) then 
        zdisp = ztot ! assuming dispersion scale = total depth 
    else 
        zdisp = 0d0
        disp_FULL_ON = .false.
    endif
    ! parameterization by Schulz-Makuch 2005 
    ! for basalt 
    c0_disp = 0.15d0
    c1_disp = 0.61d0
    ! for granite
    c0_disp = 0.21d0
    c1_disp = 0.51d0
    c_disp = c0_disp*zdisp**c1_disp

    disp = c_disp * v

    if (disp_FULL_ON) disp = disp_FULL


    ! ------------ determine calculation scheme for advection (from IMP code)
    call calcupwindscheme(  &
        up,dwn,cnr,adf & ! output 
        ,w,nz   & ! input &
        )

    ! attempting to do psd 
    if (do_psd) then 
        do ips = 1, nps
            ps(ips) = log10(ps_min) + (ips - 1d0)*(log10(ps_max) - log10(ps_min))/(nps - 1d0)
        enddo 
        dps(:) = ps(2) - ps(1)
        print *,ps
        print *,dps
        
        if (do_psd_full) then ! do psd for every mienral
            open(ipsd,file = trim(adjustl(profdir))//'/'//'psd_pr.txt',status = 'replace')
            open(ipsdv,file = trim(adjustl(profdir))//'/'//'intpsd_pr.txt',status = 'replace')
            write(ipsd,*) ' sldsp\log10(radius) ', (ps(ips),ips=1,nps), 'time'
            write(ipsdv,*) ' sldsp\diameter(um) ', (10d0**ps(ips)*1d6*2d0,ips=1,nps), 'p80(um)'
            do isps = 1, nsp_sld
                volsld = msldi(isps)*mv(isps)*1d-6
                call calc_psd_pr( &
                    & nps &! input
                    & ,pi,hrii(isps),ps_sigma_std,poroi,volsld,tol &! input
                    & ,ps,dps &! input
                    & ,msldunit &! input
                    & ,psd_pr &! output 
                    & )
                mpsd_pr(isps,:) = psd_pr(:)
                intpsd_tmp(:) = psd_pr(:)*(10d0**ps(:))**3d0
                intpsd = intpsd_tmp
                do ips = 1, nps
                    intpsd(ips) = sum(intpsd_tmp(1:ips))/sum(intpsd_tmp)
                enddo
                call calc_p80( &
                    & nps,ps,intpsd &! input 
                    & ,p80_tmp &! output
                    & )
                write(ipsd,*) chrsld(isps),(psd_pr(ips),ips=1,nps), 0d0
                write(ipsdv,*) chrsld(isps),(intpsd(ips),ips=1,nps), p80_tmp
            enddo 
            close(ipsd)
            close(ipsdv)

            ! initially particle is distributed as in parent rock 
            do isps=1,nsp_sld
                do iz = 1, nz
                    mpsd(isps,:,iz) = mpsd_pr(isps,:) 
                enddo 
            
                if (.not.incld_rough) then 
                    rough_ps(isps,:) = rough_f( 'smooth', nps, 10d0**ps(:) )
                else 
                    rough_ps(isps,:) = rough_f( roughref(isps), nps, 10d0**ps(:) )
                endif 

                do iz=1,nz
                    ssa(isps,iz) = sum( 4d0*pi*(10d0**ps(:))**2d0*rough_ps(isps,:)*mpsd(isps,:,iz)*dps(:) )
                    ssav(isps,iz) = sum( 3d0/(10d0**ps(:))*rough_ps(isps,:)*mpsd(isps,:,iz)*dps(:) )
                    ssv(isps,iz) = sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*mpsd(isps,:,iz)*dps(:))
                enddo 
            enddo 
            ! hr = ssa *(1-poro)/poro ! converting m2/sld-m3 to m2/pore-m3
            ! hr = ssa 
            do iz=1,nz
                ! hr(:,iz) = ssa(:,iz)/poro(iz)/(msldi(:)*mv(:)*1d-6) ! so that poro * hr * mv * msld becomes porosity independent
                ! hr(:,iz) = ssa(:,iz) ! so that poro * hr * mv * msld becomes porosity independent
                hr(:,iz) = ssa(:,iz)/ssv(:,iz)/poro(iz)  
                ! hr(:,iz) = ssav(:,iz)/poro(iz)  
            enddo 
        else ! do psd only for bulk 
            volsld = sum(msldi*mv*1d-6) + mblki*mvblk*1d-6
            call calc_psd_pr( &
                & nps &! input
                & ,pi,p80,ps_sigma_std,poroi,volsld,tol &! input
                & ,ps,dps &! input
                & ,msldunit &! input
                & ,psd_pr &! output 
                & )
            intpsd_tmp(:) = psd_pr(:)*(10d0**ps(:))**3d0
            intpsd = intpsd_tmp
            do ips = 1, nps
                intpsd(ips) = sum(intpsd_tmp(1:ips))/sum(intpsd_tmp)
            enddo
            call calc_p80( &
                & nps,ps,intpsd &! input 
                & ,p80_tmp &! output
                & )
            open(ipsd,file = trim(adjustl(profdir))//'/'//'psd_pr.txt',status = 'replace')
            write(ipsd,*) ' depth\log10(radius) ', (ps(ips),ips=1,nps), 'time'
            write(ipsd,*) ztot,(psd_pr(ips),ips=1,nps), 0d0
            close(ipsd)
            
            open(ipsdv,file = trim(adjustl(profdir))//'/'//'intpsd_pr.txt',status = 'replace')
            write(ipsdv,*) ' sldsp\diameter(um) ', (10d0**ps(ips)*1d6*2d0,ips=1,nps), 'p80(um)'
            write(ipsdv,*) chrsld(isps),(intpsd(ips),ips=1,nps), p80_tmp
            close(ipsdv)

            ! initially particle is distributed as in parent rock 
            do iz = 1, nz
                psd(:,iz) = psd_pr(:) 
            enddo 

            ! dM = M * [psd*dps*S(r)] * k *dt 
            ! so hr = sum (psd(:)*dps(:)*S(:) ) where S in units m2/m3 and simplest way 1/r 
            ! in this case hr = sum(  psd(:)*dps(:)*1d0/(10d0**(-ps(:))) )
            if (.not.incld_rough) then 
                rough_ps_b(:) = rough_f( roughref_b, nps, 10d0**ps(:) )
            else 
                rough_ps_b(:) = rough_f( roughref_b, nps, 10d0**ps(:) )
            endif 
            
            do iz=1,nz
                ssa(:,iz) = sum( 4d0*pi*(10d0**ps(:))**2d0 *rough_ps_b(:)*psd(:,iz)*dps(:))
                ssav(:,iz) = sum( 3d0/(10d0**ps(:)) *rough_ps_b(:)*psd(:,iz)*dps(:))
                ssv(:,iz) = sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*psd(:,iz)*dps(:))
            enddo 
            ! hr = ssa *(1-poro)/poro ! converting m2/sld-m3 to m2/pore-m3
            ! hr = ssa 
            do iz=1,nz
                hr(:,iz) = ssa(:,iz)/poro(iz) ! so that poro * hr * mv * msld becomes porosity independent
            enddo 
        endif 
    else 
        psd = 0d0
        mpsd = 0d0
    endif 

    minsld = 1d-20
    ! when minimum psd is limited, msld min is limited
    if (do_psd .and. do_psd_full .and. psd_lim_min .and. psd_vol_consv) then 
        lim_minsld_in = .true.
        do isps=1,nsp_sld
            minsld(isps) = sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*psd_th_0*dps(:)) &! m3/m3
                & /( mv(isps)*1d-6 ) !  m3/m3 / (m3 /mol ) =  mol/m3
        enddo 
    endif 


    minmaqads = 1d-20
    ! calculating threshold for psd (decided not to use)
    if (do_psd) then 
        if (do_psd_full) then
            do ips=1,nps
                do isps=1,nsp_sld
                    mpsd_th(isps,ips) = 1d-9 &! mol m-3 
                        & / ( 4d0/3d0*pi*(10d0**ps(ips))**3d0  &! m3/m3
                        & /( mv(isps)*1d-6 ) ) !  m3/m3 / (m3 /mol ) =  mol/m3
                    mpsd_th(isps,ips) = mpsd_th(isps,ips)/dps(ips) ! number/m3/log(m)
                enddo
            enddo
            ! print *,'mpsd_th'
            ! print *,mpsd_th
            ! stop
        else
            psd_th = psd_th_0
        endif 
    endif 

    ! #ifdef surfssa
    ! hri = ssa_cmn*1d6/poro
    ! mvab_save = mvab
    ! mvan_save = mvan
    ! mvcc_save = mvcc
    ! mvfo_save = mvfo
    ! mvka_save = mvka
    ! mvab = mwtab 
    ! mvan = mwtan 
    ! mvcc = mwtcc 
    ! mvfo = mwtfo 
    ! mvka = mwtka 
    ! #endif 

    dt = maxdt

    dt = 1d-20 ! for basalt exp?

    do ispa = 1, nsp_aq
        maq(ispa,:)=maqi(ispa)
    enddo 
    do ispg = 1, nsp_gas
        mgas(ispg,:)=mgasi(ispg)
    enddo 
    do isps = 1, nsp_sld
        msld(isps,:) = msldi(isps)
    enddo 

    mblk = mblki

    ! initial solid conc. modified 
    ! do iz = 1, nz
        ! msld(:,iz) = msld(:,iz)*exp( real(iz - nz) )
        ! poro(iz) = 1d0 - sum(msld(:,iz)*mv(:)*1d-6)
    ! enddo 
        
    #ifdef ksld_chk
    open (idust, file='./ksld_chk.txt', status ='unknown',action='write')
    write(chrfmt,'(i0)') nsp_sld_all
    chrfmt = '(a12,'//trim(adjustl(chrfmt))//'(1x,a5))'
    write(idust,chrfmt) 'pH\sldsp',chrsld_all
    write(chrfmt,'(i0)') nsp_sld_all
    chrfmt = '(1x,f5.2,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
    do iph = 1,nph
        pro = 10d0**(0d0 + (iph-1d0)/(nph-1d0)*(-14d0))
        
        call coefs_v2( &
            & nz,rg,rg2,25d0,sec2yr,tempk_0,pro,cec_pH_depend,mcec_all,logkhaq_all &! input
            & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
            & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
            & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,mwt_all,staq_all &!input
            & ,nsp_aq,nsp_aq_cnst,chraq,chraq_cnst,maq,maqc &!input
            & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 &! output
            & ,keqaq_oxa,keqaq_cl &! output
            & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
            & ,keqcec_all,keqiex_all &! output 
            & ) 
        ! write(idust,chrfmt) -log10(pro(1)),log10(ksld_all(:,1))
        write(idust,chrfmt) -log10(pro(1)),ksld_all(:,1)/sec2yr
    enddo 
    close(idust)

    open (idust, file='./sld_data_chk.txt', status ='unknown',action='write')
    chrfmt = '(3(1x,a5))'
    write(idust,chrfmt) 'sld','mv','mwt'
    write(chrfmt,'(i0)') nsp_sld_all
    chrfmt = '(1x,a5,2(1x,E11.3))'
    do isps = 1,nsp_sld_all
        write(idust,chrfmt) chrsld_all(isps),mv_all(isps),mwt_all(isps)
    enddo 
    close(idust)
    ! stop
    #endif 

    omega = 0d0

    pro = 1d-5

    if (allocated(kin_sld_spc)) deallocate(kin_sld_spc)
    if (allocated(chrsld_kinspc)) deallocate(chrsld_kinspc)
    nsld_kinspc = nsld_kinspc_in
    allocate(chrsld_kinspc(nsld_kinspc),kin_sld_spc(nsld_kinspc))
    chrsld_kinspc = chrsld_kinspc_in
    kin_sld_spc = kin_sld_spc_in
        
    call coefs_v2( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pro,cec_pH_depend,mcec_all,logkhaq_all &! input
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
        & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
        & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,mwt_all,staq_all &!input
        & ,nsp_aq,nsp_aq_cnst,chraq,chraq_cnst,maq,maqc &!input
        & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 &! output
        & ,keqaq_oxa,keqaq_cl &! output
        & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
        & ,keqcec_all,keqiex_all &! output 
        & ) 

    print_cb = .false. 
    print_loc = './ph.txt'

    pro = 1d0
    ios = 0d0
    if (act_ON) ios = 1d-12
    call calc_pH_v7_4( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,poro,sat,tc &! input 
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maq,maqc,mgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 &! input
        & ,keqaq_oxa,keqaq_cl &! input
        & ,print_cb,print_loc,z,act_ON &! input 
        & ,dprodmaq_all,dprodmgas_all &! output
        & ,ios,diosdmaq_all,diosdmgas_all &! output
        & ,pro,ph_error,ph_iter &! output
        & ) 
    ! print*,ios
    ! stop
    ! getting mgasx_loc & maqx_loc
    call get_maqgasx_all( &
        & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
        & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
        & ,maq,mgas,maqc,mgasc &
        & ,maqx_loc,mgasx_loc  &! output
        & )

    ! getting maqft_loc and its derivatives
    call get_maqt_all( &
    ! call get_maqt_all_v2( &
        & nz,nsp_aq_all,nsp_gas_all &
        & ,chraq_all,chrgas_all &
        & ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
        & ,mgasx_loc,maqx_loc,pro,ios,tc &
        & ,dmaqft_dpro_loc,dmaqft_dmaqf_loc,dmaqft_dmgas_loc,dmaqft_dios_loc &! output
        & ,maqft_loc  &! output
        & )
    maqft = 0d0
    do ispa=1,nsp_aq
        maqft(ispa,:)=maqft_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
    enddo 

    !!!  for adsorption 
    call get_msldx_all( &
        & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
        & ,chrsld,chrsld_all,chrsld_cnst &
        & ,msld,msldc &
        & ,msldx_loc  &! output
        & )

    ! call get_maqads_all_v3( &
    call get_maqads_all_v4( &
    ! call get_maqads_all_v4a( &
        & nz,nsp_aq_all,nsp_sld_all &
        & ,chraq_all,chrsld_all &
        & ,keqcec_all,keqiex_all,cec_pH_depend,beta_all &
        & ,msldx_loc,maqx_loc,pro &
        & ,dmaqfads_sld_dpro,dmaqfads_sld_dmaqf,dmaqfads_sld_dmsld &! output
        & ,msldf_loc,maqfads_sld_loc,beta_loc,ads_error  &! output
        & )
    maqfads_sld = 0d0
    do ispa=1,nsp_aq
        do isps=1,nsp_sld
            maqfads_sld(ispa,isps,:) &
                & =maqfads_sld_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrsld_all,chrsld(isps),dim=1),:)
        enddo
    enddo 

    maqfads = 0d0
    do ispa=1,nsp_aq
        ! maqfads(ispa,:)=maqfads_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
        do iz=1,nz
            maqfads(ispa,iz) = sum(maqfads_sld(ispa,:,iz))
        enddo
    enddo 


    ! so4fprev = so4f
    maqft_prev = maqft
    maqfads_prev = maqfads
    proi = pro(1)
    print*,proi
    ! pause

    poroprev = poro

    !  --------- read -----
    if (read_data) then 
        ! runname_save = 'test_cpl_rain-0.40E+04_pevol_sevol1_q-0.10E-01_zsat-5' ! specifiy the file where restart data is stored 
        ! runname_save = runname  ! the working folder has the restart data 
        loc_runname_save = '../'//trim(adjustl(runname_save))//'/'//trim(adjustl(profdir(3:)))
        if (trim(adjustl(runname_save)) == 'self') loc_runname_save = trim(adjustl(profdir))
        call system('cp '//trim(adjustl(loc_runname_save))//'/'//'prof_sld-save.txt '  &
            & //trim(adjustl(profdir))//'/'//'prof_sld-restart.txt')
        call system('cp '//trim(adjustl(loc_runname_save))//'/'//'prof_aq-save.txt '  &
            & //trim(adjustl(profdir))//'/'//'prof_aq-restart.txt')
        call system('cp '//trim(adjustl(loc_runname_save))//'/'//'prof_gas-save.txt '  &
            & //trim(adjustl(profdir))//'/'//'prof_gas-restart.txt')
        call system('cp '//trim(adjustl(loc_runname_save))//'/'//'bsd-save.txt '  &
            & //trim(adjustl(profdir))//'/'//'bsd-restart.txt')
        call system('cp '//trim(adjustl(loc_runname_save))//'/'//'psd-save.txt '  &
            & //trim(adjustl(profdir))//'/'//'psd-restart.txt')
        call system('cp '//trim(adjustl(loc_runname_save))//'/'//'sa-save.txt '  &
            & //trim(adjustl(profdir))//'/'//'sa-restart.txt')
            
        call get_saved_variables_num( &
            & workdir,loc_runname_save &! input
            & ,nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save,nsld_kinspc_save,nsld_sa_save &! output
            & )
        
        allocate(chraq_save(nsp_aq_save),chrsld_save(nsp_sld_save),chrgas_save(nsp_gas_save),chrrxn_ext_save(nrxn_ext_save))
        allocate(maq_save(nsp_aq_save,nz),msld_save(nsp_sld_save,nz),mgas_save(nsp_gas_save,nz))
        allocate(chrsld_kinspc_save(nsld_kinspc_save),kin_sldspc_save(nsld_kinspc_save))
        allocate(chrsld_sa_save(nsld_sa_save),hrii_save(nsld_sa_save),hr_save(nsp_sld_save,nz),mpsd_save(nsp_sld_save,nps,nz))
            
        
        call get_saved_variables( &
            & workdir,loc_runname_save &! input
            & ,nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save,nsld_kinspc_save,nsld_sa_save &! input
            & ,chraq_save,chrgas_save,chrsld_save,chrrxn_ext_save,chrsld_kinspc_save,kin_sldspc_save &! output
            & ,chrsld_sa_save,hrii_save &! output 
            & )
        
            
        open (isldprof, file=trim(adjustl(profdir))//'/'//'prof_sld-restart.txt',  &
            & status ='old',action='read')
        open (iaqprof, file=trim(adjustl(profdir))//'/'//'prof_aq-restart.txt',  &
            & status ='old',action='read')
        open (igasprof, file=trim(adjustl(profdir))//'/'//'prof_gas-restart.txt',  &
            & status ='old',action='read')
        open (ibsd, file=trim(adjustl(profdir))//'/'//'bsd-restart.txt',  &
            & status ='old',action='read')
        open (ipsd, file=trim(adjustl(profdir))//'/'//'psd-restart.txt',  &
            & status ='old',action='read')
        open (isa, file=trim(adjustl(profdir))//'/'//'sa-restart.txt',  &
            & status ='old',action='read')
        
        read (isldprof,'()')
        read (iaqprof,'()')
        read (igasprof,'()')
        read (ibsd,'()')
        read (ipsd,'()')
        read (isa,'()')
        
        do iz = 1, Nz
            ucvsld1 = 1d0
            if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
            read (isldprof,*) z(iz),(msld_save(isps,iz),isps=1,nsp_sld_save),time
            read (iaqprof,*) z(iz),(maq_save(ispa,iz),ispa=1,nsp_aq_save),pro(iz),time
            read (igasprof,*) z(iz),(mgas_save(ispg,iz),ispg=1,nsp_gas_save),time
            read (ibsd,*) z(iz),poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz),mblk(iz),cec(iz),time
            read (ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
            read (isa,*) z(iz), (hr_save(isps,iz),isps=1,nsp_sld_save), time 
            mblk(iz) = mblk(iz)/ ( mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6) )
        enddo 
        close(isldprof)
        close(iaqprof)
        close(igasprof)
        close(ibsd)
        close(ipsd)
        close(isa)
        
        if (all(psd==0d0)) then 
            no_psd_prevrun = .true.
        else 
            no_psd_prevrun = .false.
        endif 
        
        pro = 10d0**(-pro) ! read data is -log10 (pro)
        
        torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
        tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
            
        c_disp = c0_disp*zdisp**c1_disp
        disp = c_disp * v

        if (disp_FULL_ON) disp = disp_FULL
        
        do isps = 1,nsp_sld_save
            if (any(chrsld == chrsld_save(isps))) then 
                msld(findloc(chrsld,chrsld_save(isps),dim=1),:) = msld_save(isps,:)
            elseif (any(chrsld_cnst == chrsld_save(isps))) then
                msldc(findloc(chrsld_cnst,chrsld_save(isps),dim=1),:) = msld_save(isps,:)
            else 
                print *,'error in re-assignment of sld conc.'
            endif 
        enddo 
        
        do ispa = 1,nsp_aq_save
            if (any(chraq == chraq_save(ispa))) then 
                maq(findloc(chraq,chraq_save(ispa),dim=1),:) = maq_save(ispa,:)
            elseif (any(chraq_cnst == chraq_save(ispa))) then
                maqc(findloc(chraq_cnst,chraq_save(ispa),dim=1),:) = maq_save(ispa,:)
            else 
                print *,'error in re-assignment of aq conc.'
            endif 
        enddo 
        
        do ispg = 1,nsp_gas_save
            if (any(chrgas == chrgas_save(ispg))) then 
                mgas(findloc(chrgas,chrgas_save(ispg),dim=1),:) = mgas_save(ispg,:)
            elseif (any(chrgas_cnst == chrgas_save(ispg))) then
                mgasc(findloc(chrgas_cnst,chrgas_save(ispg),dim=1),:) = mgas_save(ispg,:)
            else 
                print *,'error in re-assignment of gas conc.'
            endif 
        enddo 
        
        ! counting sld species whose values are to be specificed in kinspc.save and not so yet when reading from kinspc.in
        nsld_kinspc_add = 0
        do isps_kinspc = 1,nsld_kinspc_save
            if (any(chrsld_kinspc_in == chrsld_kinspc_save(isps_kinspc))) then ! already specified 
                continue
            else 
                nsld_kinspc_add = nsld_kinspc_add + 1
            endif 
        enddo 
        
        if (nsld_kinspc_add > 0) then 
            ! deallocate 
            if (allocated(kin_sld_spc)) deallocate (kin_sld_spc)
            if (allocated(chrsld_kinspc)) deallocate (chrsld_kinspc)
            ! re-define sld species number whose rate const. is specified 
            nsld_kinspc = nsld_kinspc + nsld_kinspc_add
            ! allocate 
            allocate(kin_sld_spc(nsld_kinspc),chrsld_kinspc(nsld_kinspc))
            ! saving already specified consts. 
            chrsld_kinspc(1:nsld_kinspc_in) = chrsld_kinspc_in
            kin_sld_spc(1:nsld_kinspc_in) = kin_sld_spc_in
            ! adding previously specified rate const. 
            nsld_kinspc_add = 0
            do isps_kinspc = 1,nsld_kinspc_save
                if (any(chrsld_kinspc_in == chrsld_kinspc_save(isps_kinspc))) then 
                    continue
                else 
                    nsld_kinspc_add = nsld_kinspc_add + 1
                    chrsld_kinspc(nsld_kinspc_in + nsld_kinspc_add) = chrsld_kinspc_save(isps_kinspc)
                    kin_sld_spc(nsld_kinspc_in + nsld_kinspc_add) = kin_sldspc_save(isps_kinspc)
                endif 
            enddo 
        endif 
        
        ! overloading sa if saved 
        
        do isps = 1,nsp_sld_save
            if (any(chrsld == chrsld_save(isps))) then 
                hr(findloc(chrsld,chrsld_save(isps),dim=1),:) = hr_save(isps,:)
            endif 
        enddo     
        
        if (nsld_sa_save > 0) then 
            do isps_sa =1, nsld_sa_save 
                if (any(chrsld_sa == chrsld_sa_save(isps_sa))) then ! if SA of a sld sp. is already specified, save data does not overload 
                    continue 
                else  ! if SA is not specified some species but was specified in the previous run, saved data is loaded 
                    if (any(chrsld == chrsld_sa_save(isps_sa))) then
                        hrii(findloc(chrsld,chrsld_sa_save(isps_sa),dim=1)) = hrii_save(isps_sa)
                    endif 
                endif 
            enddo 
        endif 
        ! updating SA parameters 
        do isps = 1, nsp_sld
            hri(isps,:) = 1d0/hrii(isps)
        enddo
        
        ! calculating roughness based on newly defined reference (not reflecting old reference)
        rough = 1d0 
        if (incld_rough) then 
            do isps=1,nsp_sld
                rough(isps,:) = rough_f( roughref(isps), nz, 1d0/hri(isps,:) ) 
            enddo 
        endif 
        
        if (do_psd_full .and. .not.no_psd_prevrun) then
            ! updating parentrock psd if hrii has been loaded from a previous run 
            if (nsld_sa_save > 0) then 
                open(ipsd,file = trim(adjustl(profdir))//'/'//'psd_pr.txt',status = 'replace')
                open(ipsdv,file = trim(adjustl(profdir))//'/'//'intpsd_pr.txt',status = 'replace')
                write(ipsd,*) ' sldsp\log10(radius) ', (ps(ips),ips=1,nps), 'time'
                write(ipsdv,*) ' sldsp\diameter(um) ', (10d0**ps(ips)*1d6*2d0,ips=1,nps), 'p80(um)'
                do isps = 1, nsp_sld
                    volsld = msldi(isps)*mv(isps)*1d-6
                    call calc_psd_pr( &
                        & nps &! input
                        & ,pi,hrii(isps),ps_sigma_std,poroi,volsld,tol &! input
                        & ,ps,dps &! input
                        & ,msldunit &! input
                        & ,psd_pr &! output 
                        & )
                    mpsd_pr(isps,:) = psd_pr(:)
                    intpsd_tmp(:) = psd_pr(:)*(10d0**ps(:))**3d0
                    intpsd = intpsd_tmp
                    do ips = 1, nps
                        intpsd(ips) = sum(intpsd_tmp(1:ips))/sum(intpsd_tmp)
                    enddo
                    call calc_p80( &
                        & nps,ps,intpsd &! input 
                        & ,p80_tmp &! output
                        & )
                    write(ipsd,*) chrsld(isps),(psd_pr(ips),ips=1,nps), 0d0
                    write(ipsdv,*) chrsld(isps),(intpsd(ips),ips=1,nps), p80_tmp
                enddo 
                close(ipsd)
                close(ipsdv)
                do isps=1,nsp_sld
                    do iz = 1, nz
                        mpsd(isps,:,iz) = mpsd_pr(isps,:) 
                    enddo 
                enddo 
            endif 
            
            do isps = 1, nsp_sld_save ! loading psds 
            
                call system('cp '//trim(adjustl(loc_runname_save))//'/' &
                    & //'psd_'//trim(adjustl(chrsld_save(isps)))//'-save.txt '  &
                    & //trim(adjustl(profdir))//'/'//'psd_'//trim(adjustl(chrsld_save(isps)))//'-restart.txt')
                open (ipsd, file=trim(adjustl(profdir))//'/' &
                    & //'psd_'//trim(adjustl(chrsld_save(isps)))//'-restart.txt',  &
                    & status ='old',action='read')
                read (ipsd,'()')    
                do iz = 1, Nz
                    read (ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
                enddo 
                close(ipsd)
                mpsd_save(isps,:,:) = psd(:,:)
                if (any(chrsld == chrsld_save(isps))) then 
                    mpsd(findloc(chrsld,chrsld_save(isps),dim=1),:,:) = mpsd_save(isps,:,:)
                endif 
            enddo
            
            do isps=1,nsp_sld ! SA properties calc with updated PSDs
                do iz=1,nz
                    ssa(isps,iz) = sum( 4d0*pi*(10d0**ps(:))**2d0*rough_ps(isps,:)*mpsd(isps,:,iz)*dps(:) )
                    ssav(isps,iz) = sum( 3d0/(10d0**ps(:))*rough_ps(isps,:)*mpsd(isps,:,iz)*dps(:) )
                    ssv(isps,iz) = sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*mpsd(isps,:,iz)*dps(:) )
                enddo 
            enddo 
            
        endif 
        
        
        ! just to obtain so4f 
        print_cb = .false. 
        print_loc = './ph.txt'
        
        prox = pro
        iosx = ios
        call calc_pH_v7_4( &
            & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            & ,poro,sat,tc &! input 
            & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
            & ,maq,maqc,mgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 &! input
            & ,keqaq_oxa,keqaq_cl &! input
            & ,print_cb,print_loc,z,act_ON &! input 
            & ,dprodmaq_all,dprodmgas_all &! output
            & ,iosx,diosdmaq_all,diosdmgas_all &! output
            & ,prox,ph_error,ph_iter &! output
            & ) 
        ios = iosx
        ! stop
        ! getting mgasx_loc & maqx_loc
        call get_maqgasx_all( &
            & nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst &
            & ,maq,mgas,maqc,mgasc &
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
        maqft = 0d0
        do ispa=1,nsp_aq
            maqft(ispa,:)=maqft_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
        enddo 
        ! so4fprev = so4f
        maqft_prev = maqft
        
        !!!  for adsorption 
        call get_msldx_all( &
            & nz,nsp_sld_all,nsp_sld,nsp_sld_cnst &
            & ,chrsld,chrsld_all,chrsld_cnst &
            & ,msld,msldc &
            & ,msldx_loc  &! output
            & )

        ! call get_maqads_all_v3( &
        call get_maqads_all_v4( &
        ! call get_maqads_all_v4a( &
            & nz,nsp_aq_all,nsp_sld_all &
            & ,chraq_all,chrsld_all &
            & ,keqcec_all,keqiex_all,cec_pH_depend,beta_all &
            & ,msldx_loc,maqx_loc,prox &
            & ,dmaqfads_sld_dpro,dmaqfads_sld_dmaqf,dmaqfads_sld_dmsld &! output
            & ,msldf_loc,maqfads_sld_loc,beta_loc,ads_error  &! output
            & )
        maqfads_sld = 0d0
        do ispa=1,nsp_aq
            do isps=1,nsp_sld
                maqfads_sld(ispa,isps,:) &
                    & = maqfads_sld_loc(findloc(chraq_all,chraq(ispa),dim=1),findloc(chrsld_all,chrsld(isps),dim=1),:)
            enddo
        enddo 

        maqfads = 0d0
        do ispa=1,nsp_aq
            ! maqfads(ispa,:)=maqfads_loc(findloc(chraq_all,chraq(ispa),dim=1),:)
            do iz=1,nz
                maqfads(ispa,iz) = sum(maqfads_sld(ispa,:,iz))
            enddo
        enddo 
        
        maqfads_prev = maqfads
        
        time = 0d0
            
        if (display) then
            write(chrfmt,'(i0)') nz_disp
            chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
            
            print *
            print *,' [concs] '
            print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
            if (nsp_aq>0) then 
                print *,' < aq species >'
                do ispa = 1, nsp_aq
                    print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maq(ispa,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
            if (nsp_sld>0) then 
                print *,' < sld species >'
                do isps = 1, nsp_sld
                    print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (msld(isps,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
            if (nsp_gas>0) then 
                print *,' < gas species >'
                do ispg = 1, nsp_gas
                    print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (mgas(ispg,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
        endif      
    endif
        
    call coefs_v2( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pro,cec_pH_depend,mcec_all,logkhaq_all &! input
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
        & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
        & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,mwt_all,staq_all &!input
        & ,nsp_aq,nsp_aq_cnst,chraq,chraq_cnst,maq,maqc &!input
        & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 &! output
        & ,keqaq_oxa,keqaq_cl &! output
        & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
        & ,keqcec_all,keqiex_all &! output 
        & ) 
        
    dbl_ref = 0d0   ! diffuse boundary at the top

    save_trans = .true.
    save_trans = .false.
    zml = zml_background
    call make_transmx(  &
        & nsp_sld,imix,dz,poro,nz,z,zml,dbl_ref,tol,save_trans  &! input
        & ,trans  &! output 
        & )
        
    ! --------- loop -----
    print *, 'about to start time loop'
    it = 0
    irec_prof = 0
    irec_flx = 0

    ict = 0
    ict_prev = ict
    ict_change = .false.

    count_dtunchanged = 0

    int_flx_aq = 0d0
    int_flx_gas = 0d0
    int_flx_sld = 0d0
    int_flx_co2sp = 0d0

    int_ph = 0d0

    !! @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@

    do while (it<nt)
        ! call cpu_time(time_start)
        call system_clock(t1)
        
        if (display) then 
            print *
            print *, '-----------------------------------------'
            print '(i11,a)', it,': time iteration'
            print '(E11.3,a)',time,': time [yr]' 
            print *
        endif
        dt_prev = dt
        ! only relevant when doing non-continuous dusting    
        if (dust_step) dust_norm_prev = dust_norm
        
        if (time>rectime_prof(nrec_prof)) exit
        
        if (it == 0) then 
            maxdt = 0.2d0
        endif 
        
        if (timestep_fixed) then 
            maxdt = 0.2d0
            ! maxdt = 0.02d0 ! when calcite is included smaller time step must be assumed 
            ! maxdt = 0.005d0 ! when calcite is included smaller time step must be assumed 
            ! maxdt = 0.002d0 ! working with p80 = 10 um
            ! maxdt = 0.001d0 ! when calcite is included smaller time step must be assumed 
            ! maxdt = 0.0005d0 ! working with p80 = 1 um
            
            ! if (time<1d-2) then  
                ! maxdt = 1d-6 
            ! elseif (time>=1d-2 .and. time<1d-1) then 
            ! if (time<1d-3) then  
                ! maxdt = 1d-7 
            ! elseif (time>=1d-3 .and. time<1d-2) then  
                ! maxdt = 1d-6 
            ! elseif (time>=1d-2 .and. time<1d-1) then  
                ! maxdt = 1d-5 
            ! elseif (time>=1d-1 .and. time<1d0) then  
            if ( time<1d0) then  
                maxdt = 1d-4 
            elseif (time>=1d0 .and. time<1d1) then 
                maxdt = 1d-3 
            elseif (time>=1d1 .and. time<1d2) then 
                maxdt = 1d-2  
            ! elseif (time>=1d2 .and. time<1d3) then 
                ! maxdt = 1d-1 
            ! elseif (time>=1d3 .and. time<1d4) then 
                ! maxdt = 1d0 
            ! elseif (time>=1d4 .and. time<1d5) then 
                ! maxdt = 1d1 
            ! elseif (time>=1d5 ) then 
                ! maxdt = 1d2 
            endif 
            
            ! maxdt = maxdt * 1d-1
        endif 
        
        ! count_dtunchanged_Max = 1000
        ! count_dtunchanged_Max = 10
        ! if (sld_enforce) count_dtunchanged_Max = 10
        ! if (dt<1d-5) then 
            ! count_dtunchanged_Max = 10
        ! elseif (dt>=1d-5 .and. dt<1d0) then
            ! count_dtunchanged_Max = 100
        ! elseif (dt>=1d0 ) then 
            ! count_dtunchanged_Max = 1000
        ! endif 
        count_dtunchanged_Max_loc = count_dtunchanged_Max
        if (dt < 1d-6 ) count_dtunchanged_Max_loc = 10d0

        ! -------- modifying dt --------

        !        if ((iter <= 10).and.(dt<1d1)) then
        if (dt<maxdt) then
            ! dt = dt*1.01d0
            dt = dt*10d0
            if (dt>maxdt) dt = maxdt
        endif
        ! if (iter > 300) then
            ! dt = dt/10d0
        ! end if
        
        ! added so that saving time become more consistent
        if ( time+dt > rectime_prof(irec_prof+1) .and. time+dt > rectime_flx(irec_flx+1) ) then 
            dt = min( &
                & rectime_prof(irec_prof+1) - time + tol_step_tau &
                & ,rectime_flx(irec_flx+1) - time + tol_step_tau &
                & )
        elseif ( time+dt > rectime_prof(irec_prof+1) .and. time+dt <= rectime_flx(irec_flx+1)) then 
            dt = rectime_prof(irec_prof+1) - time + tol_step_tau
        elseif ( time+dt <= rectime_prof(irec_prof+1) .and. time+dt > rectime_flx(irec_flx+1)) then 
            dt = rectime_flx(irec_flx+1) - time + tol_step_tau
        endif 

        ! incase temperature&ph change
        
        ! if climate is changing in the model 
        if (any(climate)) then 
            ict_change = .false.
            do iclim = 1,3
                if (climate(iclim)) then
                    select case(iclim)
                        case(1)
                            if (dt > dct(iclim)/10d0) dt = dct(iclim)/10d0
                            do ict = 1, nclim(iclim)
                                if (ict /= nclim(iclim)) then 
                                    dct(iclim) = clim_T(1,ict+1) - clim_T(1,ict)
                                elseif (ict == nclim(iclim)) then 
                                    dct(iclim) = ctau(iclim) - clim_T(1,ict)
                                endif 
                                ! print *, clim_T(1,ict),mod(time,ctau(iclim)),clim_T(1,ict) + dct(iclim)
                                if ( &
                                    & clim_T(1,ict) <= mod(time,ctau(iclim)) & 
                                    & .and. clim_T(1,ict) + dct(iclim) >= mod(time,ctau(iclim)) &
                                    & ) then 
                                    ! if (  &
                                        ! & mod(time,ctau(iclim)) + dt - clim_T(1,ict) + dct(iclim) &
                                        ! & > ctau(iclim) * tol_step_tau &
                                        ! & ) then 
                                        ! dt = clim_T(1,ict) + dct(iclim) - mod(time,ctau(iclim))
                                    ! endif 
                                    ! print *, ict
                                    if (ict /= ict_prev(iclim)) ict_change(iclim) = .true.
                                    ict_prev(iclim) = ict
                                    exit 
                                endif 
                            enddo 
                            if (ict /= nclim(iclim)) then
                                tc = ( clim_T(2,ict+1) - clim_T(2,ict) ) /( clim_T(1,ict+1) - clim_T(1,ict) ) &
                                    & * ( mod(time,ctau(iclim)) - clim_T(1,ict) ) + clim_T(2,ict)
                            elseif (ict == nclim(iclim)) then 
                                tc = ( clim_T(2,1) - clim_T(2,ict) ) /( dct(iclim)  ) &
                                    & * ( mod(time,ctau(iclim)) - clim_T(1,ict) ) + clim_T(2,ict)
                            endif 
                            
                        case(2)
                            if (dt > dct(iclim)/10d0) dt = dct(iclim)/10d0
                            do ict = 1, nclim(iclim)
                                if (ict /= nclim(iclim)) then 
                                    dct(iclim) = clim_q(1,ict+1) - clim_q(1,ict)
                                elseif (ict == nclim(iclim)) then 
                                    dct(iclim) = ctau(iclim) - clim_q(1,ict)
                                endif 
                                ! print *, clim_q(1,ict),mod(time,ctau(iclim)),clim_q(1,ict) + dct(iclim)
                                if ( &
                                    & clim_q(1,ict) <= mod(time,ctau(iclim)) & 
                                    & .and. clim_q(1,ict) + dct(iclim) >= mod(time,ctau(iclim)) &
                                    & ) then 
                                    ! if (  &
                                        ! & mod(time,ctau(iclim)) + dt - clim_q(1,ict) + dct(iclim) &
                                        ! & > ctau(iclim) * tol_step_tau &
                                        ! & ) then 
                                        ! dt = clim_q(1,ict) + dct(iclim) - mod(time,ctau(iclim))
                                    ! endif 
                                    ! print *, ict
                                    if (ict /= ict_prev(iclim)) ict_change(iclim) = .true.
                                    ict_prev(iclim) = ict
                                    exit 
                                endif 
                            enddo 
                            if (ict /= nclim(iclim)) then 
                                qin = ( clim_q(2,ict+1) - clim_q(2,ict) ) /( clim_q(1,ict+1) - clim_q(1,ict) ) &
                                    & * ( mod(time,ctau(iclim)) - clim_q(1,ict) ) + clim_q(2,ict)
                            elseif (ict == nclim(iclim)) then 
                                qin = ( clim_q(2,1) - clim_q(2,ict) ) /( dct(iclim)  ) &
                                    & * ( mod(time,ctau(iclim)) - clim_q(1,ict) ) + clim_q(2,ict)
                            endif 
                            
                        case(3)
                            if (dt > dct(iclim)/10d0) dt = dct(iclim)/10d0
                            do ict = 1, nclim(iclim)
                                if (ict /= nclim(iclim)) then 
                                    dct(iclim) = clim_sat(1,ict+1) - clim_sat(1,ict)
                                elseif (ict == nclim(iclim)) then 
                                    dct(iclim) = ctau(iclim) - clim_sat(1,ict)
                                endif 
                                ! print *, clim_sat(1,ict),mod(time,ctau(iclim)),clim_sat(1,ict) + dct(iclim)
                                if ( &
                                    & clim_sat(1,ict) <= mod(time,ctau(iclim)) & 
                                    & .and. clim_sat(1,ict) + dct(iclim) >= mod(time,ctau(iclim)) &
                                    & ) then 
                                    ! if (  &
                                        ! & mod(time,ctau(iclim)) + dt - clim_sat(1,ict) + dct(iclim) &
                                        ! & > ctau(iclim) * tol_step_tau &
                                        ! & ) then 
                                        ! dt = clim_sat(1,ict) + dct(iclim) - mod(time,ctau(iclim))
                                    ! endif 
                                    ! print *, ict
                                    if (ict /= ict_prev(iclim)) ict_change(iclim) = .true.
                                    ict_prev(iclim) = ict
                                    exit 
                                endif 
                            enddo 
                            if (ict /= nclim(iclim)) then 
                                satup = ( clim_sat(2,ict+1) - clim_sat(2,ict) ) /( clim_sat(1,ict+1) - clim_sat(1,ict) ) &
                                    & * ( mod(time,ctau(iclim)) - clim_sat(1,ict) ) + clim_sat(2,ict)
                            elseif (ict == nclim(iclim)) then 
                                satup = ( clim_sat(2,1) - clim_sat(2,ict) ) /( dct(iclim)  ) &
                                    & * ( mod(time,ctau(iclim)) - clim_sat(1,ict) ) + clim_sat(2,ict)
                            endif 
                    endselect
                endif 
            enddo 
            if (dt >= minval(dct)/10d0 .or. any (ict_change) ) then 
                open(idust, file=trim(adjustl(flxdir))//'/'//'climate.txt', &
                    & status='old',action='write',position='append')
                write(idust,*) time,tc,qin,satup
                close(idust)
            endif 
            
            sat = min(1.0d0, 1d0-(1d0-satup)*(1d0-z/zsat)**2d0)
            do iz=1,nz
                if (z(iz)>=zsat) sat(iz)=1d0
            enddo 
            v = qin/poroi/sat
            torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
            tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
            
            c_disp = c0_disp*zdisp**c1_disp
            disp = c_disp * v

            if (disp_FULL_ON) disp = disp_FULL
            
        endif 
            
        call coefs_v2( &
            & nz,rg,rg2,tc,sec2yr,tempk_0,pro,cec_pH_depend,mcec_all,logkhaq_all &! input
            & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
            & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
            & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,mwt_all,staq_all &!input
            & ,nsp_aq,nsp_aq_cnst,chraq,chraq_cnst,maq,maqc &!input
            & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 &! output
            & ,keqaq_oxa,keqaq_cl &! output
            & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
            & ,keqcec_all,keqiex_all &! output 
            & ) 
        
        do isps = 1, nsp_sld
            ksld(isps,:) = ksld_all(findloc(chrsld_all,chrsld(isps),dim=1),:)
            ! print *,chrsld(isps),ksld(isps,:)
        enddo
        
        do ispa = 1, nsp_aq 
            daq(ispa) = daq_all(findloc(chraq_all,chraq(ispa),dim=1))
        enddo 
        
        do ispg = 1, nsp_gas 
            dgasa(ispg) = dgasa_all(findloc(chrgas_all,chrgas(ispg),dim=1))
            dgasg(ispg) = dgasg_all(findloc(chrgas_all,chrgas(ispg),dim=1))
        enddo 
        kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)
        kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
        knh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h0)
        kn2o = keqgas_h(findloc(chrgas_all,'pn2o',dim=1),ieqgas_h0)
        k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
        k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)
        k1nh3 = keqgas_h(findloc(chrgas_all,'pnh3',dim=1),ieqgas_h1)
        pco2i = mgasi_all(findloc(chrgas_all,'pco2',dim=1))
        pnh3i = mgasi_all(findloc(chrgas_all,'pnh3',dim=1))
        khco2i = kco2*(1d0+k1/proi + k1*k2/proi/proi)
        khnh3i = knh3*(1d0+proi/k1nh3)

        do ispg = 1, nsp_gas 
            select case(trim(adjustl(chrgas(ispg)))) 
                case('pco2')
                    khgasi(ispg) = khco2i
                case('po2')
                    khgasi(ispg) = kho
                case('pnh3')  
                    khgasi(ispg) = khnh3i
                case('pn2o')  
                    khgasi(ispg) = kn2o
            endselect 
        enddo
        
        ! print*,khgasi,knh3,k1nh3
        ! pause
        
        ! kinetic inhibition 
        if (al_inhibit) then 
            if (any(chraq == 'al')) then 
                do isps = 1, nsp_sld
                    if (staq(isps,findloc(chraq,'al',dim=1)) .ne. 0d0) then 
                        ksld(isps,:) = ksld(isps,:) &
                            & *10d0**(-4.84d0)/(10d0**(-4.84d0)+maq(findloc(chraq,'al',dim=1),:)) 
                    endif 
                enddo 
            endif 
        endif 
        
        ! time dependent kinetic inhibition
        if (anealing_dust) then
            do isps=1,nsp_sld
                if ( (rainpowder*rfrc_sld(isps)>0d0) .or. (rainpowder_2nd*rfrc_sld_2nd(isps)>0d0) ) then 
                    fkin(isps,:) = min(1d0,1d0/time)
                endif 
            enddo
        endif 
        
        save_trans = .false.
        zml = zml_background
        call make_transmx(  &
            & nsp_sld,imix,dz,poro,nz,z,zml,dbl_ref,tol,save_trans  &! input
            & ,trans  &! output 
            & )


        error = 1d4
        ! iter=0

    100 continue

        mgasx = mgas
        msldx = msld
        maqx = maq
        
        prox = pro  
        iosx = ios  
        
        ! so4f = so4fprev
        maqft = maqft_prev
        maqfads = maqfads_prev
        
        poroprev = poro
        hrprev = hr
        vprev = v
        torgprev = torg
        toraprev = tora
        wprev = w 
        
        dispprev = disp
        
        mblkx = mblk
        
        ! whether or not you are using psd
        psd_old = psd
        mpsd_old = mpsd
        psd_error_flg = .false.

        !  raining dust & OM 
        maqsupp = 0d0
        mgassupp = 0d0
        do isps = 1, nsp_sld
            if (imixtype==imixtype_nobio) then 
                ! msldsupp(isps,:) = rainpowder*rfrc_sld(isps)*exp(-z/zsupp)/zsupp &
                    ! & + rainpowder_2nd*rfrc_sld_2nd(isps)*exp(-z/zsupp)/zsupp
                ! modify to homogeneously distribute (04-10-2023)
                msldsupp(isps,:) = rainpowder*rfrc_sld(isps)/sum(dz(:)) &
                    & + rainpowder_2nd*rfrc_sld_2nd(isps)/sum(dz(:))
            else 
                msldsupp(isps,1) = rainpowder*rfrc_sld(isps)/dz(1)  &
                    & + rainpowder_2nd*rfrc_sld_2nd(isps)/dz(1)
                
            endif 
            ! checking mass balance
            if (  (rainpowder*rfrc_sld(isps) > 0d0)  & 
                & .and. ( abs(sum(msldsupp(isps,:)*dz(:))-rainpowder*rfrc_sld(isps))/rainpowder*rfrc_sld(isps) > 1d-6) ) then 
                print *, 'dust error? going to stop',chrsld(isps),sum(msldsupp(isps,:)*dz(:)),rainpowder*rfrc_sld(isps)
                stop
            endif 
        enddo 
        
        ! dust options check 
        if (dust_wave .and. dust_step) then 
            print *
            print *, 'CAUTION: options of dust_wave and dust_step are both ON'
            print * 
            stop
        endif 
        
        ! determine whether or not dust depending on pH 
        dust_off = .false.
        if (ph_limits_dust) then 
            ph_ave = 0d0
            z_ave  = 0d0
            do iz=1,nz
                if (z(iz)<=z_chk_ph) then 
                    ph_ave = ph_ave + pro(iz)*dz(iz)
                    z_ave  = z_ave  +         dz(iz)
                endif 
            enddo 
            ph_ave = -log10(ph_ave/z_ave)
            
            if (ph_ave > ph_lim) dust_off = .true.
        endif 
        
        ! when considering static system, dust is given only initial year  
        dust_off = .false.
        if (aq_close) then 
            ! homogeneously distribute
            do isps = 1, nsp_sld
                msldsupp(isps,:) = rainpowder*rfrc_sld(isps)/dz(:)/nz &
                    & + rainpowder_2nd*rfrc_sld_2nd(isps)/dz(:)/nz
                ! print*,chrsld(isps),sum(msldsupp(isps,:)*dz(:)),rainpowder*rfrc_sld(isps)
            enddo 
            ! stop
            if (dust_step .and.  time > step_tau) dust_off = .true.
        endif 
        
        ! if defined wave function is imposed on dust 
        if (dust_wave) then 
            do isps = 1, nsp_sld
                if (imixtype==imixtype_nobio) then 
                    msldsupp(isps,:) = msldsupp(isps,:)*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                else 
                    msldsupp(isps,1) = msldsupp(isps,1)*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                endif 
            enddo 
            if (time==0d0 .or. dust_norm /= merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))) then
                open(idust, file=trim(adjustl(flxdir))//'/'//'dust.txt', &
                    & status='old',action='write',position='append')
                write(idust,*) time-dt,dust_norm
                write(idust,*) time,merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                dust_norm = merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                close(idust)
            endif 
        endif 
        
        ! non continueous
        if (dust_step) then 
            
            dust_change = .false.
            
            if (.not. dust_off) then 
            
                ! dust_norm_prev = dust_norm
                
                if (dt > step_tau) then 
                    dt = step_tau
                    ! go to 100
                endif 
                
                if (step_tau + floor(time) < time  .and.  time + dt < 1d0 + floor(time)) then 
                    continue
                elseif (step_tau + floor(time) < time  .and.  time + dt >= 1d0 + floor(time)) then 
                    dt = 1d0 + floor(time) - time + step_tau * tol_step_tau 
                    dust_change = .true.
                elseif (0d0 + floor(time) <= time  .and.  time +dt <= step_tau + floor(time)) then 
                    if (dt > step_tau/10d0) then 
                        dt = step_tau/10d0
                    endif 
                elseif (0d0 + floor(time) <= time  .and.  time +dt > step_tau + floor(time)) then 
                    dt = step_tau + floor(time) - time + step_tau * tol_step_tau 
                    if (dt > step_tau/10d0) then 
                        dt = step_tau/10d0
                    endif 
                    ! checking again
                    if (0d0 + floor(time) <= time  .and.  time +dt > step_tau + floor(time)) then 
                        dust_change = .true.
                    endif  
                endif 
                
                ! if (time < step_tau + floor(time)  .and. time + dt >= step_tau + floor(time) ) then 
                    ! if ( ( step_tau + floor(time) - time  ) > step_tau * tol_step_tau ) then 
                        ! dt = step_tau - ( time - floor(time) ) !+ step_tau * tol_step_tau 
                    ! endif 
                    ! if (dt > step_tau/10d0) then 
                        ! dt = step_tau/10d0
                    ! endif 
                ! elseif (time >= step_tau + floor(time) .and. time + dt >= 1d0 + floor(time) ) then 
                    ! if ( ( 1d0 + floor(time) - time  ) > step_tau * tol_step_tau ) then 
                        ! dt = 1d0 + floor(time)  - time + step_tau * tol_step_tau 
                    ! endif
                ! endif 
                

                ! an attempt to use fickian mixing when applying rock powder
                ! when not applying the powder use the chosen mixing (can be fickian)
                
                ! if (time - floor(time) >= step_tau) then 
                ! if (time - floor(time) > step_tau) then 
                ! if (step_tau + floor(time) < time  .and.  time + dt < 1d0 + floor(time)) then 
                if (step_tau + floor(time) < time  ) then 
                    ! print *, 'no dust time', time 
                    msldsupp = 0d0
                    dust_norm = 0d0
                    ! only implement background mixing (fickian mixing as default)
                    imix    = imixtype_background
                    ! zml     = zsupp ! mixed layer depth is common to all solid sp. 
                    zml     = zml_background ! mixed layer depth is common to all solid sp. 
                ! else 
                ! elseif (0d0 + floor(time) <= time  .and.  time + dt <= step_tau + floor(time)) then 
                elseif (0d0 + floor(time) <= time  .and.  time <= step_tau + floor(time)) then 
                    ! print *, 'dust time !!', time 
                    msldsupp = msldsupp/step_tau
                    dust_norm = 1d0/step_tau
                    ! only implement chosen mixing 
                    imix    = imixtype
                    ! zml     = zml_ref ! mixed layer depth is the value specified for dust  
                    zml     = zml_dust ! mixed layer depth is the value specified for dust  
                    ! OM is mixed in fickian (implemented later)
                else 
                    print *, 'Fatale error in dusting?',time,floor(time),dt,step_tau
                    stop
                endif 
            
            ! not dusting if ph is too high
            ! if (dust_off) then 
            elseif (dust_off) then 
                msldsupp = 0d0
                dust_norm = 0d0
                ! only implement background mixing (default Fickian)
                imix    = imixtype_background 
                ! zml     = zsupp ! mixed layer depth is common to all solid sp. 
                zml     = zml_background ! mixed layer depth is common to all solid sp. 
            endif 

            ! new 5/18/2023 YK
            ! OM is mixed in fickian regardless of dust implementation
            do isps = 1, nsp_sld
                if ( rfrc_sld_plant(isps) > 0d0 ) then 
                    imix(isps)  = imixtype_OM
                    zml(isps)   = zml_OM
                endif 
            enddo
            
            ! mixing reload
            save_trans = .false.
            call make_transmx(  &
                & nsp_sld,imix,dz,poro,nz,z,zml,dbl_ref,tol,save_trans  &! input
                & ,trans  &! output 
                & )
            
            ! if ( dust_norm /= dust_norm_prev ) then
                ! open(idust, file=trim(adjustl(flxdir))//'/'//'dust.txt', &
                    ! & status='old',action='write',position='append')
                ! write(idust,*) time-dt_prev,dust_norm_prev
                ! write(idust,*) time,dust_norm
                ! close(idust)
            ! endif 
            
        endif 
        
        ! overload with OM rain 
        do isps = 1, nsp_sld
            if (imixtype==imixtype_nobio) then 
                selectcase(trim(adjustl(chrsld(isps))))
                    case('g1','g2','g3')
                        ! msldsupp(isps,:) = msldsupp(isps,:) &
                            ! & + plant_rain/12d0/((1d0-poroi)*rho_grain*1d6) &! converting g_C/g_soil/yr to mol_C/m3_soil/yr
                            ! & *1d0 &! assuming 1m depth to which plant C is supplied 
                            ! & *rfrc_sld_plant(isps) &
                            ! & *exp(-z/zsupp_plant)/zsupp_plant
                        msldsupp(isps,:) = msldsupp(isps,:) &
                            & + plant_rain/12d0*rfrc_sld_plant(isps)*exp(-z/zsupp_plant)/zsupp_plant ! when plant_
                    case('amnt','cc')
                        msldsupp(isps,:) = msldsupp(isps,:) &
                            & + plant_rain/mwt(isps)*rfrc_sld_plant(isps)*exp(-z/zsupp_plant)/zsupp_plant ! when plant_
                    case default
                        ! not doing anything
                endselect
            else 
                selectcase(trim(adjustl(chrsld(isps))))
                    case('g1','g2','g3')
                        msldsupp(isps,1) = msldsupp(isps,1) &
                            & + plant_rain/12d0*rfrc_sld_plant(isps)/dz(1) ! when plant_rain is in g_C/m2/yr
                    case('amnt','cc')
                        msldsupp(isps,1) = msldsupp(isps,1) &
                            & + plant_rain/mwt(isps)*rfrc_sld_plant(isps)/dz(1) ! when plant_rain is in g_(isps)/m2/yr
                    case default
                        ! not doing anything
                        msldsupp(isps,1) = msldsupp(isps,1) &
                            & + plant_rain/mwt(isps)*rfrc_sld_plant(isps)/dz(1) ! when plant_rain is in g_(isps)/m2/yr
                endselect
            endif 
        enddo 
        ! when enforcing solid states without previous OM spin-up
        if (sld_enforce .and. (.not.read_data)) then 
            if (any(chrgas=='pco2')) then 
                ! mgassupp(findloc(chrgas,'pco2',dim=1),:) = plant_rain/12d0*exp(-z/zsupp_plant)/zsupp_plant
                mgassupp(findloc(chrgas,'pco2',dim=1),:) = plant_rain/12d0/ztot
            endif 
        endif 
        
        if (any(isnan(msldsupp))) then 
            print *, 'error in dust'
            stop
        endif 
        
        ! do PSD for raining dust & OM 
        if (do_psd) then 

            if (dust_norm>0d0) open(ipsd,file = trim(adjustl(profdir))//'/'//'psd_rain.txt',status = 'replace')
            if (dust_norm>0d0) open(ipsdv,file = trim(adjustl(profdir))//'/'//'intpsd_rain.txt',status = 'replace')
            
            if (do_psd_full) then 
                if (dust_norm>0d0) write(ipsd,*) ' sldsp\log10(radius) ', (ps(ips),ips=1,nps), 'time'
                if (dust_norm>0d0) write(ipsdv,*) ' sldsp\diameter(um) ', (10d0**ps(ips)*1d6*2d0,ips=1,nps), 'p80(um)'
            
                do iz = 1,nz

                    ! balance for volumes
                    ! sum(msldsupp*mv*1d-6) *dt (m3/m3) must be equal to sum( 4/3(pi)r3 * psd_rain * dps) 
                    ! where psd is number / bulk m3 / log r
                    do isps = 1, nsp_sld
                    
                        if ( rfrc_sld_plant(isps) > 0d0 ) then ! those rain with OM 
                            psu_rain_list(:)        = log10(p80)
                            pssigma_rain_list(:)    = ps_sigma_std
                            psw_rain_list(:)        = 1d0
                        else  ! dust except for OM associates
                            psu_rain_list           = psu_rain_list_in
                            pssigma_rain_list       = pssigma_rain_list_in
                            psw_rain_list           = psw_rain_list_in
                        endif 
                
                        ! rained particle distribution 
                        if (read_data) then !!! only to indicate the case when dust is added (usually restart experiment)
                            psd_rain(:,iz) = 0d0
                            do ips = 1, nps_rain_char
                                psw_rain = psw_rain_list(ips)
                                psu_rain = psu_rain_list(ips)
                                pssigma_rain = pssigma_rain_list(ips)
                                psd_rain(:,iz) = psd_rain(:,iz) &
                                    & + psw_rain*1d0/pssigma_rain/sqrt(2d0*pi) &
                                    &   *exp( -0.5d0*( (ps(:) - psu_rain)/pssigma_rain )**2d0 )
                            enddo 
                        else
                            psu_rain = log10(p80)
                            pssigma_rain = ps_sigma_std
                            psd_rain(:,iz) = 1d0/pssigma_rain/sqrt(2d0*pi)*exp( -0.5d0*( (ps(:) - psu_rain)/pssigma_rain )**2d0 )
                        endif 
                        ! to ensure sum is 1
                        ! print *, sum(psd_rain*dps)
                        psd_rain(:,iz) = psd_rain(:,iz)/sum(psd_rain(:,iz)*dps(:)) 
                        ! print *, sum(psd_rain*dps)
                        ! stop
                        
                        psd_rain_tmp = psd_rain(:,iz)
                    
                    
                        volsld = msldsupp(isps,iz)*mv(isps)*1d-6
                        psd_rain(:,iz) = psd_rain_tmp * volsld *dt &
                            ! & /(1d0 - poroi)  &
                            & /sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_rain_tmp*dps(:))
                            
                        if ( abs( ( volsld *dt &
                            ! & /(1d0 - poroi) &
                            & - sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_rain(:,iz)*dps(:))) &
                            & / ( volsld*dt &
                            ! & /(1d0 - poroi) &
                            & ) ) > tol) then 
                            print *,iz, volsld*dt &
                                ! & /(1d0 - poroi) &
                                & ,sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_rain(:,iz)*dps(:))
                            stop
                        endif 
                        mpsd_rain(isps,:,iz) = psd_rain(:,iz)
                
                        if ((dust_norm>0d0) .and. (iz==1)) write(ipsd,*) chrsld(isps),(mpsd_rain(isps,ips,iz),ips=1,nps), time
                        if ((dust_norm>0d0) .and. (iz==1))  then                         
                            intpsd_tmp = psd_rain_tmp*(10d0**ps)**3d0
                            intpsd = intpsd_tmp
                            do ips = 1, nps
                                intpsd(ips) = sum(intpsd_tmp(1:ips))/sum(intpsd_tmp)
                            enddo
                            call calc_p80( &
                                & nps,ps,intpsd &! input 
                                & ,p80_tmp &! output
                                & )
                            write(ipsdv,*) chrsld(isps),(intpsd(ips),ips=1,nps), p80_tmp
                        endif 
                    enddo 
                enddo 
            else
                if (dust_norm>0d0) write(ipsd,*) ' depth\log10(radius) ', (ps(ips),ips=1,nps), 'time'
                if (dust_norm>0d0) write(ipsdv,*) ' sldsp\diameter(um) ', (10d0**ps(ips)*1d6*2d0,ips=1,nps), 'p80(um)'
            
                psu_rain_list           = psu_rain_list_in
                pssigma_rain_list       = pssigma_rain_list_in
                psw_rain_list           = psw_rain_list_in
                            
                do iz = 1,nz
                
                    ! rained particle distribution 
                    if (read_data) then 
                        psd_rain(:,iz) = 0d0
                        do ips = 1, nps_rain_char
                            psw_rain = psw_rain_list(ips)
                            psu_rain = psu_rain_list(ips)
                            pssigma_rain = pssigma_rain_list(ips)
                            psd_rain(:,iz) = psd_rain(:,iz) &
                                & + psw_rain*1d0/pssigma_rain/sqrt(2d0*pi) &
                                &   *exp( -0.5d0*( (ps(:) - psu_rain)/pssigma_rain )**2d0 )
                        enddo 
                    else
                        psu_rain = log10(p80)
                        pssigma_rain = 1d0
                        pssigma_rain = ps_sigma_std
                        psd_rain(:,iz) = 1d0/pssigma_rain/sqrt(2d0*pi)*exp( -0.5d0*( (ps(:) - psu_rain)/pssigma_rain )**2d0 )
                    endif 
                    ! to ensure sum is 1
                    ! print *, sum(psd_rain*dps)
                    psd_rain(:,iz) = psd_rain(:,iz)/sum(psd_rain(:,iz)*dps(:)) 
                    ! print *, sum(psd_rain*dps)
                    ! stop
                    
                    psd_rain_tmp = psd_rain(:,iz)

                    ! balance for volumes
                    ! sum(msldsupp*mv*1d-6) *dt (m3/m3) must be equal to sum( 4/3(pi)r3 * psd_rain * dps) 
                    ! where psd is number / bulk m3 / log r
                    volsld = sum(msldsupp(:,iz)*mv(:)*1d-6)
                    psd_rain(:,iz) = psd_rain_tmp * volsld *dt &
                        ! & /(1d0 - poroi)  &
                        & /sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_rain_tmp*dps(:))
                        
                    if ( abs( ( volsld *dt &
                        ! & /(1d0 - poroi) &
                        & - sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_rain(:,iz)*dps(:))) &
                        & / ( volsld*dt &
                        ! & /(1d0 - poroi) &
                        & ) ) > tol) then 
                        print *,iz, volsld*dt &
                            ! & /(1d0 - poroi) &
                            & ,sum(4d0/3d0*pi*(10d0**ps(:))**3d0*psd_rain(:,iz)*dps(:))
                        stop
                    endif 
                
                    if (dust_norm>0d0) write(ipsd,*) z(iz),(psd_rain(ips,iz),ips=1,nps), time
                    if (dust_norm>0d0) then
                        intpsd_tmp = psd_rain_tmp*(10d0**ps)**3d0
                        intpsd = intpsd_tmp
                        do ips = 1, nps
                            intpsd(ips) = sum(intpsd_tmp(1:ips))/sum(intpsd_tmp)
                        enddo
                        call calc_p80( &
                            & nps,ps,intpsd &! input 
                            & ,p80_tmp &! output
                            & )
                        write(ipsdv,*) chrsld(isps),(intpsd(ips),ips=1,nps), p80_tmp
                    endif 
                enddo 
                
            endif 
            
            if (dust_norm>0d0) close(ipsd)
            if (dust_norm>0d0) close(ipsdv)
        endif 

        ! if ((.not.read_data) .and. it == 0 .and. iter == 0) then 
            ! do ispa = 1, nsp_aq
                ! if (chraq(ispa)/='so4') then
                    ! maqx(ispa,1:) = 1d2
                ! endif 
            ! enddo
        ! endif 
        
        poro_iter = 0
        poro_error = 1d4
        poro_tol = 1d-6
        poro_iter_max = 50
        
    ! #ifdef poroiter
        poro_tol = 1d-9
        ! if (iwtype == iwtype_flex) poro_tol = 1d-14
        do while (poro_error > poro_tol) ! start of porosity iteration 
    ! #endif 

        porox = poro
        wx = w
        
        ads_ON_tmp = .false.
        if (ads_ON) ads_ON_tmp = .true.
            
        call alsilicate_aq_gas_1D_v3_2( &
            ! new input 
            & nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext &
            & ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  &
            & ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi &
            & ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp &
            & ,stgas_ext,stgas_dext,staq_ext,stsld_ext,staq_dext,stsld_dext &
            & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst &
            & ,chraq_cnst,chraq_all,chrgas_cnst,chrgas_all,chrsld_all &
            & ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,keqaq_s,keqaq_no3,keqaq_nh3,keqaq_oxa,keqaq_cl &
            & ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &
            & ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain,msldth_all,mv_all,staq_all,stgas_all &
            & ,trans,display,chrflx,sld_enforce &! input
            & ,nsld_kinspc,chrsld_kinspc,kin_sld_spc &! input
            & ,precstyle,solmod,fkin &! in 
            !  old inputs
            & ,hr,poro,z,dz,w_btm,sat,pro,poroprev,tora,v,tol,it,nflx,kw,maqft_prev,disp & 
            & ,ucv,torg,cplprec,rg,tc,sec2yr,tempk_0,proi,poroi,up,dwn,cnr,adf,msldunit  &
            & ,ads_ON_tmp,maqfads_prev,keqcec_all,keqiex_all,cec_pH_depend,aq_close,ios,act_ON,beta_all &
            ! old inout
            & ,dt,flgback,w &    
            ! output 
            & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec,rxnsld,flx_co2sp,maqft &
            & ,maqfads,msldf_loc,beta_loc,iosx &
            & )
        
        
        ! if (dt > 2d-2) stop
        
        
        save_trans = .false.
        call make_transmx(  &
            & nsp_sld,imix,dz,poro,nz,z,zml,dbl_ref,tol,save_trans  &! input
            & ,trans  &! output 
            & )
            
        ! sum(msldx*mv*1d-6) + mblkx*mvblk*1d-6 = 1d0 - poro
        if ( incld_blk ) then 
            do iz=1,nz
                mblkx(iz) = 1d0 - poro(iz) - sum(msldx(:,iz)*mv(:)*1d-6)
                mblkx(iz) = mblkx(iz)/(mvblk*1d-6)
            enddo 
        else 
            mblkx = 0d0
        endif 

        if (flgback) then 
            flgback = .false. 
            flgreducedt = .true.
            ! pre_calc = .true.
            dt = dt/1d1
            psd = psd_old
            mpsd = mpsd_old
            poro = poroprev
            torg = torgprev
            tora = toraprev
            disp = dispprev
            v = vprev
            hr = hrprev
            w = wprev
            call calcupwindscheme(  &
                up,dwn,cnr,adf & ! output 
                ,w,nz   & ! input &
                )
            go to 100
        endif    
        
        if (poroevol) then 
            ! poroprev = poro
    ! #ifdef surfssa
            ! mvab = mvab_save 
            ! mvan = mvan_save 
            ! mvcc = mvcc_save 
            ! mvfo = mvfo_save 
            ! mvka = mvka_save 
    ! #endif 
            ! poro = poroi + (mabi-mabx)*(mvab)*1d-6  &
                ! & +(mfoi-mfox)*(mvfo)*1d-6 &
                ! & +(mani-manx)*(mvan)*1d-6 &
                ! & +(mcci-mccx)*(mvcc)*1d-6 &
                ! & +(mkai-mkax)*(mvka)*1d-6 
            if (iwtype == iwtype_flex) then 
                poro = poroi
                ! not constant but calculated as defined (only applicable when unit of msld(x) is mol per bulk soil)
                do iz = 1,nz
                    poro(iz) = 1d0 - sum(msldx(:,iz)*mv(:)*1d-6)
                enddo 
            else 
                call calc_poro( &
                    & nz,nsp_sld,nflx,idif,irain &! in
                    & ,flx_sld,mv,poroprev,w,poroi,w_btm,dz,tol,dt &! in
                    & ,poro &! inout
                    & )
            endif 
            ! poro = poroi
            ! do isps=1,nsp_sld
                ! poro = poro + (msldi(isps)-msldx(isps,:))*mv(isps)*1d-6
            ! enddo
            ! do iz=1,nz
                ! DV(iz) = 0d0
                ! do isps = 1,nsp_sld 
                    ! DV(iz) = DV(iz) + ( flx_sld(isps, 4 + isps,iz) + flx_sld(isps, idif ,iz) + flx_sld(isps, irain ,iz) ) &
                        ! & *mv(isps)*1d-6*dt 
                ! enddo 
                ! poro(iz) = poroprev(iz) - DV(iz)
            ! enddo 
            
            if (any(poro<0d0)) then 
                print*,'negative porosity: stop'
                print*,poro
                
                flgback = .false. 
                flgreducedt = .true.
                ! pre_calc = .true.
                dt = dt/1d1
                psd = psd_old
                mpsd = mpsd_old
                poro = poroprev
                torg = torgprev
                tora = toraprev
                disp = dispprev
                v = vprev
                hr = hrprev
                w = wprev
                call calcupwindscheme(  &
                    up,dwn,cnr,adf & ! output 
                    ,w,nz   & ! input &
                    )
                go to 100
                
                ! w = w*2d0
                ! go to 100
                stop
            endif 
            if (any(poro>1d0)) then 
                print*,'porosity exceeds 1: stop'
                print*,poro
                
                flgback = .false. 
                flgreducedt = .true.
                ! pre_calc = .true.
                dt = dt/1d1
                psd = psd_old
                mpsd = mpsd_old
                poro = poroprev
                torg = torgprev
                tora = toraprev
                disp = dispprev
                v = vprev
                hr = hrprev
                w = wprev
                call calcupwindscheme(  &
                    up,dwn,cnr,adf & ! output 
                    ,w,nz   & ! input &
                    )
                go to 100
                
                ! w = w*2d0
                ! go to 100
                stop
            endif 
            
    ! #ifdef surfssa
            ! mvab = mwtab 
            ! mvan = mwtan 
            ! mvcc = mwtcc 
            ! mvfo = mwtfo 
            ! mvka = mwtka 
    ! #endif 
            v = qin/poro/sat
            torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
            tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
            
            c_disp = c0_disp*zdisp**c1_disp
            disp = c_disp * v

            if (disp_FULL_ON) disp = disp_FULL
            
    #ifndef calcw_full
            w(:) = w0 
            dwsporo = 0d0        
            wsporo = w_btm*(1d0 - poroi)
            if (noncnstw) then 
                ! w(:) = w0/(1d0-poro(:)) ! --- isovolumetric weathering?
                do iz=1,nz
                    DV(iz) = 0d0
                    do isps = 1,nsp_sld 
                        DV(iz) = DV(iz) + ( flx_sld(isps, 4 + isps,iz) + flx_sld(isps, idif ,iz) + flx_sld(isps, irain ,iz) ) &
                            & *mv(isps)*1d-6*dt 
                    enddo 
                    dwsporo(iz) = -( ( poro(iz) - poroprev(iz))/dt - DV(iz)/dt )
                enddo 
                
                do iz = nz,1,-1
                    if (iz==nz) then 
                        ! (wsporo(iz+1) - wsporo(iz))/dz(iz) = dwsporo(iz)
                        wsporo(iz) =  w_btm*(1d0 - poroi) - dwsporo(iz)*dz(iz)
                    else
                        wsporo(iz) =  wsporo(iz+1) - dwsporo(iz)*dz(iz)
                    endif 
                enddo 
                ! wsporo = w * (1d0 - poro)
                w = wsporo/(1d0-poro)
            endif 
            
            wsporo = w_btm*(1d0 - poroi)
            w = wsporo/(1d0-poro)
            
            w = w_btm
            
            call calc_uplift( &
                & nz,nsp_sld,nflx,idif,irain &! IN
                & ,iwtype &! in
                & ,flx_sld,mv,poroi,w_btm,dz,poro,poroprev,dt &! in
                & ,w &! inout
                & )
            
            ! ------------ determine calculation scheme for advection (from IMP code)
            call calcupwindscheme(  &
                up,dwn,cnr,adf & ! output 
                ,w,nz   & ! input &
                )
    #endif 
            do isps=1,nsp_sld
                hr(isps,:) = hri(isps,:)*rough(isps,:)
                if (surfevol1 ) then 
                    hr(isps,:) = hri(isps,:)*rough(isps,:)*((1d0-poro)/(1d0-poroi))**(2d0/3d0)
                endif 
                if (surfevol2 ) then 
                    hr(isps,:) = hri(isps,:)*rough(isps,:)*(poro/poroi)**(2d0/3d0)  ! SA increases with porosity 
                endif 
            enddo 
            ! if doing psd SA is calculated reflecting psd
            if (do_psd) then 
                ! hr = ssa*(1-poro)/poro ! converting m2/sld-m3 to m2/pore-m3
                ! hr = ssa
                do isps=1,nsp_sld
                    if (do_psd_full) then 
                        ! hr(isps,:) = ssa(isps,:)/poro/(msldx(isps,:)*mv(isps)*1d-6) ! so that poro * hr * mv * msld becomes porosity independent 
                        ! hr(isps,:) = ssa(isps,:)/poro/(msld(isps,:)*mv(isps)*1d-6) ! so that poro * hr * mv * msld becomes porosity independent 
                        ! hr(isps,:) = ssa(isps,:) ! so that poro * hr * mv * msld becomes porosity independent 
                        hr(isps,:) = ssa(isps,:)/ssv(isps,:)/poro 
                        ! hr(isps,:) = ssav(isps,:)/poro 
                    else
                        hr(isps,:) = ssa(isps,:)/poro ! so that poro * hr * mv * msld becomes porosity independent 
                    endif 
                enddo 
            endif 
            
            if (display .and. (.not. display_lim)) then 
                write(chrfmt,'(i0)') nz_disp
                chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
                print *
                print *,' [porosity & surface area]'
                print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
                print trim(adjustl(chrfmt)),'poro',(poro(iz),iz=1,nz, nz/nz_disp)
                print trim(adjustl(chrfmt)),'SA',(hrb(iz),iz=1,nz, nz/nz_disp)
                print *
            endif 
        endif  
        
    ! #ifdef poroiter
        if (poroiter_in) then 
            if (iwtype == iwtype_flex) then 
                ! poro_error = maxval ( abs (( w - wx )/wx ) )
                poro_error = maxval ( abs ( w - wx ) )
            else
                ! poro_error = maxval ( abs (( poro - porox )/porox ) )
                poro_error = maxval ( abs ( poro - porox ) )
            endif 
            
            print *, 'porosity iteration: ',poro_iter,poro_error
            poro_iter = poro_iter + 1
                        
            if (poro_iter > poro_iter_max) then 
                print *, 'too much porosity iteration but does not converge within assumed threshold'
                print *, 'reducing dt and move back'
                flgback = .false. 
                flgreducedt = .true.
                psd = psd_old
                mpsd = mpsd_old
                poro = poroprev
                torg = torgprev
                tora = toraprev
                disp = dispprev
                v = vprev
                hr = hrprev
                w = wprev
                call calcupwindscheme(  &
                    up,dwn,cnr,adf & ! output 
                    ,w,nz   & ! input &
                    )
                ! pre_calc = .true.
                dt = dt/1d1
                go to 100
            endif    
        else 
            poro_error = 0d0
            EXIT 
        endif 
        
        enddo ! porosity iteration end
    ! #endif 
        ! attempt to do psd
        if (do_psd) then 
            
            if (display) then 
                print *
                print *, '-- doing PSD'
            endif 
            
            ! when not doing PSDFULL, not do while loop for PSD 
            ! if (.not. do_psd_full) psd_loop = .false.
            
            
            ! print '(3(1x,a11))','>>>PSD time','dt','%done'

            dt_pbe = dt
            dt_save = dt
            time_pbe = 0d0
            ddpsd = 0d0
            do while(time_pbe < dt_save)
            ! print *
            ! print *,' ---- PSD time: ',time_pbe,' dt: ',dt, ' completeness [%]: ',100d0*time_pbe/dt_save
            dt = dt*1.05d0
            ! dt = dt*2d0
            dt_pbe = dt ! temporary recording dt 
            if (time_pbe + dt_pbe > dt_save) dt_pbe = dt_save - time_pbe ! modifying temporary dt not to exceed dt_save
            dt = dt_pbe ! dt returned 
            psd_save = psd
            mpsd_save_2 = mpsd
            
            if (.not.psd_loop) dt = dt_save
            
            if (do_psd_full) then 
                
                dmpsd = 0d0
                
                do isps = 1, nsp_sld
                
                    if ( trim(adjustl(precstyle(isps))) == 'decay') cycle ! solid species not related to SA
                    
                    if (psd_enable_skip .and. any(chrsld_nopsd == chrsld(isps)) ) cycle
                
                    DV(:) = flx_sld(isps, 4 + isps,:)*mv(isps)*1d-6*dt  

                    dpsd = 0d0
                    psd = mpsd(isps,:,:)
                    
                    if (.not. psd_impfull) then 
                        ! call psd_diss( &
                        ! call psd_diss_pbe_expall( &
                        ! call psd_diss_pbe_exp( &
                        call psd_diss_pbe( &
                            & nz,nps &! in
                            & ,z,DV,dt,pi,tol_dvd,poro &! in 
                            & ,incld_rough,roughref(isps) &! in
                            & ,psd,ps,dps,ps_min,ps_max &! in 
                            & ,chrsld(isps) &! in 
                            & ,dpsd,psd_error_flg &! inout
                            & )

                        if ( flgback .or. psd_error_flg) then 
                            print *, '*** error raised after PBE calc. for ',chrsld(isps)
                            print *, '*** escape from do-loop'
                            exit
                        endif 
                        
                        ! dt_pbe = dt
                        ! time_pbe = 0
                        ! ddpsd = 0d0
                        ! do while(time_pbe < dt)
                            ! if (time_pbe + dt_pbe > dt) dt_pbe = dt - time_pbe
                            ! DV(:) = flx_sld(isps, 4 + isps,:)*mv(isps)*1d-6*dt_pbe 
                            ! psd_save = psd
                            ! dpsd_save = dpsd
                            
                            ! call psd_diss_pbe_expall( &
                            ! call psd_diss_pbe_exp( &
                            ! call psd_diss_pbe( &
                                ! & nz,nps &! in
                                ! & ,z,DV,dt_pbe,pi,tol_dvd,poro &! in 
                                ! & ,incld_rough,rough_c0,rough_c1 &! in
                                ! & ,psd,ps,dps,ps_min,ps_max &! in 
                                ! & ,chrsld(isps) &! in 
                                ! & ,ddpsd,psd_error_flg &! inout
                                ! & )
                                
                            ! if (psd_error_flg) then 
                                ! psd_error_flg = .false.
                                ! dt_pbe = dt_pbe/10d0
                                ! psd = psd_save
                                ! dpsd = dpsd_save
                                ! cycle
                            ! endif 
                            
                            ! psd = psd + ddpsd
                            ! dpsd = dpsd + ddpsd
                            
                            ! time_pbe = time_pbe + dt_pbe
                            
                        ! enddo 
                        
                        
                        ! if (psd_error_flg) then 
                            ! psd_error_flg = .false. 
                            ! flgback = .false. 
                            ! flgreducedt = .true.
                            ! psd = psd_old
                            ! mpsd = mpsd_old
                            ! poro = poroprev
                            ! torg = torgprev
                            ! tora = toraprev
                            ! disp = dispprev
                            ! v = vprev
                            ! hr = hrprev
                            ! w = wprev
                            ! call calcupwindscheme(  &
                                ! up,dwn,cnr,adf & ! output 
                                ! ,w,nz   & ! input &
                                ! )
                            ! dt = dt/1d1
                            ! go to 100
                        ! endif 
                    else
                        do iz=1,nz
                            dpsd(:,iz) = DV(iz)/nps/dps(:)
                        enddo
                    endif 
                    
                    dmpsd(isps,:,:) = dpsd
                    
                    ! print *, chrsld(isps),DV
                
                enddo 

                if ( flgback .or. psd_error_flg) then 
                    print *,' *** because of PBE error, returning to do while loop with reduced dt'
                    if (.not.psd_loop) exit
                    flgback = .false.
                    psd_error_flg = .false.
                    dt = dt/10d0
                    psd = psd_save
                    mpsd = mpsd_save_2
                    cycle
                endif 
        
            else 

                DV = 0d0
                do isps = 1,nsp_sld 
                    do iz=1,nz
                        ! DV(iz) = DV(iz) + flx_sld(isps, 4 + isps,iz)*mv(isps)*1d-6*dt/(1d0 - poro(iz))  
                        DV(iz) = DV(iz) + flx_sld(isps, 4 + isps,iz)*mv(isps)*1d-6*dt  
                    enddo 
                enddo 

                dpsd = 0d0
                
                if (.not. psd_impfull) then 
                    ! call psd_diss( &
                        ! & nz,nps &! in
                        ! & ,z,DV,dt,pi,tol,poro &! in 
                        ! & ,incld_rough,rough_c0,rough_c1 &! in
                        ! & ,psd,ps,dps,ps_min,ps_max &! in 
                        ! & ,' blk ' &! in 
                        ! & ,dpsd,psd_error_flg &! inout
                        ! & )
                    ! call psd_diss_pbe_expall( &
                    ! call psd_diss_pbe_exp( &
                    call psd_diss_pbe( &
                        & nz,nps &! in
                        & ,z,DV,dt,pi,tol,poro &! in 
                        & ,incld_rough,roughref_b &! in
                        & ,psd,ps,dps,ps_min,ps_max &! in 
                        & ,' blk ' &! in 
                        & ,dpsd,psd_error_flg &! inout
                        & )

                    if ( flgback .or. psd_error_flg) then 
                        print *,' *** because of PBE error, returning to do while loop with reduced dt'
                        if (.not.psd_loop) exit
                        flgback = .false.
                        psd_error_flg = .false.
                        dt = dt/10d0
                        psd = psd_save
                        mpsd = mpsd_save_2
                        cycle
                    endif 
                        
                    ! if (psd_error_flg) then 
                        ! psd_error_flg = .false. 
                        ! flgback = .false. 
                        ! flgreducedt = .true.
                        ! psd = psd_old
                        ! mpsd = mpsd_old
                        ! poro = poroprev
                        ! torg = torgprev
                        ! tora = toraprev
                        ! disp = dispprev
                        ! v = vprev
                        ! hr = hrprev
                        ! w = wprev
                        ! call calcupwindscheme(  &
                            ! up,dwn,cnr,adf & ! output 
                            ! ,w,nz   & ! input &
                            ! )
                        ! dt = dt/1d1
                        ! go to 100
                    ! endif 
                else 
                    do iz=1,nz
                        dpsd(:,iz) = DV(iz)/nps/dps(:)
                    enddo
                endif 
            
            endif 
            
            ! if (psd_impfull) do_psd_norm = .false.
            
            if (do_psd_norm) then 
            
                if (do_psd_full) then 
                    
                    flx_max_max = 0d0
                    
                    do isps=1,nsp_sld
                        
                        ! if ( trim(adjustl(precstyle(isps))) == 'decay' ) then ! solid species not related to SA              
                        if ( &
                            & trim(adjustl(precstyle(isps))) == 'decay' &! solid species not related to SA              
                            & .or. ( psd_enable_skip .and. any(chrsld_nopsd == chrsld(isps)) ) &!case when not-tracking PSDs for fast reacting minerals (SA not matter?)
                            & ) then
                            flx_mpsd(isps,:,:,:) = 0d0
                            do iz=1,nz
                                mpsdx(isps,:,iz) = mpsd_pr(isps,:)
                            enddo 
                            cycle 
                        endif 
                    
                        do ips=1,nps
                            psd_norm_fact(ips) = maxval(mpsd(isps,ips,:))
                            
                            psd_norm(ips,:) = mpsd(isps,ips,:) / psd_norm_fact(ips)
                            psd_pr_norm(ips) = mpsd_pr(isps,ips) / psd_norm_fact(ips)
                            dpsd_norm(ips,:) = dmpsd(isps,ips,:) / psd_norm_fact(ips)
                            psd_rain_norm(ips,:) = mpsd_rain(isps,ips,:) / psd_norm_fact(ips) *dt /dt_save
                        enddo 
                        
                        if (.not.psd_impfull) then 
                            call psd_implicit_all_v2( &
                                & nz,nsp_sld,nps,nflx_psd &! in
                                & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in
                                & ,trans &! in
                                & ,psd_norm,psd_pr_norm,ps,dps,dpsd_norm,psd_rain_norm &! in  
                                & ,chrsld(isps) &! in 
                                & ,flgback,flx_max_max &! inout
                                & ,psdx_norm,flx_psd_norm &! out
                                & )
                        else
                            DV(:) = flx_sld(isps, 4 + isps,:)*mv(isps)*1d-6*dt 
                            call psd_implicit_all_v4( &
                                & nz,nsp_sld,nps,nflx_psd &! in
                                & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in 
                                & ,incld_rough,roughref(isps) &! in
                                & ,trans &! in
                                & ,psd_norm,psd_pr_norm,ps,dps,dpsd_norm,psd_rain_norm,DV,psd_norm_fact &! in  
                                & ,chrsld(isps) &! in 
                                & ,flgback,flx_max_max &! inout
                                & ,psdx_norm,flx_psd_norm &! out
                                & )
                        endif 

                        if ( flgback .or. psd_error_flg) then 
                            print *, '*** error raised after PSD calc. for ',chrsld(isps)
                            print *, '*** escape from do-loop'
                            exit
                        endif 
                        
                        ! if (flgback) then 
                            ! flgback = .false. 
                            ! flgreducedt = .true.
                            ! psd = psd_old
                            ! mpsd = mpsd_old
                            ! poro = poroprev
                            ! torg = torgprev
                            ! tora = toraprev
                            ! disp = dispprev
                            ! v = vprev
                            ! hr = hrprev
                            ! w = wprev
                            ! call calcupwindscheme(  &
                                ! up,dwn,cnr,adf & ! output 
                                ! ,w,nz   & ! input &
                                ! )
                            ! dt = dt/1d1
                            ! go to 100
                        ! endif 
                            
                        do ips=1,nps
                            mpsdx(isps,ips,:) = psdx_norm(ips,:)*psd_norm_fact(ips)
                            flx_mpsd(isps,ips,:,:) = flx_psd_norm(ips,:,:)*psd_norm_fact(ips)
                        enddo 
                    
                    enddo 

                    if ( flgback .or. psd_error_flg) then 
                        print *,' *** because of PSD(norm) error, returning to do while loop with reduced dt'
                        if (.not.psd_loop) exit
                        flgback = .false.
                        psd_error_flg = .false.
                        dt = dt/10d0
                        psd = psd_save
                        mpsd = mpsd_save_2
                        cycle
                    endif 
                
                else 
                    
                    do ips=1,nps
                        psd_norm_fact(ips) = maxval(psd(ips,:))
                        
                        psd_norm(ips,:) = psd(ips,:) / psd_norm_fact(ips)
                        psd_pr_norm(ips) = psd_pr(ips) / psd_norm_fact(ips)
                        dpsd_norm(ips,:) = dpsd(ips,:) / psd_norm_fact(ips)
                        psd_rain_norm(ips,:) = psd_rain(ips,:) / psd_norm_fact(ips) *dt /dt_save
                    enddo 
                    
                    flx_max_max = 0d0
                    
                    if (.not.psd_impfull) then 
                        call psd_implicit_all_v2( &
                            & nz,nsp_sld,nps,nflx_psd &! in
                            & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in
                            & ,trans &! in
                            & ,psd_norm,psd_pr_norm,ps,dps,dpsd_norm,psd_rain_norm &! in    
                            & ,' blk ' &! in 
                            & ,flgback,flx_max_max &! inout
                            & ,psdx_norm,flx_psd_norm &! out
                            & )
                    else
                        call psd_implicit_all_v4( &
                            & nz,nsp_sld,nps,nflx_psd &! in
                            & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in 
                            & ,incld_rough,roughref_b &! in
                            & ,trans &! in
                            & ,psd_norm,psd_pr_norm,ps,dps,dpsd_norm,psd_rain_norm,DV,psd_norm_fact &! in    
                            & ,' blk ' &! in 
                            & ,flgback,flx_max_max &! inout
                            & ,psdx_norm,flx_psd_norm &! out
                            & )
                    endif 

                    if ( flgback .or. psd_error_flg) then 
                        print *,' *** because of PSD(norm) error, returning to do while loop with reduced dt'
                        if (.not.psd_loop) exit
                        flgback = .false.
                        psd_error_flg = .false.
                        dt = dt/10d0
                        psd = psd_save
                        mpsd = mpsd_save_2
                        cycle
                    endif 
                    
                    do ips=1,nps
                        psdx(ips,:) = psdx_norm(ips,:)*psd_norm_fact(ips)
                        flx_psd(ips,:,:) = flx_psd_norm(ips,:,:)*psd_norm_fact(ips)
                    enddo 
                
                endif 
            else
                
                if (do_psd_full) then 
                
                    flx_max_max = 0d0
                
                    do isps = 1, nsp_sld
                        
                        if ( &
                            & trim(adjustl(precstyle(isps))) == 'decay' &! solid species not related to SA              
                            & .or. ( psd_enable_skip .and. any(chrsld_nopsd == chrsld(isps)) ) &!case when not-tracking PSDs for fast reacting minerals (SA not matter?)
                            & ) then           
                            flx_mpsd(isps,:,:,:) = 0d0
                            do iz=1,nz
                                mpsdx(isps,:,iz) = mpsd_pr(isps,:)
                            enddo 
                            cycle 
                        endif 
            
                        if (display) then 
                            print *
                            print *, '<'//trim(adjustl(chrsld(isps)))//'>'
                        endif 
                        
                        psd = mpsd(isps,:,:)
                        psd_pr = mpsd_pr(isps,:)
                        dpsd = dmpsd(isps,:,:)
                        psd_rain = mpsd_rain(isps,:,:)
                        
                        
                        if (.not.psd_impfull) then 
                            call psd_implicit_all_v2( &
                                & nz,nsp_sld,nps,nflx_psd &! in
                                & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in
                                & ,trans &! in
                                & ,psd,psd_pr,ps,dps,dpsd,psd_rain &! in    
                                & ,chrsld(isps) &! in 
                                & ,flgback,flx_max_max &! inout
                                & ,psdx,flx_psd &! out
                                & )
                        else
                            DV(:) = flx_sld(isps, 4 + isps,:)*mv(isps)*1d-6*dt  
                            psd_norm_fact = 1d0
                            call psd_implicit_all_v4( &
                                & nz,nsp_sld,nps,nflx_psd &! in
                                & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in 
                                & ,incld_rough,roughref(isps) &! in
                                & ,trans &! in
                                & ,psd,psd_pr,ps,dps,dpsd,psd_rain,DV,psd_norm_fact &! in    
                                & ,chrsld(isps) &! in
                                & ,flgback,flx_max_max &! inout
                                & ,psdx,flx_psd &! out
                                & )
                        endif 
                        
                        print *,'flx_max_max',flx_max_max
                
                        if (flgback) then 
                            flgback = .false. 
                            flgreducedt = .true.
                            psd = psd_old
                            mpsd = mpsd_old
                            poro = poroprev
                            torg = torgprev
                            tora = toraprev
                            disp = dispprev
                            v = vprev
                            hr = hrprev
                            w = wprev
                            call calcupwindscheme(  &
                                up,dwn,cnr,adf & ! output 
                                ,w,nz   & ! input &
                                )
                            dt = dt/1d1
                            go to 100
                        endif 
                        
                        mpsdx(isps,:,:) = psdx
                        flx_mpsd(isps,:,:,:) = flx_psd
                    
                    enddo
                
                else 
                
                    flx_max_max = 0d0
                    
                    if (.not.psd_impfull) then 
                        call psd_implicit_all_v2( &
                            & nz,nsp_sld,nps,nflx_psd &! in
                            & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in
                            & ,trans &! in
                            & ,psd,psd_pr,ps,dps,dpsd,psd_rain &! in    
                            & ,' blk ' &! in 
                            & ,flgback,flx_max_max &! inout
                            & ,psdx,flx_psd &! out
                            & )
                    else
                        psd_norm_fact = 1d0
                        call psd_implicit_all_v4( &
                            & nz,nsp_sld,nps,nflx_psd &! in
                            & ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev &! in 
                            & ,incld_rough,roughref_b &! in
                            & ,trans &! in
                            & ,psd,psd_pr,ps,dps,dpsd,psd_rain,DV,psd_norm_fact &! in    
                            & ,' blk ' &! in 
                            & ,flgback,flx_max_max &! inout
                            & ,psdx,flx_psd &! out
                            & )
                    endif 
                
                    if (flgback) then 
                        flgback = .false. 
                        flgreducedt = .true.
                        psd = psd_old
                        mpsd = mpsd_old
                        poro = poroprev
                        torg = torgprev
                        tora = toraprev
                        disp = dispprev
                        v = vprev
                        hr = hrprev
                        w = wprev
                        call calcupwindscheme(  &
                            up,dwn,cnr,adf & ! output 
                            ,w,nz   & ! input &
                            )
                        dt = dt/1d1
                        go to 100
                    endif    
                
                endif 
                
            endif 
            
            if (do_psd_full) then 
            
                mpsd = mpsdx 

                if (any(isnan(mpsd))) then 
                    print *, 'nan in mpsd'
                    stop
                endif 
                if (any(mpsd<0d0)) then 
                    print *, 'negative in mpsd'
                    stop
                    error_psd = 0d0
                    do isps = 1, nsp_sld
                        do iz = 1, nz
                            do ips=1,nps
                                if (mpsd(isps,ips,iz)<0d0) then 
                                    error_psd = min(error_psd,mpsd(isps,ips,iz))
                                    mpsd(isps,ips,iz) = 0d0
                                endif 
                            enddo 
                        enddo 
                    enddo 
                    if (abs(error_psd/maxval(mpsd)) > tol) then 
                        print *, 'negative mpsd'
                        flgback = .false. 
                        flgreducedt = .true.
                        psd = psd_old
                        mpsd = mpsd_old
                        poro = poroprev
                        torg = torgprev
                        tora = toraprev
                        disp = dispprev
                        v = vprev
                        hr = hrprev
                        w = wprev
                        call calcupwindscheme(  &
                            up,dwn,cnr,adf & ! output 
                            ,w,nz   & ! input &
                            )
                        dt = dt/1d1
                        go to 100
                    endif 
                endif 

                if (psd_lim_min) then    
                    where (mpsd < psd_th_0)  mpsd = psd_th_0
                endif 
            
            else 
            
                psd = psdx 

                if (any(isnan(psd))) then 
                    print *, 'nan in psd'
                    stop
                endif 
                if (any(psd<0d0)) then 
                    print *, 'negative in psd'
                    stop
                    error_psd = 0d0
                    do iz = 1, nz
                        do ips=1,nps
                            if (psd(ips,iz)<0d0) then 
                                error_psd = min(error_psd,psd(ips,iz))
                                psd(ips,iz) = 0d0
                            endif 
                        enddo 
                    enddo 
                    if (abs(error_psd/maxval(psd)) > tol) then 
                        print *, 'negative psd'
                        flgback = .false. 
                        flgreducedt = .true.
                        psd = psd_old
                        mpsd = mpsd_old
                        poro = poroprev
                        torg = torgprev
                        tora = toraprev
                        disp = dispprev
                        v = vprev
                        hr = hrprev
                        w = wprev
                        call calcupwindscheme(  &
                            up,dwn,cnr,adf & ! output 
                            ,w,nz   & ! input &
                            )
                        dt = dt/1d1
                        go to 100
                    endif 
                endif 
                if (psd_lim_min) then         
                    where (psd < psd_th_0)  psd = psd_th_0
                endif 
            endif 
            
            
            if ( flgback .or. psd_error_flg) then 
                print *, ' *** somehow error is still raised comming from PSD or PBE calc.'
                print *, ' *** retuning the initial do while loop with reduced dt'
                flgback = .false.
                psd_error_flg = .false.
                dt = dt/10d0
                psd = psd_save
                mpsd = mpsd_save_2
                cycle
            endif 
            
            time_pbe = time_pbe + dt
            
            print '(a,1x,f11.7)',' ---- PSD completeness [%]: ',100d0*time_pbe/dt_save
            
            if (time_pbe>= dt_save) then 
                print *, ' time within PSD+PBE seems to reach dt in main loop'
                print *, ' exiting the do while loop!!'
                exit
            endif 
            
            
            enddo ! end of do while loop for psd

            dt = dt_save
                    
            if (.not. psd_loop) then
                if (flgback .or. psd_error_flg) then 
                    flgback = .false. 
                    psd_error_flg = .false. 
                    flgreducedt = .true.
                    psd = psd_old
                    mpsd = mpsd_old
                    poro = poroprev
                    torg = torgprev
                    tora = toraprev
                    disp = dispprev
                    v = vprev
                    hr = hrprev
                    w = wprev
                    call calcupwindscheme(  &
                        up,dwn,cnr,adf & ! output 
                        ,w,nz   & ! input &
                        )
                    dt = dt/1d1
                    go to 100
                endif 
            endif 
            
            if (display) then 
                print *, '-- ending PSD'
                print *
            endif 
        
        else 
            psd = 0d0
            mpsd = 0d0
        endif 
            
        ! if doing psd SA is calculated reflecting psd
        if (do_psd) then 
            if (do_psd_full) then  
                do isps=1,nsp_sld
                    do iz=1,nz
                        ssa(isps,iz) = sum( 4d0*pi*(10d0**ps(:))**2d0*rough_ps(isps,:)*mpsd(isps,:,iz)*dps(:))
                        ssav(isps,iz) = sum( 3d0/(10d0**ps(:))*rough_ps(isps,:)*mpsd(isps,:,iz)*dps(:))
                        ssv(isps,iz) = sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*mpsd(isps,:,iz)*dps(:))   ! solid m3/bulk m3
                    enddo 
                enddo
            else 
                do iz=1,nz
                    ssa(:,iz) = sum( 4d0*pi*(10d0**ps(:))**2d0*rough_ps_b(:)*psd(:,iz)*dps(:))
                    ssav(:,iz) = sum( 3d0/(10d0**ps(:))*rough_ps_b(:)*psd(:,iz)*dps(:))
                    ssv(:,iz) = sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*psd(:,iz)*dps(:))
                enddo 
            endif 
            ! hr = ssa*(1-poro)/poro ! converting m2/sld-m3 to m2/pore-m3
            ! hr = ssa
            do isps=1,nsp_sld
                if (do_psd_full) then 
                    ! hr(isps,:) = ssa(isps,:)/poro/(msldx(isps,:)*mv(isps)*1d-6) ! so that poro * hr * mv * msld becomes porosity independent
                    ! hr(isps,:) = ssa(isps,:) ! so that poro * hr * mv * msld becomes porosity independent
                    hr(isps,:) = ssa(isps,:)/ssv(isps,:)/poro    !!! better for conversion?
                    ! hr(isps,:) = ssav(isps,:)/poro                ! should be this? noted in 4-24-2023
                else 
                    hr(isps,:) = ssa(isps,:)/poro ! so that poro * hr * mv * msld becomes porosity independent
                endif 
            enddo
        endif 
        
        if (any(poro < 1d-10)) then 
            print *, '***| too small porosity: going to end sim as likely ending up crogging '
            print *, '***| ... and no more reasonable simulation ... ! '
            stop
        endif 
        
        ! calculating specific surface area (m2/g)
        ssas = 0d0
        do isps=1,nsp_sld
            if (do_psd_full) then 
                ssas(isps,:) = hr(isps,:)*poro*ssv(isps,:)/ ( msldx(isps,:)*mwt(isps) )
                ! ssas(isps,:) = hr(isps,:)*poro/( msldx(isps,:)*mwt(isps) )
            else
                ssas(isps,:) = hr(isps,:)*poro/( msldx(isps,:)*mwt(isps) )
            endif 
        enddo
        
        ! activity coefficient for H+
        gamma = 1d0
        if (act_ON) then         
            rcharge = 1d0
            call calc_gamma_davies(  &
                & nz,iosx,tc,rcharge &
                & ,gamma_tmp,dgamma_dios_tmp &
                & )
            gamma = gamma_tmp
        endif 

        if (display  .and. (.not. display_lim)) then 
            write(chrfmt,'(i0)') nz_disp
            chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
            
            print *
            print *,' [concs] '
            print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
            if (nsp_aq>0) then 
                print *,' < aq species >'
                do ispa = 1, nsp_aq
                    ! print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/nz_disp)
                    print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maqx(ispa,iz)*maqft(ispa,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
            if (nsp_sld>0) then 
                print *,' < sld species >'
                do isps = 1, nsp_sld
                    print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
            if (nsp_gas>0) then 
                print *,' < gas species >'
                do ispg = 1, nsp_gas
                    print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
            
            print *
            print *,' [saturation & pH] '
            if (nsp_sld>0) then 
                print *,' < sld species omega >'
                do isps = 1, nsp_sld
                    print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/nz_disp)
                enddo 
            endif 
            print *,' < pH >'
            print trim(adjustl(chrfmt)), 'ph', (-log10(gamma(iz)*prox(iz)),iz=1,nz, nz/nz_disp)
            
            
            
            write(chrfmt,'(i0)') nflx
            chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'
            
            print *
            print *,' [fluxes] '
            print trim(adjustl(chrfmt)),' ',(chrflx(iflx),iflx=1,nflx)
            
            write(chrfmt,'(i0)') nflx
            chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
            if (nsp_aq>0) then 
                print *,' < aq species >'
                do ispa = 1, nsp_aq
                    print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
                enddo 
            endif 
            if (nsp_sld>0) then 
                print *,' < sld species >'
                do isps = 1, nsp_sld
                    print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
                enddo 
            endif 
            if (nsp_gas>0) then 
                print *,' < gas species >'
                do ispg = 1, nsp_gas
                    print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
                enddo 
            endif 
            
            if (do_psd) then 
                if (do_psd_full) then
                    write(chrfmt,'(i0)') nflx_psd
                    chrfmt = '(a5,1x,a5,'//trim(adjustl(chrfmt))//'(1x,a11))'

                    print *
                    print *,' [fluxes -- PSD] '
                    print trim(adjustl(chrfmt)),'sld','rad','tflx','adv','dif','rain','rxn','res'
                    
                    do isps = 1, nsp_sld
                        write(chrfmt,'(i0)') nflx_psd
                        chrfmt = '(a5,1x,f5.2,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
                        do ips = 1, nps
                            print trim(adjustl(chrfmt)), chrsld(isps), ps(ips) &
                                & ,(sum(flx_mpsd(isps,ips,iflx,:)*dz(:)),iflx=1,nflx_psd)
                        enddo 
                        print *
                    enddo 
                else 
                    write(chrfmt,'(i0)') nflx_psd
                    chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'

                    print *
                    print *,' [fluxes -- PSD] '
                    print trim(adjustl(chrfmt)),'rad','tflx','adv','dif','rain','rxn','res'

                    write(chrfmt,'(i0)') nflx_psd
                    chrfmt = '(f5.2,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
                    do ips = 1, nps
                        print trim(adjustl(chrfmt)), ps(ips), (sum(flx_psd(ips,iflx,:)*dz(:)),iflx=1,nflx_psd)
                    enddo 
                endif 
            endif
            
    ! #ifdef disp_lim
            if (display_lim_in) display_lim = .true.
    ! #endif 
            
        endif 

        ! stop
    ! #ifdef lim_minsld
        if (lim_minsld_in) then 
            ! where (msldx < 1d-20)  msldx = 1d-20
            do isps = 1, nsp_sld
                ! if ( trim(adjustl(precstyle(isps))) /= 'decay' ) cycle 
                ! if ( trim(adjustl(precstyle(isps))) /= 'decay' ) cycle 
                ! if ( trim(adjustl(chrsld(isps))) /= 'cc' ) cycle 
                ! if ( trim(adjustl(chrsld(isps))) /= 'dlm' ) cycle 
                ! if ( trim(adjustl(chrsld(isps))) /= 'arg' ) cycle 
                ! if ( trim(adjustl(chrsld(isps))) /= 'gps' ) cycle 
                do iz=1,nz
                    if (msldx(isps,iz) < minsld(isps)) msldx(isps,iz) = minsld(isps)
                enddo 
            enddo 
            if (ads_ON) then
                do ispa=1,nsp_aq
                    selectcase(trim(adjustl(chraq(ispa))))
                        case('ca','mg','k','na')
                            do iz=1,nz
                                if (maqfads(ispa,iz) < minmaqads(ispa)) maqfads(ispa,iz) = minmaqads(ispa)
                            enddo 
                        case default
                    endselect
                enddo 
            endif 
        endif 
    ! #endif 
        
        mgas = mgasx
        maq = maqx
        msld = msldx
        
        pro = prox
        ios = iosx
        ! so4fprev = so4f
        maqft_prev = maqft
        maqfads_prev = maqfads
        
        mblk = mblkx

        it = it + 1
        time = time + dt
        count_dtunchanged = count_dtunchanged + 1
        
        do iz = 1, nz
            ! rho_grain_z(iz) = sum(msldx(:,iz)*mwt(:)*1d-6)
            ! sldvolfrac(iz) = sum(msldx(:,iz)*mv(:)*1d-6)
            ! accounting for blk soil
            rho_grain_z(iz) = sum(msldx(:,iz)*mwt(:)*1d-6) + mblkx(iz)*mwtblk*1d-6
            sldvolfrac(iz) = sum(msldx(:,iz)*mv(:)*1d-6) + mblkx(iz)*mvblk*1d-6
            
            if (msldunit=='blk') then 
                rho_grain_z(iz) = rho_grain_z(iz) / ( 1d0 - poro(iz) )
                sldvolfrac(iz) = sldvolfrac(iz) / ( 1d0 - poro(iz) )
            endif 
                
        enddo 
        
        ! calculating cec
        cec     = 0d0
        proxads = 0d0
        bs      = 0d0
        do iz=1,nz
            ucvsld1 = 1d0
            if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
            
            do isps=1,nsp_sld_all    
                if (keqcec_all(isps) == 0d0) cycle
                
                if (any(chrsld_all(isps) == chrsld)) then 
                    proxads(iz) = proxads(iz) + ( &
                        & + keqcec_all(isps) &
                        & * msldx(findloc(chrsld,chrsld_all(isps),dim=1),iz) &
                        & * msldf_loc(isps,iz) &
                        & * beta_loc(isps,iz) &
                        & )
                        
                    cec(iz) = cec(iz) + ( &
                        & + keqcec_all(isps) &
                        & * msldx(findloc(chrsld,chrsld_all(isps),dim=1),iz) &
                        & )
                endif 
                
            enddo 
            
            proxads(iz)  =  proxads(iz) *1d5/ucvsld1/(rho_grain_z(iz)*1d6)
            cec(iz)      =  cec(iz)     *1d5/ucvsld1/(rho_grain_z(iz)*1d6)
            bs(iz)       =  proxads(iz)/cec(iz)
            
        enddo 
        
        do iz=1,nz
            ucvsld1 = 1d0
            if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
            
            do ispa=1,nsp_aq
                cecaq(ispa,iz)  = maqx(ispa,iz)*maqfads(ispa,iz)*1d5/ucvsld1/(rho_grain_z(iz)*1d6) ! converting mol/m3 to mol/g then to cmol/kg
                cecaqr(ispa,iz) = cecaq(ispa,iz)*base_charge(ispa) / cec(iz) 
                cecaqwt(ispa,iz)  = mwtaq(ispa)*10d0*cecaq(ispa,iz) ! converting cmol/kg to ppm
            enddo 
        enddo 
        
        ! calculating volume weighted average surface area 
        do iz=1,nz
            hrb(iz) = sum( hr(:,iz)* msldx(:,iz)*mv(:)*1d-6) / sum(msldx(:,iz)*mv(:)*1d-6)
            ssab(iz) = sum( ssa(:,iz)* msldx(:,iz)*mv(:)*1d-6) / sum(msldx(:,iz)*mv(:)*1d-6)
        enddo 
        
        
        do iflx=1,nflx
            do isps=1,nsp_sld 
                int_flx_sld(isps,iflx) = int_flx_sld(isps,iflx) + sum(flx_sld(isps,iflx,:)*dz(:))*dt
            enddo 
            
            do ispa=1,nsp_aq 
                int_flx_aq(ispa,iflx) = int_flx_aq(ispa,iflx) + sum(flx_aq(ispa,iflx,:)*dz(:))*dt
            enddo 
            
            do ispg=1,nsp_gas 
                int_flx_gas(ispg,iflx) = int_flx_gas(ispg,iflx) + sum(flx_gas(ispg,iflx,:)*dz(:))*dt
            enddo 
            
            do ico2=1,6 
                if (ico2 .le. 4) then 
                    int_flx_co2sp(ico2,iflx) = int_flx_co2sp(ico2,iflx) + sum(flx_co2sp(ico2,iflx,:)*dz(:))*dt
                elseif (ico2 .eq. 5) then 
                    int_flx_co2sp(ico2,iflx) = int_flx_co2sp(ico2,iflx) + ( &
                        & sum(flx_co2sp(2,iflx,:)*dz(:))+sum(flx_co2sp(3,iflx,:)*dz(:))+sum(flx_co2sp(4,iflx,:)*dz(:)) &
                        & ) * dt
                elseif (ico2 .eq. 6) then 
                    int_flx_co2sp(ico2,iflx) = int_flx_co2sp(ico2,iflx) + ( &
                        & sum(flx_co2sp(3,iflx,:)*dz(:))+2d0*sum(flx_co2sp(4,iflx,:)*dz(:)) &
                        & ) *dt
                endif 
            enddo 
        enddo 
        
        do iz=1,nz
            int_ph(iz) = int_ph(iz) + sum(gamma(1:iz)*prox(1:iz)*dz(1:iz))/z(iz) * dt
        enddo 

        if (time >= savetime) then 
            
            open(isldprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_sld-save.txt', status='replace')
            open(igasprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_gas-save.txt', status='replace')
            open(iaqprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq-save.txt', status='replace')
            open(ibsd, file=trim(adjustl(profdir))//'/'  &
                & //'bsd-save.txt', status='replace')
            open(isa,file=trim(adjustl(profdir))//'/' &
                & //'sa-save.txt', status='replace')
                
            write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
            write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
            write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
            write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' , ' w[m/yr] '  &
                & , ' vol[m3/m3] ',' dens[g/cm3] ', ' blk[wt%] ',' cec[cmol/kg] ',' time '
            write(isa,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '

            do iz = 1, Nz
                ucvsld1 = 1d0
                if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
                write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
                write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
                write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
                write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz) &
                    & ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6), cec(iz), time
                write(isa,*) z(iz),(hr(isps,iz),isps = 1, nsp_sld),time
            end do

            close(isldprof)
            close(iaqprof)
            close(igasprof)
            close(ibsd)
            close(isa)
            
            if (do_psd) then 
                if (do_psd_full) then 
                    do isps=1,nsp_sld
                        open(ipsd, file=trim(adjustl(profdir))//'/'  &
                            & //'psd_'//trim(adjustl(chrsld(isps)))//'-save.txt', status='replace')
                        write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
                        do iz = 1, Nz
                            write(ipsd,*) z(iz), (mpsd(isps,ips,iz),ips=1,nps), time 
                        end do
                        close(ipsd)
                    enddo 
                else 
                    open(ipsd, file=trim(adjustl(profdir))//'/'  &
                        & //'psd-save.txt', status='replace')
                    write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
                    do iz = 1, Nz
                        write(ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
                    end do
                    close(ipsd)
                endif 
            endif 
            
            savetime = savetime + dsavetime
            
        endif 
        
        flx_recorded = .false.
        
        if (time>=rectime_prof(irec_prof+1)) then
            write(chr,'(i3.3)') irec_prof+1
            
            
            print_cb = .true. 
            print_loc = trim(adjustl(profdir))//'/' &
                & //'charge_balance-'//chr//'.txt'

                
            call calc_pH_v7_4( &
                & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                & ,poro,sat,tc &! input 
                & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
                & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 &! input
                & ,keqaq_oxa,keqaq_cl &! input
                & ,print_cb,print_loc,z,act_ON &! input 
                & ,dprodmaq_all,dprodmgas_all &! output
                & ,iosx,diosdmaq_all,diosdmgas_all &! output
                & ,prox,ph_error,ph_iter &! output
                & ) 

            ! activity coefficient for H+
            gamma = 1d0
            if (act_ON) then         
                rcharge = 1d0
                call calc_gamma_davies(  &
                    & nz,iosx,tc,rcharge &
                    & ,gamma_tmp,dgamma_dios_tmp &
                    & )
                gamma = gamma_tmp
            endif 
            
            open(isldprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_sld-'//chr//'.txt', status='replace')
            open(isldprof2,file=trim(adjustl(profdir))//'/' &
                & //'prof_sld(wt%)-'//chr//'.txt', status='replace')
            open(isldprof3,file=trim(adjustl(profdir))//'/' &
                & //'prof_sld(v%)-'//chr//'.txt', status='replace')
            open(isldsat,file=trim(adjustl(profdir))//'/' &
                & //'sat_sld-'//chr//'.txt', status='replace')
            open(igasprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_gas-'//chr//'.txt', status='replace')
            open(iaqprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq-'//chr//'.txt', status='replace')
            open(iaqprof2,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq(tot)-'//chr//'.txt', status='replace')
            open(iaqprof3,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq(ads)-'//chr//'.txt', status='replace')
            open(iaqprof4,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq(ads%cec)-'//chr//'.txt', status='replace')
            open(iaqprof5,file=trim(adjustl(profdir))//'/' &
                & //'prof_ex(tot)-'//chr//'.txt', status='replace')
            open(iaqprof6,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq(adsppm)-'//chr//'.txt', status='replace')
            open(ibsd, file=trim(adjustl(profdir))//'/'  &
                & //'bsd-'//chr//'.txt', status='replace')
            open(irate, file=trim(adjustl(profdir))//'/'  &
                & //'rate-'//chr//'.txt', status='replace')
            open(isa,file=trim(adjustl(profdir))//'/' &
                & //'sa-'//chr//'.txt', status='replace')
            open(isa2,file=trim(adjustl(profdir))//'/' &
                & //'ssa-'//chr//'.txt', status='replace')
                
            write(chrfmt,'(i0)') nsp_sld+2
            chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
            write(isldprof,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
            write(isldprof2,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
            write(isldprof3,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
            write(isldsat,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
            write(chrfmt,'(i0)') nsp_aq+3
            chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
            write(iaqprof,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'ph','time'
            write(iaqprof2,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'aph','time'
            write(iaqprof3,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'h','time'
            write(iaqprof4,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'h','time'
            write(iaqprof5,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'h','time'
            write(iaqprof6,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'h','time'
            write(chrfmt,'(i0)') nsp_gas+2
            chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
            write(igasprof,trim(adjustl(chrfmt))) 'z',(chrgas(isps),isps=1,nsp_gas),'time'
            write(chrfmt,'(i0)') 11
            chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a12))'
            write(ibsd,trim(adjustl(chrfmt))) 'z','poro', 'sat', 'v[m/yr]', 'm2/m3' , 'w[m/yr]' &
                & , 'vol[m3/m3]','dens[g/cm3]','blk[wt%]','cec[cmol/kg]','time'
            write(chrfmt,'(i0)') 2 + nsp_sld + nrxn_ext
            chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
            write(irate,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),(chrrxn_ext(irxn),irxn=1,nrxn_ext),'time'
            write(chrfmt,'(i0)') nsp_sld+2
            chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
            write(isa,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
            write(isa2,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'

            do iz = 1, Nz
                ucvsld1 = 1d0
                if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
                
                write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
                write(isldprof2,*) z(iz),(msldx(isps,iz)*mwt(isps)*1d2/ucvsld1/(rho_grain_z(iz)*1d6),isps = 1, nsp_sld),time
                write(isldprof3,*) z(iz),(msldx(isps,iz)*mv(isps)/ucvsld1*1d-6*1d2,isps = 1, nsp_sld),time
                write(isldsat,*) z(iz),(omega(isps,iz),isps = 1, nsp_sld),time
                write(igasprof,*) z(iz),(mgasx(ispg,iz),ispg = 1, nsp_gas),time
                write(iaqprof,*) z(iz),(maqx(ispa,iz),ispa = 1, nsp_aq),-log10(prox(iz)),time
                write(iaqprof2,*) z(iz),(maqx(ispa,iz)*maqft(ispa,iz),ispa = 1, nsp_aq),-log10(gamma(iz)*prox(iz)),time
                write(iaqprof3,*) z(iz),(cecaq(ispa,iz),ispa = 1, nsp_aq) ,proxads(iz),time
                write(iaqprof4,*) z(iz),(cecaqr(ispa,iz)*1d2,ispa = 1, nsp_aq) ,bs(iz)*1d2,time
                write(iaqprof5,*) z(iz),(poro(iz)*sat(iz)*1d3*maqx(ispa,iz)*maqft(ispa,iz) &
                    & + maqx(ispa,iz)*maqfads(ispa,iz),ispa = 1, nsp_aq) &
                    & ,poro(iz)*sat(iz)*1d3*prox(iz) + proxads(iz) / (1d5/ucvsld1/(rho_grain_z(iz)*1d6)) ,time
                write(iaqprof6,*) z(iz),(cecaqwt(ispa,iz),ispa = 1, nsp_aq) ,1d0*10d0*proxads(iz),time
                write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz)  &
                    & ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6),cec(iz),time
                write(irate,*) z(iz), (rxnsld(isps,iz),isps=1,nsp_sld),(rxnext(irxn,iz),irxn=1,nrxn_ext), time 
                write(isa,*) z(iz),(hr(isps,iz),isps = 1, nsp_sld),time
                write(isa2,*) z(iz),(ssas(isps,iz),isps = 1, nsp_sld),time
            end do

            close(isldprof)
            close(isldprof2)
            close(isldprof3)
            close(isldsat)
            close(iaqprof)
            close(iaqprof2)
            close(iaqprof3)
            close(iaqprof4)
            close(iaqprof5)
            close(iaqprof6)
            close(igasprof)
            close(ibsd)
            close(irate)
            close(isa)
            close(isa2)
            
            if (do_psd) then 
                
                if (do_psd_full) then 
                    
                    do isps = 1, nsp_sld 
                        open(ipsd, file=trim(adjustl(profdir))//'/'  &
                            & //'psd_'//trim(adjustl(chrsld(isps)))//'-'//chr//'.txt', status='replace')
                        open(ipsdv, file=trim(adjustl(profdir))//'/'  &
                            & //'psd_'//trim(adjustl(chrsld(isps)))//'(v%)-'//chr//'.txt', status='replace')
                        open(ipsds, file=trim(adjustl(profdir))//'/'  &
                            & //'psd_'//trim(adjustl(chrsld(isps)))//'(SA%)-'//chr//'.txt', status='replace')
                        open(ipsdflx, file=trim(adjustl(flxdir))//'/'  &
                            & //'flx_psd_'//trim(adjustl(chrsld(isps)))//'-'//chr//'.txt', status='replace')
                        
                        write(chrfmt,'(i0)') nps
                        chrfmt = '(1x,a16,'//trim(adjustl(chrfmt))//'(1x,f11.6),1x,a5)'
                        write(ipsd,trim(adjustl(chrfmt))) 'z[m]\log10(r[m])',(ps(ips),ips=1,nps),'time'
                        write(ipsdv,trim(adjustl(chrfmt))) 'z[m]\log10(r[m])',(ps(ips),ips=1,nps),'time'
                        write(ipsds,trim(adjustl(chrfmt))) 'z[m]\log10(r[m])',(ps(ips),ips=1,nps),'time'
                        write(chrfmt,'(i0)') nflx_psd
                        chrfmt = '(1x,a5,1x,a16,'//trim(adjustl(chrfmt))//'(1x,a11))'
                        write(ipsdflx,trim(adjustl(chrfmt))) 'time','log10(r[m])\flx','tflx','adv','dif','rain','rxn','res'
                        
                        do iz = 1, Nz
                            ucvsld2 = 1d0 - poro(iz)
                            if (msldunit == 'blk') ucvsld2 = 1d0
                            
                            write(ipsd,*) z(iz), (mpsd(isps,ips,iz),ips=1,nps), time 
                            write(ipsdv,*) z(iz), (4d0/3d0*pi*(10d0**ps(ips))**3d0*mpsd(isps,ips,iz)*dps(ips) &
                                ! & /sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*psd(:,iz)*dps(:))  * 1d2 &
                                & / (  msld(isps,iz)*mv(isps)*1d-6  ) * 1d2 &
                                & /ucvsld2  &
                                & ,ips=1,nps), time 
                            write(ipsds,*) z(iz), (4d0*pi*(10d0**ps(ips))**2d0 &
                                & *rough_ps(isps,ips)*mpsd(isps,ips,iz)*dps(ips) &
                                ! & /sum( 4d0*pi*(10d0**ps(:))**2d0*rough_c0(isps)*(10d0**ps(:))**rough_c1(isps)*psd(:,iz)*dps(:))  * 1d2 &
                                & / ssa(isps,iz)  * 1d2 &
                                ! & /( poro(iz)/(1d0 - poro(iz) ))  &
                                & ,ips=1,nps), time 
                        end do
                        
                        do ips=1,nps
                            write(ipsdflx,*) time,ps(ips), (sum(flx_mpsd(isps,ips,iflx,:)*dz(:)),iflx=1,nflx_psd)
                        enddo 
                        
                        close(ipsd)
                        close(ipsdv)
                        close(ipsds)
                        close(ipsdflx)
                    enddo 
                else 
                    
                    open(ipsd, file=trim(adjustl(profdir))//'/'  &
                        & //'psd-'//chr//'.txt', status='replace')
                    open(ipsdv, file=trim(adjustl(profdir))//'/'  &
                        & //'psd(v%)-'//chr//'.txt', status='replace')
                    open(ipsds, file=trim(adjustl(profdir))//'/'  &
                        & //'psd(SA%)-'//chr//'.txt', status='replace')
                    open(ipsdflx, file=trim(adjustl(flxdir))//'/'  &
                        & //'flx_psd-'//chr//'.txt', status='replace')
                    
                    write(chrfmt,'(i0)') nps
                    chrfmt = '(1x,a16,'//trim(adjustl(chrfmt))//'(1x,f11.6),1x,a5)'
                    write(ipsd,trim(adjustl(chrfmt))) 'z[m]\log10(r[m])',(ps(ips),ips=1,nps),'time'
                    write(ipsdv,trim(adjustl(chrfmt))) 'z[m]\log10(r[m])',(ps(ips),ips=1,nps),'time'
                    write(ipsds,trim(adjustl(chrfmt))) 'z[m]\log10(r[m])',(ps(ips),ips=1,nps),'time'
                    write(chrfmt,'(i0)') nflx_psd
                    chrfmt = '(1x,a5,1x,a16,'//trim(adjustl(chrfmt))//'(1x,a11))'
                    write(ipsdflx,trim(adjustl(chrfmt))) 'time','log10(r[m])\flx','tflx','adv','dif','rain','rxn','res'
                    
                    do iz = 1, Nz
                        ucvsld2 = 1d0 - poro(iz)
                        if (msldunit == 'blk') ucvsld2 = 1d0
                        
                        write(ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
                        write(ipsdv,*) z(iz), (4d0/3d0*pi*(10d0**ps(ips))**3d0*psd(ips,iz)*dps(ips) &
                            ! & /sum( 4d0/3d0*pi*(10d0**ps(:))**3d0*psd(:,iz)*dps(:))  * 1d2 &
                            & / ( sum( msld(:,iz)*mv(:)*1d-6) + mblk(iz)*mvblk*1d-6 ) * 1d2 &
                            & /ucvsld2  &
                            & ,ips=1,nps), time 
                        write(ipsds,*) z(iz), (4d0*pi*(10d0**ps(ips))**2d0 &
                            & *rough_ps_b(ips)*psd(ips,iz)*dps(ips) &
                            ! & /sum( 4d0*pi*(10d0**ps(:))**2d0*rough_c0_b*(10d0**ps(:))**rough_c1*psd(:,iz)*dps(:))  * 1d2 &
                            & / ssab(iz)  * 1d2 &
                            ! & /( poro(iz)/(1d0 - poro(iz) ))  &
                            & ,ips=1,nps), time 
                    end do
                    
                    do ips=1,nps
                        write(ipsdflx,*) time,ps(ips), (sum(flx_psd(ips,iflx,:)*dz(:)),iflx=1,nflx_psd)
                    enddo 
                    
                    close(ipsd)
                    close(ipsdv)
                    close(ipsds)
                    close(ipsdflx)
                
                endif 
            
            endif 
            
            irec_prof=irec_prof+1
            
    #ifdef full_flux_report
            do isps=1,nsp_sld 
                do iz=1,nz
                    write(chriz,'(i3.3)') iz
                    open(isldflx(isps,iz), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'-'//trim(adjustl(chriz))//'.txt' &
                        & , action='write',status='old',position='append')
                    write(isldflx(isps,iz),*) time,z(iz),(sum(flx_sld(isps,iflx,1:iz)*dz(1:iz)),iflx=1,nflx)
                    close(isldflx(isps,iz))
                enddo 
            enddo 
            
            do ispa=1,nsp_aq 
                do iz= 1,nz
                    write(chriz,'(i3.3)') iz
                    open(iaqflx(ispa,iz), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'-'//trim(adjustl(chriz))//'.txt' &
                        & , action='write',status='old',position='append')
                    write(iaqflx(ispa,iz),*) time,z(iz),(sum(flx_aq(ispa,iflx,1:iz)*dz(1:iz)),iflx=1,nflx)
                    close(iaqflx(ispa,iz))
                enddo 
            enddo 
            
            do ispg=1,nsp_gas 
                do iz=1,nz
                    write(chriz,'(i3.3)') iz
                    open(igasflx(ispg,iz), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'-'//trim(adjustl(chriz))//'.txt'  &
                        & , action='write',status='old',position='append')
                    write(igasflx(ispg,iz),*) time,z(iz),(sum(flx_gas(ispg,iflx,1:iz)*dz(1:iz)),iflx=1,nflx)
                    close(igasflx(ispg,iz))
                enddo 
            enddo 
            
            do ico2=1,6 
                do iz=1,nz
                    write(chriz,'(i3.3)') iz
                    open(ico2flx(ico2,iz), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'-'//trim(adjustl(chriz))//'.txt' &
                        & , action='write',status='old',position='append')
                    if (ico2 .le. 4) then 
                        write(ico2flx(ico2,iz),*) time,z(iz),(sum(flx_co2sp(ico2,iflx,1:iz)*dz(1:iz)),iflx=1,nflx)
                    elseif (ico2 .eq. 5) then 
                        write(ico2flx(ico2,iz),*) time,z(iz) &
                            & ,(sum(flx_co2sp(2,iflx,1:iz)*dz(1:iz))+sum(flx_co2sp(3,iflx,1:iz)*dz(1:iz)) &
                            &       +sum(flx_co2sp(4,iflx,1:iz)*dz(1:iz)) &
                            & ,iflx=1,nflx)
                    elseif (ico2 .eq. 6) then 
                        write(ico2flx(ico2,iz),*) time,z(iz) &
                            & ,(sum(flx_co2sp(3,iflx,1:iz)*dz(1:iz))+2d0*sum(flx_co2sp(4,iflx,1:iz)*dz(1:iz)) &
                            & ,iflx=1,nflx)
                    endif 
                    close(ico2flx(ico2,iz))
                enddo 
            enddo 
    #else
            do isps=1,nsp_sld 
                open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
                    & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', action='write',status='old',position='append')
                write(isldflx(isps),*) time,(sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
                close(isldflx(isps))
                
                open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
                    & //'int_flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', action='write',status='old',position='append')
                write(isldflx(isps),*) time,(int_flx_sld(isps,iflx)/time,iflx=1,nflx)
                close(isldflx(isps))
            enddo 
            
            do ispa=1,nsp_aq 
                open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
                    & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', action='write',status='old',position='append')
                write(iaqflx(ispa),*) time,(sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
                close(iaqflx(ispa))
                
                open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
                    & //'int_flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', action='write',status='old',position='append')
                write(iaqflx(ispa),*) time,(int_flx_aq(ispa,iflx)/time,iflx=1,nflx)
                close(iaqflx(ispa))
            enddo 
            
            do ispg=1,nsp_gas 
                open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
                    & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', action='write',status='old',position='append')
                write(igasflx(ispg),*) time,(sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
                close(igasflx(ispg))
                
                open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
                    & //'int_flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', action='write',status='old',position='append')
                write(igasflx(ispg),*) time,(int_flx_gas(ispg,iflx)/time,iflx=1,nflx)
                close(igasflx(ispg))
            enddo 
            
            do ico2=1,6 
                open(ico2flx(ico2), file=trim(adjustl(flxdir))//'/' &
                    & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', action='write',status='old',position='append')
                if (ico2 .le. 4) then 
                    write(ico2flx(ico2),*) time,(sum(flx_co2sp(ico2,iflx,:)*dz(:)),iflx=1,nflx)
                elseif (ico2 .eq. 5) then 
                    write(ico2flx(ico2),*) time &
                        & ,(sum(flx_co2sp(2,iflx,:)*dz(:))+sum(flx_co2sp(3,iflx,:)*dz(:))+sum(flx_co2sp(4,iflx,:)*dz(:)) &
                        & ,iflx=1,nflx)
                elseif (ico2 .eq. 6) then 
                    write(ico2flx(ico2),*) time &
                        & ,(sum(flx_co2sp(3,iflx,:)*dz(:))+2d0*sum(flx_co2sp(4,iflx,:)*dz(:)) &
                        & ,iflx=1,nflx)
                endif 
                close(ico2flx(ico2))
                
                open(ico2flx(ico2), file=trim(adjustl(flxdir))//'/' &
                    & //'int_flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', action='write',status='old',position='append')
                write(ico2flx(ico2),*) time,(int_flx_co2sp(ico2,iflx)/time,iflx=1,nflx)
                close(ico2flx(ico2))
            enddo 
            
            
            open(iphint, file=trim(adjustl(flxdir))//'/'//'int_ph.txt', action='write',status='old',position='append')
            write(iphint,*) time,(-log10(int_ph(iz)/time),iz=1,nz)
            close(iphint)
            
            open(iphint2, file=trim(adjustl(flxdir))//'/'//'ph.txt', action='write',status='old',position='append')
            write(iphint2,*) time,(-log10(gamma(iz)*prox(iz)),iz=1,nz)
            close(iphint2)
            
    #endif 
            flx_recorded = .true.
            
            open(isldprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_sld-save.txt', status='replace')
            open(igasprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_gas-save.txt', status='replace')
            open(iaqprof,file=trim(adjustl(profdir))//'/' &
                & //'prof_aq-save.txt', status='replace')
            open(ibsd, file=trim(adjustl(profdir))//'/'  &
                & //'bsd-save.txt', status='replace')
            open(ipsd, file=trim(adjustl(profdir))//'/'  &
                & //'psd-save.txt', status='replace')
            open(isa,file=trim(adjustl(profdir))//'/' &
                & //'sa-save.txt', status='replace')
                
            write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
            write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
            write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
            write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' w[m/yr] '  &
                & , ' vol[m3/m3] ',' dens[g/cm3] ', ' blk[wt%] ', ' cec[cmol/kg] ' , ' time '
            write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
            write(isa,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '

            do iz = 1, Nz
                ucvsld1 = 1d0
                if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
                write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
                write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
                write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
                write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz) &
                    & ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6),cec(iz),time
                write(ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
                write(isa,*) z(iz),(hr(isps,iz),isps = 1, nsp_sld),time
            end do

            close(isldprof)
            close(iaqprof)
            close(igasprof)
            close(ibsd)
            close(ipsd)
            close(isa)
            
            if (do_psd) then 
                if (do_psd_full) then 
                    do isps=1,nsp_sld
                        open(ipsd, file=trim(adjustl(profdir))//'/'  &
                            & //'psd_'//trim(adjustl(chrsld(isps)))//'-save.txt', status='replace')
                        write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
                        do iz = 1, Nz
                            write(ipsd,*) z(iz), (mpsd(isps,ips,iz),ips=1,nps), time 
                        end do
                        close(ipsd)
                    enddo 
                else 
                    open(ipsd, file=trim(adjustl(profdir))//'/'  &
                        & //'psd-save.txt', status='replace')
                    write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
                    do iz = 1, Nz
                        write(ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
                    end do
                    close(ipsd)
                endif 
            endif 
            
    ! #ifdef disp_lim
            if (display_lim_in) display_lim = .false.
    ! #endif 
            
        end if    
        
        ! saving flx when climate is changed within model 
        if ( (any(climate) .and. any (ict_change))  & 
            ! or when time to record flx 
            & .or.(time>=rectime_flx(irec_flx+1)) &
            ! or when definining flx_save_alltime
            & .or. flx_save_alltime &
            ) then 
            
            if (.not. flx_save_alltime .and. time>=rectime_flx(irec_flx+1)) then 
                irec_flx = irec_flx + 1
            endif 
            
            if (.not. flx_recorded) then 
                
                do isps=1,nsp_sld 
                    open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', action='write',status='old',position='append')
                    write(isldflx(isps),*) time,(sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
                    close(isldflx(isps))
                
                    open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
                        & //'int_flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', action='write',status='old',position='append')
                    write(isldflx(isps),*) time,(int_flx_sld(isps,iflx)/time,iflx=1,nflx)
                    close(isldflx(isps))
                enddo 
                
                do ispa=1,nsp_aq 
                    open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', action='write',status='old',position='append')
                    write(iaqflx(ispa),*) time,(sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
                    close(iaqflx(ispa))
                
                    open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
                        & //'int_flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', action='write',status='old',position='append')
                    write(iaqflx(ispa),*) time,(int_flx_aq(ispa,iflx)/time,iflx=1,nflx)
                    close(iaqflx(ispa))
                enddo 
                
                do ispg=1,nsp_gas 
                    open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', action='write',status='old',position='append')
                    write(igasflx(ispg),*) time,(sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
                    close(igasflx(ispg))
                
                    open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
                        & //'int_flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', action='write',status='old',position='append')
                    write(igasflx(ispg),*) time,(int_flx_gas(ispg,iflx)/time,iflx=1,nflx)
                    close(igasflx(ispg))
                enddo 
                
                do ico2=1,6 
                    open(ico2flx(ico2), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', action='write',status='old',position='append')
                    if (ico2 .le. 4) then 
                        write(ico2flx(ico2),*) time,(sum(flx_co2sp(ico2,iflx,:)*dz(:)),iflx=1,nflx)
                    elseif (ico2 .eq. 5) then 
                        write(ico2flx(ico2),*) time &
                            & ,(sum(flx_co2sp(2,iflx,:)*dz(:))+sum(flx_co2sp(3,iflx,:)*dz(:))+sum(flx_co2sp(4,iflx,:)*dz(:)) &
                            & ,iflx=1,nflx)
                    elseif (ico2 .eq. 6) then 
                        write(ico2flx(ico2),*) time &
                            & ,(sum(flx_co2sp(3,iflx,:)*dz(:))+2d0*sum(flx_co2sp(4,iflx,:)*dz(:)) &
                            & ,iflx=1,nflx)
                    endif 
                    close(ico2flx(ico2))
                    
                    open(ico2flx(ico2), file=trim(adjustl(flxdir))//'/' &
                        & //'int_flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', action='write',status='old',position='append')
                    write(ico2flx(ico2),*) time,(int_flx_co2sp(ico2,iflx)/time,iflx=1,nflx)
                    close(ico2flx(ico2))
                enddo 
                
                
                open(iphint, file=trim(adjustl(flxdir))//'/'//'int_ph.txt', action='write',status='old',position='append')
                write(iphint,*) time,(-log10(int_ph(iz)/time),iz=1,nz)
                close(iphint)
            
                open(iphint2, file=trim(adjustl(flxdir))//'/'//'ph.txt', action='write',status='old',position='append')
                write(iphint2,*) time,(-log10(gamma(iz)*prox(iz)),iz=1,nz)
                close(iphint2)
                
                
            endif 
        endif 
        
        ! save flux all time for specified species
        if ( .not. flx_save_alltime  ) then 
            
            do isps=1,nsp_sld 
                if ( any ( chrsp_saveall == chrsld(isps)) ) then
                    open(isldflx(isps), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_all_sld-'//trim(adjustl(chrsld(isps)))//'.txt', action='write',status='old',position='append')
                    write(isldflx(isps),*) time,(sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
                    close(isldflx(isps))
                endif 
            enddo 
            
            do ispa=1,nsp_aq 
                if ( any ( chrsp_saveall == chraq(ispa) ) ) then
                    open(iaqflx(ispa), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_all_aq-'//trim(adjustl(chraq(ispa)))//'.txt', action='write',status='old',position='append')
                    write(iaqflx(ispa),*) time,(sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
                    close(iaqflx(ispa))
                endif 
            enddo 
            
            do ispg=1,nsp_gas 
                if ( any ( chrsp_saveall == chrgas(ispg) ) ) then
                    open(igasflx(ispg), file=trim(adjustl(flxdir))//'/' &
                        & //'flx_all_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', action='write',status='old',position='append')
                    write(igasflx(ispg),*) time,(sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
                    close(igasflx(ispg))
                endif 
            enddo 
            
        endif 
        
        ! save only if everything goes well
        if (dust_step) then
            if ( dust_norm /= dust_norm_prev .or. dust_change) then
            ! if ( dust_change ) then
                open(idust, file=trim(adjustl(flxdir))//'/'//'dust.txt', &
                    & status='old',action='write',position='append')
                ! write(idust,*) time-dt_prev,dust_norm_prev
                ! integration is now from time - dt to time with dust_norm
                write(idust,*) time-dt,dust_norm
                write(idust,*) time,dust_norm
                close(idust)
            endif 
            ! print *, time-dt,dust_norm
            ! print *, time,dust_norm
        endif
        
        progress_rate_prev = progress_rate
        
        ! call cpu_time(time_fin)
        call system_clock(t2,t_rate,t_max)
        if ( t2 < t1 ) then
            diff = (t_max - t1) + t2 + 1
        else
            diff = t2 - t1
        endif
        
        ! progress_rate = dt/(time_fin-time_start)*sec2yr ! (model yr)/(computer yr)
        ! progress_rate = (time_fin-time_start) ! (computer sec)
        progress_rate = diff/dble(t_rate) ! (computer sec)
        
        if (.not.timestep_fixed) then 
            if (it/=1) then 
                if (flgreducedt) then 
                    maxdt = maxdt/10d0
                    flgreducedt = .false.
                    count_dtunchanged = 0
                else
                    ! maxdt = maxdt* (progress_rate/progress_rate_prev)**0.33d0
                    maxdt = maxdt* (progress_rate/progress_rate_prev)**(-0.33d0)
                    if (maxdt > maxdt_max) maxdt = maxdt_max
                    if (dt < maxdt) count_dtunchanged = 0
                    ! if (dt > maxdt) dt = maxdt
                endif 
                
                if (count_dtunchanged > count_dtunchanged_Max_loc) then 
                    maxdt = maxdt*10d0
                    count_dtunchanged = 0
                endif 
            endif 
        endif 
        
        ! if (progress_rate ==0d0 .or. progress_rate_prev ==0d0) maxdt = 1d2
        ! print *,progress_rate,progress_rate_prev,maxdt,time_fin,time_start
        if (isnan(maxdt).or.maxdt ==0d0) then 
        ! if (.true.) then 
            print *
            print *, 'maxdt is nan or zero',progress_rate,progress_rate_prev,maxdt,time_fin,time_start
            stop
        endif 
        
        if (display  .and. (.not. display_lim)) then 
            print *
            print '(E11.3,a)',progress_rate,': computation time per iteration [sec]'
            print '(E11.3,a)',maxdt, ': maxdt [yr]'
            print '(i11,a)',count_dtunchanged,': count_dtunchanged'
            print *, '-----------------------------------------'
            print *
        endif 
        
    end do
            
            
    open(isldprof,file=trim(adjustl(profdir))//'/' &
        & //'prof_sld-save.txt', status='replace')
    open(igasprof,file=trim(adjustl(profdir))//'/' &
        & //'prof_gas-save.txt', status='replace')
    open(iaqprof,file=trim(adjustl(profdir))//'/' &
        & //'prof_aq-save.txt', status='replace')
    open(ibsd, file=trim(adjustl(profdir))//'/'  &
        & //'bsd-save.txt', status='replace')
    open(ipsd, file=trim(adjustl(profdir))//'/'  &
        & //'psd-save.txt', status='replace')
    open(isa,file=trim(adjustl(profdir))//'/' &
        & //'sa-save.txt', status='replace')
                
    write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
    write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
    write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
    write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' w[m/yr] ' &
        & , ' vol[m3/m3] ',' dens[g/cm3] ', ' blk[wt%] ',' cec[cmol/kg] ', ' time '
    write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
    write(isa,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '

    do iz = 1, Nz
        ucvsld1 = 1d0
        if (msldunit == 'blk') ucvsld1 = 1d0 - poro(iz)
        write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
        write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
        write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
        write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz)  &
            & ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6),cec(iz),time
        write(ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
        write(isa,*) z(iz),(hr(isps,iz),isps = 1, nsp_sld),time
    end do

    close(isldprof)
    close(iaqprof)
    close(igasprof)
    close(ibsd)
    close(ipsd)
    close(isa)
            
    if (do_psd) then 
        if (do_psd_full) then 
            do isps=1,nsp_sld
                open(ipsd, file=trim(adjustl(profdir))//'/'  &
                    & //'psd_'//trim(adjustl(chrsld(isps)))//'-save.txt', status='replace')
                write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
                do iz = 1, Nz
                    write(ipsd,*) z(iz), (mpsd(isps,ips,iz),ips=1,nps), time 
                end do
                close(ipsd)
            enddo 
        else 
            open(ipsd, file=trim(adjustl(profdir))//'/'  &
                & //'psd-save.txt', status='replace')
            write(ipsd,*) ' z[m]\log10(r[m]) ',(ps(ips),ips=1,nps),' time '
            do iz = 1, Nz
                write(ipsd,*) z(iz), (psd(ips,iz),ips=1,nps), time 
            end do
            close(ipsd)
        endif 
    endif 

    call system ('cp gases.in '//trim(adjustl(profdir))//'/gases.save')
    call system ('cp solutes.in '//trim(adjustl(profdir))//'/solutes.save')
    call system ('cp slds.in '//trim(adjustl(profdir))//'/slds.save')
    call system ('cp extrxns.in '//trim(adjustl(profdir))//'/extrxns.save')
    call system ('cp kinspc.in '//trim(adjustl(profdir))//'/kinspc.save')
    call system ('cp sa.in '//trim(adjustl(profdir))//'/sa.save')

endsubroutine weathering_main

endmodule scepter_weathering_main