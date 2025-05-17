module scepter_variables
    use scepter_constants
    implicit none

    public

    integer nsp_gas,nrxn_ext
    character(5),dimension(:),allocatable::chraq,chrsld,chrgas,chrrxn_ext,chrsld_kinspc 
    real(kind=8),dimension(:),allocatable::kin_sld_spc
    character(500) sim_name,runname_save,cwd,path,path2,cmd
    real(kind=8) ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,plant_rain,zml_ref,tc,rainpowder_2nd &
        & ,step_tau
    integer count_dtunchanged_Max

    real(kind=8) dt  ! yr 
    real(kind=8) time
    real(kind=8) pco2i,pnh3i,proi
    real(kind=8) :: rho_grain = 2.7d0 ! g/cm3 as soil grain density 
    real(kind=8) :: rho_grain_calc,rho_grain_calcx != 2.7d0 ! g/cm3 as soil grain density 
    real(kind=8) :: rho_error,rho_tol, poroi_calc 
    logical(kind=8) :: incld_blk
    real(kind=8)::zsupp_plant = 0.3d0 !  e-folding decrease
    real(kind=8) :: rcharge
    real(kind=8) :: rough_c0_b = 10d0**(3.3d0)
    real(kind=8) :: rough_c1_b = 0.33d0
    real(kind=8) :: c_disp,c0_disp,c1_disp,zdisp  ! dispersion coefficients 
    real(kind=8) kho,ucv,kco2,k1,kw,k2,khco2i,knh3,k1nh3,khnh3i,kn2o

    !-----------------------------
    ! Loop indices and counters
    !-----------------------------
    integer iz,it,ispa,ispg,isps,irxn,ispa2,ispg2,isps2,ico2,ph_iter,isps_kinspc,isps_sa

    !-----------------------------
    ! Error and tolerance parameters
    !-----------------------------
    real(kind=8) error 
    real(kind=8) :: tol = 1d-6

    !-----------------------------
    ! Output and file parameters
    !-----------------------------
    real(kind=8) rectime_prof(nrec_prof)
    real(kind=8) rectime_flx(nrec_flx)
    character(3) chr
    character(256) runname,workdir, chrz(3), chrq(3),base,fname, chrrain, flxdir, profdir
    character(500) loc_runname_save
    integer irec_prof, irec_flx, iter
    logical flx_recorded

    !-----------------------------
    ! Time step parameters
    !-----------------------------
    integer  iflx
    real(kind=8) :: maxdt = 0.2d0 ! for basalt exp?
    real(kind=8) :: maxdt_max = 1d2  ! default   
    ! real(kind=8) :: maxdt_max = 1d0   ! when time step matters a reduced value might work 

    !-----------------------------
    ! Model switches and flags
    !-----------------------------
    logical :: read_data = .false.
    logical :: incld_rough = .true.
    logical :: cplprec = .true.
    logical :: dust_wave = .false.
    logical :: al_inhibit = .false.
    logical :: timestep_fixed = .false.
    logical :: display = .true.
    logical :: regular_grid = .true.
    logical :: sld_enforce = .false.
    logical :: poroevol = .false.
    logical :: surfevol1 = .false.
    logical :: surfevol2 = .false.
    logical :: noncnstw = .true.  ! varied with porosity
    logical :: display_lim = .false. ! limiting display fluxes and concs. 
    logical :: dust_step = .true.
    logical,dimension(3) :: climate != .false.
    logical :: season = .false.
    logical :: disp_ON = .false.
    logical :: disp_FULL_ON = .false.
    logical :: ads_ON = .true.
    logical :: ph_limits_dust = .false.
    logical :: aq_close = .false.
    logical :: act_ON = .false.
    logical ads_ON_tmp,dust_Off

    !-----------------------------
    ! pH and chemical parameters
    !-----------------------------
    real(kind=8)::z_chk_ph = 0.5d0
    real(kind=8)::ph_lim = 6.8d0 
    real(kind=8) ph_ave,z_ave

    !-----------------------------
    ! Dust and wave parameters
    !-----------------------------
#ifdef def_flx_save_alltime
    logical :: flx_save_alltime = .true.
#else
    logical :: flx_save_alltime = .false.
#endif

    real(kind=8) :: tol_step_tau = 1d-6 ! yr time duration during which dust is added
    real(kind=8) :: wave_tau = 2d0 ! yr periodic time for wave 
    real(kind=8) :: dust_norm = 0d0
    real(kind=8) :: dust_norm_prev = 0d0
    logical :: dust_change 

    !-----------------------------
    ! Climate parameters
    !-----------------------------
    real(kind=8),dimension(:,:),allocatable :: clim_T,clim_q,clim_sat
    real(kind=8),dimension(3) :: dct,ctau
    integer iclim,ict
    integer,dimension(3)::nclim,ict_prev
    logical,dimension(3)::ict_change
    character(50),dimension(3) :: clim_file

    !-----------------------------
    ! Model type parameters
    !-----------------------------
    integer iwtype 
    integer imixtype 
    integer imixtype_background
    integer imixtype_OM
    integer iroughtype 

    !-----------------------------
    ! Input parameters
    !-----------------------------
    logical display_lim_in !  defining whether limiting display or not  (input from input file swtiches.in)
    logical poroiter_in !  true if porosity (or w) is iteratively checked  (input from input file swtiches.in)
    logical lim_minsld_in !  true if minimum sld conc. is enforced  (input from input file swtiches.in)

    !-----------------------------
    ! Recording time parameters
    !-----------------------------
#ifndef nrec_prof_in
    data rectime_prof /1d1,3d1,1d2,3d2,1d3,3d3,1d4,3d4 &
        & ,1d5,2d5,3d5,4d5,5d5,6d5,7d5,8d5,9d5,1d6,1.1d6,1.2d6/
#endif 
    real(kind=8) :: savetime = 1d3
    real(kind=8) :: dsavetime = 1d3
    logical :: rectime_scheme_old = .false.

    integer poro_iter , poro_iter_max
    real(kind=8) beta 

    !-----------------------------
    ! Time step control parameters
    !-----------------------------
    logical :: flgback = .false.
    logical :: flgreducedt = .false.
    logical :: flgreducedt_prev = .false.
    real(kind=8) time_start, time_fin, progress_rate, progress_rate_prev
    integer count_dtunchanged,count_dtunchanged_Max_loc  

    integer ::nsp_sld_cnst != nsp_sld_all - nsp_sld
    integer ::nsp_aq_cnst != nsp_aq_all - nsp_aq
    integer ::nsp_gas_cnst != nsp_gas_all - nsp_gas
    integer ::nsp3 != nsp_sld + nsp_aq + nsp_gas
    integer :: nflx ! = 5 + nrxn_ext + nsp_sld  
    integer :: nsld_kinspc,nsld_kinspc_add

    character(5),dimension(:),allocatable::chrsld_2
    character(5),dimension(nsp_sld_all)::chrsld_all
    character(5),dimension(nsp_aq_ph)::chraq_ph
    character(5),dimension(nsp_aq_all)::chraq_all
    character(5),dimension(nsp_gas_ph)::chrgas_ph
    character(5),dimension(nsp_gas_all)::chrgas_all
    character(5),dimension(nrxn_ext_all)::chrrxn_ext_all

    real(kind=8) psu_pr,pssigma_pr,psu_rain,psw_rain,pssigma_rain,ps_new,ps_newp,dvd_res,error_psd,volsld,flx_max_max,psd_th_flex
    real(kind=8) p80_tmp
    real(kind=8) :: ps_sigma_std = 1d0
    integer ips,iips,ips_new
    logical psd_error_flg,no_psd_prevrun

    !-----------------------------
    ! Particle size distribution flux parameters
    !-----------------------------
    logical :: do_psd = .true.
    logical :: do_psd_norm = .true.
    logical :: do_psd_full = .true.
    logical :: psd_lim_min = .true.
    logical :: psd_vol_consv = .false.
    logical :: psd_impfull = .false.
    logical :: psd_loop = .true.
    logical :: psd_enable_skip = .true.

    real(kind=8),dimension(:),allocatable::pssigma_rain_list,psu_rain_list,psw_rain_list
    real(kind=8),dimension(:),allocatable::pssigma_rain_list_in,psu_rain_list_in,psw_rain_list_in
 
    real(kind=8):: rough_c0 != 10d0**(3.3d0)
    real(kind=8):: rough_c1 != 0.33d0
    character(10)::roughref_b
    integer nsld_sa
    character(5),dimension(:),allocatable::chrsld_sa
    real(kind=8) time_pbe,dt_pbe,dt_save
    integer nsld_nopsd
    character(5),dimension(:),allocatable::chrsld_nopsd ! minerals whose PSDs tracking is not conducted for some reasons (e.g., too fast; mostly precipitating etc.)

    !-----------------------------
    ! CEC parameters
    !-----------------------------
    integer nsld_cec
    character(5),dimension(:),allocatable::chrsld_cec
    real(kind=8),dimension(nsp_sld_all):: mcec_all,mcec_all_def
    real(kind=8),dimension(nsp_sld_all,nsp_aq_all):: logkhaq_all,logkhaq_all_def
    real(kind=8),dimension(nsp_sld_all):: beta_all,beta_all_def


    logical,dimension(nsp_sld_all)::cec_pH_depend

    logical:: anealing_dust = .false.

    integer ieqgas_h0,ieqgas_h1,ieqgas_h2
    data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

    integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
    data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

    integer ieqaq_co3,ieqaq_hco3
    data ieqaq_co3,ieqaq_hco3/1,2/

    integer ieqaq_so4,ieqaq_so42
    data ieqaq_so4,ieqaq_so42/1,2/



    integer isldprof,isldprof2,isldprof3,iaqprof,igasprof,isldsat,ibsd,irate,ipsd,ipsdv,ipsds,ipsdflx  &
        & ,isa,isa2,iaqprof2,iaqprof3,iaqprof4,iaqprof5,iaqprof6

    real(kind=8) zml_background,zml_OM,zml_dust
    real(kind=8) dbl_ref
    integer :: nz_disp = 10
    real(kind=8) dt_prev

    logical print_cb,ph_error,save_trans,ads_error
    character(500) print_loc

    real(kind=8) def_dust,def_rain,def_pr,def_OM_frc
    character(3) chriz
    character(50) chrfmt
    integer::itflx,iadv,idif,irain,ires
    data itflx,iadv,idif,irain/1,2,3,4/


    !-----------------------------
    ! Save parameters
    !-----------------------------
    integer nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save,nsld_kinspc_save,nsld_sa_save 
    character(5),dimension(:),allocatable::chraq_save,chrsld_save,chrgas_save,chrrxn_ext_save &
        & ,chrsld_kinspc_save,chrsld_sa_save
    real(kind=8),dimension(:,:),allocatable::msld_save,mgas_save,maq_save,hr_save
    real(kind=8),dimension(:,:,:),allocatable::mpsd_save
    real(kind=8),dimension(:),allocatable::kin_sldspc_save,hrii_save

    !-----------------------------
    ! Timing parameters
    !-----------------------------
    integer t1,t2,t_rate,t_max,diff
    character(3):: msldunit = 'blk'
    real(kind=8) ucvsld1,ucvsld2

    integer iph

endmodule scepter_variables 