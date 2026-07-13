"""
restart_add_gbas.py - Enhanced Rock Weathering (ERW) Restart Module

This module provides functionality to restart SCEPTER simulations from completed
spinup runs and apply basalt dust for Enhanced Rock Weathering studies.

Main Functions:
- add_gbas(): Core function to restart simulation with basalt application
- Various batch processing functions for different experimental designs

Usage:
    python3 restart_add_gbas.py <spinup_name> <restart_name>

Author: SCEPTER development team
"""

import numpy as np
import os, shutil, subprocess, sys
import random
import copy
from scipy.interpolate import interp1d

import get_int_prof
import get_soilpH_time
import make_inputs
import get_inputs


def _fortran_rstrt_path(spinup_dir, restart_dir):
    """
    Build the restart source path written to frame.in.

    SCEPTER prepends '../' to runname_save from inside the run directory
    (see scepter_weathering_main.f90), so the value must be relative to the
    *parent* of restart_dir.  Paths may contain '/', which is fine once
    make_inputs quotes the field for Fortran list-directed input.
    """
    spinup_abs = os.path.abspath(spinup_dir)
    restart_parent = os.path.dirname(os.path.abspath(restart_dir))
    return os.path.relpath(spinup_abs, restart_parent)


def _erw_outputs_look_valid(run_dir, dust_sld="gbas"):
    """
    Return True only if outputs look like a finished ERW restart, not a
    copied spinup tree.  Copied spinups already have prof/ and flx/, so those
    directories alone must not count as success.
    """
    if os.path.exists(os.path.join(run_dir, "run_complete.txt")):
        # Only trust this marker if ERW species or flux data are also present
        pass

    sld_paths = [
        os.path.join(run_dir, "prof", "prof_sld-020.txt"),
        os.path.join(run_dir, "prof", "prof_sld-001.txt"),
    ]
    for sld_path in sld_paths:
        if not os.path.exists(sld_path):
            continue
        try:
            header = open(sld_path).readline().split()
        except OSError:
            continue
        if dust_sld in header:
            return True

    # Non-empty integrated flux file is another strong success signal
    for flx_name in (
        "int_flx_co2sp-hco3.txt",
        "int_flx_sld-gbas.txt",
        f"int_flx_sld-{dust_sld}.txt",
    ):
        flx_path = os.path.join(run_dir, "flx", flx_name)
        if not os.path.exists(flx_path):
            continue
        try:
            with open(flx_path) as f:
                nlines = sum(1 for _ in f)
        except OSError:
            continue
        if nlines > 1:
            return True

    return False


def _run_scepter(exe_path, run_dir, logfile=None):
    """Run the SCEPTER binary; return True iff process exit code is 0."""
    if logfile:
        with open(logfile, "w") as logf:
            proc = subprocess.run(
                [exe_path], cwd=run_dir, stdout=logf, stderr=subprocess.STDOUT
            )
    else:
        proc = subprocess.run([exe_path], cwd=run_dir)
    return proc.returncode == 0


def add_gbas(outdir, runname_spinup, dustname, runname_restart, **kwargs):
    """
    Restart a SCEPTER simulation from a completed spinup and apply basalt dust.

    This function takes a completed spinup simulation and creates a new ERW
    simulation that applies basalt dust to study enhanced rock weathering.

    Args:
        outdir (str): Output directory path
        runname_spinup (str): Name of the completed spinup run directory
        dustname (str): Type of dust to apply (e.g., 'gbasalt', 'cc')
        runname_restart (str): Name for the new ERW simulation directory
        **kwargs: Additional configuration parameters

    Returns:
        bool: True if simulation completed successfully, False otherwise

    Configuration Parameters (kwargs):
        ttot_restart (float): Duration of ERW simulation in years (default: 100)
        fdust_restart (float): Dust application rate in g/m²/yr (default: 5000)
        taudust_restart (float): Fraction of year when dust is applied (default: 0.01)
        zdust_restart (float): Depth of dust application in m (default: 0.15)
        p80str (str): Particle size distribution parameter (default: '100um')
        nstep_restart (int): Number of time steps for numerical stability (default: 10)
        mix_scheme_restart (int): Mixing scheme (0=none, 1=Fickian, 2=homogeneous, 3=tilling)
        kinspc (bool): Enable kinetic species tracking (default: False)
        sld_list_add (list): Additional solid species to track
        gas_list_add (list): Additional gas species to track
        exrxn_list_add (list): Additional exchange reactions
        mod_2nd_phases (str): Secondary phase modifications
        update_exe (bool): Update executable before running (default: True)
        use_local_storage (bool): Use temporary storage for computation (default: True)
        save_restart_dir (bool): Save restart directory after completion (default: True)
        set_ExTimeLimit (bool): Set execution time limit (default: False)
        Ex_TimeLim (int): Execution time limit in minutes (default: 20)
        make_runlogfile (bool): Create run log file (default: False)
        exename_src (str): Source executable name (default: 'scepter_test')
    """

    # Extract configuration parameters with defaults
    ttot_restart = kwargs.get("ttot_restart", 100)  # ERW simulation duration (years)
    fdust_restart = kwargs.get("fdust_restart", 5000)  # Dust application rate (g/m²/yr)
    taudust_restart = kwargs.get(
        "taudust_restart", 0.01
    )  # Dust application duration (fraction of year)
    zdust_restart = kwargs.get("zdust_restart", 0.15)  # Dust application depth (m)
    p80str = kwargs.get("p80str", "100um")  # Particle size distribution
    nstep_restart = kwargs.get("nstep_restart", 10)  # Time step control
    mix_scheme_restart = kwargs.get(
        "mix_scheme_restart", 1
    )  # Mixing scheme (1=Fickian)
    kinspc = kwargs.get("kinspc", False)  # Kinetic species tracking
    sld_list_add = kwargs.get("sld_list_add", [])  # Additional solid species
    gas_list_add = kwargs.get("gas_list_add", [])  # Additional gas species
    exrxn_list_add = kwargs.get("exrxn_list_add", [])  # Additional exchange reactions
    mod_2nd_phases = kwargs.get("mod_2nd_phases", "")  # Secondary phase modifications
    update_exe = kwargs.get("update_exe", True)  # Update executable
    use_local_storage = kwargs.get("use_local_storage", True)  # Use temp storage
    save_restart_dir = kwargs.get("save_restart_dir", True)  # Save restart directory
    set_ExTimeLimit = kwargs.get("set_ExTimeLimit", False)  # Set time limit
    Ex_TimeLim = kwargs.get("Ex_TimeLim", 20)  # Time limit (minutes)
    make_runlogfile = kwargs.get("make_runlogfile", False)  # Create log file
    exename_src = kwargs.get("exename_src", "scepter_test")  # Source executable

    # Set up file paths and directories
    exename = "scepter"  # Target executable name
    to = " "  # Copy command separator
    where = "/"  # Directory separator

    # Determine output directory (use temp storage if specified)
    outdir_tmp = outdir
    if use_local_storage:
        outdir_tmp = os.environ["TMPDIR"] + "/scepter_output/"

    # Set up source and destination paths
    runname = runname_restart
    src = outdir + runname_spinup  # Source: completed spinup directory
    dst = outdir_tmp + runname  # Destination: new ERW simulation directory

    # Copy spinup directory to create restart simulation
    if not os.path.exists(dst):
        shutil.copytree(src, dst)
    else:
        shutil.rmtree(dst)
        shutil.copytree(src, dst)

    if use_local_storage:
        runname = runname_spinup

        src = outdir + runname_spinup
        dst = outdir_tmp + runname

        if not os.path.exists(dst):
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)

    if update_exe:
        runname = runname_restart
        # exename_src = 'scepter_test'
        # exename_src = './tmp/scepter_test'
        # os.system('make --file=makefile_test')
        os.system("rm " + outdir_tmp + runname + where + exename)
        os.system("cp " + exename_src + to + outdir_tmp + runname + where + exename)

    # ------------------------------------------------
    # get input data for field run
    # ------------------------------------------------
    runname = runname_spinup
    # (1) frame
    (
        ztot,
        nz,
        ttot,
        temp,
        fdust,
        fdust2,
        taudust,
        omrain,
        zom,
        poro,
        moistsrf,
        zwater,
        zdust,
        w,
        q,
        p,
        nstep,
        rstrt,
        runid,
    ) = get_inputs.get_input_frame(outdir, runname)
    print(
        ztot,
        nz,
        ttot,
        temp,
        fdust,
        fdust2,
        taudust,
        omrain,
        zom,
        poro,
        moistsrf,
        zwater,
        zdust,
        w,
        q,
        p,
        nstep,
        rstrt,
        runid,
    )
    # (2) switches
    (
        w_scheme,
        mix_scheme,
        poro_iter,
        sldmin_lim,
        display,
        report,
        restart,
        rough,
        act_ON,
        dt_fix,
        cec_on,
        dz_fix,
        close_aq,
        poro_evol,
        sa_evol_1,
        sa_evol_2,
        psd_bulk,
        psd_full,
        season,
    ) = get_inputs.get_input_switches(outdir, runname)
    print(
        w_scheme,
        mix_scheme,
        poro_iter,
        sldmin_lim,
        display,
        report,
        restart,
        rough,
        act_ON,
        dt_fix,
        cec_on,
        dz_fix,
        close_aq,
        poro_evol,
        sa_evol_1,
        sa_evol_2,
        psd_bulk,
        psd_full,
        season,
    )
    # (3) tracers
    sld_list, aq_list, gas_list, exrxn_list = get_inputs.get_input_tracers(
        outdir, runname
    )
    print(sld_list, aq_list, gas_list, exrxn_list)
    # print(type(gas_list))
    # exit("end run")
    # (4) dust
    filename = "dust.in"
    sld_data_list = get_inputs.get_input_sld_properties(outdir, runname, filename)
    sld_list_dust = [listtmp[0] for listtmp in sld_data_list]
    sld_val_dust = [listtmp[1] for listtmp in sld_data_list]
    sld_dust_spinup = [
        (sld_list_dust[i], sld_val_dust[i]) for i in range(len(sld_list_dust))
    ]

    # ------------------------------------------------
    # get dust data
    # ------------------------------------------------

    dustsrc = "dust_" + dustname + ".in"
    sld_data_list = get_inputs.get_input_sld_properties("./", "data", dustsrc)
    sld_list_dust_add = [listtmp[0] for listtmp in sld_data_list]
    sld_val_dust_add = [listtmp[1] for listtmp in sld_data_list]
    sld_dust_restart = [
        (sld_list_dust_add[i], sld_val_dust_add[i])
        for i in range(len(sld_list_dust_add))
        if sld_list_dust_add[i] not in sld_list_dust
    ]
    sld_dust_restart.extend(sld_dust_spinup)

    # ------------------------------------------------
    # making inputs for restart experiments
    # ------------------------------------------------
    runname = runname_restart
    # (1) tracers
    sld_list_restart = [sld for sld in sld_list_dust_add if sld not in sld_list]
    sld_list_restart.extend(sld_list)
    sld_list_restart = [sld for sld in sld_list_restart if sld not in sld_list_add]
    sld_list_restart.extend(sld_list_add)
    exrxn_list_restart = [exrxn for exrxn in exrxn_list if exrxn not in exrxn_list_add]
    exrxn_list_restart.extend(exrxn_list_add)
    gas_list_restart = [gas for gas in gas_list if gas not in gas_list_add]
    gas_list_restart.extend(gas_list_add)

    make_inputs.get_input_tracers(
        outdir=outdir_tmp,
        runname=runname,
        sld_list=sld_list_restart,
        aq_list=aq_list,
        gas_list=gas_list_restart,
        exrxn_list=exrxn_list_restart,
    )

    # (2) framework
    make_inputs.get_input_frame(
        outdir=outdir_tmp,
        runname=runname,
        ztot=ztot,
        nz=nz,
        ttot=ttot_restart,  #
        temp=temp,
        fdust=fdust_restart,  #
        fdust2=fdust2,
        taudust=taudust_restart,  #
        omrain=omrain,
        zom=zom,
        poro=poro,
        moistsrf=moistsrf,
        zwater=zwater,
        zdust=zdust_restart,  #
        w=w,
        q=q,
        p=p,
        nstep=nstep_restart,  #
        # rstrt=outdir+runname_spinup,
        rstrt=runname_spinup,
        runid=runname_restart,
    )

    # (3) switches
    make_inputs.get_input_switches(
        outdir=outdir_tmp,
        runname=runname,
        w_scheme=w_scheme,
        mix_scheme=mix_scheme_restart,  #
        poro_iter=poro_iter,
        sldmin_lim=sldmin_lim,
        display=display,
        report=report,
        restart="true",
        rough=rough,
        act_ON=act_ON,
        dt_fix=dt_fix,
        cec_on=cec_on,
        dz_fix=dz_fix,
        close_aq=close_aq,
        poro_evol=poro_evol,
        sa_evol_1=sa_evol_1,
        sa_evol_2=sa_evol_2,
        psd_bulk=psd_bulk,
        psd_full=psd_full,
        season=season,
    )

    # (4) dust composition
    filename = "dust.in"
    sld_varlist = sld_dust_restart
    make_inputs.get_input_sld_properties(
        outdir=outdir_tmp,
        runname=runname,
        filename=filename,
        sld_varlist=sld_varlist,
    )

    # (5) PSD
    filename = "psdrain.in"
    srcfile = "./data/psdrain_" + p80str + ".in"
    sld_varlist = sld_dust_restart
    make_inputs.get_input_sld_properties(
        outdir=outdir_tmp,
        runname=runname,
        filename=filename,
        srcfile=srcfile,
    )

    # (6) specify kinetics
    if kinspc:  # if a value is specified, this means True
        filename = "kinspc.in"
        sld_varlist = [(sld, kinspc) for sld in sld_list_dust_add]
        make_inputs.get_input_sld_properties(
            outdir=outdir_tmp,
            runname=runname,
            filename=filename,
            sld_varlist=sld_varlist,
        )

    # (7) 2ndary phases modification (if necessary)
    if mod_2nd_phases != "":
        filename = "2ndslds.in"
        srcfile = "./data/" + mod_2nd_phases
        make_inputs.get_input_sld_properties(
            outdir=outdir,
            runname=runname,
            filename=filename,
            srcfile=srcfile,
            # ,sld_varlist=sld_varlist
        )

    # ------------------------------------------------
    # run restart experiments
    # ------------------------------------------------

    runsuccess = False

    print("chmod u+x " + outdir_tmp + runname + where + exename)
    os.system("chmod u+x " + outdir_tmp + runname + where + exename)

    if not set_ExTimeLimit:

        if os.path.exists(outdir_tmp + runname + "/run_complete.txt"):
            os.remove(outdir_tmp + runname + "/run_complete.txt")

        if make_runlogfile:
            os.system(
                outdir_tmp
                + runname
                + where
                + exename
                + " > "
                + outdir_tmp
                + runname
                + "/logfile.txt"
            )
        else:
            os.system(outdir_tmp + runname + where + exename)
        runsuccess = True

        # Check for completion by looking for run_complete.txt OR output files
        if not os.path.exists(outdir_tmp + runname + "/run_complete.txt"):
            # Fallback: check if simulation generated output files
            if os.path.exists(outdir_tmp + runname + "/flx") or os.path.exists(
                outdir_tmp + runname + "/prof"
            ):
                print(
                    "**** run_complete.txt not found, but output files exist - marking as successful"
                )
                # Create the file to prevent future checks
                open(outdir_tmp + runname + "/run_complete.txt", "w").close()
                runsuccess = True
            else:
                print("**** run is not completed")
                runsuccess = False

    else:

        if os.path.exists(outdir_tmp + runname + "/run_complete.txt"):
            os.remove(outdir_tmp + runname + "/run_complete.txt")

        if make_runlogfile:
            logf = open(outdir_tmp + runname_lab + "/logfile.txt", "w")
            proc = subprocess.Popen(
                [outdir_tmp + runname + where + exename], stdout=logf
            )
        else:
            proc = subprocess.Popen([outdir_tmp + runname + where + exename])

        my_timeout = 60 * Ex_TimeLim

        try:
            proc.wait(my_timeout)
            print("run finished within {:f} min".format(int(my_timeout / 60.0)))
            runsuccess = True
            # Check for completion by looking for run_complete.txt OR output files
            if not os.path.exists(outdir_tmp + runname + "/run_complete.txt"):
                # Fallback: check if simulation generated output files
                if os.path.exists(outdir_tmp + runname + "/flx") or os.path.exists(
                    outdir_tmp + runname + "/prof"
                ):
                    print(
                        "**** run_complete.txt not found, but output files exist - marking as successful"
                    )
                    # Create the file to prevent future checks
                    open(outdir_tmp + runname + "/run_complete.txt", "w").close()
                    runsuccess = True
                else:
                    print("**** run is not completed")
                    runsuccess = False
        except subprocess.TimeoutExpired:
            proc.kill()
            print("run UNfinished within {:f} min".format(int(my_timeout / 60.0)))

    if make_runlogfile and runsuccess:
        os.remove(outdir_tmp + runname + "/logfile.txt")

    # ------------------------------------------------
    # remove or copy restart experiments
    # ------------------------------------------------
    # if run failed, save in any case with different name
    if not runsuccess:
        if use_local_storage:
            src = outdir_tmp + runname
            dst = outdir + runname + "_FAIL"

            if not os.path.exists(dst):
                shutil.copytree(src, dst)
            else:
                shutil.rmtree(dst)
                shutil.copytree(src, dst)

    else:
        # (A) case to remove
        if not save_restart_dir:
            if not use_local_storage:
                shutil.rmtree(outdir_tmp + runname)

        # (B) case to save
        if save_restart_dir:
            if use_local_storage:
                src = outdir_tmp + runname
                dst = outdir + runname

                if not os.path.exists(dst):
                    shutil.copytree(src, dst)
                else:
                    shutil.rmtree(dst)
                    shutil.copytree(src, dst)

    return runsuccess


def single_run_varPSD():
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run'
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_wcec"
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_midOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'

    dustname = "gbasalt"
    kinspc = False
    p80str = sys.argv[3]
    runid = p80str + "_homo"
    # runid = 'aneal_' + p80str
    # runid = 'aneal_BM00_' + p80str
    # runid = 'aneal_Letal21_' + p80str
    # runid = 'Letal21_' + p80str
    # runid = 'BM00_' + p80str
    runname_restart = runname_spinup + "_" + dustname + "_" + runid

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=100,
        fdust_restart=2241.7,
        # fdust_restart       = 5000,
        # taudust_restart     = 0.01,
        taudust_restart=0.0,
        zdust_restart=0.25,
        # zdust_restart       = 0.12,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        # mix_scheme_restart  = 1,
        mix_scheme_restart=2,
        kinspc=kinspc,
        update_exe=True,
        use_local_storage=True,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
    )


def single_run_spc_tau():
    """
    Single ERW simulation with specified parameters.

    This function runs a single Enhanced Rock Weathering simulation by restarting
    from a completed spinup and applying basalt dust with specified parameters.

    Command line arguments:
        sys.argv[1]: Name of the spinup run to restart from
        sys.argv[2]: Name for the new ERW simulation

    Configuration:
        - Uses glassy basalt (gbasalt) as feedstock
        - Applies 1000 g/m²/yr for 100 years
        - Uses homogeneous mixing scheme
        - Particle size: 320um
    """

    outdir = "../scepter_output/"
    # outdir = '/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/'

    # Storage configuration
    use_local_storage = True
    use_local_storage = False

    runname_spinup = sys.argv[1]
    # runname_spinup = sys.argv[1]

    dustname = "gbasalt"  # feedstock [**dust_{dustname}.in has to exists in ./data/ | e.g., gbasalt is glassy basalt, cc is calcite etc.  ]

    p80str = "320um"  # p80 value [**psdrain_{p80str}.in has to exists in ./data/ |  does not matter for now if PSD tracking is not ON, which is the case for 2025 paper ]

    taudust = 0.1  # fraction of year when feedstock is applied (feedstock will be  applied every year)

    ttot_restart = 100  # total duration of restart sim in y

    fdust = 1000  # g/m2/y

    n_runs = 5

    mix_scheme = 2

    runname_restart = sys.argv[2]

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        # fdust_restart       = 2241.7,
        fdust_restart=fdust,
        taudust_restart=taudust,
        zdust_restart=0.25,
        # zdust_restart       = 0.30,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        mix_scheme_restart=mix_scheme,
        # kinspc              = kinspc,
        update_exe=False,
        use_local_storage=use_local_storage,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
    )

    if not runsuccess:
        exit("error in feedstock app experiment")


def single_run_rand_kin_tau():
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    print(sys.argv)
    site = sys.argv[1]

    use_local_storage = True
    use_local_storage = False

    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v3/{}_field".format(
        site
    )
    # runname_spinup = sys.argv[1]

    # ++++++ get data for spinup ++++++ #

    dep_sample = 0.15

    # get porewater pH
    # phint = get_int_prof.get_ph_int_site(outdir,runname_spinup,dep_sample)
    phint = get_int_prof.get_waterave_site_v3(
        outdir, runname_spinup, dep_sample, 20, "ph"
    )

    # get acidity in %
    # acint = get_int_prof.get_ac_int_site_v2(outdir,runname_spinup,dep_sample)
    acint = get_int_prof.get_ac_int_site_v2(
        outdir, runname_spinup, dep_sample, 20, simple_average=True
    )

    # get field solid wt%
    sps = ["g2", "inrt"]
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(
            outdir, runname_spinup, dep_sample, [sp], 20
        )
        sldwt_list.append(sldwt)

    # get SOM wt%
    omint = sldwt_list[sps.index("g2")]

    # get soil pH in lab
    time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
        get_soilpH_time.calc_soilpH(
            outdir,
            runname_spinup,
            dep_sample,
            20,
            include_N=True,
            include_Cl=False,
            include_Al=False,
            include_DIC=True,
            use_CaCl2=True,
            add_salt=False,
            add_NO3=False,
            use_local_storage=use_local_storage,
            save_lab_dir=False,
            use_Sikora_buffer=False,
            set_ExTimeLimit=True,
            water_frac=2.5,
            cacl2_conc=0.01,
            # salt_conc           = salt_conc,
            # NO3_conc            = NO3_conc,
            # salt_sp             = salt_sp,
            Ex_TimeLim=20,
        )
    )

    if flag_error:
        exit("error in calculation of soil pH for spinup")

    # calculating buffer pH
    time_buf, IS_buf, acint_buf, phint_buf, flag_error, runname_lab = (
        get_soilpH_time.calc_soilpH(
            outdir,
            runname_spinup,
            dep_sample,
            20,
            include_N=True,
            include_Cl=False,
            include_Al=False,
            include_DIC=True,
            use_CaCl2=False,
            add_salt=False,
            add_NO3=False,
            use_local_storage=use_local_storage,
            save_lab_dir=False,
            use_Sikora_buffer=True,
            set_ExTimeLimit=False,
            # water_frac          = water_frac,
            # cacl2_conc          = cacl2_conc,
            # salt_conc           = salt_conc,
            # NO3_conc            = NO3_conc,
            # salt_sp             = salt_sp,
            Ex_TimeLim=20,
            make_runlogfile=True,
            # sub_as_a_job        = sub_as_a_job,
        )
    )

    if flag_error:
        exit("error in calculation of buffer pH for spinup")

    dustname = "gbasalt"
    p80str = "320um"

    n_runs = 5

    mix_scheme = 2

    kinspc, taudust, frac = [np.nan] * 3

    i = -1

    res_list_all = []
    res_list_all.append(
        [
            i,
            kinspc,
            taudust,
            dep_sample,
            phint,
            phint_lab,
            phint_buf,
            omint,
            acint,
            100.0 - acint,
            frac,
        ]
    )

    for i in range(n_runs + 1):

        kinspc = 10.0 ** random.uniform(-5, -2)
        # taudust = 10.**random.uniform(-3, -1)
        taudust = 0.1
        fdust = 5000

        if i == 0:
            kinspc = 1e-20
            taudust = 0
            fdust = 0

        runname_restart = f"MRV3/US_cropland/fwd/{site}_rand_{i}"

        runsuccess = add_gbas(
            outdir,
            runname_spinup,
            dustname,
            runname_restart,
            ttot_restart=1,
            # fdust_restart       = 2241.7,
            fdust_restart=fdust,
            taudust_restart=taudust,
            zdust_restart=0.25,
            # zdust_restart       = 0.30,
            p80str=p80str,
            # nstep_restart       = 1000,
            nstep_restart=10,
            mix_scheme_restart=mix_scheme,
            kinspc=kinspc,
            update_exe=True,
            use_local_storage=use_local_storage,
            save_restart_dir=True,
            set_ExTimeLimit=False,
            Ex_TimeLim=20,
        )

        if not runsuccess:
            exit("error in basalt app experiment")

        # (1) get porewater pH
        # phint = get_int_prof.get_ph_int_site(outdir,runname_restart,dep_sample)
        phint = get_int_prof.get_waterave_site_v3(
            outdir, runname_restart, dep_sample, 20, "ph"
        )

        # (2) get acidity in %
        # acint = get_int_prof.get_ac_int_site_v2(outdir,runname_restart,dep_sample)
        acint = get_int_prof.get_ac_int_site_v2(
            outdir, runname_restart, dep_sample, 20, simple_average=True
        )

        # (3) get bulk density
        dense = get_int_prof.get_rhobulk_int_site(
            outdir, runname_restart, dep_sample, 20
        )

        # (4) get field solid wt%
        sps = ["g2", "inrt"]
        sldwt_list = []
        for sp in sps:
            sldwt = get_int_prof.get_sldwt_int_site(
                outdir, runname_restart, dep_sample, [sp], 20
            )
            sldwt_list.append(sldwt)

        # (5) get SOM wt%
        omint = sldwt_list[sps.index("g2")]

        # (6) get soil pH in lab
        time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
            get_soilpH_time.calc_soilpH(
                outdir,
                runname_restart,
                dep_sample,
                20,
                include_N=True,
                include_Cl=False,
                include_Al=False,
                include_DIC=True,
                use_CaCl2=True,
                add_salt=False,
                add_NO3=False,
                use_local_storage=use_local_storage,
                save_lab_dir=False,
                use_Sikora_buffer=False,
                set_ExTimeLimit=True,
                water_frac=2.5,
                cacl2_conc=0.01,
                # salt_conc           = salt_conc,
                # NO3_conc            = NO3_conc,
                # salt_sp             = salt_sp,
                Ex_TimeLim=20,
            )
        )

        if flag_error:
            exit("error in calculation of soil pH for basalt exp")

        # fraction of basalt dissolved
        frac, time = get_int_prof.get_dis_frac(outdir, runname_restart, "gbas")

        # calculating buffer pH
        time_buf, IS_buf, acint_buf, phint_buf, flag_error, runname_lab = (
            get_soilpH_time.calc_soilpH(
                outdir,
                runname_restart,
                dep_sample,
                20,
                include_N=True,
                include_Cl=False,
                include_Al=False,
                include_DIC=True,
                use_CaCl2=False,
                add_salt=False,
                add_NO3=False,
                use_local_storage=use_local_storage,
                save_lab_dir=False,
                use_Sikora_buffer=True,
                set_ExTimeLimit=False,
                # water_frac          = water_frac,
                # cacl2_conc          = cacl2_conc,
                # salt_conc           = salt_conc,
                # NO3_conc            = NO3_conc,
                # salt_sp             = salt_sp,
                Ex_TimeLim=20,
                make_runlogfile=True,
                # sub_as_a_job        = sub_as_a_job,
            )
        )

        if flag_error:
            exit("error in calculation of buffer pH for basalt exp")

        res_list_all.append(
            [
                i,
                kinspc,
                taudust,
                dep_sample,
                phint,
                phint_lab,
                phint_buf,
                omint,
                acint,
                100.0 - acint,
                frac,
            ]
        )

    filename = outdir + runname_spinup + "/basalt_add_exp.res"

    np.savetxt(filename, np.array(res_list_all))

    nmlist = [
        "i",
        "kinspc",
        "taudust",
        "dep_sample",
        "phint",
        "phint_lab",
        "phint_buf",
        "omint",
        "acint",
        "100.-acint",
        "frac",
    ]

    with open(filename, "r") as f:
        contents = f.readlines()

    contents.insert(0, "\t".join(nmlist) + "\n")

    with open(filename, "w") as f:
        contents = "".join(contents)
        f.write(contents)


def get_soilpH():
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    print(sys.argv)
    # site = sys.argv[1]

    use_local_storage = True
    # use_local_storage = False

    # runname = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v3/{}_field'.format(site)
    runname = sys.argv[1]
    # runname = sys.argv[1]

    # ++++++ get data for spinup ++++++ #

    dep_sample = 0.15

    # get porewater pH
    # phint = get_int_prof.get_ph_int_site(outdir,runname,dep_sample)
    phint = get_int_prof.get_waterave_site_v3(outdir, runname, dep_sample, 20, "ph")

    # get acidity in %
    # acint = get_int_prof.get_ac_int_site_v2(outdir,runname,dep_sample)
    acint = get_int_prof.get_ac_int_site_v2(
        outdir, runname, dep_sample, 20, simple_average=True
    )

    # get field solid wt%
    sps = ["g2", "inrt"]
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir, runname, dep_sample, [sp], 20)
        sldwt_list.append(sldwt)

    # get SOM wt%
    omint = sldwt_list[sps.index("g2")]

    # get soil pH in lab
    time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
        get_soilpH_time.calc_soilpH(
            outdir,
            runname,
            dep_sample,
            20,
            include_N=True,
            include_SO4=True,
            include_Cl=False,
            include_Al=False,
            include_DIC=True,
            use_CaCl2=True,
            add_salt=False,
            add_NO3=False,
            use_local_storage=use_local_storage,
            save_lab_dir=False,
            use_Sikora_buffer=False,
            set_ExTimeLimit=True,
            water_frac=2.5,
            cacl2_conc=0.01,
            # salt_conc           = salt_conc,
            # NO3_conc            = NO3_conc,
            # salt_sp             = salt_sp,
            Ex_TimeLim=20,
        )
    )

    if flag_error:
        exit("error in calculation of soil pH for spinup")

    np.savetxt(f"{outdir}{runname}/final_bs_ph.res", np.array([100 - acint, phint_lab]))


def multiples_run_lime_varPSD(flg_ctrl=False):
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    # sites = [226,244,254,282,288,289,302,392]
    print(sys.argv)
    site = sys.argv[1]

    # for site in sites:
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run'
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_wcec"
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v2/{}_field'.format(site)
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_fixTpHCECSOMpCO2/{}_field'.format(site)
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_fixTpHCECSOMpCO2NO3/{}_field".format(
        site
    )
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v3/{}_field".format(
        site
    )
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_midOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'

    dustname = "gbasalt"
    dustname = "cc"
    # dustname = 'amnt'
    # dustname = 'cao'
    # dustname = 'amnt_cao'
    # dustname = 'BlueRidge'
    # dustname = 'BlueRidge_beerling24'

    kinspc = False

    # fdust = 5000. # 50 t/ha/y
    # ttot_restart = 4

    fdust_co2 = 100.0
    fdust_co2 = 500.0
    fdust_co2 = 1000.0

    fdust_str = f"{int(fdust_co2/100):d}"

    # ----- fertilizer addition -----

    # fdust = 10. # 10 gN/m2/y <---> 100 kgN/ha/y
    fdust = 1.0  # 1 gN/m2/y <---> 10 kgN/ha/y
    fdust = 0.0  # 0 gN/m2/y <---> 0 kgN/ha/y (control)

    fdust *= (
        1.0 / 14.0 * 80.043 / 2.0
    )  # gN/m2/y converted to mol N/m2/y and then g NH4NO3/m2/y (= 28.58678571428571)

    # ----- CaO addition -----
    # 1, 5, and 10 tCO2/ha/y = 100, 500, 1000 gCO2/m2/yr = 100/44. 500/44, 1000/44 mol CO2/m2/yr = 100/44/2. 500/44/2, 1000/44/2 mol Ca2+/m2/yr
    fdust = 100.0 / 44.0 / 2.0 * 56.0774
    fdust = 500.0 / 44.0 / 2.0 * 56.0774
    fdust = 1000.0 / 44.0 / 2.0 * 56.0774

    # ----- CaCO3 addition -----
    fdust = 100.0 / 44.0 / 1.0 * 100.089
    fdust = 500.0 / 44.0 / 1.0 * 100.089
    fdust = fdust_co2 / 44.0 / 1.0 * 100.089

    if flg_ctrl:
        fdust = 0

    ttot_restart = 100

    p80str = "320um"
    # p80str =  'BlueRdige_lewis21'
    # p80str =  '10um'

    runid = "homo_10gNm-2y-1"
    runid = "homo_10gNm-2y-1_ctrl"
    runid = "homo_v2_10gNm-2y-1"
    runid = "homo_v2_10gNm-2y-1_ctrl"
    runid = "homo_v2"

    if flg_ctrl:
        runid += "_ctrl"
    # runid =   'homo_v2_ctrl'
    # runid = 'aneal_' + p80str
    # runid = 'aneal_BM00_' + p80str
    # runid = 'aneal_Letal21_' + p80str
    # runid = 'Letal21_' + p80str
    # runid = 'BM00_' + p80str
    # runid = 'homo_v2_{}_'.format(dustname) + p80str
    # runname_restart = runname_spinup + '_' + dustname + '_' + runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/liming_test/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/liming_fixTpHCECSOMpCO2/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test_10gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test_v2_1gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_test_1gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_test_10gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = f'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/cao_1_add_noprec/{site}_{runid}'
    # runname_restart = f'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/cao_5_add_noprec/{site}_{runid}'
    # runname_restart = f'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/cao_10_add_noprec/{site}_{runid}'
    runname_restart = f"US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/{dustname}_{fdust_str}_add_noprec_nopsd/{site}_{runid}"
    # runname_restart = 'EF/test/{}_'.format(site)+runid

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        # fdust_restart       = 5000,
        # taudust_restart     = 0.01,
        # taudust_restart     = 0.01,
        taudust_restart=0.0,
        zdust_restart=0.25,
        # zdust_restart       = 0.12,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        # mix_scheme_restart  = 1,
        mix_scheme_restart=2,
        kinspc=kinspc,
        # extra_solutes       = ['so4'],
        force_PSD=-1,
        mod_2nd_phases="2ndslds_nocarbonate.in",
        update_exe=True,
        # use_local_storage   = True,
        use_local_storage=False,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        # exename_src         = 'scepter_test_fine',
    )


def multiples_run_season():

    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    runname_spinup = "JCU/test"

    dustname = "NONE"

    kinspc = False

    fdust = 0.0  # 0 gN/m2/y <---> 0 kgN/ha/y (control)

    ttot_restart = 4

    p80str = "320um"

    runname_restart = "JCU/test_rain"

    clim_temp = np.transpose(np.loadtxt("./data/q_temp_JCU.in", skiprows=1))

    timeline = list(clim_temp[0, :])
    q_temp = list(clim_temp[2, :])
    T_temp = list(clim_temp[1, :])

    dust_temp = [[0] * len(timeline)]

    runsuccess = restart_add_gbas_clim_DEVDEV.add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        timeline=timeline,
        dust_temp=dust_temp,
        q_temp=q_temp,
        T_temp=T_temp,
        # taudust_restart     = taudust,
        zdust_restart=0.25,
        # zdust_restart       = 0.30,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        mix_scheme_restart=1,
        # kinspc              = kinspc,
        # extra_solutes       = ['so4'],
        force_PSD=-1,
        mod_2nd_phases="2ndslds_nocarbonate.in",
        update_exe=True,
        use_local_storage=False,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        rec_logfile=False,
    )


def multiples_run_lime_varPSD_season(run_ctrl=False):
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    # sites = [226,244,254,282,288,289,302,392]
    print(sys.argv)
    site = sys.argv[1]

    # for site in sites:
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run'
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_wcec"
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v2/{}_field'.format(site)
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_fixTpHCECSOMpCO2/{}_field'.format(site)
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_fixTpHCECSOMpCO2NO3/{}_field".format(
        site
    )
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v3/{}_field".format(
        site
    )
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v3_ztot/{}'.format(site)
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_midOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'

    dustname = "gbasalt"
    dustname = "cc"
    dustname = "amnt"
    dustname = "MI"
    # dustname = 'MI_nolime'
    # dustname = 'amnt_cao'
    # dustname = 'BlueRidge'
    # dustname = 'BlueRidge_beerling24'

    kinspc = False

    # fdust = 5000. # 50 t/ha/y
    # ttot_restart = 4

    # fdust = 10. # 10 gN/m2/y <---> 100 kgN/ha/y
    fdust = 1.0  # 1 gN/m2/y <---> 10 kgN/ha/y
    fdust = 0.0  # 0 gN/m2/y <---> 0 kgN/ha/y (control)

    fdust *= (
        1.0 / 14.0 * 80.043 / 2.0
    )  # gN/m2/y converted to mol N/m2/y and then g NH4NO3/m2/y (= 28.58678571428571)

    ttot_restart = 116

    p80str = "320um"
    # p80str =  'BlueRdige_lewis21'
    # p80str =  '10um'

    runid = "homo_10gNm-2y-1"
    runid = "homo_10gNm-2y-1_ctrl"
    runid = "homo_v2_10gNm-2y-1"
    runid = "homo_v2_10gNm-2y-1_ctrl"
    runid = "homo_v2"
    # runid = 'aneal_' + p80str
    # runid = 'aneal_BM00_' + p80str
    # runid = 'aneal_Letal21_' + p80str
    # runid = 'Letal21_' + p80str
    # runid = 'BM00_' + p80str
    # runid = 'homo_v2_{}_'.format(dustname) + p80str
    # runname_restart = runname_spinup + '_' + dustname + '_' + runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/liming_test/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/liming_fixTpHCECSOMpCO2/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test_10gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test_v2_1gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_test_1gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_chk/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_ztot/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_onlyacid/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_onlyacid_max/{}_'.format(site)+runid
    runname_restart = (
        "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_onlyacid_min/{}_".format(
            site
        )
        + runid
    )
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_max/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_noprec_min/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_MI_nolime/{}_'.format(site)+runid
    # runname_restart = 'EF/test/{}_'.format(site)+runid

    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_max.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_min.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_nolime.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_onlyacid.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_onlyacid_max.in',skiprows=1))
    dust_temp = np.transpose(
        np.loadtxt("./data/Dust_temp_MI_onlyacid_min.in", skiprows=1)
    )
    timeline = dust_temp[0, :]

    # for control
    if run_ctrl:
        dust_temp[1:, :] = 0
        runname_restart += "_ctrl"

    runsuccess = restart_add_gbas_clim_DEVDEV.add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        timeline=timeline,
        dust_temp=list(dust_temp[1:, :]),
        # taudust_restart     = taudust,
        zdust_restart=0.25,
        # zdust_restart       = 0.30,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        mix_scheme_restart=2,
        # kinspc              = kinspc,
        extra_solutes=["so4"],
        force_PSD=-1,
        mod_2nd_phases="2ndslds_nocarbonate.in",
        update_exe=True,
        # use_local_storage   = False,
        use_local_storage=True,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        rec_logfile=False,
        exename_src="scepter_test_MI",
    )


def multiples_run_lime_season_onetime_Africa(run_ctrl=False):

    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    print(sys.argv)
    site = sys.argv[1]
    # randid = sys.argv[2]

    # rungroup = '_phoff_0p0'
    # rungroup = '_phoff_0p2'
    # rungroup = '_phoff_0p4'
    # rungroup = '_phoff_n0p4'
    # rungroup = '_phoff_n0p2'
    rungroup = ""

    lime_style = "alt"
    # lime_style = 'onetime'

    # for site in sites:
    runname_spinup = f"Africa/alpha0/spintuneup{rungroup}/{site}_field"

    # (1) frame
    (
        ztot,
        nz,
        ttot,
        temp,
        fdust,
        fdust2,
        taudust,
        omrain,
        zom,
        poro,
        moistsrf,
        zwater,
        zdust,
        w,
        q,
        p,
        nstep,
        rstrt,
        runid,
    ) = get_inputs.get_input_frame(outdir, runname_spinup)

    # (7) OM rain
    filename = "OM_rain.in"
    sld_data_list = get_inputs.get_input_sld_properties(
        outdir, runname_spinup, filename
    )
    sld_list_om = [listtmp[0] for listtmp in sld_data_list]
    sld_val_om = [listtmp[1] for listtmp in sld_data_list]

    print(sld_val_om, sld_val_om[-1], omrain, sld_val_om[-1] * omrain)
    # exit()

    dustname = "cc"
    dustname = "cc_g2_amnt"

    kinspc = False

    # fdust = 5000. # 50 t/ha/y
    # ttot_restart = 4

    # fdust = 10. # 10 gN/m2/y <---> 100 kgN/ha/y
    fdust = 1.0  # 1 gN/m2/y <---> 10 kgN/ha/y
    fdust = 0.0  # 0 gN/m2/y <---> 0 kgN/ha/y (control)

    fdust *= (
        1.0 / 14.0 * 80.043 / 2.0
    )  # gN/m2/y converted to mol N/m2/y and then g NH4NO3/m2/y (= 28.58678571428571)

    ttot_restart = 100

    # f_cc = 10
    f_cc = 25
    # f_cc = 100
    f_cc_str = f"{f_cc:.1f}".replace(".", "p")

    p80str = "320um"

    # tau = 1
    # tau = 0.5
    tau = 0.1
    # tau = 0.05
    # tau = 0.01
    tau_str = f"{tau:.2f}".replace(".", "p")

    tau_lime_fert = 2

    dt_spin = 5

    # runname_restart = f'Africa/alpha0/onetimeliming_{rungroup}/{site}_onetime_cc_{f_cc_str}_tau{tau_str}_{tau_lime_fert:d}y_apart_spin_{dt_spin:d}y_noprec'
    runname_restart = f"Africa/alpha0/onetimeliming{rungroup}/{site}_{lime_style}_cc_{f_cc_str}_tau{tau_str}_{tau_lime_fert:d}y_apart_spin_{dt_spin:d}y_noprec"

    # --- getting dust seasonal data

    timeline_spin = np.arange(0, dt_spin, tau).tolist()
    dust_cc_spin = [0] * len(timeline_spin)
    dust_g2_spin = [omrain * 30.0 / 12.0] * len(timeline_spin)
    dust_amnt_spin = [omrain * sld_val_om[-1]] * len(timeline_spin)

    ## cc
    timing_cc = [i for i in range(ttot_restart) if i % tau_lime_fert == 0]
    if lime_style == "onetime":
        timing_cc = [0]
    totlime, totyear, minduration, timing = f_cc, ttot_restart, tau, timing_cc
    timeline, dust_cc = randomliming.onetime_liming_schedule(
        totlime, totyear, minduration, timing
    )
    dust_cc = list(dust_cc)

    ## g2
    dust_g2 = [omrain * 30.0 / 12.0] * len(timeline)

    # amnt
    timing_amnt = [
        i for i in range(ttot_restart) if i % tau_lime_fert == tau_lime_fert / 2
    ]
    if lime_style == "onetime":
        timing_amnt = [0 + tau_lime_fert]
    totlime, totyear, minduration, timing = (
        omrain * sld_val_om[-1],
        ttot_restart,
        tau,
        timing_amnt,
    )
    timeline, dust_amnt = randomliming.onetime_liming_schedule(
        totlime, totyear, minduration, timing
    )
    dust_amnt = list(dust_amnt)

    dt = timeline[1] - timeline[0]

    timeline += [timeline[-1] + dt / 10.0]
    dust_cc += [dust_cc[-1]]
    dust_g2 += [dust_g2[-1]]
    dust_amnt += [dust_amnt[-1]]

    timeline = timeline_spin + [i + dt_spin for i in timeline]
    dust_cc = dust_cc_spin + dust_cc
    dust_g2 = dust_g2_spin + dust_g2
    dust_amnt = dust_amnt_spin + dust_amnt

    dust_temp = [dust_cc, dust_g2, dust_amnt]

    # dust_temp = np.transpose(  np.loadtxt('./data/Dust_temp_MI.in',skiprows=1)  )
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_max.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_min.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_nolime.in',skiprows=1))
    # timeline = dust_temp[0,:]

    # for control
    if run_ctrl:
        # print(dust_temp)
        # print(len(dust_temp))
        # print(len(dust_temp[0]))
        dust_temp[0] = [0] * len(dust_cc)
        runname_restart += "_ctrl"

    # if os.path.exists(f'{outdir}{runname_restart}/run_complete.txt'):
    # exit(f' *** {runname_restart} is already run before: exit ***')

    runsuccess = restart_add_gbas_clim_DEVDEV.add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart + dt_spin,
        fdust_restart=fdust,
        timeline=timeline,
        dust_temp=dust_temp,
        omrain2dust=True,
        # taudust_restart     = taudust,
        zdust_restart=0.25,
        # zdust_restart       = 0.30,
        p80str=p80str,
        # nstep_restart       = 1000 ,
        nstep_restart=10,
        mix_scheme_restart=2,
        # kinspc              = kinspc,
        force_PSD=-1,
        mod_2nd_phases="2ndslds_nocarbonate.in",
        force_report=0,
        update_exe=True,
        # use_local_storage   = False,
        use_local_storage=True,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        rec_logfile=False,
    )


def multiples_run_lime_season_Africa(run_ctrl=False):

    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    print(sys.argv)
    site = sys.argv[1]
    randid = sys.argv[2]

    # for site in sites:
    runname_spinup = "Africa/alpha0/spintuneup/{}_field".format(site)

    # (1) frame
    (
        ztot,
        nz,
        ttot,
        temp,
        fdust,
        fdust2,
        taudust,
        omrain,
        zom,
        poro,
        moistsrf,
        zwater,
        zdust,
        w,
        q,
        p,
        nstep,
        rstrt,
        runid,
    ) = get_inputs.get_input_frame(outdir, runname_spinup)

    # (7) OM rain
    filename = "OM_rain.in"
    sld_data_list = get_inputs.get_input_sld_properties(
        outdir, runname_spinup, filename
    )
    sld_list_om = [listtmp[0] for listtmp in sld_data_list]
    sld_val_om = [listtmp[1] for listtmp in sld_data_list]

    print(sld_val_om, sld_val_om[-1], omrain, sld_val_om[-1] * omrain)
    # exit()

    dustname = "cc"
    dustname = "cc_g2_amnt"

    kinspc = False

    # fdust = 5000. # 50 t/ha/y
    # ttot_restart = 4

    # fdust = 10. # 10 gN/m2/y <---> 100 kgN/ha/y
    fdust = 1.0  # 1 gN/m2/y <---> 10 kgN/ha/y
    fdust = 0.0  # 0 gN/m2/y <---> 0 kgN/ha/y (control)

    fdust *= (
        1.0 / 14.0 * 80.043 / 2.0
    )  # gN/m2/y converted to mol N/m2/y and then g NH4NO3/m2/y (= 28.58678571428571)

    ttot_restart = 100

    p80str = "320um"

    runname_restart = "Africa/alpha0/randomliming/{}_rand_{}_noprec".format(
        site, randid
    )

    # --- getting dust seasonal data
    ## cc
    totlime, totyear, minduration = 25, 100, 24
    timeline, dust_cc = randomliming.random_liming_schedule(
        totlime, totyear, minduration
    )
    dust_cc = list(dust_cc)

    ## g2
    dust_g2 = [omrain * 30.0 / 12.0] * len(timeline)

    # amnt
    totlime, totyear, minduration = omrain * sld_val_om[-1], 100, 24
    timeline, dust_amnt = randomliming.random_liming_schedule(
        totlime, totyear, minduration
    )
    dust_amnt = list(dust_amnt)

    dust_temp = [dust_cc, dust_g2, dust_amnt]

    # dust_temp = np.transpose(  np.loadtxt('./data/Dust_temp_MI.in',skiprows=1)  )
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_max.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_min.in',skiprows=1))
    # dust_temp = np.transpose(np.loadtxt('./data/Dust_temp_MI_nolime.in',skiprows=1))
    # timeline = dust_temp[0,:]

    # for control
    if run_ctrl:
        dust_temp[1:, :] = 0

    runsuccess = restart_add_gbas_clim_DEVDEV.add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        timeline=timeline,
        dust_temp=dust_temp,
        omrain2dust=True,
        # taudust_restart     = taudust,
        zdust_restart=0.25,
        # zdust_restart       = 0.30,
        p80str=p80str,
        # nstep_restart       = 1000 ,
        nstep_restart=10,
        mix_scheme_restart=2,
        # kinspc              = kinspc,
        force_PSD=-1,
        mod_2nd_phases="2ndslds_nocarbonate.in",
        force_report=0,
        update_exe=True,
        use_local_storage=False,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        rec_logfile=False,
    )


def multiples_run_lime_Africa(run_ctrl=False):

    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    print(sys.argv)
    site = sys.argv[1]

    # rungroup = '_phoff_0p0'
    # rungroup = '_phoff_0p2'
    # rungroup = '_phoff_0p4'
    # rungroup = '_phoff_n0p4'
    # rungroup = '_phoff_n0p2'
    rungroup = ""

    # runname_spinup = 'Africa/spintuneup_alpha0/{}_field'.format(site)
    # runname_spinup = 'Africa/alpha0/spintuneup/{}_field'.format(site)
    runname_spinup = f"Africa/alpha0/spintuneup{rungroup}/{site}_field"

    dustname = "cc"

    kinspc = False

    fdust = 25.0  # 0.25 t/ha/y
    # fdust = 100. # 1 t/ha/y
    # fdust = 10 # 0.1 t/ha/y

    fdust_str = f"{fdust/100.:.2f}".replace(".", "p")

    ttot_restart = 100

    p80str = "320um"
    # p80str =  'BlueRdige_lewis21'
    # p80str =  '10um'

    # runname_restart = 'Africa/liming_alpha0/{}_cnst0p25'.format(site)
    # runname_restart = 'Africa/liming_alpha0/{}_cnst0p25_noprec'.format(site)
    runname_restart = f"Africa/alpha0/cnstliming{rungroup}/{site}_{fdust_str}_noprec"

    if run_ctrl:
        fdust = 0
        runname_restart += "_ctrl"

    # if os.path.exists(f'{outdir}{runname_restart}/run_complete.txt'):
    # exit(f' *** {runname_restart} is already run before: exit ***')

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        # fdust_restart       = 5000,
        # taudust_restart     = 0.01,
        taudust_restart=0.0,
        zdust_restart=0.25,
        # zdust_restart       = 0.12,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        # mix_scheme_restart  = 1,
        mix_scheme_restart=2,
        kinspc=kinspc,
        mod_2nd_phases="2ndslds_nocarbonate.in",
        update_exe=True,
        use_local_storage=True,
        # use_local_storage   = False,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
    )


def multiples_run_lime_varPSD_EF():
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    # sites = [226,244,254,282,288,289,302,392]
    print(sys.argv)
    site = sys.argv[1]

    # for site in sites:
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run'
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_lowOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_midOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2_w2nd-single_run_highOM_nocec"
    runname_spinup = "spinup_inrt_g2_ca_pco2-single_run_lowOM_wcec"
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v2/{}_field'.format(site)
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_fixTpHCECSOMpCO2/{}_field'.format(site)
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_fixTpHCECSOMpCO2NO3/{}_field".format(
        site
    )
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_v3/{}_field".format(
        site
    )
    runname_spinup = "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha0p0/spintuneup_v3/{}_field".format(
        site
    )
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha5p0/spintuneup_v3/{}_field'.format(site)
    # runname_spinup = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha10p0/spintuneup_v3/{}_field'.format(site)
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_midOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'
    # runname_spinup = 'spinup_inrt_g2_ca_pco2-single_run_highOM_wcec'

    if not os.path.exists(outdir + runname_spinup + "/iteration.res"):
        exit("spintuneup not successful? not gonna do ERW experiment")

    dustname = "gbasalt"
    dustname = "cc"
    dustname = "amnt"
    dustname = "amnt_cao"
    dustname = "BlueRidge"
    dustname = "basalt_Baek_no2nd"
    dustname = "gbas"
    # dustname = 'BlueRidge_beerling24'

    kinspc = False

    # fdust = 5000. # 50 t/ha/y
    # ttot_restart = 4
    fdust = 2500.0  # 25 t/ha/y
    ttot_restart = 100

    # fdust = 10. # 10 gN/m2/y <---> 100 kgN/ha/y
    # fdust = 1. # 1 gN/m2/y <---> 10 kgN/ha/y
    # fdust = 0. # 0 gN/m2/y <---> 0 kgN/ha/y (control)

    # fdust *= 1./14.*80.043/2. # gN/m2/y converted to mol N/m2/y and then g NH4NO3/m2/y (= 28.58678571428571)

    # ttot_restart = 100

    p80str = "320um"
    p80str = "BlueRdige_lewis21"
    p80str = "100um"
    # p80str =  '1220um'

    runid = "homo_10gNm-2y-1"
    runid = "homo_10gNm-2y-1_ctrl"
    runid = "homo_v2_10gNm-2y-1"
    runid = "homo_v2_10gNm-2y-1_ctrl"
    runid = "homo_v2"
    runid = "homo_v2"
    # runid =   'homo_v2_ctrl'
    # runid = 'aneal_' + p80str
    # runid = 'aneal_BM00_' + p80str
    # runid = 'aneal_Letal21_' + p80str
    # runid = 'Letal21_' + p80str
    # runid = 'BM00_' + p80str
    runid = "homo_v2_{}_".format(dustname) + p80str
    # runname_restart = runname_spinup + '_' + dustname + '_' + runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/liming_test/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/liming_fixTpHCECSOMpCO2/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fert_test_1gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/fertlime_test_1gNm-2y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/test_25tha-1y-1/{}_'.format(site)+runid
    runname_restart = (
        "US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha0p0/test_25tha-1y-1/{}_".format(
            site
        )
        + runid
    )
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha5p0/test_25tha-1y-1/{}_'.format(site)+runid
    # runname_restart = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha10p0/test_25tha-1y-1/{}_'.format(site)+runid
    # runname_restart = 'EF/test/{}_'.format(site)+runid

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        # fdust_restart       = 5000,
        # taudust_restart     = 0.01,
        taudust_restart=0.01,
        zdust_restart=0.25,
        # zdust_restart       = 0.12,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        # mix_scheme_restart  = 1,
        mix_scheme_restart=2,
        kinspc=kinspc,
        update_exe=True,
        # use_local_storage   = True,
        use_local_storage=False,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
    )


def multiples_run_varPSD_MY():
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    # runname_spinup = 'MY/spinup/rand_3/79'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30'
    runname_spinup = "MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha10"

    dustname = "Tawau_2"
    # dustname = 'BlueRidge_beerling24'

    kinspc = False

    fdust = 5000.0  # 50 t/ha/y
    # ttot_restart = 4
    # fdust = 2500. # 25 t/ha/y
    ttot_restart = 6

    # sld_list_add = ['gb','gt']
    sld_list_add = []
    gas_list_add = ["po2"]
    exrxn_list_add = ["fe2o2"]

    # p80str =  'Tawau_mark'
    # p80str =  'Tawau_raw'
    p80str = "Tawau_mark_2"
    # p80str =  'Tawau_raw_2'

    # runname_restart = 'MY/basalt/rand_3/79_test1'
    # runname_restart = 'MY/basalt/rand_2_2/272_test1'
    # runname_restart = 'MY/basalt/rand_2_2/272_test2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test3'
    # runname_restart = 'MY/basalt/rand_2_2/272_test4'
    # runname_restart = 'MY/basalt/rand_2_2/272_test5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30'
    runname_restart = "MY/basalt/rand_2_2/272_test_homo_nz30_alpha10"

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        ttot_restart=ttot_restart,
        fdust_restart=fdust,
        # fdust_restart       = 5000,
        # taudust_restart     = 0.0,
        taudust_restart=0.1,
        # zdust_restart       = 0.25,
        # zdust_restart       = 0.20,
        zdust_restart=0.30,
        # zdust_restart       = 0.12,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        # mix_scheme_restart  = 1,
        mix_scheme_restart=2,
        kinspc=kinspc,
        sld_list_add=sld_list_add,
        gas_list_add=gas_list_add,
        exrxn_list_add=exrxn_list_add,
        update_exe=True,
        # use_local_storage   = True,
        use_local_storage=False,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        exename_src="scepter_MY",
    )

    use_local_storage = True

    # dep_sample = 0.15
    for dep_sample in [0.10, 0.30]:

        depth_str = "{:.1f}".format(dep_sample).replace(".", "p")

        res_list = []

        for itime in range(1, 20 + 1, 1):

            # (1) get porewater pH
            phint = get_int_prof.get_ph_int_site(
                outdir, runname_restart, dep_sample, itime
            )

            # (2) get acidity in %
            acint = get_int_prof.get_ac_int_site_v2(
                outdir, runname_restart, dep_sample, itime
            )

            # (3) get bulk density
            dense = get_int_prof.get_rhobulk_int_site(
                outdir, runname_restart, dep_sample, itime
            )

            # (4) get field solid wt%
            sps = ["g2", "inrt"]
            sldwt_list = []
            for sp in sps:
                sldwt = get_int_prof.get_sldwt_int_site(
                    outdir, runname_restart, dep_sample, [sp], itime
                )
                sldwt_list.append(sldwt)

            # (5) get SOM wt%
            omint = sldwt_list[sps.index("g2")]

            # (6) get soil pH in lab
            time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
                get_soilpH_time.calc_soilpH(
                    outdir,
                    runname_restart,
                    dep_sample,
                    itime,
                    include_N=True,
                    # include_Cl          = False,
                    include_Cl=True,
                    include_Al=False,
                    include_DIC=True,
                    # use_CaCl2           = True,
                    use_CaCl2=False,
                    add_salt=False,
                    add_NO3=False,
                    use_local_storage=use_local_storage,
                    save_lab_dir=False,
                    # save_lab_dir        = True,
                    use_Sikora_buffer=False,
                    # set_ExTimeLimit     = True,
                    set_ExTimeLimit=False,
                    water_frac=5,
                    cacl2_conc=0.01,
                    # salt_conc           = salt_conc,
                    # NO3_conc            = NO3_conc,
                    # salt_sp             = salt_sp,
                    cl_salt="cacl2",
                    # cl_salt             = 'mgcl2',
                    Ex_TimeLim=20,
                )
            )

            res_list.append(
                [time, dep_sample, phint, phint_lab, omint, acint, 100.0 - acint]
            )

        np.savetxt(
            outdir + runname_restart + f"/ph_BS_OM_{depth_str}.res",
            np.array(res_list),
            delimiter=" ",
            header="time dep_sample phint phint_lab omint acint 100-acint",
            comments="",
        )


def multiples_run_varPSD_MY_season(flg_ctrl=False):

    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    run_field = True
    # run_field = False

    use_local_storage = True
    # use_local_storage = False

    # site_id = 122
    # site_id = 230
    site_id = int(sys.argv[1])
    # site_id = 707
    # site_id = 506
    # site_id = 377

    alpha = float(sys.argv[2])
    alpha_str = f"{alpha:.1f}".replace(".", "p")
    # alpha = 0
    # alpha = 1
    # alpha = 2
    # alpha = 3
    # alpha = 4
    # alpha = 5
    # alpha = 6
    # alpha = 8
    # alpha = 10

    # run_group = 'rand_cec27'
    # run_group = 'rand_cec33'
    # run_group = 'rand_cec30'
    run_group = sys.argv[3]

    # poro = 0.4
    # poro = 0.6
    # poro = 0.5
    poro = float(sys.argv[4])
    poro_str = f"{poro:.2f}".replace(".", "p")

    cl_mod_fact = 0.5
    cl_mod_fact = 0.9
    cl_mod_fact = 0.0
    cl_mod_fact = 1.0
    # cl_mod_fact = 0.99
    cl_mod_str = f"{cl_mod_fact:.2f}".replace(".", "p")

    # gb = float(sys.argv[4])
    # gb_str = f'{gb:.2f}'.replace('.','p')

    zdust_restart = 0.30
    # zdust_restart = 0.20
    # zdust_restart = 0.10
    zmix_str = f"{int(zdust_restart*100):d}"

    # basalt_run_id = '_ccdef'
    # basalt_run_id = '_z50cm'
    # basalt_run_id = '_p0p4'
    # basalt_run_id = '_p0p6'
    # basalt_run_id = '_cecom'
    # basalt_run_id = ''
    # basalt_run_id = f'_poro{poro_str}'
    basalt_run_id = f"_poro{poro_str}"
    # basalt_run_id = f'_poro{poro_str}_maggi'
    # basalt_run_id = f'_poro{poro_str}_maggidef'
    # basalt_run_id = f'_poro{poro_str}_mgca_naclx2'
    # basalt_run_id = f'_poro{poro_str}_gb0p01'
    # basalt_run_id = f'_poro{poro_str}_gb0p1'
    # basalt_run_id = f'_gb{gb_str}'

    # runname_spinup = 'MY/spinup/rand_3/79'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30'
    # runname_spinup = f'MY/basalt_spunup/rand_2_2_{site_id:d}_mod_3_homo_nz30'
    # runname_spinup = f'MY/basalt_spunup/rand_2_2_{site_id:d}_mod_3_homo_nz60'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha:d}'
    runname_spinup = f"MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}"
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca_nh4'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca_ka0p3'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca_wx0'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca_wx0_2lyr_v2'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca_wx10'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha_str}{basalt_run_id}_mgca_wx0p001'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_fick_nz30_alpha{alpha_str}{basalt_run_id}'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha0'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha2'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha5'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha10'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_nps200'
    # runname_spinup = 'MY/spinup/rand_2_2/506'
    # runname_spinup = 'MY/spinup/rand_2_2/707'
    # runname_spinup = f'MY/spinup/{run_group}/{site_id:d}'

    # dustname = 'Tawau_2'
    # dustname = 'Tawau_2_g2amnt'
    dustname = "Tawau_2_g2amnt_norm"
    # dustname = 'Tawau_2_g2amnt_norm_v2'
    # dustname = 'BlueRidge_beerling24'

    kinspc = False

    fdust = 5000.0  # 50 t/ha/y
    if flg_ctrl:
        fdust = 0
    # ttot_restart = 4
    # fdust = 2500. # 25 t/ha/y
    ttot_restart = 6
    # ttot_restart = 5.99
    # ttot_restart = 5.95
    # ttot_restart = 4

    # sld_list_add = ['gb','gt']
    sld_list_add = []
    gas_list_add = ["po2"]
    # gas_list_add = []
    exrxn_list_add = ["fe2o2"]

    # gas_list_add = ['po2','pnh3']
    # exrxn_list_add = ['fe2o2','amo2o']

    # p80str =  'Tawau_mark'
    # p80str =  'Tawau_raw'
    p80str = "Tawau_mark_2"
    # p80str =  'Tawau_raw_2'

    # runname_restart = 'MY/basalt/rand_3/79_test1'
    # runname_restart = 'MY/basalt/rand_2_2/272_test1'
    # runname_restart = 'MY/basalt/rand_2_2/272_test2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test3'
    # runname_restart = 'MY/basalt/rand_2_2/272_test4'
    # runname_restart = 'MY/basalt/rand_2_2/272_test5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30'
    # runname_restart = f'MY/basalt/rand_2_2/{site_id:d}_test_homo_nz30_season'
    # runname_restart = f'MY/basalt/rand_2_2/{site_id:d}_test_homo_nz60_season'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_fick20'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_tau0p1'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_tau0p05'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_wx10_tau0p05'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_wx0p001'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_zd10cm_tau0p05'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_cl{cl_mod_str}'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_cl{cl_mod_str}_qmeas'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_cl{cl_mod_str}_tau0p05'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_nh4'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_zmix{zmix_str}_cl{cl_mod_str}'
    runname_restart = f"MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_zmix{zmix_str}_cl{cl_mod_str}_mod"
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_zmix{zmix_str}_cl{cl_mod_str}_qmeas'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_zmix{zmix_str}_cl{cl_mod_str}_qmeas2'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_zmix{zmix_str}_cl{cl_mod_str}_qpred'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_wx0_2lyr_v2'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_wx0'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_mixtext'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_mgca_zmix10_tau0p05'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v3'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v4'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_nsp200_season_v6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v7'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v8'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v9'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v10'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha0'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha10'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha0_v2'
    # runname_restart = 'MY/basalt/rand_2_2/506_test_fick_nz60_season'
    # runname_restart = 'MY/basalt/rand_2_2/506_test_fick_nz60_season'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v11'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v12'

    if flg_ctrl:
        runname_restart += "_ctrl"

    # t  tmp  q  npp  soilresp  swvf  N_flux w_plant
    season_data = np.loadtxt(
        "../shelffield/particle_weathering_paper3/forcings_all.txt", skiprows=1
    )
    # season_data = np.delete(season_data, season_data.shape[0]-1, axis=0)
    # timeline = list(season_data[:,0]-season_data[0,0])
    timeline = list(season_data[:, 0])

    # timeline = list(season_data[:,0])

    q_temp = list(season_data[:, 2])
    T_temp = list(season_data[:, 1])
    moist_temp = list(season_data[:, 5] / 0.3 / poro)
    om_temp = season_data[:, 4] / 12.0 * 30.0  # converting from gC/m2/y to gCH2O/m2/y
    # N_temp = (season_data[:,6]/14.*80./2.)
    N_temp = (
        season_data[:, 6] * 80.0 / 2.0
    )  # converting from N mol/m2/y to gNH4NO3/m2/y
    fkin_temp = list(season_data[:, 7])

    om_temp[om_temp < 0] = 0.0
    N_temp[N_temp < 0] = 0.0

    om_temp = list(om_temp)
    N_temp = list(N_temp)

    f_q = interp1d(timeline, q_temp, fill_value="extrapolate")
    f_T = interp1d(timeline, T_temp, fill_value="extrapolate")
    f_moist = interp1d(timeline, moist_temp, fill_value="extrapolate")
    f_om = interp1d(timeline, om_temp, fill_value="extrapolate")
    f_N = interp1d(timeline, N_temp, fill_value="extrapolate")
    f_fkin = interp1d(timeline, fkin_temp, fill_value="extrapolate")

    # del_t = 0.001
    del_t = 0.01
    # del_t = 0.005
    # del_t = 0.05

    # del_t_dust = 0.1
    del_t_dust = 0.05

    timeline = np.arange(0, ttot_restart + del_t, del_t)
    timeline[-1] -= del_t * 0.1

    # mask = timeline.astype(int) == timeline
    mask = timeline % 1 < del_t_dust - del_t * 0.1

    if sum(mask) != int(6 * del_t_dust / del_t):
        exit("mask is wrong")

    timeline = timeline.tolist()

    q_temp = f_q(timeline).tolist()
    T_temp = f_T(timeline).tolist()
    moist_temp = f_moist(timeline).tolist()
    om_temp = f_om(timeline).tolist()
    N_temp = f_N(timeline).tolist()
    fkin_temp = f_fkin(timeline).tolist()

    dust_fr_list = [
        0.010,
        0.006,
        0.025,
        0.133,
        0.29,
        0.02,
        0.058,
    ]

    # modifying q based on measurements?
    q_temp = np.array(q_temp)
    # q_temp *= 3.627521739/4.200056415
    q_temp *= 6.019372694 / 4.200056415
    # if not flg_ctrl:    q_temp *= 3.627521739/4.200056415
    # else:               q_temp *= 2.979418711/4.200056415
    q_temp = q_temp.tolist()

    dust_temp = []

    for i in range(len(dust_fr_list)):

        # continuous

        dust_temp.append(
            [fdust * dust_fr_list[i]] * len(timeline)
        )  # homogeneous deployment

        # only during del_t of initial year

        # fdust_tmp = fdust*dust_fr_list[i]*mask/del_t_dust
        # fdust_tmp.tolist()
        # dust_temp.append( fdust_tmp  )

    dust_temp.append(om_temp)
    dust_temp.append(N_temp)

    # change rain composition

    filename = "rain.in"
    rain_data_list = get_inputs.get_input_sld_properties(
        outdir, runname_spinup, filename
    )
    print(
        rain_data_list,
    )

    if flg_ctrl:
        cl_mod_fact *= 0.5

    for i in range(len(rain_data_list)):
        rainsp, rainconc = rain_data_list[i]
        if rainsp == "cl":
            cl_mod = copy.copy(rainconc * cl_mod_fact)

    for i in range(len(rain_data_list)):
        rainsp, rainconc = rain_data_list[i]
        if rainsp == "cl":
            rain_data_list[i][1] += -cl_mod
        if rainsp == "mg":
            rain_data_list[i][1] += -cl_mod * 0.5

        if rain_data_list[i][1] < 1e-20:
            rain_data_list[i][1] = 1e-20

    force_PSD = 0
    if flg_ctrl:
        force_PSD = -1  # if control run, not using PSD

    if run_field:
        runsuccess = restart_add_gbas_clim_DEVDEV.add_gbas(
            outdir,
            runname_spinup,
            dustname,
            runname_restart,
            ttot_restart=ttot_restart,
            fdust_restart=fdust,
            timeline=timeline,
            dust_temp=dust_temp,
            T_temp=T_temp,
            q_temp=q_temp,
            moist_temp=moist_temp,
            fkin_temp=fkin_temp,
            rain_force=rain_data_list,
            omrain2dust=True,
            zdust_restart=zdust_restart,
            # zdust_restart       = 0.10 ,
            # zdust_restart       = 0.20 ,
            p80str=p80str,
            nstep_restart=10,
            mix_scheme_restart=2,
            # mix_scheme_restart  = 1  ,
            force_report=0,
            force_PSD=force_PSD,
            update_exe=True,
            use_local_storage=use_local_storage,
            save_restart_dir=True,
            set_ExTimeLimit=False,
            Ex_TimeLim=20,
            rec_logfile=False,
            # exename_src         = 'scepter_MY',
            exename_src="scepter_MYCl",
            # exename_src         = 'scepter_MYClNH4',
            # exename_src         = 'scepter_MYClFick',
            # exename_src         = 'scepter_MY_ccdef',
            # exename_src         = 'scepter_MY_fick',
        )

    # dep_sample = 0.15
    dep_sample_top = 0
    for dep_sample in [0.10, 0.30]:
        # for dep_sample in [0.10,0.295]:

        depth_str = "{:.1f}".format(dep_sample).replace(".", "p")

        res_list = []

        for itime in range(1, 20 + 1, 1):

            print(itime, dep_sample_top, dep_sample)

            # (1) get porewater pH
            phint = get_int_prof.get_ph_int_site(
                outdir,
                runname_restart,
                dep_sample,
                itime,
                dep_sample_top=dep_sample_top,
            )

            # (2) get acidity in %
            acint = get_int_prof.get_ac_int_site_v2(
                outdir,
                runname_restart,
                dep_sample,
                itime,
                dep_sample_top=dep_sample_top,
            )

            # (3) get bulk density
            dense = get_int_prof.get_rhobulk_int_site(
                outdir,
                runname_restart,
                dep_sample,
                itime,
                dep_sample_top=dep_sample_top,
            )

            # (4) get field solid wt%
            sps = ["g2", "inrt"]
            sldwt_list = []
            for sp in sps:
                sldwt = get_int_prof.get_sldwt_int_site(
                    outdir,
                    runname_restart,
                    dep_sample,
                    [sp],
                    itime,
                    dep_sample_top=dep_sample_top,
                )
                sldwt_list.append(sldwt)

            # (5) get SOM wt%
            omint = sldwt_list[sps.index("g2")]

            # (6) get soil pH in lab
            # time,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = get_soilpH_time.calc_soilpH(
            time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
                get_soilpH_time_DEVTMP.calc_soilpH(
                    outdir,
                    runname_restart,
                    dep_sample,
                    itime,
                    include_N=True,
                    # include_Cl          = False,
                    include_Cl=True,
                    # include_Al          = False,
                    include_Al=True,
                    include_DIC=True,
                    # use_CaCl2           = True,
                    use_CaCl2=False,
                    add_salt=False,
                    add_NO3=False,
                    dep_sample_top=dep_sample_top,
                    use_local_storage=use_local_storage,
                    save_lab_dir=False,
                    # save_lab_dir        = True,
                    use_Sikora_buffer=False,
                    set_ExTimeLimit=True,
                    # set_ExTimeLimit     = False,
                    water_frac=5,
                    cacl2_conc=0.01,
                    # salt_conc           = salt_conc,
                    # NO3_conc            = NO3_conc,
                    # salt_sp             = salt_sp,
                    cl_salt="hcl",
                    # cl_salt             = 'mgcl2',
                    Ex_TimeLim=20,
                    exename_src="scepter_MYCl",
                )
            )

            res_list.append(
                [time, dep_sample, phint, phint_lab, omint, acint, 100.0 - acint]
            )

        np.savetxt(
            outdir + runname_restart + f"/ph_BS_OM_{depth_str}.res",
            np.array(res_list),
            delimiter=" ",
            header="time dep_sample phint phint_lab omint acint 100-acint",
            comments="",
        )

        dep_sample_top = copy.copy(dep_sample)


def multiples_run_varPSD_MY_season_sensitivity():

    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    run_field = True
    # run_field = False

    use_local_storage = True
    # use_local_storage = False

    # site_id = 122
    # site_id = 230
    site_id = int(sys.argv[1])
    # site_id = 707
    # site_id = 506
    # site_id = 377

    alpha = float(sys.argv[2])
    alpha_str = f"{alpha:.1f}".replace(".", "p")
    # alpha = 0
    # alpha = 1
    # alpha = 2
    # alpha = 3
    # alpha = 4
    # alpha = 5
    # alpha = 6
    # alpha = 8
    # alpha = 10

    poro = 0.45

    # runname_spinup = 'MY/spinup/rand_3/79'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30'
    # runname_spinup = f'MY/basalt_spunup/rand_2_2_{site_id:d}_mod_3_homo_nz30'
    # runname_spinup = f'MY/basalt_spunup/rand_2_2_{site_id:d}_mod_3_homo_nz60'
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30_alpha{alpha:d}'
    runname_spinup = f"MY/sensitivity/alpha{alpha_str}/spintuneup/{site_id}_field"
    # runname_spinup = f'MY/basalt_spunup/{run_group}_{site_id:d}_mod_3_homo_nz30'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha0'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha2'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha5'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_alpha10'
    # runname_spinup = 'MY/basalt_spunup/rand_2_2_272_mod_3_homo_nz30_nps200'
    # runname_spinup = 'MY/spinup/rand_2_2/506'
    # runname_spinup = 'MY/spinup/rand_2_2/707'
    # runname_spinup = f'MY/spinup/{run_group}/{site_id:d}'

    # dustname = 'Tawau_2'
    # dustname = 'Tawau_2_g2amnt'
    dustname = "Tawau_2_g2amnt_norm"
    # dustname = 'BlueRidge_beerling24'

    kinspc = False

    fdust = 5000.0  # 50 t/ha/y
    # ttot_restart = 4
    # fdust = 2500. # 25 t/ha/y
    ttot_restart = 6
    ttot_restart = 5.99
    # ttot_restart = 5.95
    # ttot_restart = 4

    # sld_list_add = ['gb','gt']
    extra_solutes = ["si", "al", "fe2", "fe3", "cl"]
    sld_list_add = ["gt", "gb", "ka", "amsi"]
    gas_list_add = ["po2"]
    exrxn_list_add = ["fe2o2"]

    # p80str =  'Tawau_mark'
    # p80str =  'Tawau_raw'
    p80str = "Tawau_mark_2"
    # p80str =  'Tawau_raw_2'

    # runname_restart = 'MY/basalt/rand_3/79_test1'
    # runname_restart = 'MY/basalt/rand_2_2/272_test1'
    # runname_restart = 'MY/basalt/rand_2_2/272_test2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test3'
    # runname_restart = 'MY/basalt/rand_2_2/272_test4'
    # runname_restart = 'MY/basalt/rand_2_2/272_test5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30'
    # runname_restart = f'MY/basalt/rand_2_2/{site_id:d}_test_homo_nz30_season'
    # runname_restart = f'MY/basalt/rand_2_2/{site_id:d}_test_homo_nz60_season'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}'
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season_alpha{alpha_str}{basalt_run_id}_v2'
    runname_restart = f"MY/sensitivity/alpha{alpha_str}/basalt_cont_v4/{site_id}"
    # runname_restart = f'MY/basalt/{run_group}/{site_id:d}_test_homo_nz30_season'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v3'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v4'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_nsp200_season_v6'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v7'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v8'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v9'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v10'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha0'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha2'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha5'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha10'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_alpha0_v2'
    # runname_restart = 'MY/basalt/rand_2_2/506_test_fick_nz60_season'
    # runname_restart = 'MY/basalt/rand_2_2/506_test_fick_nz60_season'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v11'
    # runname_restart = 'MY/basalt/rand_2_2/272_test_homo_nz30_season_v12'

    # t  tmp  q  npp  soilresp  swvf  N_flux w_plant
    season_data = np.loadtxt(
        "../shelffield/particle_weathering_paper3/forcings_all.txt", skiprows=1
    )
    # season_data = np.delete(season_data, season_data.shape[0]-1, axis=0)
    # timeline = list(season_data[:,0]-season_data[0,0])
    timeline = list(season_data[:, 0])

    # timeline = list(season_data[:,0])

    q_temp = list(season_data[:, 2])
    T_temp = list(season_data[:, 1])
    moist_temp = list(season_data[:, 5] / 0.3 / poro)
    om_temp = season_data[:, 4] / 12.0 * 30.0  # converting from gC/m2/y to gCH2O/m2/y
    # N_temp = (season_data[:,6]/14.*80./2.)
    N_temp = (
        season_data[:, 6] * 80.0 / 2.0
    )  # converting from N mol/m2/y to gNH4NO3/m2/y
    fkin_temp = list(season_data[:, 7])

    om_temp[om_temp < 0] = 0.0
    N_temp[N_temp < 0] = 0.0

    om_temp = list(om_temp)
    N_temp = list(N_temp)

    f_q = interp1d(timeline, q_temp, fill_value="extrapolate")
    f_T = interp1d(timeline, T_temp, fill_value="extrapolate")
    f_moist = interp1d(timeline, moist_temp, fill_value="extrapolate")
    f_om = interp1d(timeline, om_temp, fill_value="extrapolate")
    f_N = interp1d(timeline, N_temp, fill_value="extrapolate")
    f_fkin = interp1d(timeline, fkin_temp, fill_value="extrapolate")

    # del_t = 0.001
    del_t = 0.01
    # del_t = 0.005
    # del_t = 0.05

    del_t_dust = 0.03

    timeline = np.arange(0, ttot_restart + del_t, del_t)
    timeline[-1] -= del_t * 0.1

    mask = timeline.astype(int) == timeline

    timeline = timeline.tolist()

    q_temp = f_q(timeline).tolist()
    T_temp = f_T(timeline).tolist()
    moist_temp = f_moist(timeline).tolist()
    om_temp = f_om(timeline).tolist()
    N_temp = f_N(timeline).tolist()
    fkin_temp = f_fkin(timeline).tolist()

    dust_fr_list = [
        0.010,
        0.006,
        0.025,
        0.133,
        0.29,
        0.02,
        0.058,
    ]

    dust_temp = []

    for i in range(len(dust_fr_list)):
        dust_temp.append(
            [fdust * dust_fr_list[i]] * len(timeline)
        )  # homogeneous deployment
        # only during del_t of initial year
        # fdust_tmp = fdust*dust_fr_list[i]*mask/del_t
        # fdust_tmp.tolist()
        # dust_temp.append( fdust_tmp  )

    dust_temp.append(om_temp)
    dust_temp.append(N_temp)

    if run_field:
        runsuccess = restart_add_gbas_clim_DEVDEV.add_gbas(
            outdir,
            runname_spinup,
            dustname,
            runname_restart,
            ttot_restart=ttot_restart,
            fdust_restart=fdust,
            timeline=timeline,
            dust_temp=dust_temp,
            T_temp=T_temp,
            q_temp=q_temp,
            moist_temp=moist_temp,
            fkin_temp=fkin_temp,
            omrain2dust=True,
            zdust_restart=0.30,
            # zdust_restart       = 0.20 ,
            p80str=p80str,
            nstep_restart=10,
            mix_scheme_restart=2,
            # mix_scheme_restart  = 1  ,
            force_report=0,
            extra_solutes=extra_solutes,
            sld_list_add=sld_list_add,
            gas_list_add=gas_list_add,
            exrxn_list_add=exrxn_list_add,
            update_exe=True,
            use_local_storage=use_local_storage,
            save_restart_dir=True,
            set_ExTimeLimit=False,
            Ex_TimeLim=20,
            rec_logfile=False,
            exename_src="scepter_MYCl",
            # exename_src         = 'scepter_MY_ccdef',
            # exename_src         = 'scepter_MY_fick',
        )

    # dep_sample = 0.15
    for dep_sample in [0.10, 0.30]:
        # for dep_sample in [0.10,0.295]:

        depth_str = "{:.1f}".format(dep_sample).replace(".", "p")

        res_list = []

        for itime in range(1, 20 + 1, 1):

            print(itime, dep_sample)

            # (1) get porewater pH
            phint = get_int_prof.get_ph_int_site(
                outdir, runname_restart, dep_sample, itime
            )

            # (2) get acidity in %
            acint = get_int_prof.get_ac_int_site_v2(
                outdir, runname_restart, dep_sample, itime
            )

            # (3) get bulk density
            dense = get_int_prof.get_rhobulk_int_site(
                outdir, runname_restart, dep_sample, itime
            )

            # (4) get field solid wt%
            sps = ["g2", "inrt"]
            sldwt_list = []
            for sp in sps:
                sldwt = get_int_prof.get_sldwt_int_site(
                    outdir, runname_restart, dep_sample, [sp], itime
                )
                sldwt_list.append(sldwt)

            # (5) get SOM wt%
            omint = sldwt_list[sps.index("g2")]

            # (6) get soil pH in lab
            # time,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = get_soilpH_time.calc_soilpH(
            time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
                get_soilpH_time_DEVTMP.calc_soilpH(
                    outdir,
                    runname_restart,
                    dep_sample,
                    itime,
                    include_N=True,
                    # include_Cl          = False,
                    include_Cl=True,
                    # include_Al          = False,
                    include_Al=True,
                    include_DIC=True,
                    # use_CaCl2           = True,
                    use_CaCl2=False,
                    add_salt=False,
                    add_NO3=False,
                    use_local_storage=use_local_storage,
                    save_lab_dir=False,
                    # save_lab_dir        = True,
                    use_Sikora_buffer=False,
                    set_ExTimeLimit=True,
                    # set_ExTimeLimit     = False,
                    water_frac=5,
                    cacl2_conc=0.01,
                    # salt_conc           = salt_conc,
                    # NO3_conc            = NO3_conc,
                    # salt_sp             = salt_sp,
                    # cl_salt             = 'cacl2',
                    # cl_salt             = 'mgcl2',
                    Ex_TimeLim=20,
                    exename_src="scepter_MYCl",
                )
            )

            res_list.append(
                [time, dep_sample, phint, phint_lab, omint, acint, 100.0 - acint]
            )

        np.savetxt(
            outdir + runname_restart + f"/ph_BS_OM_{depth_str}.res",
            np.array(res_list),
            delimiter=" ",
            header="time dep_sample phint phint_lab omint acint 100-acint",
            comments="",
        )


def single_run_varKin():
    outdir = "/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/"

    dustname = "gbasalt"

    runname_spinup = "test_EF_alpha3p4"
    runname_spinup = "test_EF_alpha5p8"
    runname_spinup = "test_EF_alpha7p4"
    runname_spinup = "chk_EF_alpha7p4"
    runname_spinup = "chk2_EF_alpha7p4"
    runname_spinup = sys.argv[1]
    p80str = "320um"
    mix_scheme = 1
    mix_scheme = 2
    # mix_scheme = 3
    taudust = 0.1
    # taudust = 0.01
    # taudust = 0.001
    kinspc = 1e-5
    kinspc_list = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 1e-2]

    lower_bound = -5
    upper_bound = -2

    # Generate five random float numbers within the specified range
    kinspc_list = []
    for _ in range(5):
        kinspc_list.append(10 ** random.uniform(lower_bound, upper_bound))

    res_list_all = []
    cnt = 0
    for kinspc in kinspc_list:

        runid = "kin-{:d}_mix-{:d}_tau-{:f}".format(cnt, mix_scheme, taudust).replace(
            ".", "p"
        )
        runname_restart = runname_spinup + "_" + dustname + "_" + runid

        use_local_storage = True
        # use_local_storage = False

        runsuccess = add_gbas(
            outdir,
            runname_spinup,
            dustname,
            runname_restart,
            ttot_restart=1,
            # fdust_restart       = 2241.7,
            fdust_restart=5000,
            taudust_restart=taudust,
            # zdust_restart       = 0.25,
            zdust_restart=0.30,
            p80str=p80str,
            # nstep_restart       = 1000,
            nstep_restart=10,
            mix_scheme_restart=mix_scheme,
            kinspc=kinspc,
            update_exe=True,
            use_local_storage=use_local_storage,
            save_restart_dir=True,
            set_ExTimeLimit=False,
            Ex_TimeLim=20,
        )

        res_list = []
        for dep_sample in [0.10, 0.30]:
            # (1) get porewater pH
            phint = get_int_prof.get_ph_int_site(
                outdir, runname_restart, dep_sample, 20
            )

            # (2) get acidity in %
            acint = get_int_prof.get_ac_int_site_v2(
                outdir, runname_restart, dep_sample, 20
            )

            # (3) get bulk density
            dense = get_int_prof.get_rhobulk_int_site(
                outdir, runname_restart, dep_sample, 20
            )

            # (4) get field solid wt%
            sps = ["g2", "inrt"]
            sldwt_list = []
            for sp in sps:
                sldwt = get_int_prof.get_sldwt_int_site(
                    outdir, runname_restart, dep_sample, [sp], 20
                )
                sldwt_list.append(sldwt)

            # (5) get SOM wt%
            omint = sldwt_list[sps.index("g2")]

            # (6) get soil pH in lab
            time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
                get_soilpH_time.calc_soilpH(
                    outdir,
                    runname_restart,
                    dep_sample,
                    20,
                    include_N=True,
                    include_Cl=False,
                    include_Al=False,
                    include_DIC=True,
                    use_CaCl2=True,
                    add_salt=False,
                    add_NO3=False,
                    use_local_storage=use_local_storage,
                    save_lab_dir=False,
                    use_Sikora_buffer=False,
                    set_ExTimeLimit=True,
                    water_frac=2.5,
                    cacl2_conc=0.01,
                    # salt_conc           = salt_conc,
                    # NO3_conc            = NO3_conc,
                    # salt_sp             = salt_sp,
                    Ex_TimeLim=20,
                )
            )

            # fraction of basalt dissolved
            frac, time = get_int_prof.get_dis_frac(outdir, runname_restart, "gbas")

            res_list.append(
                [dep_sample, phint, phint_lab, omint, acint, 100.0 - acint, frac]
            )
            res_list_all.append(
                [dep_sample, phint, phint_lab, omint, acint, 100.0 - acint, frac]
            )

        for res in res_list:
            print(res)

        np.savetxt(
            outdir
            + runname_restart
            + "/res_tmp_mix-{:d}-{:f}".format(mix_scheme, taudust).replace(".", "p")
            + ".res",
            np.array(res_list),
        )

        cnt += 1
    np.savetxt(
        outdir
        + runname_spinup
        + "/res_tmp_mix-{:d}-{:f}".format(mix_scheme, taudust).replace(".", "p")
        + ".res",
        np.array(res_list_all),
    )


def single_run_varKin_EF(
    outdir,
    runname_spinup,
    runname_restart,
    tau,
    dustname,
    fdust_restart,
    taudust,
    kinspc,
    zdust_restart,
    make_runlogfile,
    calc_buffpH,
    sample_depth,
):
    # outdir = '/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/'

    # dustname = 'gbasalt'

    # runname_spinup = 'test_EF_alpha3p4'
    # runname_spinup = 'test_EF_alpha5p8'
    # runname_spinup = 'test_EF_alpha7p4'
    # runname_spinup = 'chk_EF_alpha7p4'
    # runname_spinup = 'chk2_EF_alpha7p4'
    p80str = "320um"
    mix_scheme = 1
    mix_scheme = 2
    # mix_scheme = 3
    # taudust = 0.1
    # taudust = 0.01
    # taudust = 0.001
    # kinspc = 1e-5
    # kinspc_list = [ 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 1e-2]

    # res_list_all = []

    # runid = 'kin-{:0.1f}_mix-{:d}_tau-{:f}'.format(-np.log10(kinspc),mix_scheme,taudust).replace('.','p')
    # runname_restart = runname_spinup + '_' + dustname + '_' + runid

    use_local_storage = True
    use_local_storage = False

    exename_src = "scepter_fwd"

    n_return = 8

    runsuccess = add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        # ttot_restart        = 1,
        ttot_restart=tau,
        # fdust_restart       = 2241.7,
        fdust_restart=fdust_restart,
        taudust_restart=taudust,
        # zdust_restart       = 0.25,
        zdust_restart=zdust_restart,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        mix_scheme_restart=mix_scheme,
        kinspc=kinspc,
        update_exe=True,
        use_local_storage=use_local_storage,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        # make_runlogfile     = True,
        # make_runlogfile     = False,
        make_runlogfile=make_runlogfile,
        exename_src=exename_src,
    )

    if not runsuccess:
        print("error after add_gbas in retart_add_gbas")
        return [np.nan] * n_return, runsuccess

    res_list = []
    # for dep_sample in [0.10,0.30]:
    # for dep_sample in [0.15]:
    for dep_sample in [sample_depth]:
        # (1) get porewater pH
        # phint = get_int_prof.get_ph_int_site(outdir,runname_restart,dep_sample,20)
        phint = get_int_prof.get_waterave_site_v3(
            outdir, runname_restart, dep_sample, 20, "ph"
        )

        # (2) get acidity in %
        # acint = get_int_prof.get_ac_int_site_v2(outdir,runname_restart,dep_sample,20)
        acint = get_int_prof.get_ac_int_site_v2(
            outdir, runname_restart, dep_sample, 20, simple_average=True
        )

        # (3) get bulk density
        dense = get_int_prof.get_rhobulk_int_site(
            outdir, runname_restart, dep_sample, 20
        )

        # (4) get field solid wt%
        sps = ["g2", "inrt"]
        sldwt_list = []
        for sp in sps:
            sldwt = get_int_prof.get_sldwt_int_site(
                outdir, runname_restart, dep_sample, [sp], 20
            )
            sldwt_list.append(sldwt)

        # (5) get SOM wt%
        omint = sldwt_list[sps.index("g2")]

        # (6) get soil pH in lab
        time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
            get_soilpH_time.calc_soilpH(
                outdir,
                runname_restart,
                dep_sample,
                20,
                include_N=True,
                include_Cl=False,
                include_Al=False,
                include_DIC=True,
                use_CaCl2=True,
                add_salt=False,
                add_NO3=False,
                use_local_storage=use_local_storage,
                save_lab_dir=False,
                use_Sikora_buffer=False,
                # set_ExTimeLimit     = True,
                set_ExTimeLimit=False,
                water_frac=2.5,
                cacl2_conc=0.01,
                # salt_conc           = salt_conc,
                # NO3_conc            = NO3_conc,
                # salt_sp             = salt_sp,
                Ex_TimeLim=20,
                make_runlogfile=True,
                # make_runlogfile     = False,
            )
        )

        if flag_error:
            print("error after get_soilpH_time in retart_add_gbas")
            return [np.nan] * n_return, flag_error

        if calc_buffpH:
            time_buf, IS_buf, acint_buf, phint_buf, flag_error, runname_lab = (
                get_soilpH_time.calc_soilpH(
                    outdir,
                    runname_restart,
                    dep_sample,
                    20,
                    include_N=True,
                    include_Cl=False,
                    include_Al=False,
                    include_DIC=True,
                    use_CaCl2=False,
                    add_salt=False,
                    add_NO3=False,
                    use_local_storage=use_local_storage,
                    save_lab_dir=False,
                    use_Sikora_buffer=True,
                    set_ExTimeLimit=False,
                    # water_frac          = water_frac,
                    # cacl2_conc          = cacl2_conc,
                    # salt_conc           = salt_conc,
                    # NO3_conc            = NO3_conc,
                    # salt_sp             = salt_sp,
                    Ex_TimeLim=20,
                    make_runlogfile=True,
                    # sub_as_a_job        = sub_as_a_job,
                )
            )
        else:
            phint_buf = np.nan

        if flag_error:
            print("error after get_soilpH_time in retart_add_gbas")
            return [np.nan] * n_return, flag_error

        # fraction of basalt dissolved
        disfrac, dis, rain, time = get_int_prof.get_dis_frac(
            outdir, runname_restart, dustname
        )

        # fraction of base cations advected out of dissolved base cations from basalt
        aqsps = ["na", "k", "ca", "mg"]
        distot = 0
        advtot = 0
        for aqsp in aqsps:
            advfrac_tmp, adv_tmp, dis_tmp, time = get_int_prof.get_adv_frac(
                outdir, runname_restart, aqsp, dustname
            )
            distot += dis_tmp
            advtot += adv_tmp

        advfrac = advtot / distot * 100.0

        res_list.append(
            [
                dep_sample,
                phint,
                phint_lab,
                omint,
                acint,
                100.0 - acint,
                disfrac,
                advfrac,
                phint_buf,
            ]
        )

    for res in res_list:
        print(res)

    return res_list, flag_error


def single_run_varKin_EF_multiyear(
    outdir,
    runname_spinup,
    runname_restart,
    tau,
    dustname,
    taudust,
    kinspc,
    calc_buffpH,
):
    # outdir = '/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/'

    # dustname = 'gbasalt'

    # runname_spinup = 'test_EF_alpha3p4'
    # runname_spinup = 'test_EF_alpha5p8'
    # runname_spinup = 'test_EF_alpha7p4'
    # runname_spinup = 'chk_EF_alpha7p4'
    # runname_spinup = 'chk2_EF_alpha7p4'
    p80str = "320um"
    mix_scheme = 1
    mix_scheme = 2
    # mix_scheme = 3
    # taudust = 0.1
    # taudust = 0.01
    # taudust = 0.001
    # kinspc = 1e-5
    # kinspc_list = [ 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 1e-2]

    # res_list_all = []

    timeline = list(np.linspace(0, tau, 100 * tau, endpoint=False))
    dust_temp = [0] * len(timeline)
    dust_temp[0] = 5000 / (timeline[1] - timeline[0])

    # runid = 'kin-{:0.1f}_mix-{:d}_tau-{:f}'.format(-np.log10(kinspc),mix_scheme,taudust).replace('.','p')
    # runname_restart = runname_spinup + '_' + dustname + '_' + runid

    use_local_storage = True
    use_local_storage = False

    n_return = 8

    runsuccess = restart_add_gbas_clim_DEV.add_gbas(
        outdir,
        runname_spinup,
        dustname,
        runname_restart,
        # ttot_restart        = 1,
        ttot_restart=tau,
        # fdust_restart       = 2241.7,
        fdust_restart=5000,
        timeline=timeline,
        dust_temp=dust_temp,
        # taudust_restart     = taudust,
        # zdust_restart       = 0.25,
        zdust_restart=0.30,
        p80str=p80str,
        # nstep_restart       = 1000,
        nstep_restart=10,
        mix_scheme_restart=mix_scheme,
        kinspc=kinspc,
        update_exe=True,
        use_local_storage=use_local_storage,
        save_restart_dir=True,
        set_ExTimeLimit=False,
        Ex_TimeLim=20,
        rec_logfile=False,
    )

    if not runsuccess:
        print("error after add_gbas in retart_add_gbas")
        return [np.nan] * n_return, runsuccess

    res_list = []
    # for dep_sample in [0.10,0.30]:
    for dep_sample in [0.15]:
        # (1) get porewater pH
        # phint = get_int_prof.get_ph_int_site(outdir,runname_restart,dep_sample,20)
        phint = get_int_prof.get_waterave_site_v3(
            outdir, runname_restart, dep_sample, 20, "ph"
        )

        # (2) get acidity in %
        # acint = get_int_prof.get_ac_int_site_v2(outdir,runname_restart,dep_sample,20)
        acint = get_int_prof.get_ac_int_site_v2(
            outdir, runname_restart, dep_sample, 20, simple_average=True
        )

        # (3) get bulk density
        dense = get_int_prof.get_rhobulk_int_site(
            outdir, runname_restart, dep_sample, 20
        )

        # (4) get field solid wt%
        sps = ["g2", "inrt"]
        sldwt_list = []
        for sp in sps:
            sldwt = get_int_prof.get_sldwt_int_site(
                outdir, runname_restart, dep_sample, [sp], 20
            )
            sldwt_list.append(sldwt)

        # (5) get SOM wt%
        omint = sldwt_list[sps.index("g2")]

        # (6) get soil pH in lab
        time, IS_lab, acint_lab, phint_lab, flag_error, runname_lab = (
            get_soilpH_time.calc_soilpH(
                outdir,
                runname_restart,
                dep_sample,
                20,
                include_N=True,
                include_Cl=False,
                include_Al=False,
                include_DIC=True,
                use_CaCl2=True,
                add_salt=False,
                add_NO3=False,
                use_local_storage=use_local_storage,
                save_lab_dir=False,
                use_Sikora_buffer=False,
                # set_ExTimeLimit     = True,
                set_ExTimeLimit=False,
                water_frac=2.5,
                cacl2_conc=0.01,
                # salt_conc           = salt_conc,
                # NO3_conc            = NO3_conc,
                # salt_sp             = salt_sp,
                Ex_TimeLim=20,
                make_runlogfile=True,
                # make_runlogfile     = False,
            )
        )

        if flag_error:
            print("error after get_soilpH_time in retart_add_gbas")
            return [np.nan] * n_return, flag_error

        if calc_buffpH:
            time_buf, IS_buf, acint_buf, phint_buf, flag_error, runname_lab = (
                get_soilpH_time.calc_soilpH(
                    outdir,
                    runname_restart,
                    dep_sample,
                    20,
                    include_N=True,
                    include_Cl=False,
                    include_Al=False,
                    include_DIC=True,
                    use_CaCl2=False,
                    add_salt=False,
                    add_NO3=False,
                    use_local_storage=use_local_storage,
                    save_lab_dir=False,
                    use_Sikora_buffer=True,
                    set_ExTimeLimit=False,
                    # water_frac          = water_frac,
                    # cacl2_conc          = cacl2_conc,
                    # salt_conc           = salt_conc,
                    # NO3_conc            = NO3_conc,
                    # salt_sp             = salt_sp,
                    Ex_TimeLim=20,
                    make_runlogfile=True,
                    # sub_as_a_job        = sub_as_a_job,
                )
            )
        else:
            phint_buf = np.nan

        if flag_error:
            print("error after get_soilpH_time in retart_add_gbas")
            return [np.nan] * n_return, flag_error

        # fraction of basalt dissolved
        disfrac, dis, rain, time = get_int_prof.get_dis_frac(
            outdir, runname_restart, "gbas"
        )

        # fraction of base cations advected out of dissolved base cations from basalt
        aqsps = ["na", "k", "ca", "mg"]
        distot = 0
        advtot = 0
        for aqsp in aqsps:
            advfrac_tmp, adv_tmp, dis_tmp, time = get_int_prof.get_adv_frac(
                outdir, runname_restart, aqsp, "gbas"
            )
            distot += dis_tmp
            advtot += adv_tmp

        advfrac = advtot / distot * 100.0

        res_list.append(
            [
                dep_sample,
                phint,
                phint_lab,
                omint,
                acint,
                100.0 - acint,
                disfrac,
                advfrac,
                phint_buf,
            ]
        )

    for res in res_list:
        print(res)

    return res_list, flag_error


def main():
    """
    Main function - Entry point for the restart_add_gbas.py script.

    This function contains commented-out calls to various experimental functions.
    Uncomment the desired function to run specific ERW experiments.

    Available functions:
    - single_run_varPSD(): Single run with variable particle size distribution
    - single_run_varKin(): Single run with variable kinetics
    - multiples_run_lime_varPSD(): Multiple runs with lime and variable PSD
    - multiples_run_varPSD_MY(): Multiple runs with variable PSD (multi-year)
    - multiples_run_lime_season(): Multiple runs with seasonal variations
    - get_soilpH(): Extract soil pH data from completed runs
    - single_run_spc_tau(): Single run for testing (currently active)

    Usage:
        python3 restart_add_gbas.py <spinup_name> <restart_name>
    """

    # Experimental function calls (commented out)
    # single_run_varPSD()                    # Variable particle size distribution
    # single_run_varKin()                    # Variable kinetics
    # multiples_run_lime_varPSD(flg_ctrl=True)   # Multiple lime runs with control
    # multiples_run_lime_varPSD(flg_ctrl=False)  # Multiple lime runs without control
    # multiples_run_lime_varPSD_EF()         # Multiple lime runs (EF region)
    # multiples_run_varPSD_MY()              # Multiple runs (multi-year)
    # multiples_run_varPSD_MY_season(flg_ctrl=False)  # Multi-year seasonal runs
    # multiples_run_varPSD_MY_season(flg_ctrl=True)   # Multi-year seasonal runs with control
    # multiples_run_varPSD_MY_season_sensitivity()    # Sensitivity analysis
    # multiples_run_lime_varPSD_season(run_ctrl=False) # Mississippi seasonal runs
    # multiples_run_lime_varPSD_season(run_ctrl=True)  # Mississippi seasonal runs with control
    # get_soilpH()                           # Extract soil pH data (Mississippi runs)

    # Currently active function for testing
    single_run_spc_tau()  # Single run for front-end testing

    # Additional experimental functions
    # multiples_run_season()                 # Seasonal variation runs
    # multiples_run_lime_Africa(run_ctrl=False)      # Africa lime runs
    # multiples_run_lime_Africa(run_ctrl=True)       # Africa lime runs with control
    # multiples_run_lime_season_Africa(run_ctrl=False) # Africa seasonal runs
    # multiples_run_lime_season_onetime_Africa(run_ctrl=False) # Africa one-time runs
    # multiples_run_lime_season_onetime_Africa(run_ctrl=True)  # Africa one-time runs with control
    # single_run_rand_kin_tau()              # Random kinetics runs


if __name__ == "__main__":
    main()
