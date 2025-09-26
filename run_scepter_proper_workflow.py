#!/usr/bin/env python3
"""
Proper SCEPTER Workflow: Spinup First, Then ERW
Based on developer guidance using tunespin_3_newton_inert_buff_v2_clean.py and basalt_buff_tunespin_bisec_v2.py

This script implements the correct two-step process:
1. Run spinup simulation to establish initial soil conditions
2. Use spinup results as starting point for ERW simulation
"""

import pandas as pd
import numpy as np
import os
import shutil
import subprocess
import sys
import json
import time
from datetime import datetime
from spinup import run_a_scepter_run
import math


def load_csv_data(csv_file="./data/inputdata_depres.csv"):
    """Load and process the CSV data"""
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} data points from {csv_file}")
    return df


def calculate_distance(lat1, lon1, lat2, lon2):
    """Calculate distance between two points using Haversine formula"""
    R = 6371  # Earth's radius in km

    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)
    delta_lat = math.radians(lat2 - lat1)
    delta_lon = math.radians(lon2 - lon1)

    a = (
        math.sin(delta_lat / 2) ** 2
        + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2) ** 2
    )
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return R * c


def find_nearest_csv_site(target_lat, target_lon, df=None, max_distance_km=None):
    """Find the nearest CSV data point to target coordinates"""
    if df is None:
        df = load_csv_data()

    if max_distance_km is None:
        max_distance_km = 100

    distances = []
    for idx, row in df.iterrows():
        dist = calculate_distance(target_lat, target_lon, row["lat"], row["lon"])
        distances.append(dist)

    df["distance"] = distances
    nearest_idx = df["distance"].idxmin()
    nearest_row = df.loc[nearest_idx]

    if nearest_row["distance"] > max_distance_km:
        return {
            "found": False,
            "message": f"No CSV data within {max_distance_km} km of {target_lat:.3f}°N, {target_lon:.3f}°W",
        }

    return {
        "found": True,
        "row": nearest_row,
        "distance_km": nearest_row["distance"],
        "csv_lat": nearest_row["lat"],
        "csv_lon": nearest_row["lon"],
        "target_lat": target_lat,
        "target_lon": target_lon,
    }


def csv_row_to_spinup_params(row):
    """Convert CSV row to SCEPTER spinup parameters"""
    # Extract key values from CSV row
    lat = row["lat"]
    lon = row["lon"]
    temp = row["temperature [oC]"]
    soil_moisture = row["soil moisture [m3/m3]"]
    ph_0_30cm = row["pH_H2O_0-30cm [-]"]
    soc_0_30cm = row["SOC_0-30cm [wt%]"]
    cec_0_30cm = row["CEC_0-30cm [cmol/kg]"]
    npp = row["NPP [gC/m/yr]"]
    porosity = row["porosity [-]"]
    runoff = row["runoff [m/yr]"]
    erosion_mm_yr = row["erosion [mm/yr]"]
    cropland_pct = row["cropland [%]"]

    # Additional CSV values
    bs_0_20cm = row["BS_0-20cm [%]"]
    nitrification_rate = row["nitrification rate [kgN/ha/day]"]

    # Calculate initial conditions
    h_initial = 10 ** (-ph_0_30cm)
    soc_initial = soc_0_30cm / 100.0
    base_cations_total = (bs_0_20cm / 100.0) * cec_0_30cm
    ca_initial = base_cations_total * 0.6
    mg_initial = base_cations_total * 0.25
    k_initial = base_cations_total * 0.1
    na_initial = base_cations_total * 0.05

    # Create runname from coordinates
    runname = f"site_{lat:.1f}_{lon:.1f}_spinup".replace("-", "neg").replace(".", "p")

    # Map CSV data to SCEPTER parameters for SPINUP
    # Using same defaults as spinup.py, only overriding with CSV data
    params = {
        # Basic simulation setup
        "runname": runname,
        "lat": lat,
        "lon": lon,
        # ---- frame.in parameters (from spinup.py defaults) ----
        "ztot": 0.5,
        "nz": 30,
        "ttot": 1e5,
        "temp": temp,  # From CSV
        "fdust": 0,  # NO dust in spinup
        "fdust2": 0,
        "taudust": 0,
        "omrain": npp * 1.5,  # Use NPP with developer's conversion factor
        "zom": 0.25,
        "poro": porosity,  # Use CSV porosity
        "moistsrf": soil_moisture,  # Use CSV soil moisture
        "zwater": 1000,
        "zdust": 0.25,
        "w": erosion_mm_yr / 1000.0,  # Use CSV erosion (convert mm/yr to m/yr)
        "q": runoff,  # Use CSV runoff
        "p": 1e-5,
        "nstep": 10,  # Default from spinup.py
        "rstrt": "self",
        "runid": runname,
        # ---- switches.in parameters (from spinup.py defaults) ----
        "w_scheme": 0,  # Default from spinup.py
        "mix_scheme": 0,  # Default from spinup.py
        "poro_iter": "true",  # Default from spinup.py
        "sldmin_lim": "true",
        "display": "true",
        "disp_lim": "true",
        "restart": "false",
        "rough": "true",
        "act_ON": "true",
        "dt_fix": "false",
        "cec_on": "true",
        "dz_fix": "true",
        "close_aq": "false",
        "poro_evol": "true",
        "sa_evol_1": "true",
        "sa_evol_2": "false",
        "psd_bulk": "true",
        "psd_full": "true",
        "season": "false",
        # ---- tracers (from spinup.py defaults) ----
        "sld_list": ["inrt", "g2"],
        "aq_list": ["ca", "k", "mg", "na"],
        "gas_list": ["pco2"],
        "exrxn_list": [],
        # ---- boundary values (from spinup.py defaults) ----
        "pr_list": [("inrt", 1.0)],
        "rain_list": [("ca", 5.0e-6)],  # Default from spinup.py
        "atm_list": [
            ("pco2", 3.16e-4),
            ("po2", 0.21),
            ("pnh3", 1e-50),
            ("pn2o", 1e-50),
        ],
        # ---- solid phase properties (from spinup.py defaults) ----
        "sld_varlist_dust": [],
        "sld_varlist_cec": [
            ("inrt", 10, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4)
        ],  # Default from spinup.py
        "sld_varlist_omrain": [("g2", 1.0)],
        "sld_varlist_kinspc": [],
        "sld_varlist_2ndslds": [],
        "srcfile_dust": None,
        "srcfile_omrain": None,
        "srcfile_cec": None,
        "srcfile_kinspc": None,
        "srcfile_2ndslds": "./data/2ndslds_def.in",
        # ---- python stuff (from spinup.py defaults) ----
        "use_local_storage": True,
        "lim_calc_time": False,
        "max_calc_time": 20,
        # Store original CSV data for reference
        "csv_data": {
            "ph_0_30cm": ph_0_30cm,
            "soc_0_30cm": soc_0_30cm,
            "cropland_pct": cropland_pct,
            "bs_0_20cm": bs_0_20cm,
            "nitrification_rate": nitrification_rate,
            "h_initial": h_initial,
            "soc_initial": soc_initial,
            "base_cations_total": base_cations_total,
        },
    }

    return params


def csv_row_to_erw_params(row, erw_rate, spinup_runname):
    """Convert CSV row to SCEPTER ERW parameters using spinup results"""
    # Extract key values from CSV row
    lat = row["lat"]
    lon = row["lon"]
    temp = row["temperature [oC]"]
    soil_moisture = row["soil moisture [m3/m3]"]
    ph_0_30cm = row["pH_H2O_0-30cm [-]"]
    soc_0_30cm = row["SOC_0-30cm [wt%]"]
    cec_0_30cm = row["CEC_0-30cm [cmol/kg]"]
    npp = row["NPP [gC/m/yr]"]
    porosity = row["porosity [-]"]
    runoff = row["runoff [m/yr]"]
    erosion_mm_yr = row["erosion [mm/yr]"]
    cropland_pct = row["cropland [%]"]

    # Additional CSV values
    bs_0_20cm = row["BS_0-20cm [%]"]
    nitrification_rate = row["nitrification rate [kgN/ha/day]"]

    # Calculate initial conditions
    h_initial = 10 ** (-ph_0_30cm)
    soc_initial = soc_0_30cm / 100.0
    base_cations_total = (bs_0_20cm / 100.0) * cec_0_30cm
    ca_initial = base_cations_total * 0.6
    mg_initial = base_cations_total * 0.25
    k_initial = base_cations_total * 0.1
    na_initial = base_cations_total * 0.05

    # Create runname from coordinates
    runname = f"site_{lat:.1f}_{lon:.1f}_erw_{erw_rate}g".replace("-", "neg").replace(
        ".", "p"
    )

    # Determine if site is agricultural or natural
    is_agricultural = cropland_pct >= 50

    # ERW parameters - with safeguards for PBE stability
    safe_erw_rate = min(erw_rate, 3000.0)
    if erw_rate > 3000.0:
        print(
            f"⚠️  Warning: ERW rate {erw_rate} g/m²/yr capped to {safe_erw_rate} g/m²/yr for PBE stability"
        )

    # Map CSV data to SCEPTER parameters for ERW
    # Using same defaults as basalt_buff_tunespin_bisec_v2.py, only overriding with CSV data
    params = {
        # Basic simulation setup
        "runname": runname,
        "lat": lat,
        "lon": lon,
        "spinup_runname": spinup_runname,  # Reference to spinup results
        # ---- frame.in parameters (from basalt_buff_tunespin_bisec_v2.py) ----
        "ztot": 0.5,
        "nz": 30,
        "ttot": 1e5,
        "temp": temp,  # From CSV
        "fdust": safe_erw_rate,  # ERW application rate
        "fdust2": 0,  # Default from basalt script
        "taudust": 0.005,  # Default from basalt script (very short duration)
        "omrain": npp * 1.5,  # Use NPP with developer's conversion factor
        "zom": 0.25,
        "poro": porosity,  # Use CSV porosity
        "moistsrf": soil_moisture,  # Use CSV soil moisture
        "zwater": 1000,
        "zdust": 0.25,
        "w": erosion_mm_yr / 1000.0,  # Use CSV erosion (convert mm/yr to m/yr)
        "q": runoff,  # Use CSV runoff
        "p": 1e-5,
        "nstep": 10,  # Default from spinup.py
        "rstrt": spinup_runname,  # Restart from spinup results
        "runid": runname,
        # ---- switches.in parameters (from basalt_buff_tunespin_bisec_v2.py) ----
        "w_scheme": 0,  # Default from spinup.py
        "mix_scheme": 2,  # Default from basalt script (homogeneous mixing)
        "poro_iter": "true",  # Default from spinup.py
        "sldmin_lim": "true",
        "display": "true",
        "disp_lim": "true",
        "restart": "true",  # Restart from spinup
        "rough": "true",
        "act_ON": "true",
        "dt_fix": "false",
        "cec_on": "true",
        "dz_fix": "true",
        "close_aq": "false",
        "poro_evol": "true",
        "sa_evol_1": "true",
        "sa_evol_2": "false",
        "psd_bulk": "true",
        "psd_full": "true",
        "season": "false",
        # ---- tracers (from basalt_buff_tunespin_bisec_v2.py) ----
        "sld_list": ["inrt", "g2", "gbas"],  # Include glass basalt (gbas)
        "aq_list": ["ca", "k", "mg", "na"],  # Basic aqueous species
        "gas_list": ["pco2"],
        "exrxn_list": [],
        # ---- boundary values (from spinup.py defaults) ----
        "pr_list": [("inrt", 1.0)],
        "rain_list": [("ca", 5.0e-6)],  # Default from spinup.py
        "atm_list": [
            ("pco2", 3.16e-4),
            ("po2", 0.21),
            ("pnh3", 1e-50),
            ("pn2o", 1e-50),
        ],
        # ---- solid phase properties (from spinup.py defaults) ----
        "sld_varlist_dust": [],
        "sld_varlist_cec": [
            ("inrt", 10, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4)
        ],  # Default from spinup.py
        "sld_varlist_omrain": [("g2", 1.0)],
        "sld_varlist_kinspc": [],
        "sld_varlist_2ndslds": [],
        "srcfile_dust": "./data/dust_gbasalt.in",  # Use glass basalt composition
        "srcfile_omrain": None,
        "srcfile_cec": None,
        "srcfile_kinspc": None,
        "srcfile_2ndslds": "./data/2ndslds_def.in",
        # ---- python stuff (from spinup.py defaults) ----
        "use_local_storage": True,
        "lim_calc_time": False,
        "max_calc_time": 20,
        # Store original CSV data for reference
        "csv_data": {
            "ph_0_30cm": ph_0_30cm,
            "soc_0_30cm": soc_0_30cm,
            "cropland_pct": cropland_pct,
            "bs_0_20cm": bs_0_20cm,
            "nitrification_rate": nitrification_rate,
            "h_initial": h_initial,
            "soc_initial": soc_initial,
            "base_cations_total": base_cations_total,
            "erw_rate": safe_erw_rate,
        },
    }

    return params


def run_spinup_simulation(
    target_lat, target_lon, max_distance_km=100, outdir_src="../scepter_output/"
):
    """Run spinup simulation to establish initial soil conditions"""
    print(f"🔄 Running SPINUP simulation for {target_lat:.3f}°N, {target_lon:.3f}°W")

    # Load CSV data
    df = load_csv_data()

    # Find nearest site
    nearest = find_nearest_csv_site(target_lat, target_lon, df, max_distance_km)

    if not nearest["found"]:
        print(f"❌ {nearest['message']}")
        return {"success": False, "error": nearest["message"]}

    print(f"✓ Found nearest CSV data point:")
    print(f"  CSV coordinates: {nearest['csv_lat']:.3f}°N, {nearest['csv_lon']:.3f}°W")
    print(f"  Distance: {nearest['distance_km']:.1f} km")

    # Convert CSV row to spinup parameters
    params = csv_row_to_spinup_params(nearest["row"])

    print(f"\nRunning SPINUP simulation: {params['runname']}")
    print(f"Target location: {target_lat:.3f}°N, {target_lon:.3f}°W")
    print(
        f"Using CSV data from: {nearest['csv_lat']:.3f}°N, {nearest['csv_lon']:.3f}°W"
    )
    print(f"Simulation time: {params['ttot']:.0e} years")
    print(f"No dust application (spinup only)")

    # Run spinup simulation
    result = run_a_scepter_run(
        runname=params["runname"],
        outdir_src=outdir_src,
        ztot=params["ztot"],
        nz=params["nz"],
        ttot=params["ttot"],
        temp=params["temp"],
        fdust=params["fdust"],
        fdust2=params["fdust2"],
        taudust=params["taudust"],
        omrain=params["omrain"],
        zom=params["zom"],
        poro=params["poro"],
        moistsrf=params["moistsrf"],
        zwater=params["zwater"],
        zdust=params["zdust"],
        w=params["w"],
        q=params["q"],
        p=params["p"],
        nstep=params["nstep"],
        rstrt=params["rstrt"],
        w_scheme=params["w_scheme"],
        mix_scheme=params["mix_scheme"],
        poro_iter=params["poro_iter"],
        sldmin_lim=params["sldmin_lim"],
        display=params["display"],
        disp_lim=params["disp_lim"],
        restart=params["restart"],
        rough=params["rough"],
        act_ON=params["act_ON"],
        dt_fix=params["dt_fix"],
        cec_on=params["cec_on"],
        dz_fix=params["dz_fix"],
        close_aq=params["close_aq"],
        poro_evol=params["poro_evol"],
        sa_evol_1=params["sa_evol_1"],
        sa_evol_2=params["sa_evol_2"],
        psd_bulk=params["psd_bulk"],
        psd_full=params["psd_full"],
        season=params["season"],
        sld_list=params["sld_list"],
        aq_list=params["aq_list"],
        gas_list=params["gas_list"],
        exrxn_list=params["exrxn_list"],
        pr_list=params["pr_list"],
        rain_list=params["rain_list"],
        atm_list=params["atm_list"],
        sld_varlist_dust=params["sld_varlist_dust"],
        sld_varlist_cec=params["sld_varlist_cec"],
        sld_varlist_omrain=params["sld_varlist_omrain"],
        sld_varlist_kinspc=params["sld_varlist_kinspc"],
        sld_varlist_2ndslds=params["sld_varlist_2ndslds"],
        srcfile_dust=params["srcfile_dust"],
        srcfile_omrain=params["srcfile_omrain"],
        srcfile_cec=params["srcfile_cec"],
        srcfile_kinspc=params["srcfile_kinspc"],
        srcfile_2ndslds=params["srcfile_2ndslds"],
        use_local_storage=params["use_local_storage"],
        lim_calc_time=params["lim_calc_time"],
        max_calc_time=params["max_calc_time"],
        runid=params["runid"],
    )

    if result["success"]:
        print(f"SPINUP completed successfully: {params['runname']}")
    else:
        print(f"SPINUP failed: {result.get('error', 'Unknown error')}")

    return result


def run_erw_simulation(
    target_lat,
    target_lon,
    erw_rate,
    spinup_runname,
    max_distance_km=100,
    outdir_src="../scepter_output/",
):
    """Run ERW simulation using spinup results as starting point"""
    print(f"🌋 Running ERW simulation for {target_lat:.3f}°N, {target_lon:.3f}°W")
    print(f"   ERW rate: {erw_rate} g/m²/yr")
    print(f"   Starting from spinup: {spinup_runname}")

    # Load CSV data
    df = load_csv_data()

    # Find nearest site
    nearest = find_nearest_csv_site(target_lat, target_lon, df, max_distance_km)

    if not nearest["found"]:
        print(f"❌ {nearest['message']}")
        return {"success": False, "error": nearest["message"]}

    # Convert CSV row to ERW parameters
    params = csv_row_to_erw_params(nearest["row"], erw_rate, spinup_runname)

    print(f"\nRunning ERW simulation: {params['runname']}")
    print(f"Target location: {target_lat:.3f}°N, {target_lon:.3f}°W")
    print(
        f"Using CSV data from: {nearest['csv_lat']:.3f}°N, {nearest['csv_lon']:.3f}°W"
    )
    print(f"Simulation time: {params['ttot']:.0e} years")
    print(f"ERW application rate: {params['fdust']} g/m²/yr")
    print(f"Application duration: {params['taudust']} years")
    print(f"Restarting from: {params['rstrt']}")

    # Run ERW simulation
    result = run_a_scepter_run(
        runname=params["runname"],
        outdir_src=outdir_src,
        ztot=params["ztot"],
        nz=params["nz"],
        ttot=params["ttot"],
        temp=params["temp"],
        fdust=params["fdust"],
        fdust2=params["fdust2"],
        taudust=params["taudust"],
        omrain=params["omrain"],
        zom=params["zom"],
        poro=params["poro"],
        moistsrf=params["moistsrf"],
        zwater=params["zwater"],
        zdust=params["zdust"],
        w=params["w"],
        q=params["q"],
        p=params["p"],
        nstep=params["nstep"],
        rstrt=params["rstrt"],
        w_scheme=params["w_scheme"],
        mix_scheme=params["mix_scheme"],
        poro_iter=params["poro_iter"],
        sldmin_lim=params["sldmin_lim"],
        display=params["display"],
        disp_lim=params["disp_lim"],
        restart=params["restart"],
        rough=params["rough"],
        act_ON=params["act_ON"],
        dt_fix=params["dt_fix"],
        cec_on=params["cec_on"],
        dz_fix=params["dz_fix"],
        close_aq=params["close_aq"],
        poro_evol=params["poro_evol"],
        sa_evol_1=params["sa_evol_1"],
        sa_evol_2=params["sa_evol_2"],
        psd_bulk=params["psd_bulk"],
        psd_full=params["psd_full"],
        season=params["season"],
        sld_list=params["sld_list"],
        aq_list=params["aq_list"],
        gas_list=params["gas_list"],
        exrxn_list=params["exrxn_list"],
        pr_list=params["pr_list"],
        rain_list=params["rain_list"],
        atm_list=params["atm_list"],
        sld_varlist_dust=params["sld_varlist_dust"],
        sld_varlist_cec=params["sld_varlist_cec"],
        sld_varlist_omrain=params["sld_varlist_omrain"],
        sld_varlist_kinspc=params["sld_varlist_kinspc"],
        sld_varlist_2ndslds=params["sld_varlist_2ndslds"],
        srcfile_dust=params["srcfile_dust"],
        srcfile_omrain=params["srcfile_omrain"],
        srcfile_cec=params["srcfile_cec"],
        srcfile_kinspc=params["srcfile_kinspc"],
        srcfile_2ndslds=params["srcfile_2ndslds"],
        use_local_storage=params["use_local_storage"],
        lim_calc_time=params["lim_calc_time"],
        max_calc_time=params["max_calc_time"],
        runid=params["runid"],
    )

    if result["success"]:
        print(f"ERW completed successfully: {params['runname']}")
    else:
        print(f"ERW failed: {result.get('error', 'Unknown error')}")

    return result


def run_proper_workflow(
    target_lat,
    target_lon,
    erw_rate,
    max_distance_km=100,
    outdir_src="../scepter_output/",
):
    """Run the proper SCEPTER workflow: spinup first, then ERW"""
    print("=" * 80)
    print("🌱 SCEPTER PROPER WORKFLOW: SPINUP → ERW")
    print("=" * 80)

    # Step 1: Run spinup simulation
    print("\n📋 STEP 1: SPINUP SIMULATION")
    print("-" * 40)
    spinup_result = run_spinup_simulation(
        target_lat, target_lon, max_distance_km, outdir_src
    )

    if not spinup_result["success"]:
        print(f"❌ Spinup failed, cannot proceed with ERW")
        return {
            "success": False,
            "error": "Spinup failed",
            "spinup_result": spinup_result,
        }

    # Extract spinup runname
    spinup_runname = spinup_result.get("runname", "unknown")

    # Step 2: Run ERW simulation
    print("\n📋 STEP 2: ERW SIMULATION")
    print("-" * 40)
    erw_result = run_erw_simulation(
        target_lat, target_lon, erw_rate, spinup_runname, max_distance_km, outdir_src
    )

    # Return combined results
    return {
        "success": erw_result["success"],
        "spinup_result": spinup_result,
        "erw_result": erw_result,
        "spinup_runname": spinup_runname,
        "erw_runname": erw_result.get("runname", "unknown"),
    }


def load_sites_config(json_file):
    """Load sites configuration from JSON file"""
    with open(json_file, "r") as f:
        config = json.load(f)
    return config


def run_single_site_proper_workflow(lat, lon, erw_rate, site_name, max_distance_km=100):
    """Run a single site through the proper workflow"""
    print(f"\n{'='*80}")
    print(f"Starting proper workflow for: {site_name}")
    print(f"Location: {lat:.3f}°N, {lon:.3f}°W")
    print(f"ERW Rate: {erw_rate} g/m²/yr")
    print(f"{'='*80}")

    # Build command
    cmd = [
        "python3",
        "run_scepter_proper_workflow.py",
        str(lat),
        str(lon),
        str(erw_rate),
        str(max_distance_km),
    ]

    print(f"Command: {' '.join(cmd)}")

    # Run the command
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )  # 1 hour timeout
        end_time = time.time()
        duration = end_time - start_time

        if result.returncode == 0:
            print(f"SUCCESS: {site_name} completed in {duration:.1f} seconds")
            print("STDOUT:")
            print(result.stdout)
            return True, result.stdout
        else:
            print(f"FAILED: {site_name} failed after {duration:.1f} seconds")
            print("STDERR:")
            print(result.stderr)
            print("STDOUT:")
            print(result.stdout)
            return False, result.stderr

    except subprocess.TimeoutExpired:
        print(f"TIMEOUT: {site_name} exceeded 1 hour limit")
        return False, "Timeout exceeded"
    except Exception as e:
        print(f"ERROR: {site_name} - {str(e)}")
        return False, str(e)


def run_batch_sites(json_file):
    """Run multiple sites from JSON configuration"""
    if not os.path.exists(json_file):
        print(f"Error: JSON file '{json_file}' not found")
        sys.exit(1)

    # Load configuration
    config = load_sites_config(json_file)
    sites = config["sites"]
    settings = config.get("settings", {})

    max_distance_km = settings.get("max_distance_km", 100)

    print(f"USGS Sites Proper Workflow")
    print(f"Config file: {json_file}")
    print(f"Max distance: {max_distance_km} km")
    print(f"Total sites: {len(sites)}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Track results
    results = []
    successful = 0
    failed = 0

    # Run each site
    for i, site in enumerate(sites, 1):
        print(f"\n{'#'*80}")
        print(f"Site {i}/{len(sites)}: {site['name']}")
        print(f"{'#'*80}")

        lat = site["lat"]
        lon = site["lon"]
        erw_rate = site["erw_rate"]
        site_name = site["name"]

        success, output = run_single_site_proper_workflow(
            lat, lon, erw_rate, site_name, max_distance_km
        )

        results.append(
            {
                "site": site_name,
                "lat": lat,
                "lon": lon,
                "erw_rate": erw_rate,
                "success": success,
                "output": output,
            }
        )

        if success:
            successful += 1
        else:
            failed += 1

        # Brief pause between sites
        if i < len(sites):
            print(f"\nWaiting 5 seconds before next site...")
            time.sleep(5)

    # Final summary
    print(f"\n{'='*80}")
    print(f"Sites Proper Workflow - COMPLETED")
    print(f"{'='*80}")
    print(f"Successful: {successful}/{len(sites)}")
    print(f"Failed: {failed}/{len(sites)}")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Detailed results
    print(f"\nDETAILED RESULTS:")
    print(f"{'='*80}")
    for result in results:
        status = "SUCCESS" if result["success"] else "FAILED"
        print(
            f"{status}: {result['site']} ({result['lat']:.3f}°N, {result['lon']:.3f}°W)"
        )
        if not result["success"]:
            print(f"    Error: {result['output'][:200]}...")

    # Save results to file
    results_file = (
        f"sites_proper_workflow_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    )
    with open(results_file, "w") as f:
        json.dump(
            {
                "summary": {
                    "total_sites": len(sites),
                    "successful": successful,
                    "failed": failed,
                    "start_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "config_file": json_file,
                },
                "results": results,
            },
            f,
            indent=2,
        )

    print(f"\nResults saved to: {results_file}")

    # Exit with error code if any failed
    if failed > 0:
        sys.exit(1)
    else:
        print(f"\nAll sites completed successfully!")
        sys.exit(0)


def main():
    """Main function for proper SCEPTER workflow"""
    if len(sys.argv) < 2:
        print("Usage:")
        print("  # Run single site")
        print(
            "  python3 run_scepter_proper_workflow.py <lat> <lon> <erw_rate> [max_distance_km]"
        )
        print("  # Run batch sites from JSON")
        print("  python3 run_scepter_proper_workflow.py --batch <json_file>")
        print("")
        print("Examples:")
        print("  # Run proper workflow for a site")
        print("  python3 run_scepter_proper_workflow.py 39.34 -82.97 2000")
        print("  # Run with custom distance limit")
        print("  python3 run_scepter_proper_workflow.py 39.34 -82.97 2000 50")
        print("  # Run batch of sites")
        print("  python3 run_scepter_proper_workflow.py --batch usgs_12_sites_erw.json")
        sys.exit(1)

    # Check if batch mode
    if sys.argv[1] == "--batch":
        if len(sys.argv) < 3:
            print("Error: JSON file required for batch mode")
            print("Usage: python3 run_scepter_proper_workflow.py --batch <json_file>")
            sys.exit(1)
        json_file = sys.argv[2]
        run_batch_sites(json_file)
        return

    # Single site mode
    if len(sys.argv) < 4:
        print("Error: Insufficient arguments for single site mode")
        print(
            "Usage: python3 run_scepter_proper_workflow.py <lat> <lon> <erw_rate> [max_distance_km]"
        )
        sys.exit(1)

    target_lat = float(sys.argv[1])
    target_lon = float(sys.argv[2])
    erw_rate = float(sys.argv[3])
    max_distance_km = float(sys.argv[4]) if len(sys.argv) > 4 else 100

    # Validate ERW rate
    if erw_rate <= 0:
        print("Error: ERW rate must be > 0")
        print("Typical rates: 500-3000 g/m²/yr")
        sys.exit(1)

    # Run proper workflow
    result = run_proper_workflow(target_lat, target_lon, erw_rate, max_distance_km)

    if result["success"]:
        print("\n" + "=" * 80)
        print("PROPER WORKFLOW COMPLETED SUCCESSFULLY!")
        print("=" * 80)
        print(f"Spinup simulation: {result['spinup_runname']}")
        print(f"ERW simulation: {result['erw_runname']}")
        print(f"Location: {target_lat:.3f}°N, {target_lon:.3f}°W")
        print(f"ERW rate: {erw_rate} g/m²/yr")
    else:
        print("\n" + "=" * 80)
        print("PROPER WORKFLOW FAILED")
        print("=" * 80)
        print(f"Error: {result.get('error', 'Unknown error')}")


if __name__ == "__main__":
    main()
