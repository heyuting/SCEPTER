#!/usr/bin/env python3
"""
Run SCEPTER from Coordinates
Main entry point to run SCEPTER simulations anywhere in the US
Uses real soil and climate data from ./data/inputdata_depres.csv (822 sites)
Based on the structure from spinup.py
"""

import pandas as pd
import numpy as np
from spinup import run_a_scepter_run
import math
import json


def load_csv_data(csv_file="./data/inputdata_depres.csv"):
    """
    Load and process the CSV data
    """
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} data points from {csv_file}")
    print(f"Latitude range: {df['lat'].min():.1f} to {df['lat'].max():.1f}")
    print(f"Longitude range: {df['lon'].min():.1f} to {df['lon'].max():.1f}")
    print(
        f"Temperature range: {df['temperature [oC]'].min():.1f} to {df['temperature [oC]'].max():.1f}°C"
    )
    return df


def calculate_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.asin(math.sqrt(a))

    # Radius of earth in kilometers
    r = 6371
    return c * r


def find_nearest_csv_site(target_lat, target_lon, df=None, max_distance_km=None):
    """
    Find the nearest CSV data point to the given coordinates

    Parameters:
    -----------
    target_lat : float
        Target latitude in decimal degrees
    target_lon : float
        Target longitude in decimal degrees
    df : pandas.DataFrame, optional
        CSV dataframe. If None, loads from file
    max_distance_km : float, optional
        Maximum distance in km to search within

    Returns:
    --------
    dict : Contains the nearest site data and distance
    """
    if df is None:
        df = load_csv_data()

    # Calculate distances to all points
    distances = []
    for idx, row in df.iterrows():
        dist = calculate_distance(target_lat, target_lon, row["lat"], row["lon"])
        distances.append((idx, dist, row))

    # Sort by distance
    distances.sort(key=lambda x: x[1])

    nearest_idx, nearest_dist, nearest_row = distances[0]

    # Check if within maximum distance
    if max_distance_km is not None and nearest_dist > max_distance_km:
        return {
            "found": False,
            "distance_km": nearest_dist,
            "message": f"Nearest site is {nearest_dist:.1f} km away, exceeds maximum distance of {max_distance_km} km",
        }

    return {
        "found": True,
        "index": nearest_idx,
        "distance_km": nearest_dist,
        "row": nearest_row,
        "csv_lat": nearest_row["lat"],
        "csv_lon": nearest_row["lon"],
        "target_lat": target_lat,
        "target_lon": target_lon,
    }


def run_scepter_from_coordinates(
    target_lat,
    target_lon,
    runname=None,
    simulation_type="spinup",
    erw_rate=0,
    max_distance_km=100,
    outdir_src="../scepter_output/",
):
    """
    Run SCEPTER simulation using CSV data nearest to the given coordinates

    Parameters:
    -----------
    target_lat : float
        Target latitude in decimal degrees
    target_lon : float
        Target longitude in decimal degrees
    runname : str, optional
        Custom run name. If None, generates from coordinates
    simulation_type : str
        'spinup' for control simulation or 'erw' for Enhanced Rock Weathering
    erw_rate : float
        ERW application rate in g/m²/yr (only used if simulation_type='erw')
    max_distance_km : float
        Maximum distance to search for CSV data (default: 100 km)
    outdir_src : str
        Output directory

    Returns:
    --------
    dict : Results of the simulation
    """
    print(f"Searching for CSV data near {target_lat:.3f}°N, {target_lon:.3f}°W")

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

    # Generate runname if not provided
    if runname is None:
        coord_str = f"site_{target_lat:.3f}_{target_lon:.3f}".replace(
            "-", "neg"
        ).replace(".", "p")
        if simulation_type == "erw":
            runname = f"{coord_str}_erw_{erw_rate}g"
        else:
            runname = f"{coord_str}_spinup"

    # Convert CSV row to SCEPTER parameters
    params = csv_row_to_scepter_params(nearest["row"])

    # Apply ERW modifications if requested
    if simulation_type == "erw":
        params = apply_erw_modifications(params, erw_rate, nearest["row"])

    # Override runname and coordinates with user-provided values
    params["runname"] = runname
    params["target_lat"] = target_lat
    params["target_lon"] = target_lon

    print(f"\nRunning SCEPTER simulation: {runname}")
    print(f"Simulation type: {simulation_type.upper()}")
    if simulation_type == "erw":
        print(f"ERW application rate: {erw_rate} g/m²/yr")
    print(f"Target location: {target_lat:.3f}°N, {target_lon:.3f}°W")
    print(
        f"Using CSV data from: {nearest['csv_lat']:.3f}°N, {nearest['csv_lon']:.3f}°W"
    )
    print(f"Temperature: {params['temp']:.1f}°C")
    print(f"CEC: {params['sld_varlist_cec'][0][1]:.1f} cmol/kg")
    print(f"Porosity: {params['poro']:.2f}")
    print(f"Cropland: {params['csv_data']['cropland_pct']:.0f}%")

    # Run SCEPTER simulation
    success = run_a_scepter_run(
        params["runname"],
        outdir_src,
        # ---- frame.in ----
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
        runid=params["runname"],
        # ---- switches.in ----
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
        # ---- tracers ----
        sld_list=params["sld_list"],
        aq_list=params["aq_list"],
        gas_list=params["gas_list"],
        exrxn_list=params["exrxn_list"],
        # ---- boundary values ----
        pr_list=params["pr_list"],
        rain_list=params["rain_list"],
        atm_list=params["atm_list"],
        # ---- solid phase properties ----
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
        # ---- python stuff ----
        use_local_storage=params["use_local_storage"],
    )

    result = {
        "success": success,
        "runname": runname,
        "target_coordinates": (target_lat, target_lon),
        "csv_coordinates": (nearest["csv_lat"], nearest["csv_lon"]),
        "distance_km": nearest["distance_km"],
        "parameters": params,
    }

    if success:
        print(f"✓ Simulation completed successfully: {runname}")
    else:
        print(f"❌ Simulation failed: {runname}")

    return result


def apply_erw_modifications(params, erw_rate, csv_row):
    """
    Apply Enhanced Rock Weathering modifications to SCEPTER parameters

    Parameters:
    -----------
    params : dict
        Base SCEPTER parameters
    erw_rate : float
        ERW application rate in g/m²/yr
    csv_row : pandas.Series
        Original CSV data row

    Returns:
    --------
    dict : Modified parameters for ERW simulation
    """

    # Determine if site is agricultural or natural based on cropland percentage
    cropland_pct = csv_row["cropland [%]"]
    is_agricultural = cropland_pct >= 50

    # ERW parameters
    params["fdust"] = erw_rate  # Basalt application rate

    # Set application duration based on site type
    if is_agricultural:
        params["taudust"] = 10  # 10 years for agricultural sites (annual applications)
        params["fdust2"] = erw_rate * 0.1  # Small amount of fertilizer dust
    else:
        params["taudust"] = 20  # 20 years for natural sites (longer term)
        params["fdust2"] = 0

    # Add basalt minerals to species lists
    if "fo" not in params["sld_list"]:
        params["sld_list"].extend(["fo", "ab", "an", "fa", "hb"])  # Basalt minerals

    if "si" not in params["aq_list"]:
        params["aq_list"].extend(["si", "fe2", "fe3"])  # Additional aqueous species

    # Enhanced CEC due to basalt weathering
    original_cec = params["sld_varlist_cec"][0][1]
    enhanced_cec = original_cec + (erw_rate / 1000.0) * 10  # Rough scaling

    # Update CEC for all solid species
    params["sld_varlist_cec"] = [
        ("inrt", enhanced_cec, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
        ("g2", enhanced_cec, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
        ("fo", enhanced_cec * 0.1, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
        ("ab", enhanced_cec * 0.1, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
        ("an", enhanced_cec * 0.1, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
    ]

    # Set dust source file
    params["srcfile_dust"] = "./data/dust_basalt.in"
    params["sld_varlist_dust"] = []  # Use source file instead

    # Enhanced organic matter due to improved plant growth
    # Enhance organic matter input for ERW (basalt weathering increases nutrient availability)
    original_om = params["omrain"]
    params["omrain"] = original_om * (1.0 + erw_rate / 10000.0)
    
    # Enhance SOC content for ERW (basalt improves soil structure and organic matter retention)
    if "soc_initial" in params:
        original_soc = params["soc_initial"]
        params["soc_initial"] = original_soc * (1.0 + erw_rate / 20000.0)  # Slight OM enhancement

    # Adjust mixing if agricultural site with ERW
    if is_agricultural and erw_rate > 1000:  # High ERW rate
        params["mix_scheme"] = 3  # Tilling mixing for incorporation

    return params


def csv_row_to_scepter_params(row):
    """
    Convert a CSV row to SCEPTER parameters
    Maps CSV columns to spinup.py hardcoded values
    """

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
    
    # Additional CSV values for more realistic simulations
    bs_0_20cm = row["BS_0-20cm [%]"]  # Base saturation percentage
    nitrification_rate = row["nitrification rate [kgN/ha/day]"]  # Nitrogen cycling
    
    # Calculate initial pH from soil pH (convert to H+ concentration)
    h_initial = 10**(-ph_0_30cm)  # H+ concentration from pH
    
    # Calculate initial SOC content (convert wt% to concentration)
    soc_initial = soc_0_30cm / 100.0  # Convert percentage to fraction
    
    # Calculate base cation concentrations from base saturation and CEC
    # Base saturation = (Ca + Mg + K + Na) / CEC * 100%
    base_cations_total = (bs_0_20cm / 100.0) * cec_0_30cm  # cmol/kg
    # Distribute among cations (typical ratios)
    ca_initial = base_cations_total * 0.6  # ~60% Ca
    mg_initial = base_cations_total * 0.25  # ~25% Mg  
    k_initial = base_cations_total * 0.1   # ~10% K
    na_initial = base_cations_total * 0.05  # ~5% Na

    # Create runname from coordinates
    runname = f"site_{lat:.1f}_{lon:.1f}".replace("-", "neg").replace(".", "p")

    # Map CSV data to SCEPTER parameters (following spinup.py structure)
    params = {
        # Basic simulation setup
        "runname": runname,
        "lat": lat,
        "lon": lon,
        # ---- frame.in parameters ----
        "ztot": 0.5,  # Keep standard depth
        "nz": 30,  # Keep standard grid
        "ttot": 1e5,  # Keep standard time
        "temp": temp,  # Use CSV temperature
        "fdust": 0,  # No dust initially
        "fdust2": 0,
        "taudust": 0,
        "omrain": npp * 0.3,  # Convert NPP to organic matter input (rough conversion)
        "soc_initial": soc_initial,  # Initial SOC content from CSV
        "zom": 0.25,
        "poro": porosity,  # Use CSV porosity
        "moistsrf": soil_moisture,  # Use CSV soil moisture
        "zwater": 1000,  # Keep default
        "zdust": 0.25,
        "w": erosion_mm_yr / 1000.0,  # Convert mm/yr to m/yr for erosion rate
        "q": runoff,  # Use CSV runoff
        "p": 1e-5,  # Keep default particle size
        "nstep": 10,
        "rstrt": "self",
        # ---- switches.in parameters ----
        "w_scheme": 1,
        "mix_scheme": (
            1 if cropland_pct < 50 else 3
        ),  # Fickian for natural, tilling for cropland
        "poro_iter": "false",
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
        # ---- tracers ----
        "sld_list": ["inrt", "g2"],
        "aq_list": ["ca", "k", "mg", "na"],
        "gas_list": ["pco2"],
        "exrxn_list": [],
        # ---- boundary values ----
        "pr_list": [("inrt", 1.0)],
        "rain_list": [
            ("ca", ca_initial * 1e-6),  # Convert cmol/kg to mol/L (rough conversion)
            ("mg", mg_initial * 1e-6),
            ("k", k_initial * 1e-6), 
            ("na", na_initial * 1e-6),
        ],
        "atm_list": [
            ("pco2", 3.16e-4),
            ("po2", 0.21),
            ("pnh3", 1e-50),
            ("pn2o", 1e-50),
        ],
        # ---- solid phase properties ----
        "sld_varlist_dust": [],
        "sld_varlist_cec": [
            ("inrt", cec_0_30cm, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
            ("g2", cec_0_30cm, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4),
        ],
        "sld_varlist_omrain": [("g2", 1.0)],
        "sld_varlist_kinspc": [],
        "sld_varlist_2ndslds": [],
        "srcfile_dust": None,
        "srcfile_omrain": None,
        "srcfile_cec": None,
        "srcfile_kinspc": None,
        "srcfile_2ndslds": "./data/2ndslds_def.in",
        # ---- python stuff ----
        "use_local_storage": True,
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


def run_csv_site(row_index, df, outdir_src="../scepter_output/"):
    """
    Run SCEPTER simulation for a specific CSV row
    """
    row = df.iloc[row_index]
    params = csv_row_to_scepter_params(row)

    print(f"\nRunning site {row_index + 1}/{len(df)}: {params['runname']}")
    print(f"Location: {params['lat']:.1f}°N, {params['lon']:.1f}°W")
    print(f"Temperature: {params['temp']:.1f}°C")
    print(f"CEC: {params['sld_varlist_cec'][0][1]:.1f} cmol/kg")
    print(f"Porosity: {params['poro']:.2f}")
    print(f"Cropland: {params['csv_data']['cropland_pct']:.0f}%")

    # Run SCEPTER simulation
    success = run_a_scepter_run(
        params["runname"],
        outdir_src,
        # ---- frame.in ----
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
        runid=params["runname"],
        # ---- switches.in ----
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
        # ---- tracers ----
        sld_list=params["sld_list"],
        aq_list=params["aq_list"],
        gas_list=params["gas_list"],
        exrxn_list=params["exrxn_list"],
        # ---- boundary values ----
        pr_list=params["pr_list"],
        rain_list=params["rain_list"],
        atm_list=params["atm_list"],
        # ---- solid phase properties ----
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
        # ---- python stuff ----
        use_local_storage=params["use_local_storage"],
    )

    return success, params


def run_csv_batch(start_index=0, end_index=None, max_sites=10):
    """
    Run SCEPTER simulations for a batch of CSV sites
    """
    df = load_csv_data()

    if end_index is None:
        end_index = min(start_index + max_sites, len(df))

    print(f"\nRunning batch simulation for sites {start_index} to {end_index-1}")
    print(f"Total sites: {end_index - start_index}")

    results = []
    successful = 0

    for i in range(start_index, end_index):
        try:
            success, params = run_csv_site(i, df)
            results.append(
                {
                    "index": i,
                    "runname": params["runname"],
                    "success": success,
                    "lat": params["lat"],
                    "lon": params["lon"],
                    "temp": params["temp"],
                    "cec": params["sld_varlist_cec"][0][1],
                    "cropland_pct": params["csv_data"]["cropland_pct"],
                }
            )
            if success:
                successful += 1
                print(f"✓ Site {i+1} completed successfully")
            else:
                print(f"✗ Site {i+1} failed")
        except Exception as e:
            print(f"✗ Site {i+1} error: {e}")
            results.append({"index": i, "success": False, "error": str(e)})

    print(f"\nBatch completed: {successful}/{end_index - start_index} sites successful")
    return results


def load_sites_from_json(json_file):
    """
    Load site definitions from JSON file
    """
    try:
        with open(json_file, "r") as f:
            data = json.load(f)

        sites = data.get("sites", [])
        settings = data.get("settings", {})

        print(f"Loaded {len(sites)} sites from {json_file}")
        return sites, settings

    except FileNotFoundError:
        print(f"Error: JSON file {json_file} not found")
        return [], {}
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON format in {json_file}: {e}")
        return [], {}


def run_sites_from_json(json_file):
    """
    Run SCEPTER simulations for all sites defined in JSON file
    """
    sites, settings = load_sites_from_json(json_file)

    if not sites:
        return {}

    # Get settings
    max_distance_km = settings.get("max_distance_km", 100)
    outdir = settings.get("outdir", "../scepter_output/")
    timeout_minutes = settings.get("timeout_minutes", 30)

    print(f"\nRunning {len(sites)} sites from JSON file")
    print(f"Settings: max_distance={max_distance_km}km, timeout={timeout_minutes}min")
    print("=" * 60)

    results = {}
    successful = 0

    for i, site in enumerate(sites):
        try:
            # Extract site parameters
            name = site.get("name", f"site_{i}")
            lat = site["lat"]
            lon = site["lon"]
            sim_type = site.get("simulation_type", "spinup")
            erw_rate = site.get("erw_rate", 0)
            description = site.get("description", "")

            print(f"\n[{i+1}/{len(sites)}] {name}")
            print(f"  Location: {lat:.3f}°N, {lon:.3f}°W")
            print(f"  Type: {sim_type.upper()}")
            if sim_type == "erw":
                print(f"  ERW rate: {erw_rate} g/m²/yr")
            if description:
                print(f"  Description: {description}")

            # Run simulation
            result = run_scepter_from_coordinates(
                lat, lon, name, sim_type, erw_rate, max_distance_km, outdir
            )

            results[name] = result
            if result["success"]:
                successful += 1
                print(f"  ✓ Completed successfully")
            else:
                print(f"  ❌ Failed")

        except Exception as e:
            print(f"  ❌ Error: {e}")
            results[name] = {"success": False, "error": str(e)}

    print(f"\n" + "=" * 60)
    print(f"Batch completed: {successful}/{len(sites)} sites successful")
    print("=" * 60)

    # Summary by simulation type
    spinup_sites = [s for s in sites if s.get("simulation_type", "spinup") == "spinup"]
    erw_sites = [s for s in sites if s.get("simulation_type", "spinup") == "erw"]

    if spinup_sites:
        spinup_success = sum(
            1
            for s in spinup_sites
            if results.get(s.get("name", ""), {}).get("success", False)
        )
        print(f"Spinup simulations: {spinup_success}/{len(spinup_sites)} successful")

    if erw_sites:
        erw_success = sum(
            1
            for s in erw_sites
            if results.get(s.get("name", ""), {}).get("success", False)
        )
        print(f"ERW simulations: {erw_success}/{len(erw_sites)} successful")

        # ERW rate summary
        erw_rates = [s.get("erw_rate", 0) for s in erw_sites]
        if erw_rates:
            print(f"ERW rates used: {min(erw_rates)} - {max(erw_rates)} g/m²/yr")

    return results


def analyze_csv_data():
    """
    Analyze the CSV data to understand parameter ranges
    """
    df = load_csv_data()

    print("\nCSV Data Analysis:")
    print("=" * 50)

    key_columns = [
        "temperature [oC]",
        "soil moisture [m3/m3]",
        "pH_H2O_0-30cm [-]",
        "SOC_0-30cm [wt%]",
        "CEC_0-30cm [cmol/kg]",
        "NPP [gC/m/yr]",
        "porosity [-]",
        "runoff [m/yr]",
        "erosion [mm/yr]",
        "cropland [%]",
    ]

    for col in key_columns:
        if col in df.columns:
            print(f"{col}:")
            print(f"  Range: {df[col].min():.2f} - {df[col].max():.2f}")
            print(f"  Mean: {df[col].mean():.2f}")
            print(f"  Std: {df[col].std():.2f}")
            print()

    # Analyze geographic distribution
    print("Geographic Distribution:")
    print(f"  Latitude: {df['lat'].min():.1f}° to {df['lat'].max():.1f}°N")
    print(f"  Longitude: {df['lon'].min():.1f}° to {df['lon'].max():.1f}°W")

    # Analyze land use
    print(f"\nLand Use Analysis:")
    print(f"  Natural sites (cropland < 10%): {len(df[df['cropland [%]'] < 10])}")
    print(
        f"  Mixed sites (10-50% cropland): {len(df[(df['cropland [%]'] >= 10) & (df['cropland [%]'] < 50)])}"
    )
    print(f"  Agricultural sites (>50% cropland): {len(df[df['cropland [%]'] >= 50])}")


def main():
    """
    Main function for CSV to SCEPTER conversion
    """
    import sys

    if len(sys.argv) < 2:
        print("Usage:")
        print(
            "  python run_scepter_from_coordinates.py analyze                    # Analyze CSV data"
        )
        print(
            "  python run_scepter_from_coordinates.py run <index>                # Run single site by index"
        )
        print(
            "  python run_scepter_from_coordinates.py batch <start> <end>        # Run batch by index"
        )
        print(
            "  python run_scepter_from_coordinates.py test                       # Run first 3 sites"
        )
        print(
            "  python run_scepter_from_coordinates.py coords <lat> <lon> [runname] [spinup|erw] [erw_rate]"
        )
        print(
            "  python run_scepter_from_coordinates.py find <lat> <lon>           # Find nearest CSV site"
        )
        print(
            "  python run_scepter_from_coordinates.py json <json_file>           # Run multiple sites from JSON"
        )
        print("")
        print("Examples:")
        print("  # Spinup (control) simulation")
        print(
            "  python run_scepter_from_coordinates.py coords 40.5 -85.2 Indiana_site spinup"
        )
        print("  # ERW simulation with 2000 g/m²/yr basalt")
        print(
            "  python run_scepter_from_coordinates.py coords 40.5 -85.2 Indiana_ERW erw 2000"
        )
        print("  # Run multiple sites from JSON file")
        print("  python run_scepter_from_coordinates.py json sites_template.json")
        print(
            "  python run_scepter_from_coordinates.py find 40.5 -85.2            # Find nearest data point"
        )
        sys.exit(1)

    command = sys.argv[1]

    if command == "analyze":
        analyze_csv_data()
    elif command == "run":
        if len(sys.argv) < 3:
            print("Please specify site index")
            sys.exit(1)
        index = int(sys.argv[2])
        df = load_csv_data()
        run_csv_site(index, df)
    elif command == "batch":
        start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        end = int(sys.argv[3]) if len(sys.argv) > 3 else start + 10
        run_csv_batch(start, end)
    elif command == "test":
        run_csv_batch(0, 3, 3)
    elif command == "coords":
        if len(sys.argv) < 4:
            print("Please specify latitude and longitude")
            print("Examples:")
            print("  python run_scepter_from_coordinates.py coords 40.5 -85.2")
            print(
                "  python run_scepter_from_coordinates.py coords 40.5 -85.2 my_site spinup"
            )
            print(
                "  python run_scepter_from_coordinates.py coords 40.5 -85.2 my_site erw 2000"
            )
            sys.exit(1)

        lat = float(sys.argv[2])
        lon = float(sys.argv[3])
        runname = sys.argv[4] if len(sys.argv) > 4 else None
        simulation_type = sys.argv[5] if len(sys.argv) > 5 else "spinup"
        erw_rate = float(sys.argv[6]) if len(sys.argv) > 6 else 0

        # Validate simulation type
        if simulation_type not in ["spinup", "erw"]:
            print(
                f"Error: simulation_type must be 'spinup' or 'erw', got '{simulation_type}'"
            )
            sys.exit(1)

        # Validate ERW rate
        if simulation_type == "erw" and erw_rate <= 0:
            print("Error: ERW rate must be > 0 for ERW simulations")
            print("Typical rates: 500-5000 g/m²/yr")
            sys.exit(1)

        run_scepter_from_coordinates(lat, lon, runname, simulation_type, erw_rate)
    elif command == "json":
        if len(sys.argv) < 3:
            print("Please specify JSON file")
            print(
                "Example: python run_scepter_from_coordinates.py json sites_template.json"
            )
            sys.exit(1)
        json_file = sys.argv[2]
        run_sites_from_json(json_file)
    elif command == "find":
        if len(sys.argv) < 4:
            print("Please specify latitude and longitude")
            print("Example: python run_scepter_from_coordinates.py find 40.5 -85.2")
            sys.exit(1)
        lat = float(sys.argv[2])
        lon = float(sys.argv[3])
        df = load_csv_data()
        nearest = find_nearest_csv_site(lat, lon, df)
        if nearest["found"]:
            print(f"\n✓ Nearest CSV data point:")
            print(f"  Target: {lat:.3f}°N, {lon:.3f}°W")
            print(f"  CSV: {nearest['csv_lat']:.3f}°N, {nearest['csv_lon']:.3f}°W")
            print(f"  Distance: {nearest['distance_km']:.1f} km")
            print(f"  CSV Index: {nearest['index']}")

            # Show key parameters
            row = nearest["row"]
            print(f"\n  Key Parameters:")
            print(f"    Temperature: {row['temperature [oC]']:.1f}°C")
            print(f"    CEC: {row['CEC_0-30cm [cmol/kg]']:.1f} cmol/kg")
            print(f"    Porosity: {row['porosity [-]']:.2f}")
            print(f"    Cropland: {row['cropland [%]']:.0f}%")
            print(f"    pH: {row['pH_H2O_0-30cm [-]']:.1f}")
        else:
            print(f"❌ {nearest['message']}")
    else:
        print("Invalid command")


if __name__ == "__main__":
    main()
