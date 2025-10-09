#!/usr/bin/env python3
"""
CSV Reader: Extract CSV data based on coordinates
This script provides functions to load CSV data and find the nearest data point
to target coordinates, returning the relevant parameters for SCEPTER.
"""

import pandas as pd
import numpy as np
from math import radians, cos, sin, asin, sqrt


def load_csv_data(csv_file="./data/inputdata_with_tunedpars.csv"):
    """Load CSV data with soil and climate parameters"""
    try:
        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} data points from {csv_file}")
        return df
    except FileNotFoundError:
        print(f"Error: CSV file {csv_file} not found")
        return None
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return None


def calculate_distance(lat1, lon1, lat2, lon2):
    """Calculate distance between two points using Haversine formula"""
    # Convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))

    # Radius of earth in kilometers
    r = 6371
    return c * r


def find_nearest_csv_site(target_lat, target_lon, df, max_distance_km=100):
    """Find the nearest CSV data point to target coordinates"""
    min_distance = float("inf")
    nearest_row = None

    for idx, row in df.iterrows():
        distance = calculate_distance(target_lat, target_lon, row["lat"], row["lon"])
        if distance < min_distance:
            min_distance = distance
            nearest_row = row

    if min_distance > max_distance_km:
        print(
            f"Warning: Nearest data point is {min_distance:.1f} km away (max: {max_distance_km} km)"
        )

    return nearest_row, min_distance


def get_csv_parameters(
    target_lat,
    target_lon,
    csv_file="./data/inputdata_with_tunedpars.csv",
    max_distance_km=100,
):
    """
    Main function to get CSV parameters for a target location

    Args:
        target_lat: Target latitude
        target_lon: Target longitude
        csv_file: Path to CSV file
        max_distance_km: Maximum distance to search for nearest data point

    Returns:
        dict: Dictionary containing all relevant parameters for SCEPTER
    """

    # Load CSV data
    df = load_csv_data(csv_file)
    if df is None:
        return None

    # Find nearest site
    nearest_row, distance = find_nearest_csv_site(
        target_lat, target_lon, df, max_distance_km
    )

    print(
        f"Using CSV data from: {nearest_row['lat']:.3f}°N, {nearest_row['lon']:.3f}°W"
    )
    print(f"Distance: {distance:.1f} km")

    # Print CSV data summary
    print(f"CSV data summary:")
    print(f"  Temperature: {nearest_row['temperature [oC]']:.1f}°C")
    print(f"  Soil moisture: {nearest_row['soil moisture [m3/m3]']:.3f} m³/m³")
    print(f"  CEC (0-30cm): {nearest_row['CEC_0-30cm [cmol/kg]']:.1f} cmol/kg")
    print(f"  Porosity: {nearest_row['porosity [-]']:.3f}")
    print(f"  Runoff: {nearest_row['runoff [m/yr]']:.3f} m/yr")
    print(f"  Erosion: {nearest_row['erosion [mm/yr]']:.1f} mm/yr")

    # Additional CSV parameters (not used by spinup.py)
    print(f"  SOC (0-30cm): {nearest_row['SOC_0-30cm [wt%]']:.2f} wt%")
    print(f"  Base saturation: {nearest_row['BS_0-20cm [%]']:.1f}%")
    print(f"  Cropland: {nearest_row['cropland [%]']:.1f}%")
    print(
        f"  Nitrification rate: {nearest_row['nitrification rate [kgN/ha/day]']:.3f} kgN/ha/day"
    )
    print(f"  NPP: {nearest_row['NPP [gC/m/yr]']:.1f} gC/m²/yr")
    print(f"  pH (0-30cm): {nearest_row['pH_H2O_0-30cm [-]']:.2f}")

    # Tuned parameters (if available)
    if "OC_input_gC_per_m2_per_yr" in nearest_row and pd.notna(
        nearest_row["OC_input_gC_per_m2_per_yr"]
    ):
        print(
            f"  [TUNED] OC input: {nearest_row['OC_input_gC_per_m2_per_yr']:.1f} gC/m²/yr"
        )
    if "log_tau_OC_yr" in nearest_row and pd.notna(nearest_row["log_tau_OC_yr"]):
        print(f"  [TUNED] log(tau_OC): {nearest_row['log_tau_OC_yr']:.2f}")
    if "log_KH_per_Na" in nearest_row and pd.notna(nearest_row["log_KH_per_Na"]):
        print(f"  [TUNED] log(KH/Na): {nearest_row['log_KH_per_Na']:.2f}")
    if "log_Ca_M" in nearest_row and pd.notna(nearest_row["log_Ca_M"]):
        print(f"  [TUNED] log(Ca): {nearest_row['log_Ca_M']:.2f}")

    # Extract CSV parameters
    temp_csv = nearest_row["temperature [oC]"]
    soil_moisture = nearest_row["soil moisture [m3/m3]"]
    cec_csv = nearest_row["CEC_0-30cm [cmol/kg]"]
    porosity = nearest_row["porosity [-]"]
    runoff = nearest_row["runoff [m/yr]"]
    erosion_mm_yr = nearest_row["erosion [mm/yr]"]

    # Additional CSV parameters (not used by spinup.py)
    soc_0_30cm = nearest_row["SOC_0-30cm [wt%]"]
    bs_0_20cm = nearest_row["BS_0-20cm [%]"]
    cropland_pct = nearest_row["cropland [%]"]
    nitrification_rate = nearest_row["nitrification rate [kgN/ha/day]"]
    npp = nearest_row["NPP [gC/m/yr]"]
    ph_0_30cm = nearest_row["pH_H2O_0-30cm [-]"]

    # Extract tuned parameters (if available)
    oc_input = None
    log_tau_oc = None
    log_kh_na = None
    log_ca = None

    if "OC_input_gC_per_m2_per_yr" in nearest_row and pd.notna(
        nearest_row["OC_input_gC_per_m2_per_yr"]
    ):
        oc_input = nearest_row["OC_input_gC_per_m2_per_yr"]
    if "log_tau_OC_yr" in nearest_row and pd.notna(nearest_row["log_tau_OC_yr"]):
        log_tau_oc = nearest_row["log_tau_OC_yr"]
    if "log_KH_per_Na" in nearest_row and pd.notna(nearest_row["log_KH_per_Na"]):
        log_kh_na = nearest_row["log_KH_per_Na"]
    if "log_Ca_M" in nearest_row and pd.notna(nearest_row["log_Ca_M"]):
        log_ca = nearest_row["log_Ca_M"]

    # Return parameters in format expected by spinup.py
    params = {
        # Location
        "lat": nearest_row["lat"],
        "lon": nearest_row["lon"],
        "distance_km": distance,
        "temp": temp_csv,
        "moistsrf": soil_moisture
        / porosity,  # moistsrf needs normalization by poro (moistsrf=moistsrf/poro)
        "poro": porosity,
        "q": runoff,
        "w": erosion_mm_yr / 1000.0,  # Convert mm/yr to m/yr
        "cec": cec_csv,
        # Tuned parameters (if available)
        "omrain": oc_input,  # OC input (tuned) - replaces default omrain calculation
        "log_tau_oc": log_tau_oc,  # log(tau_OC) - OC turnover time
        "log_kh_na": log_kh_na,  # log(KH/Na) - cation exchange coefficient
        "log_ca": log_ca,  # log(Ca) - calcium concentration
        # Additional CSV parameters (not directly used by spinup.py)
        # "soc_0_30cm": soc_0_30cm,
        # "bs_0_20cm": bs_0_20cm,
        # "cropland_pct": cropland_pct,
        # "nitrification_rate": nitrification_rate,
        # "npp": npp,
        # "ph_0_30cm": ph_0_30cm,
    }

    return params


def main():
    """Test function"""
    # Example usage
    target_lat = 39.34  # Scioto River, Ohio
    target_lon = -82.97

    params = get_csv_parameters(target_lat, target_lon)

    if params:
        print(f"\nExtracted parameters for spinup.py:")
        print(f"  temp = {params['temp']:.1f}")
        print(f"  moistsrf = {params['moistsrf']:.3f}")
        print(f"  poro = {params['poro']:.3f}")
        print(f"  q = {params['q']:.3f}")
        print(f"  w = {params['w']:.6f}")
        print(f"  cec = {params['cec']:.1f}")


if __name__ == "__main__":
    main()
