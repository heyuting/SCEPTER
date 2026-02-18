# SCEPTER Workflow Summary

## Overview

Complete workflow for running SCEPTER (Soil Chemical Evolution and Physical Transport of Elements in Regolith) simulations from site selection to final output analysis.

---

## Phase 1: Site Selection and Configuration

### 1.1 Define Sites (JSON Configuration)

**File**: `usgs_12_sites_control.json`

Each site entry contains:

- **Location**: `lat`, `lon` (decimal degrees)
- **Metadata**: `name`, `station_id`, `elevation`
- **Site Type**: `site_type` (forested/agricultural), `forest_pct`, `agriculture_pct`
- **Simulation Settings**: `simulation_type`, `erw_rate`, `description`

**Example**:

```json
{
  "name": "Scioto_River_Chillicothe_OH_control",
  "station_id": "03231500",
  "lat": 39.34,
  "lon": -82.97,
  "site_type": "agricultural",
  "forest_pct": 12.0,
  "agriculture_pct": 80.0
}
```

### 1.2 Parameter Database

**File**: `data/inputdata_with_tunedpars.csv`

Contains site-specific tuned parameters:

- **Organic Carbon Input (omrain)**: Site-specific OC flux (g C/m²/yr)
- **OC Turnover Time (log_tau_oc)**: Decomposition rate (log₁₀ years)
- **Cation Exchange (log_kh_na)**: CEC coefficients (log₁₀)
- **Calcium Concentration (log_ca)**: Dissolved Ca in rainwater (log₁₀ mol/L)
- **Physical Properties**: `temp`, `poro`, `moistsrf`, `w`, `q`, `cec`

**Module**: `csv_reader.py`

- `load_csv_data()`: Loads CSV into pandas DataFrame
- `find_nearest_csv_site()`: Finds nearest CSV entry to target coordinates (Haversine distance)
- `get_csv_parameters()`: Returns parameter dictionary for given lat/lon

---

## Phase 2: Spinup Simulation (Baseline Weathering)

### 2.1 Generate SLURM Job Scripts

**Script**: `create_spinup_slurm_jobs.py`

**Command**:

```bash
python3 create_spinup_slurm_jobs.py [json_file]
# Default: usgs_12_sites_control.json
```

**What it does**:

1. Reads site configuration from JSON
2. Creates individual SLURM script for each site in `slurm_scripts/`
3. Each script runs `spinup.run_single_site(site_name, target_lat, target_lon)`

**SLURM Configuration** (for Yale Grace):

- `--user=yhs5`
- `--partition=normal`
- `--time=24:00:00`
- `--nodes=1`, `--ntasks=1`, `--cpus-per-task=1`
- Logs: `logs/spinup_<site_name>_<job_id>.log`

### 2.2 Submit Jobs

**Script**: `submit_all_spinup.sh`

**Command**:

```bash
./submit_all_spinup.sh
```

**What it does**:

- Submits all spinup jobs to SLURM queue
- Each site runs in parallel on separate nodes

### 2.3 Monitor Jobs

**Commands**:

```bash
# Check queue status
squeue -u yhs5

# Check status script
./check_spinup_status.sh

# Watch logs
tail -f logs/spinup_<site_name>_*.log
```

---

## Phase 3: Spinup Execution

### 3.1 Core Function: `spinup.run_single_site()`

**File**: `spinup.py`

**Function Signature**:

```python
def run_single_site(site_name, target_lat, target_lon, outdir_src="../scepter_output/")
```

**Workflow**:

1. **Load Tuned Parameters** (Line 237 in `spinup.py`):

   ```python
   csv_params = csv_reader.get_csv_parameters(target_lat, target_lon)
   ```

   **This is where parameters are READ and STORED:**

   **Reading Process** (`csv_reader.py`, lines 62-189):
   - **Step 1**: `load_csv_data(csv_file)` loads entire CSV into pandas DataFrame (line 82)
   - **Step 2**: `find_nearest_csv_site()` finds nearest CSV entry to target coordinates using Haversine distance (lines 87-89)
   - **Step 3**: Extracts parameters from nearest CSV row (lines 130-160):
     - Basic: `temp`, `poro`, `moistsrf`, `q`, `w`, `cec`
     - Tuned: `omrain` (OC_input_gC_per_m2_per_yr), `log_tau_oc`, `log_kh_na`, `log_ca`
   - **Step 4**: Creates parameter dictionary (lines 163-187):
     ```python
     params = {
         "temp": temp_csv,
         "moistsrf": soil_moisture / porosity,
         "poro": porosity,
         "q": runoff,
         "w": erosion_mm_yr / 1000.0,  # Convert mm/yr to m/yr
         "cec": cec_csv,
         "omrain": oc_input,  # Tuned OC input
         "log_tau_oc": log_tau_oc,  # Tuned OC turnover time
         "log_kh_na": log_kh_na,  # Tuned CEC coefficient
         "log_ca": log_ca,  # Tuned Ca concentration
         ...
     }
     ```
   - **Step 5**: Returns dictionary `csv_params` (line 189)

   **Storage**: Parameters are stored in the `csv_params` dictionary variable in `run_single_site()` function (line 237)

2. **Extract and Use Parameters** (Lines 240-283):
   Parameters are extracted from `csv_params` dictionary and assigned to simulation variables:

   ```python
   cec = csv_params["cec"] if csv_params else 10.0
   temp = csv_params["temp"] if csv_params else 15
   poro = csv_params["poro"] if csv_params else 0.5
   moistsrf = csv_params["moistsrf"] if csv_params else 0.5
   w = csv_params["w"] if csv_params else 1e-3
   q = csv_params["q"] if csv_params else 0.3
   omrain = csv_params.get("omrain") if csv_params else 300
   ```

3. **Calculate CEC Coefficients** (Lines 243-252):
   - `logkhna = csv_params.get("log_kh_na", 5.9)`
   - `logkhk = logkhna - 1.1`
   - `logkhca = (logkhna - 0.665) * 2.0`
   - `logkhmg = (logkhna - 0.507) * 2.0`
   - `logkhal = (logkhna - 0.41) * 3.0`

4. **Set Simulation Parameters**:
   - **Grid**: `ztot = 0.5 m`, `nz = 30` layers
   - **Time**: `ttot = 1e5` years
   - **Mixing**: `mix_scheme = 2` (homogeneous, 2025 paper)
   - **Tracers**: `sld_list = ["inrt", "g2"]`, `aq_list = ["ca", "k", "mg", "na"]`, `gas_list = ["pco2"]`
   - **No Dust**: `fdust = 0` (baseline simulation)

5. **Call Core Runner** (Line 351):
   All parameters (from CSV + calculated) are passed to `run_a_scepter_run()`:
   ```python
   run_a_scepter_run(runname, outdir_src,
                     ztot=ztot, nz=nz, temp=temp, poro=poro,
                     moistsrf=moistsrf, w=w, q=q, omrain=omrain,
                     sld_varlist_cec=sld_varlist_cec, ...)
   ```

### 3.2 Core Runner: `run_a_scepter_run()`

**File**: `spinup.py`

**What it does**:

1. **Compile SCEPTER**:

   ```python
   os.system("make")  # Compiles Fortran code
   ```

2. **Create Output Directory**:

   ```python
   os.makedirs(f"{outdir}/{runname}")
   ```

3. **Generate Input Files** (via `make_inputs.py`):
   - `frame.in`: Grid, time, climate parameters
   - `switches.in`: Model configuration flags
   - `tracers.in`: Chemical species lists
   - `boundary.in`: Boundary conditions (rain, atmosphere, parent rock)
   - `cec.in`: Cation exchange coefficients
   - `kinspc.in`: Kinetic species parameters (OC turnover time)
   - `2ndslds.in`: Secondary solid phases

4. **Copy Executable**:

   ```python
   os.system(f"cp scepter {outdir}/{runname}/scepter")
   ```

5. **Run SCEPTER**:
   ```python
   subprocess.run([f"{outdir}/{runname}/scepter"], cwd=f"{outdir}/{runname}")
   ```

### 3.3 SCEPTER Fortran Execution

**Main Program**: `scepter.f90`

**Execution Flow**:

1. **Read Input Files**: `scepter_input.f90`
2. **Initialize Grid**: `scepter_makegrid.f90` (30 layers, 0.5 m total)
3. **Time Loop**: `scepter_weathering_main.f90`
   - Chemical equilibrium: `scepter_equilibrium.f90`
   - Transport: `scepter_transport.f90`
   - Kinetics: `scepter_kinetics.f90`
   - PSD evolution: `scepter_psd.f90`
4. **Write Outputs**: `scepter_IO.f90`

### 3.4 Output Files

**Location**: `../scepter_output/<site_name>/`

**Directory Structure**:

```
<site_name>/
├── flx/                    # Flux output files
│   ├── flx_sld-inrt.txt   # Solid phase fluxes
│   ├── flx_aq-ca.txt      # Aqueous species fluxes
│   ├── flx_gas-pco2.txt   # Gas phase fluxes
│   └── int_flx_co2sp-hco3.txt  # Integrated CO₂/HCO₃⁻ fluxes
├── prof/                   # Profile output files
│   ├── prof_sld-inrt-020.txt  # Solid profiles at time step 20
│   ├── prof_aq-ca-020.txt     # Aqueous profiles
│   ├── prof_gas-pco2-020.txt  # Gas profiles
│   └── charge_balance-020.txt # Charge balance details
├── frame.in                # Input files (for reference)
├── switches.in
├── tracers.in
├── boundary.in
├── cec.in
├── kinspc.in
├── 2ndslds.in
├── extrxns.in              # Extracted reactions (if any)
├── run_complete.txt        # Completion flag
└── scepter                 # Compiled executable
```

**Key Output Variables**:

- **Fluxes** (`flx/`): `tflx`, `adv`, `dif`, `rain`, `gbas`, `inrt`, `g2`, `res` (units: mol/m²/yr)
- **Profiles** (`prof/`): Depth-resolved concentrations (mol/L pore water)
- **Charge Balance** (`prof/charge_balance-*.txt`): All species concentrations for charge balance

---

## Phase 4: ERW Restart Simulation (Enhanced Rock Weathering)

### 4.1 Generate ERW SLURM Jobs

**Script**: `create_restart_add_gbas_batch.py` (or `run_restart_add_gbas_batch.py`)

**Command**:

```bash
python3 create_restart_add_gbas_batch.py [json_file]
```

**What it does**:

1. Reads same JSON configuration
2. Creates SLURM scripts for ERW simulations
3. Each script runs `restart_add_gbas.py` with basalt application

**SLURM Configuration**:

- `--user=yhs5`
- `--partition=normal`
- `--time=02:00:00`
- `--mem=8G`

### 4.2 Submit ERW Jobs

**Script**: `submit_all_restart_add_gbas.sh`

**Command**:

```bash
./submit_all_restart_add_gbas.sh
```

### 4.3 ERW Execution

**File**: `restart_add_gbas.py`

**Key Parameters**:

- **Basalt Application**: `fdust = 1000` g/m²/yr (10 t/ha/yr)
- **Particle Size**: `p80str = "320um"`
- **Simulation Time**: `ttot_restart = 100` years
- **Time Steps**: `nstep_restart = 5` (reduced for stability)
- **Mixing**: `mix_scheme = 2` (homogeneous)

**Workflow**:

1. **Load Spinup Output**: Reads final state from `../scepter_output/<site_name>_control/`
2. **Apply Basalt**: Adds basalt dust at specified rate
3. **Run Simulation**: Continues from spinup state for 100 years
4. **Output**: Saves to `../scepter_output/<site_name>_erw/`

---

## Phase 5: Output Analysis

### 5.1 Flux Analysis

**Files**: `flx/int_flx_co2sp-hco3.txt`

**Columns**:

- `time`: Simulation time (years)
- `tflx`: Total flux (mol/m²/yr)
- `adv`: Advective flux
- `dif`: Diffusive flux
- `rain`: Rainwater input
- `gbas`: Basalt contribution (ERW only)
- `inrt`, `g2`: Solid phase contributions
- `res`: Residual/other sources

**Units**: All fluxes in mol/m²/yr

### 5.2 Profile Analysis

**Files**: `prof/charge_balance-*.txt`

**Key Columns for Bicarbonate Transport**:

- `z`: Depth (m, center of layer)
- `hco3`: Free HCO₃⁻ (mol/L pore water)
- `mg(hco3)`, `na(hco3)`, `ca(hco3)`, `fe2(hco3)`: Complexed bicarbonate (mol/L pore water)
- `tot_charge`: Total charge balance

**Units**: All concentrations in mol/L of pore water

### 5.3 Analysis Scripts

- `get_int_prof.py`: Extract integrated profiles
- `get_int_prof_time.py`: Time-dependent profiles
- `get_soilpH_time.py`: pH evolution over time

---

## Complete Workflow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│ Phase 1: Site Selection & Configuration                     │
├─────────────────────────────────────────────────────────────┤
│ 1. Define sites in JSON (lat, lon, metadata)                │
│ 2. CSV contains tuned parameters (omrain, tau_oc, CEC, etc)│
│ 3. csv_reader.py matches coordinates to CSV entries         │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Phase 2: Generate SLURM Jobs                               │
├─────────────────────────────────────────────────────────────┤
│ 1. create_spinup_slurm_jobs.py → slurm_scripts/*.sh       │
│ 2. Each script calls spinup.run_single_site()               │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Phase 3: Submit & Monitor                                   │
├─────────────────────────────────────────────────────────────┤
│ 1. submit_all_spinup.sh → sbatch all jobs                  │
│ 2. squeue -u yhs5 → monitor queue                          │
│ 3. Check logs/ for progress                                │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Phase 4: Spinup Execution                                  │
├─────────────────────────────────────────────────────────────┤
│ 1. spinup.run_single_site()                                │
│    ├─ Load CSV parameters                                   │
│    ├─ Calculate CEC coefficients                           │
│    └─ Call run_a_scepter_run()                             │
│ 2. run_a_scepter_run()                                      │
│    ├─ Compile SCEPTER (make)                               │
│    ├─ Generate input files (make_inputs.py)                 │
│    └─ Execute Fortran (scepter.f90)                        │
│ 3. SCEPTER Fortran                                          │
│    ├─ Read inputs                                           │
│    ├─ Initialize grid (30 layers, 0.5 m)                   │
│    ├─ Time loop (1e5 years)                                 │
│    └─ Write outputs (flx/, prof/)                          │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Phase 5: ERW Restart (Optional)                            │
├─────────────────────────────────────────────────────────────┤
│ 1. create_restart_add_gbas_batch.py → ERW SLURM scripts   │
│ 2. submit_all_restart_add_gbas.sh → Submit ERW jobs        │
│ 3. restart_add_gbas.py                                     │
│    ├─ Load spinup final state                              │
│    ├─ Apply basalt (1000 g/m²/yr)                          │
│    └─ Run 100 years                                         │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Phase 6: Analysis                                           │
├─────────────────────────────────────────────────────────────┤
│ 1. Extract fluxes: flx/int_flx_co2sp-hco3.txt             │
│ 2. Extract profiles: prof/charge_balance-*.txt             │
│ 3. Compare forested vs agricultural sites                   │
│ 4. Analyze bicarbonate transport down profile               │
└─────────────────────────────────────────────────────────────┘
```

---

## Key Files Reference

| File                                | Purpose                                 |
| ----------------------------------- | --------------------------------------- |
| `usgs_12_sites_control.json`        | Site configuration (lat, lon, metadata) |
| `data/inputdata_with_tunedpars.csv` | Tuned parameters database               |
| `csv_reader.py`                     | Match coordinates to CSV parameters     |
| `spinup.py`                         | Spinup simulation logic                 |
| `create_spinup_slurm_jobs.py`       | Generate spinup SLURM scripts           |
| `submit_all_spinup.sh`              | Submit all spinup jobs                  |
| `restart_add_gbas.py`               | ERW restart simulation                  |
| `create_restart_add_gbas_batch.py`  | Generate ERW SLURM scripts              |
| `make_inputs.py`                    | Generate SCEPTER input files            |
| `scepter.f90`                       | Main Fortran program                    |
| `scepter_weathering_main.f90`       | Core simulation loop                    |

---

## Quick Start Commands

```bash
# 1. Generate spinup jobs
python3 create_spinup_slurm_jobs.py usgs_12_sites_control.json

# 2. Submit all spinup jobs
./submit_all_spinup.sh

# 3. Monitor jobs
squeue -u yhs5
./check_spinup_status.sh

# 4. After spinup completes, generate ERW jobs
python3 create_restart_add_gbas_batch.py usgs_12_sites_control.json

# 5. Submit ERW jobs
./submit_all_restart_add_gbas.sh

# 6. Check outputs
ls -lh ../scepter_output/*/flx/
ls -lh ../scepter_output/*/prof/
```

---

## Notes

- **Depth**: Default simulation depth is 0.5 m (50 cm) with 30 layers
- **Time**: Spinup runs for 1e5 years; ERW runs for 100 years
- **Parallelization**: Each site runs on separate SLURM node
- **Parameters**: Site-specific parameters automatically loaded from CSV based on coordinates
- **Output Location**: `../scepter_output/<site_name>/` (or `_erw` for ERW simulations)
