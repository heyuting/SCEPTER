# SCEPTER SLURM Parallel Workflow

## Overview

This workflow runs SCEPTER spinup simulations for multiple sites in parallel on SLURM, with each site running on its own node.

## Files Created

- `create_slurm_jobs.py` - Generates individual SLURM job scripts for each site
- `submit_all_spinup.sh` - Submits all jobs at once
- `check_spinup_status.sh` - Checks status of running/completed jobs
- `slurm_scripts/` - Directory containing individual SLURM job scripts (one per site)
- `logs/` - Directory for SLURM output logs

## Workflow Steps

### 1. Generate SLURM Job Scripts

```bash
python3 create_slurm_jobs.py [json_file]
```

Default: `usgs_12_sites_control.json`

This creates 12 individual SLURM scripts in `slurm_scripts/`, one for each site.

### 2. Submit All Jobs

```bash
./submit_all_spinup.sh
```

Or manually:

```bash
for script in slurm_scripts/spinup_*.sh; do sbatch $script; done
```

### 3. Monitor Jobs

```bash
# Check job queue
squeue -u $USER

# Check status script
./check_spinup_status.sh

# Watch logs in real-time
tail -f logs/spinup_Scioto_River_Chillicothe_OH_control_*.log
```

### 4. Check Results

```bash
# List completed simulations
ls -lh ../scepter_output/*control*/

# Check for completion flags
ls ../scepter_output/*control*/run_complete.txt
```

## Job Configuration

Each job is configured with:

- **Account**: smeglin
- **QOS**: regular
- **Constraint**: cpu
- **Nodes**: 1 per site
- **CPUs**: 1 per site
- **Time limit**: 24 hours
- **Logs**: `logs/spinup_<site_name>_<job_id>.log`

## Tuned Parameters

Each site automatically uses tuned parameters from `data/inputdata_with_tunedpars.csv`:

1. **Organic Carbon Input (Jorg)** - Site-specific OC flux
2. **OC Turnover Time (τorg)** - Site-specific decomposition rate
3. **Cation Exchange Coefficients (KH/Na)** - Site-specific CEC parameters
4. **Calcium Concentration ([Ca])** - Site-specific dissolved Ca
5. **Mixing Scheme** - Homogeneous mixing (2025 paper)

## Sites Processed

1. Scioto River, Chillicothe, OH
2. Scioto River, Higby, OH
3. Great Miami River, Hamilton, OH
4. Vermilion River, Danville, IL
5. Embarras River, Diona, IL
6. White River, Centerton, IN
7. Allegheny River, Salamanca, NY
8. Little Kanawha River, Palestine, WV
9. Levisa Fork, Paintsville, KY
10. Big Sandy River, Louisa, KY
11. Cumberland River, Williamsburg, KY
12. Cumberland River, Carthage, TN

## Output Structure

```
../scepter_output/
├── Scioto_River_Chillicothe_OH_control/
│   ├── flx/           # Flux output files
│   ├── prof/          # Profile output files
│   ├── frame.in       # Input files
│   ├── cec.in
│   ├── kinspc.in
│   └── ...
├── Scioto_River_Higby_OH_control/
└── ...
```

## Troubleshooting

### Jobs not submitting

- Check SLURM account: `sacctmgr show assoc user=$USER`
- Verify scripts exist: `ls slurm_scripts/`

### Jobs failing

- Check error logs: `cat logs/spinup_*_*.err`
- Check output logs: `cat logs/spinup_*_*.log`

### Simulations not completing

- Check if output files exist: `ls ../scepter_output/*/flx/`
- Manually create completion flag: `touch ../scepter_output/<site>/run_complete.txt`

## Cancel Jobs

```bash
# Cancel all your jobs
scancel -u $USER

# Cancel specific job
scancel <job_id>

# Cancel all spinup jobs
scancel -u $USER -n "spinup_*"
```

## Re-run Failed Sites

If some sites fail, you can re-submit individual jobs:

```bash
sbatch slurm_scripts/spinup_<site_name>.sh
```
