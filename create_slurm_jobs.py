#!/usr/bin/env python3
"""
Create SLURM job scripts for parallel SCEPTER spinup simulations
Each site gets its own SLURM job script
"""

import json
import os
import sys


def create_slurm_jobs(json_file="usgs_12_sites_control.json"):
    """Create individual SLURM job scripts for each site"""

    # Read JSON file
    print(f"Reading sites from: {json_file}")
    with open(json_file, "r") as f:
        config = json.load(f)

    sites = config.get("sites", [])
    print(f"Found {len(sites)} sites\n")

    # Create logs directory if it doesn't exist
    os.makedirs("logs", exist_ok=True)
    os.makedirs("slurm_scripts", exist_ok=True)

    # Create individual SLURM scripts
    for idx, site in enumerate(sites):
        site_name = site.get("name", f"site_{idx}")
        target_lat = site.get("lat")
        target_lon = site.get("lon")

        slurm_script = f"""#!/bin/bash
#SBATCH --job-name=spinup_{site_name}
#SBATCH --account=m4259
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=logs/spinup_{site_name}_%j.log
#SBATCH --error=logs/spinup_{site_name}_%j.err

# Print job info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Site: {site_name}"
echo "Coordinates: {target_lat}°N, {target_lon}°W"
echo ""

# Load environment
module load python

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Run spinup for this site
python3 << 'PYEOF'
import spinup

site_name = "{site_name}"
target_lat = {target_lat}
target_lon = {target_lon}

print(f"Running spinup for: {{site_name}}")
print(f"Coordinates: {{target_lat}}°N, {{target_lon}}°W\\n")

result = spinup.run_single_site(site_name, target_lat, target_lon)

if result.get("success"):
    print(f"\\nSpinup completed successfully for {{site_name}}")
    exit(0)
else:
    print(f"\\nSpinup failed for {{site_name}}")
    exit(1)
PYEOF

echo ""
echo "Job finished at: $(date)"
"""

        # Write SLURM script
        script_file = f"slurm_scripts/spinup_{site_name}.sh"
        with open(script_file, "w") as f:
            f.write(slurm_script)

        # Make executable
        os.chmod(script_file, 0o755)

        print(f"Created: {script_file}")

    print(f"\n{'='*80}")
    print(f"Created {len(sites)} SLURM job scripts in slurm_scripts/")
    print(f"{'='*80}")
    print(f"\nTo submit all jobs:")
    print(f"  for script in slurm_scripts/spinup_*.sh; do sbatch $script; done")
    print(f"\nTo submit individual jobs:")
    print(f"  sbatch slurm_scripts/spinup_<site_name>.sh")
    print(f"\nMonitor jobs:")
    print(f"  squeue -u $USER")
    print(f"\nCheck logs:")
    print(f"  ls -lh logs/")


if __name__ == "__main__":
    json_file = sys.argv[1] if len(sys.argv) > 1 else "usgs_12_sites_control.json"
    create_slurm_jobs(json_file)
