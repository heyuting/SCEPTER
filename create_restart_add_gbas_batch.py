#!/usr/bin/env python3
"""
Run restart_add_gbas.py for all 12 USGS sites
This script creates SLURM jobs for parallel execution of ERW simulations
"""

import json
import os
import sys
import subprocess
from datetime import datetime


def create_restart_add_gbas_jobs(json_file="usgs_12_sites_control.json"):
    """Create SLURM job scripts for ERW restart simulations"""

    # Read sites configuration
    print(f"Reading sites from: {json_file}")
    with open(json_file, "r") as f:
        config = json.load(f)

    sites = config.get("sites", [])
    print(f"Found {len(sites)} sites\n")

    # Create directories
    os.makedirs("restart_logs", exist_ok=True)
    os.makedirs("restart_slurm_scripts", exist_ok=True)

    # Create job scripts
    job_scripts = []

    for idx, site in enumerate(sites):
        site_name = site.get("name", f"site_{idx}")
        station_id = site.get("station_id", f"station_{idx}")

        # Define spinup and restart names
        spinup_name = site_name  # Spinup output directory name
        restart_name = f"{site_name}_erw"  # ERW simulation name

        # Create SLURM script filename
        script_name = f"restart_add_gbas_{site_name}.sh"
        script_path = f"restart_slurm_scripts/{script_name}"

        # Create the SLURM script content
        slurm_content = f"""#!/bin/bash
#SBATCH --job-name=erw_{site_name[:20]}
#SBATCH --output=restart_logs/restart_add_gbas_{site_name}_%j.log
#SBATCH --error=restart_logs/restart_add_gbas_{site_name}_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=regular_m

# Load modules if needed
# module load python/3.9

# Set working directory
cd {os.getcwd()}

# Print job info
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Site: {site_name}"
echo "Station ID: {station_id}"
echo "Spinup name: {spinup_name}"
echo "Restart name: {restart_name}"
echo "Working directory: $(pwd)"
echo ""

# Check if spinup output exists
echo "Checking for spinup output..."
if [ ! -d "../scepter_output/{spinup_name}" ]; then
    echo "ERROR: Spinup output directory not found: ../scepter_output/{spinup_name}"
    echo "Please run spinup simulation first!"
    exit 1
fi

echo "Spinup output found: ../scepter_output/{spinup_name}"
echo ""

# Run ERW restart simulation
echo "Starting ERW restart simulation..."
python3 restart_add_gbas.py {spinup_name} {restart_name}

# Check exit status
if [ $? -eq 0 ]; then
    echo "ERW restart simulation completed successfully"
else
    echo "ERW restart simulation failed with exit code $?"
    exit 1
fi

echo "Job completed at: $(date)"
"""

        # Write SLURM script
        with open(script_path, "w") as f:
            f.write(slurm_content)

        # Make script executable
        os.chmod(script_path, 0o755)

        job_scripts.append(script_path)
        print(f"Created: {script_path}")

    return job_scripts


def submit_jobs(job_scripts):
    """Submit all SLURM jobs"""
    print(f"\nSubmitting {len(job_scripts)} jobs...")

    job_ids = []
    for script in job_scripts:
        try:
            result = subprocess.run(
                ["sbatch", script], capture_output=True, text=True, check=True
            )
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            print(f"Submitted {script}: Job ID {job_id}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit {script}: {e}")
            print(f"Error output: {e.stderr}")

    return job_ids


def create_submit_script(job_scripts):
    """Create a script to submit all jobs"""
    submit_script = "submit_all_restart_add_gbas.sh"

    content = f"""#!/bin/bash
# Submit all ERW restart jobs
# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

echo "Submitting {len(job_scripts)} ERW restart jobs..."

"""

    for script in job_scripts:
        content += f"sbatch {script}\n"

    content += f"""
echo "All jobs submitted!"
echo "Monitor with: squeue -u $USER"
echo "Check logs in: restart_logs/"
"""

    with open(submit_script, "w") as f:
        f.write(content)

    os.chmod(submit_script, 0o755)
    print(f"Created submit script: {submit_script}")

    return submit_script


def check_spinup_outputs(sites):
    """Check if spinup outputs exist for all sites"""
    print("Checking for spinup outputs...")
    missing_outputs = []

    for site in sites:
        site_name = site.get("name")
        spinup_path = f"../scepter_output/{site_name}"

        if os.path.exists(spinup_path):
            print(f"✓ Found: {spinup_path}")
        else:
            print(f"✗ Missing: {spinup_path}")
            missing_outputs.append(site_name)

    if missing_outputs:
        print(f"\n⚠️  WARNING: {len(missing_outputs)} spinup outputs are missing:")
        for site in missing_outputs:
            print(f"   - {site}")
        print("\nPlease run spinup simulations first before running ERW restarts!")
        return False
    else:
        print(f"\n✅ All {len(sites)} spinup outputs found!")
        return True


def main():
    """Main function"""
    print("=" * 60)
    print("ERW RESTART BATCH JOB CREATOR")
    print("=" * 60)

    # Read sites configuration
    with open("usgs_12_sites_control.json", "r") as f:
        config = json.load(f)

    sites = config.get("sites", [])

    # Check if spinup outputs exist
    if not check_spinup_outputs(sites):
        print("\n❌ Cannot proceed without spinup outputs.")
        print("Please run spinup simulations first:")
        print("  ./submit_all_spinup.sh")
        return

    print()

    # Create job scripts
    job_scripts = create_restart_add_gbas_jobs(json_file="usgs_12_sites_control.json")

    # Create submit script
    submit_script = create_submit_script(job_scripts)

    print("\n" + "=" * 60)
    print("SETUP COMPLETE!")
    print("=" * 60)
    print(f"Created {len(job_scripts)} SLURM job scripts")
    print(f"Created submit script: {submit_script}")
    print()
    print("To submit all jobs:")
    print(f"  ./{submit_script}")
    print()
    print("To submit individual jobs:")
    for script in job_scripts[:3]:  # Show first 3 as examples
        print(f"  sbatch {script}")
    if len(job_scripts) > 3:
        print(f"  ... and {len(job_scripts) - 3} more")
    print()
    print("Monitor jobs:")
    print("  squeue -u $USER")
    print("  watch -n 30 'squeue -u $USER'")
    print()
    print("Check logs:")
    print("  ls -lh restart_logs/")
    print("  tail -f restart_logs/restart_add_gbas_*_*.log")
    print()
    print("Expected output directories:")
    for site in sites[:3]:  # Show first 3 as examples
        site_name = site.get("name")
        print(f"  ../scepter_output/{site_name}_erw/")
    if len(sites) > 3:
        print(f"  ... and {len(sites) - 3} more")


if __name__ == "__main__":
    main()
