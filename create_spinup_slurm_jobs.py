#!/usr/bin/env python3
"""
Create SLURM job scripts for parallel SCEPTER spinup simulations
Each job gets its own SLURM job script based on parameters.json in jobs/{job_id}/
"""

import json
import os
import subprocess
import sys


def _resolve_job_path(job_input, jobs_base):
    """
    Resolve job_input to (job_id, params_path).
    job_input can be a short job_id (e.g. baseline_19579) or full path
    (e.g. /home/yhs5/project/SCEPTER/jobs/baseline_19579).
    """
    if os.path.sep in str(job_input) and os.path.isdir(job_input):
        # Full path to job folder
        job_folder = job_input.rstrip(os.path.sep)
        job_id = os.path.basename(job_folder)
        params_path = os.path.join(job_folder, "parameters.json")
    else:
        # Short job_id
        job_id = job_input
        params_path = os.path.join(jobs_base, job_id, "parameters.json")
    return job_id, params_path


def create_slurm_jobs(job_inputs, jobs_base="jobs"):
    """
    Create individual SLURM job scripts for each job.

    For each job:
    - Reads parameters from {jobs_base}/{job_id}/parameters.json (or full path)
    - Output directory: jobs/{job_id}/output
    - JSON structure: coordinate [lat, lon], location_name, job_folder, output_dir, job_id
    """
    # Create logs directory if it doesn't exist (for sbatch output)
    os.makedirs("logs", exist_ok=True)

    created = 0
    for job_input in job_inputs:
        job_id, params_path = _resolve_job_path(job_input, jobs_base)
        if not os.path.exists(params_path):
            print(f"Warning: Skipping {job_id} - {params_path} not found")
            continue

        with open(params_path, "r") as f:
            params = json.load(f)

        # Prefer job_id from JSON (canonical short form)
        job_id = params.get("job_id", job_id)

        # Extract parameters from JSON
        coordinate = params.get("coordinate", [None, None])
        target_lat = coordinate[0]
        target_lon = coordinate[1]
        site_name = params.get("location_name", job_id)
        # outdir: jobs/{job_id}/output
        outdir = os.path.join("jobs", job_id, "output")

        if target_lat is None or target_lon is None:
            print(f"Warning: Skipping {job_id} - missing coordinate in parameters.json")
            continue

        # Project root = SCEPTER folder (two layers above job_id). Path: SCEPTER/jobs/{job_id}/parameters.json
        project_root = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(params_path)))
        )

        outdir_repr = repr(outdir)
        slurm_script = f"""#!/bin/bash
#SBATCH --job-name=spinup_{job_id}
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=logs/spinup_{job_id}_%j.log
#SBATCH --error=logs/spinup_{job_id}_%j.err

# Print job info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Site: {site_name}"
echo "Coordinates: {target_lat}°N, {target_lon}°W"
echo ""


# Run from SCEPTER project root (path embedded when this script was generated)
SCEPTER_ROOT="{project_root}"
cd "$SCEPTER_ROOT" || {{ echo "ERROR: cannot cd to $SCEPTER_ROOT"; exit 1; }}
echo "Working directory: $(pwd)"

# Run spinup for this job
python3 << 'PYEOF'
import time
import spinup

site_name = "{site_name}"
target_lat = {target_lat}
target_lon = {target_lon}
outdir_src = {outdir_repr}

print(f"Running spinup for: {{site_name}}")
print(f"Coordinates: {{target_lat}}°N, {{target_lon}}°W")
print(f"Output directory: {{outdir_src}}\\n")

start_time = time.time()
result = spinup.run_single_site(site_name, target_lat, target_lon, outdir_src=outdir_src)
end_time = time.time()
elapsed = (end_time - start_time) / 60  # minutes
print(f"Spinup completed in {{elapsed:.2f}} minutes")

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

        # Write SLURM script to jobs/{job_id}/
        job_folder = os.path.dirname(params_path)
        os.makedirs(job_folder, exist_ok=True)
        script_file = os.path.join(job_folder, f"spinup_{job_id}.sh")
        with open(script_file, "w") as f:
            f.write(slurm_script)

        # Make executable
        os.chmod(script_file, 0o755)

        print(f"Created: {script_file}")

        # Run the script with sbatch
        result = subprocess.run(["sbatch", script_file], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  Submitted: {result.stdout.strip()}")
        else:
            print(f"  sbatch failed: {result.stderr.strip()}")

        created += 1

    print(f"\n{'='*80}")
    print(f"Created and submitted {created} SLURM job(s)")
    print(f"Scripts saved to jobs/{{job_id}}/spinup_{{job_id}}.sh")
    print(f"{'='*80}")
    print(f"\nMonitor jobs: squeue -u $USER")
    print(f"Check logs:  ls -lh logs/")


def discover_job_ids(jobs_base="jobs"):
    """Discover job_ids by scanning jobs/ directory for subdirectories with parameters.json"""
    if not os.path.isdir(jobs_base):
        return []
    job_ids = []
    for name in os.listdir(jobs_base):
        path = os.path.join(jobs_base, name)
        if os.path.isdir(path) and os.path.exists(
            os.path.join(path, "parameters.json")
        ):
            job_ids.append(name)
    return sorted(job_ids)


if __name__ == "__main__":
    jobs_base = "jobs"
    if len(sys.argv) > 1:
        # Job IDs passed as command line arguments
        job_ids = sys.argv[1:]
        print(f"Creating SLURM jobs for {len(job_ids)} job(s)")
        print(f"Parameters from: {jobs_base}/{{job_id}}/parameters.json")
        print(f"Output: jobs/{{job_id}}/output\n")
    else:
        # Discover job IDs from jobs/ directory
        job_ids = discover_job_ids(jobs_base)
        if not job_ids:
            print("Usage: python create_spinup_slurm_jobs.py [job_id1 job_id2 ...]")
            print("   Or: place job folders with parameters.json in jobs/ directory")
            print("\nExample: python create_spinup_slurm_jobs.py baseline_18419")
            sys.exit(1)
        print(f"Discovered {len(job_ids)} job(s) in {jobs_base}/")
        print(f"Parameters from: {jobs_base}/{{job_id}}/parameters.json")
        print(f"Output: jobs/{{job_id}}/output\n")

    create_slurm_jobs(job_ids, jobs_base=jobs_base)
