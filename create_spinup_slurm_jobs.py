#!/usr/bin/env python3
"""
Create SLURM job scripts for parallel SCEPTER spinup simulations.

Modes:
1. From sites JSON (frontend): python create_spinup_slurm_jobs.py sites.json
   - Reads sites[] with name, lat, lon; creates job folders and SLURM scripts for each.
2. From job folders: python create_spinup_slurm_jobs.py [job_id1 job_id2 ...]
   - Uses existing jobs/{job_id}/parameters.json
3. Auto-discover: python create_spinup_slurm_jobs.py
   - Finds all job folders in jobs/ with parameters.json
"""

import json
import os
import re
import subprocess
import sys


def _sanitize_job_id(name):
    """Convert site name to valid job_id (no spaces, slashes, etc.)."""
    s = str(name).strip().replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
    return re.sub(r"[^\w\-.]", "_", s) or "site"


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


def _build_and_submit_slurm_script(job_id, site_name, target_lat, target_lon, outdir, project_root, submit=True):
    """Build SLURM script for one location and optionally submit. Returns path to script."""
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

module purge
module load GCC/12.2.0
module load OpenBLAS/0.3.21-GCC-12.2.0

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
result = spinup.run_single_site(site_name, target_lat, target_lon, outdir_src=outdir_src, output_in_place=True)
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
    script_file = os.path.join(outdir, f"spinup_{job_id}.sh")
    os.makedirs(os.path.dirname(script_file), exist_ok=True)
    with open(script_file, "w") as f:
        f.write(slurm_script)
    os.chmod(script_file, 0o755)

    if submit:
        result = subprocess.run(["sbatch", script_file], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  Submitted: {result.stdout.strip()}")
        else:
            print(f"  sbatch failed: {result.stderr.strip()}")
    return script_file


def create_slurm_jobs_from_sites(sites, jobs_base="jobs", submit=True):
    """
    Create job folders, parameters.json, and SLURM scripts for each location from frontend sites JSON.

    sites: list of dicts with keys: name, lat, lon (station_id, etc. optional)
    """
    os.makedirs("logs", exist_ok=True)
    os.makedirs(jobs_base, exist_ok=True)
    project_root = os.path.dirname(os.path.abspath(os.path.join(jobs_base, ".")))

    created = 0
    for idx, site in enumerate(sites):
        site_name = site.get("name", f"site_{idx}")
        target_lat = site.get("lat") if site.get("lat") is not None else site.get("latitude")
        target_lon = site.get("lon") if site.get("lon") is not None else site.get("longitude")

        if target_lat is None or target_lon is None:
            print(f"Warning: Skipping {site_name} - missing lat/lon")
            continue

        job_id = site.get("job_id") or _sanitize_job_id(site_name)
        # Ensure unique job_id if names collide
        job_folder = os.path.join(jobs_base, job_id)
        if os.path.exists(job_folder) and not site.get("job_id"):
            job_id = f"{_sanitize_job_id(site_name)}_{idx}"
            job_folder = os.path.join(jobs_base, job_id)

        params = {
            "coordinate": [float(target_lat), float(target_lon)],
            "location_name": site_name,
            "job_id": job_id,
            "job_folder": job_folder,
        }
        params_path = os.path.join(job_folder, "parameters.json")
        os.makedirs(job_folder, exist_ok=True)
        with open(params_path, "w") as f:
            json.dump(params, f, indent=2)

        script_file = _build_and_submit_slurm_script(
            job_id, site_name, target_lat, target_lon, job_folder, project_root, submit=submit
        )
        print(f"Created: {script_file}")
        created += 1

    print(f"\n{'='*80}")
    print(f"Created and submitted {created} SLURM job(s) for {created} location(s)")
    print(f"Output: jobs/{{job_id}}/")
    print(f"{'='*80}")
    print(f"\nMonitor jobs: squeue -u $USER")
    print(f"Check logs:  ls -lh logs/")
    return created


def create_slurm_jobs(job_inputs, jobs_base="jobs", submit=True):
    """
    Create individual SLURM job scripts for each job.

    For each job:
    - Reads parameters from {jobs_base}/{job_id}/parameters.json (or full path)
    - Output directory: jobs/{job_id}
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
        # outdir: jobs/{job_id}
        outdir = os.path.join("jobs", job_id)

        if target_lat is None or target_lon is None:
            print(f"Warning: Skipping {job_id} - missing coordinate in parameters.json")
            continue

        project_root = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(params_path)))
        )
        job_folder = os.path.dirname(params_path)

        script_file = _build_and_submit_slurm_script(
            job_id, site_name, target_lat, target_lon, job_folder, project_root, submit=submit
        )
        print(f"Created: {script_file}")
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
    submit = "--no-submit" not in sys.argv
    args = [a for a in sys.argv[1:] if a != "--no-submit"]

    if len(args) >= 1 and args[0].endswith(".json") and os.path.isfile(args[0]):
        # Sites JSON from frontend
        json_path = args[0]
        with open(json_path, "r") as f:
            config = json.load(f)
        if isinstance(config, list):
            sites = config
        else:
            sites = config.get("sites", [])
            if not isinstance(sites, list):
                sites = [config] if isinstance(config, dict) else []
        print(f"Creating SLURM jobs for {len(sites)} location(s) from {json_path}")
        print(f"Output: jobs/{{job_id}}/\n")
        create_slurm_jobs_from_sites(sites, jobs_base=jobs_base, submit=submit)
    elif len(args) >= 1:
        # Job IDs passed as command line arguments
        job_ids = args
        print(f"Creating SLURM jobs for {len(job_ids)} job(s)")
        print(f"Parameters from: {jobs_base}/{{job_id}}/parameters.json")
        print(f"Output: jobs/{{job_id}}\n")
        create_slurm_jobs(job_ids, jobs_base=jobs_base, submit=submit)
    else:
        # Discover job IDs from jobs/ directory
        job_ids = discover_job_ids(jobs_base)
        if not job_ids:
            print("Usage: python create_spinup_slurm_jobs.py <sites.json>")
            print("       python create_spinup_slurm_jobs.py [job_id1 job_id2 ...]")
            print("   Or: place job folders with parameters.json in jobs/ directory")
            print("\nSites JSON format: {{\"sites\": [{{\"name\": \"...\", \"lat\": 39.3, \"lon\": -82.9}}]}}")
            sys.exit(1)
        print(f"Discovered {len(job_ids)} job(s) in {jobs_base}/")
        print(f"Parameters from: {jobs_base}/{{job_id}}/parameters.json")
        print(f"Output: jobs/{{job_id}}\n")
        create_slurm_jobs(job_ids, jobs_base=jobs_base, submit=submit)
