#!/bin/bash
#
# Submit all SCEPTER spinup jobs to SLURM
#

echo "Submitting all spinup jobs..."
echo ""

# Create logs directory if it doesn't exist
mkdir -p logs

# Counter for submitted jobs
count=0

# Submit each job
for script in slurm_scripts/spinup_*.sh; do
    if [ -f "$script" ]; then
        job_output=$(sbatch "$script")
        job_id=$(echo "$job_output" | awk '{print $NF}')
        site_name=$(basename "$script" .sh | sed 's/spinup_//')
        echo "Submitted: $site_name (Job ID: $job_id)"
        ((count++))
    fi
done

echo ""
echo "=================================="
echo "Submitted $count jobs"
echo "=================================="
echo ""
echo "Monitor jobs:    squeue -u \$USER"
echo "Cancel all jobs: scancel -u \$USER"
echo "Check logs:      ls -lh logs/"
echo ""

