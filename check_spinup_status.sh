#!/bin/bash
#
# Check status of SCEPTER spinup jobs
#

echo "=================================="
echo "SCEPTER Spinup Job Status"
echo "=================================="
echo ""

# Check SLURM queue
echo "Running jobs:"
squeue -u $USER -n "spinup_*" --format="%.18i %.30j %.8T %.10M %.9l %.6D %R" 2>/dev/null || echo "  No jobs in queue"

echo ""
echo "=================================="
echo "Completed Simulations"
echo "=================================="
echo ""

# Check for completed output directories
output_dir="../scepter_output"
if [ -d "$output_dir" ]; then
    for site_dir in $output_dir/*/; do
        if [ -d "$site_dir" ]; then
            site_name=$(basename "$site_dir")
            if [[ "$site_name" == *"control"* ]]; then
                if [ -f "$site_dir/run_complete.txt" ]; then
                    echo "✓ $site_name (COMPLETE)"
                elif [ -d "$site_dir/flx" ] || [ -d "$site_dir/prof" ]; then
                    echo "⚠ $site_name (OUTPUT EXISTS, no run_complete.txt)"
                else
                    echo "✗ $site_name (INCOMPLETE)"
                fi
            fi
        fi
    done
else
    echo "Output directory not found: $output_dir"
fi

echo ""
echo "=================================="
echo "Recent Log Files"
echo "=================================="
echo ""
ls -lht logs/spinup_*control*.log 2>/dev/null | head -5 || echo "No log files found"

echo ""

