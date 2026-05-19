#!/bin/bash
# ============================================================
# MSE/MAE Batch Analysis Script
# Loop through specified directory and silently analyze test_result.h5 files
# ============================================================

SCRIPT_NAME="process.sh"

while pgrep -f "$SCRIPT_NAME" > /dev/null; do
    sleep 3600
done

source /opt/conda/etc/profile.d/conda.sh
conda activate deeph
work_dir=/root/local-disk/structure
SCRIPT_PATH="/root/local-disk/DeepH-E3/deephe3/error_analysis.py"

sources=("bilayer" "twisted")
systems=("bond" "d" "strain")

for source in "${sources[@]}"; do
    for system in "${systems[@]}"; do
        dir=$work_dir/${source}_data/carbon_eval/$system/output_dir
        struct_dir=$work_dir/${source}_data/carbon_dft/$system
        if [ ! -d "$dir" ]; then
            continue
        fi
        if [ ! -f "$dir/test_result.h5" ]; then
            echo "no h5: $dir"
            continue
        fi
        if [ -d "$struct_dir" ]; then
            python $SCRIPT_PATH -i "$dir/test_result.h5" -o "$dir/error_analyze.csv" --json "$dir/error_analyze.json" --d-cc 1.42 --d-vdw 3.35 --struct-dir "$struct_dir" -q
        else
            python $SCRIPT_PATH -i "$dir/test_result.h5" -o "$dir/error_analyze.csv" --json "$dir/error_analyze.json" --d-cc 1.42 --d-vdw 3.35 -q
        fi
        if [ ! -f "$dir/error_analyze.csv" ]; then
            echo "check: $dir"
        else
            echo "done: $dir"
        fi
    done
done