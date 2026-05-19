#!/bin/bash
# ============================================================
# Error Analysis Visualization Script Invocation
# ============================================================

source /opt/conda/etc/profile.d/conda.sh
conda activate deeph

sources=("twisted_data" "bilayer_data")
systems=("bond" "d" "strain")

for source in "${sources[@]}"; do
    mkdir /root/local-disk/structure/pearson/$source
    for system in "${systems[@]}"; do
        dir=/root/local-disk/structure/$source/carbon_eval/$system
        mkdir /root/local-disk/structure/pearson/$source/$system
        #if [ "$system" = "bond" ]; then
        #    delta="${3:-Delta_cc}"
        #elif [ "$system" = "d" ]; then
        #    delta="${4:-Delta_vdw}"
        #else
        #    delta="${5:-Delta_strain}"
        #fi
        #python /root/local-disk/code/plot_error_analysis.py -i "$dir/error_analyze.csv" -o "$dir/error_analyze.png" --x-delta $delta
        python /root/local-disk/structure/script/plot_mae_std_correlation.py -i $dir/error_analyze.csv --figsize 10,8 --dpi 600
        cp $dir/error_analyze.csv /root/local-disk/structure/pearson/$source/$system
        cp $dir/mae_std_correlation.png /root/local-disk/structure/pearson/$source/$system
    done
done