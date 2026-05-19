#!/bin/bash

SCRIPT_NAME="deeph-band.sh"

while pgrep -f "$SCRIPT_NAME" > /dev/null; do
    sleep 3600
done

source /opt/conda/etc/profile.d/conda.sh
conda activate deeph

work_dir=/root/local-disk
plot_py=$work_dir/code/band_plotter.py
ymin=-0.2
ymax=0.2

dirs=("structure_test" "crystal_test" "wrinkle_test")
tasks=("cnt" "cnt_big" "fullerene" "graphite_bilayer" "graphite_twisted" "wrinkle_1D_amp" "wrinkle_1D_lambda" "wrinkle_2D_amp" "wrinkle_2D_lambda")
datas=("archive_test1")
xs=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
ys=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
for dir in "${dirs[@]}"; do
    for data in "${datas[@]}"; do
        for task in "${tasks[@]}"; do
            for x in "${xs[@]}"; do
                for y in "${ys[@]}"; do
                    if [ "$task" = "graphite_bilayer" ]; then
                        ymin=-0.6
                        ymax=0.6
                    fi
                    if [ "$task" = "graphite_twisted" ]; then
                        ymin=-0.2
                        ymax=0.2
                    fi
                    if [ "$task" = "cnt" ]; then
                        ymin=-0.15
                        ymax=0.15
                    fi
                    if [ "$task" = "cnt_big" ]; then
                        ymin=-0.075
                        ymax=0.075
                    fi
                    if [ "$dir" = "wrinkle_test" ]; then
                        ymin=-0.05
                        ymax=0.05
                    fi
                    if [ "$task" = "wrinkle_1D_amp" ]; then
                        ymin=-0.02
                        ymax=0.02
                    fi
                    if [ "$task" = "wrinkle_2D_amp" ]; then
                        ymin=-0.02
                        ymax=0.02
                    fi
                    output_dir=/root/local-disk/structure/$dir/$data/carbon_eval/$task/output_dir/t-$x-$y
                    if [ ! -f "$output_dir/openmx_dft.Band" ]; then
                        continue
                    fi
                    python $plot_py -i $output_dir -o $output_dir/band.png --ymin $ymin --ymax $ymax --figsize 6.5,8 #--align sort
                    ymin=-0.2
                    ymax=0.2
                done
            done
        done
    done
done