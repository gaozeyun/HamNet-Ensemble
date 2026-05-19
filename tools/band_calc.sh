#!/bin/bash
#PBS -N deeph_ensemble_calc
#PBS -l nodes=1:ppn=24
#PBS -l Qlist=n24

SCRIPT_NAME="run_error_analysis.sh"

while pgrep -f "$SCRIPT_NAME" > /dev/null; do
    sleep 600
done

work_dir=/root/local-disk
sparse_calc=$work_dir/code/sparse_calc.jl
julia=$work_dir/julia-1.6.6/bin/julia
outfile=$work_dir/structure/band_calc.log
num_band=60
max_iter=600
num_threads=0

if [ $num_threads -eq 0 ]; then
    export JULIA_NUM_THREADS=$(nproc)
    echo "Setting JULIA_NUM_THREADS to $JULIA_NUM_THREADS"
else
    export JULIA_NUM_THREADS=$num_threads
fi

echo "Starting ensemble post-processing at: $(date)" > $outfile
source /opt/conda/etc/profile.d/conda.sh
conda activate deeph-pack

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
                    input_dir="/root/local-disk/structure/$dir/$data/carbon_dataset/$task/t-$x-$y"
                    output_dir="/root/local-disk/structure/$dir/$data/carbon_eval/$task/output_dir/t-$x-$y"
                    if [ ! -f "$output_dir/hamiltonians_pred.h5" ]; then
                        continue
                    fi
                    if [ "$task" = "graphite_bilayer" ]; then
                        num_band=16
                    fi

                    if [ "$dir" = "structure_test" ]; then
                        if [ "$y" != "1" ]; then
                            continue
                        fi
                    fi

                    echo "--------------------------------------------------" >> $outfile
                    echo ">>> Processing dir: $dir ..." >> $outfile

                    # calculate the kmesh
                    # element_file="$input_dir/element.dat"
                    # if [ -f "$element_file" ]; then
                    #     natoms=$(wc -l < "$element_file")
                    #     k_val=$(python3 -c "import math; print(max(1, round(150.0 / math.sqrt($natoms))))")  #let ka>150
                    #     echo "Set kmesh to [$k_val, $k_val, 1]." >> $outfile
                    # else
                    #     echo "Warning: element.dat not found in $dir! Defaulting kmesh to [1, 1, 1] for safety." >> $outfile
                    #     k_val=1
                    # fi

                    #read the json
                    info_file="$input_dir/info.json"
                    if [ -f "$info_file" ]; then
                        fermi_level=$(python3 -c "import json; f=open('$info_file'); data=json.load(f); print(data.get('fermi_level', 0.0)); f.close()")
                        echo "Successfully extracted fermi_level: $fermi_level" >> $outfile
                    else
                        echo "Warning: info.json not found! Defaulting fermi_level to -3.825771247967394" >> $outfile
                        fermi_level=-3.825771247967394
                    fi

                    targets=("0" "1" "2" "3" "dft" "mean")
                    for target in "${targets[@]}"
                    do
                        echo "--------------------------------------------------" >> $outfile
                        echo ">>> Processing target: $target ..." >> $outfile

                        #prepare the input and output
                        tmp_dir="$output_dir/tmp_calc_$target"
                        mkdir -p $tmp_dir
                        ln -sf $input_dir/rlat.dat $tmp_dir/
                        ln -sf $input_dir/orbital_types.dat $tmp_dir/
                        ln -sf $input_dir/site_positions.dat $tmp_dir/
                        [ -f "$info_file" ] && ln -sf $info_file $tmp_dir/
                        if [ "$target" == "dft" ]; then
                            ln -sf $input_dir/hamiltonians.h5 $tmp_dir/hamiltonians_pred.h5
                        elif [ "$target" == "mean" ]; then
                            ln -sf $output_dir/hamiltonians_pred.h5 $tmp_dir/hamiltonians_pred.h5
                        else
                            ln -sf $output_dir/hamiltonians_pred_$target.h5 $tmp_dir/hamiltonians_pred.h5
                        fi
                        if [ -f "$input_dir/overlaps.h5" ]; then
                            ln -sf $input_dir/overlaps.h5 $tmp_dir/
                        else
                            echo "    No overlaps.h5 found. Exit" >> $outfile
                            rm -r $tmp_dir
                            continue
                        fi
                        if [[ "$dir" == "orth-graphene" || "$dir" == "cnt" ]]; then
                        cat > $tmp_dir/band_config.json <<!
{
    "calc_job": "band",
    "which_k": 0,
    "fermi_level": $fermi_level,
    "max_iter": $max_iter,
    "num_band": $num_band,
    "k_data": ["13 0.0 0.0 0.0 0.5 0.0 0.0 Gamma X", "10 0.5 0.0 0.0 0.5 0.5 0.0 X S", "10 0.5 0.5 0.0 0.0 0.5 0.0 S Y", "17 0.0 0.5 0.0 0.0 0.0 0.0 Y Gamma"]
}
!
                        else
                        cat > $tmp_dir/band_config.json <<!
{
    "calc_job": "band",
    "which_k": 0,
    "fermi_level": $fermi_level,
    "max_iter": $max_iter,
    "num_band": $num_band,
    "k_data": ["20 0.0 0.0 0.0 0.5 0.0 0.0 Gamma M", "15 0.5 0.0 0.0 0.33333333 0.33333333 0.0 M K", "10 0.33333333 0.33333333 0.0 0.0 0.0 0.0 K Gamma"]
}
!
                        fi
# ========== DOS 计算已注释 ==========
                # cat > $tmp_dir/dos_config.json <<!
# {
#     "calc_job": "dos",
#     "fermi_level": $fermi_level,
#     "max_iter": $max_iter,
#     "num_band": $num_band,
#     "kmesh": [$k_val, $k_val, 1],
#     "epsilon": 0.005,
#     "omegas": [-0.5, 0.5, 1000]
# }
# !
# =====================================

                        #calculate the band
                        if [ ! -f "$output_dir/openmx_$target.Band" ]; then
                            echo "    Calculating Band Structure (parallel k-point mode)..." >> $outfile
                            # 使用优化后的 Julia 脚本，添加 --num_threads 参数
                            # num_threads=0 表示使用所有可用线程
                            $julia $sparse_calc \
                                --input_dir $tmp_dir \
                                --output_dir $tmp_dir \
                                --config $tmp_dir/band_config.json \
                                --ill_project true > $tmp_dir/band.log 2>&1
                            cp $tmp_dir/openmx.Band $output_dir/openmx_$target.Band
                            cp $tmp_dir/band.log $output_dir/band_$target.log
                            cp $tmp_dir/band_config.json $output_dir/band.json
                        fi

                        # ========== DOS 计算已注释 ==========
                        # if [ ! -f "$output_dir/dos_$target.dat" ]; then
                        #     echo "    Calculating DOS..." >> $outfile
                        #     $julia $sparse_calc --input_dir $tmp_dir --output_dir $tmp_dir --config $tmp_dir/dos_config.json --ill_project true > $tmp_dir/dos.log 2>&1
                        #     cp $tmp_dir/dos.dat $output_dir/dos_$target.dat
                        #     cp $tmp_dir/dos.log $output_dir/dos_$target.log
                        #     cp $tmp_dir/dos_config.json $output_dir/dos.json
                        # fi
                        # =====================================

                        #check if the calculation is successful
                        if [ -f "$output_dir/openmx_$target.Band" ]; then
                            rm -r $output_dir/band_$target.log
                            echo "    Band data saved." >> $outfile
                        else
                            echo "    Failed to calculate band." >> $outfile
                            exit 1
                        fi

                        # ========== DOS 检查已注释 ==========
                        # if [ -f "$output_dir/dos_$target.dat" ]; then
                        #     rm -r $output_dir/dos_$target.log
                        #     echo "    Dos data saved." >> $outfile
                        # else
                        #     echo "    Failed to calculate dos." >> $outfile
                        # fi
                        # =====================================

                        echo "    Target $target finished." >> $outfile

                        # ========== 额外检查已注释 ==========
                        # if [ -f "$output_dir/band_$target.log" ] && [ -f "$output_dir/dos_$target.log" ]; then
                        #     echo "Wrong in $dir, please check"
                        #     break
                        # fi
                        # =====================================

                        rm -rf $tmp_dir
                    done
                done
            done
        done
    done
done
echo "All tasks completed" >> $outfile
date >> $outfile
bash /root/local-disk/structure/band_plot.sh