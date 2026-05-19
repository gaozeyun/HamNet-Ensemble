#!/bin/bash

work_dir=/root/local-disk/structure
rm -r $work_dir/input
mkdir $work_dir/input
source /opt/conda/etc/profile.d/conda.sh
conda activate deeph

tasks=("bilayer" "twisted")
tests=("bond" "d" "strain") # "bond" "d" "strain" "strain_xy"
stackings=("AB" "AA")
m_values=(1)
n_values=(2)
e_xxs=(0.001 0.002 0.003 0.004 0.005)
e_yys=(0.001 0.002 0.003 0.004 0.005)
w_amps=(0.01 0.02 0.03 0.04 0.05)
w_lambdas=(10 20 30 40 50)
w_types=("1D" "2D")

task_type=data           #wrinkle, twisting, topo, data
crystal="graphene"
m=1
n=2

if [ "$task_type" = "data" ]; then
    for task in "${tasks[@]}"; do
        mkdir -p $work_dir/input/$task
        for test in "${tests[@]}"; do
            mkdir -p $work_dir/input/$task/$test
            idx1=0
            idx2=0
            for stacking in "${stackings[@]}"; do
                for ((i=-20; i<21; i++)); do
                    outdir=$work_dir/input/$task/$test/t-$idx1-$idx2
                    mkdir -p $outdir
                    if [ "$test" = "bond" ]; then
                        dcc=$(echo "scale=3; 1.42 * (1 + 0.004 * $i)" | bc)
                        dvdw=3.35
                        xx=0
                    elif [ "$test" = "d" ]; then
                        dcc=1.42
                        dvdw=$(echo "scale=3; 3.35 * (1 + 0.015 * $i)" | bc)
                        xx=0
                    elif [ "$test" = "strain" ]; then
                        dcc=1.42
                        dvdw=3.35
                        xx=$(echo "scale=3; 0.0025 * $i" | bc)
                    fi
                    python3 ./graphene_builder.py \
                        -t $task \
                        -o $outdir \
                        -s $stacking \
                        --d_cc $dcc \
                        --d_vdw $dvdw \
                        -m $m \
                        -n $n \
                        --eps_xx $xx
                    if [ -f "$outdir/openmx.dat" ]; then
                        echo "[Success] $task: $outdir"
                    else
                        echo "[Failed] $outdir"
                        rm -r $outdir
                    fi
                    if [ $idx2 -lt 9 ]; then
                        idx2=$((idx2 + 1))
                    else
                        idx2=0
                        idx1=$((idx1 + 1))
                    fi
                done
            done
        done
    done
else
    continue
fi

for task in "${tasks[@]}"; do
    if [ "$task_type" = "topo" ]; then
        mkdir -p $work_dir/input/graphite_$task
    else
        continue
    fi
    idx1=0
    idx2=0
    for stacking in "${stackings[@]}"; do
        outdir=$work_dir/input/graphite_$task/t-$idx1-$idx2
        if [ "$task_type" = "topo" ]; then
            mkdir -p $outdir
            python3 ./graphene_builder.py \
                -c $crystal \
                -t $task \
                -o $outdir \
                -s $stacking \
                --d_cc 1.42 \
                --d_vdw 3.35 \
                -m $m \
                -n $n
            if [ -f "$outdir/openmx.dat" ]; then
                echo "[Success] $task: $outdir"
            else
                echo "[Failed] $outdir"
                rm -r $outdir
            fi
            if [ $idx2 -lt 9 ]; then
                idx2=$((idx2 + 1))
            else
                idx2=0
                idx1=$((idx1 + 1))
            fi
        else
            continue
        fi
    done
done

task="twisted"
stacking="AB"
m=1
n=2
for test in "${tests[@]}"; do
    if [ "$task_type" = "wrinkle" ]; then
        mkdir -p $work_dir/input/$test
    else
        continue
    fi
    idx1=0
    idx2=0
    d_cc_0=1.42
    d_vdw_0=3.35
    for ((i=-10; i<11; i++))
    do
        outdir=$work_dir/input/$test/t-$idx1-$idx2
        d_cc=$d_cc_0
        d_vdw=$d_vdw_0
        e_xx=0
        e_yy=0
        w_amp=0
        w_lambda=50
        w_type="1D"
        if [ "$test" = "bond" ]; then
            d_cc=$(echo "scale=3; $d_cc_0 * (1 + (0.004 * $i))" | bc)
        elif [ "$test" = "d" ]; then
            d_vdw=$(echo "scale=3; $d_vdw_0 * (1 + (0.015 * $i))" | bc)
        elif [ "$test" = "strain" ]; then
            e_xx=$(echo "scale=2; 0.005 * $i" | bc)
        elif [ "$test" = "strain_xy" ]; then
            e_xx=$(echo "scale=2; 0.005 * $i" | bc)
            e_yy=$(echo "scale=2; 0.005 * $i" | bc)
        elif [ "$test" = "wrinkle_1D_amp" ]; then
            w_amp=$(echo "scale=2; 0.5 + 0.45 * ($i + 5)" | bc)
        elif [ "$test" = "wrinkle_2D_amp" ]; then
            w_amp=$(echo "scale=2; 0.5 + 0.45 * ($i + 5)" | bc)
            w_type="2D"
        elif [ "$test" = "wrinkle_1D_lambda" ]; then
            w_amp=0.75
            w_lambda=$(echo "scale=2; 20 + 6 * ($i + 5)" | bc)
        elif [ "$test" = "wrinkle_2D_lambda" ]; then
            w_amp=0.75
            w_lambda=$(echo "scale=2; 20 + 6 * ($i + 5)" | bc)
            w_type="2D"
        fi
        if [ "$task_type" = "wrinkle" ]; then
            mkdir -p $outdir
            python3 ./graphene_builder.py \
                -t $task \
                -o $outdir \
                -s $stacking \
                --d_cc $d_cc \
                --d_vdw $d_vdw \
                --eps_xx $e_xx \
                --eps_yy $e_yy \
                --w_amp $w_amp \
                --w_lambda $w_lambda \
                --w_type $w_type \
                -m $m \
                -n $n
            if [ -f "$outdir/openmx.dat" ]; then
                echo "[Success] $task: $outdir"
            else
                echo "[Failed] $outdir"
                rm -r $outdir
            fi
            if [ $idx2 -lt 9 ]; then
                idx2=$((idx2 + 1))
            else
                idx2=0
                idx1=$((idx1 + 1))
            fi
        else
            continue
        fi
    done
done

task="twisted"
if [ "$task_type" = "twisting" ]; then
    mkdir -p $work_dir/input/$task
else
    continue
fi
for stacking in "${stackings[@]}"; do
    idx1=0
    idx2=0
    for m in "${m_values[@]}"; do
        for n in "${n_values[@]}"; do
            outdir=$work_dir/input/$task/t-$idx1-$idx2
            if [ "$task_type" = "twisting" ]; then
                mkdir -p $outdir
                python3 ./graphene_builder.py \
                    -t $task \
                    -o $outdir \
                    -s $stacking \
                    --d_cc 1.42 \
                    --d_vdw 3.35 \
                    -m $m \
                    -n $n
                if [ -f "$outdir/openmx.dat" ]; then
                    echo "[Success] $task: $outdir"
                else
                    echo "[Failed] $outdir"
                    rm -r $outdir
                fi
                if [ $idx2 -lt 9 ]; then
                    idx2=$((idx2 + 1))
                else
                    idx2=0
                    idx1=$((idx1 + 1))
                fi
            else
                continue
            fi
        done
    done
done