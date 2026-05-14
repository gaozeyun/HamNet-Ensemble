import numpy as np
import argparse
import os
import json

def convert_to_vaspkit_fermi_aligned(input_band, input_json, output_dat, output_kpath):
    if not os.path.exists(input_band) or not os.path.exists(input_json):
        print("Error: 输入文件缺失。")
        return

    with open(input_json, 'r') as jf:
        config = json.load(jf)
    with open(input_band, 'r') as f:
        lines = f.readlines()

    header = lines[0].split()
    e_fermi_raw = float(header[2])
    rlv = np.array([float(x) for x in lines[1].split()]).reshape(3, 3)

    num_paths_in_file = int(lines[2])
    line_idx = 3 + num_paths_in_file
    
    kpts_frac = []
    
    # 储存分别归属于价带和导带的能量
    all_valence = []
    all_conduction = []

    print("启动【费米基准排序】引擎...")
    while line_idx < len(lines):
        line = lines[line_idx].strip()
        if not line:
            line_idx += 1
            continue
        
        parts = line.split()
        if len(parts) != 4: break
        
        kpts_frac.append([float(x) for x in parts[1:4]])
        line_idx += 1
        
        e_k_raw = []
        while line_idx < len(lines):
            next_line = lines[line_idx].strip()
            if not next_line:
                line_idx += 1
                continue
            
            next_parts = next_line.split()
            if len(next_parts) == 4 and next_parts[0] == header[0]:
                break 
                
            e_k_raw.extend([float(x) for x in next_parts])
            line_idx += 1
            
        # 1. 转换为 eV 并对齐到费米面 (删除错误乘数)
        e_k = np.array(e_k_raw) - e_fermi_raw
        
        # 2. 强制从小到大排序
        e_k.sort()
        
        # 3. 【核心创新】：以费米面(0.0 eV)为刀，将能带劈成两半
        val = e_k[e_k <= 0.0]
        con = e_k[e_k > 0.0]
        
        all_valence.append(val)
        all_conduction.append(con)

    nkpts = len(kpts_frac)
    kpts_frac = np.array(kpts_frac)
    kpts_cart = kpts_frac @ rlv
    distances = np.insert(np.cumsum(np.linalg.norm(np.diff(kpts_cart, axis=0), axis=1)), 0, 0.0)
    
    # 4. 计算全局最大价带数和导带数
    max_v = max(len(v) for v in all_valence)
    max_c = max(len(c) for c in all_conduction)
    total_aligned_bands = max_v + max_c
    
    # 构建填充有 NaN 的统一矩阵
    energies_aligned = np.full((nkpts, total_aligned_bands), np.nan)
    
    for k in range(nkpts):
        v = all_valence[k]
        c = all_conduction[k]
        # 将价带贴着费米面往下填，缺失的深层核心态留空(NaN)
        energies_aligned[k, max_v - len(v) : max_v] = v
        # 将导带贴着费米面往上填，缺失的高能态留空(NaN)
        energies_aligned[k, max_v : max_v + len(c)] = c

    # 输出 K-Path 标签
    tick_pos = [distances[0]]
    tick_labels = []
    current_idx = 0
    for i, path_str in enumerate(config['k_data']):
        parts = path_str.split()
        if i == 0: tick_labels.append(parts[7])
        current_idx += int(parts[0])
        if current_idx <= nkpts:
            tick_pos.append(distances[current_idx - 1])
            tick_labels.append(parts[8])

    # 5. 输出对齐后的 band.dat
    with open(output_dat, 'w') as f:
        f.write("#K-Path(1/A) Energy-Level(eV)\n")
        f.write(f"# NKPTS & NBANDS:  {nkpts:d} {total_aligned_bands:d}\n")
        
        for b in range(total_aligned_bands):
            f.write(f"# Band-Index    {b + 1:d}\n")
            for k in range(nkpts):
                val_e = energies_aligned[k, b]
                if np.isnan(val_e):
                    # 如果该 K 点此处无能带(被投影掉)，输出 NaN，绘图软件会自动断线
                    f.write(f"   {distances[k]:.5f}    NaN\n")
                else:
                    f.write(f"   {distances[k]:.5f}    {val_e:.6f}\n")
            f.write("\n")

    with open(output_kpath, 'w') as f:
        for label, pos in zip(tick_labels, tick_pos):
            f.write(f"{label:15s} {pos:15.8f}\n")
            
    print(f"--- 转换成功：费米对齐完成 ---")
    print(f"统一价带数: {max_v}, 统一导带数: {max_c}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-j", "--json", required=True)
    parser.add_argument("-o", "--output_dat", default="band.dat")
    parser.add_argument("-k", "--output_kpath", default="kpath.dat")
    args = parser.parse_args()
    convert_to_vaspkit_fermi_aligned(args.input, args.json, args.output_dat, args.output_kpath)