#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import sys
import json
import argparse
import math
from pathlib import Path
from ase import Atoms
from ase.build import make_supercell
from ase.io import write

GLOBAL_TOLERANCE = 1e-8

def align_to_x(atoms_obj):
    vec_a = atoms_obj.cell[0]
    angle_rad = np.arctan2(vec_a[1], vec_a[0])
    angle_deg = np.degrees(angle_rad)
    atoms_obj.rotate(a=-angle_deg, v='z', center=(0, 0, 0), rotate_cell=True)
    
    return atoms_obj

def calc_iso_arclength_strain(L_flat, amplitude, n_waves):
    if amplitude == 0.0 or n_waves == 0:
        return 0.0
    l_flat_single = L_flat / n_waves
    max_amp = l_flat_single / 2.0
    if amplitude >= max_amp:
        print(f"[!] 警告: 振幅 {amplitude:.2f} Å 超过临界值 {max_amp:.2f} Å")
        return -0.5
    l_proj_single = np.sqrt(l_flat_single**2 - 4.0 * amplitude**2)
    eps = (l_proj_single - l_flat_single) / l_flat_single
    return eps

def build_graphene_primitive(d_cc=1.42, crystal='graphene', d_vdw=3.35):
    """
    构建石墨烯/石墨单层原胞

    参数:
        d_cc: 碳-碳键长 (Å)
        crystal: 'graphene' 或 'graphite'
        d_vdw: 层间距 (Å)，仅对 graphite 有效

    石墨烯: c轴=30Å, 原子位于z=0.5 (晶胞中央)
    石墨: c轴=2*d_vdw, 原子位于z=0.25 (便于第二层叠加时落在晶胞内)
    """
    a = np.sqrt(3) * d_cc
    is_graphite = crystal.lower() == 'graphite'

    if is_graphite:
        c = 2 * d_vdw
        z_frac = 0.25
    else:
        c = 30.0
        z_frac = 0.5

    lattice_vectors = [[a, 0.0, 0.0], [0.5 * a, 0.5 * np.sqrt(3) * a, 0.0], [0.0, 0.0, c]]
    scaled_positions = [[0.0, 0.0, z_frac], [1/3, 1/3, z_frac]]
    primitive_cell = Atoms(
        symbols='C2',
        scaled_positions=scaled_positions,
        cell=lattice_vectors,
        pbc=[True, True, True]
    )
    primitive_cell.wrap()

    return primitive_cell


def build_bilayer_graphene(stacking='AA', d_cc=1.42, d_vdw=3.35, crystal='graphene'):
    """
    构建双层石墨烯/石墨结构

    参数:
        stacking: 堆垛方式 ('AA' 或 'AB')
        d_cc: 碳-碳键长 (Å)
        d_vdw: 层间距 (Å)
        crystal: 'graphene' 或 'graphite'

    石墨烯: c轴=30Å, 第一层在中央
    石墨: c轴=2*d_vdw, 第一层在z=0.25处，第二层+d_vdw后落在晶胞内
    """
    layer1 = build_graphene_primitive(d_cc=d_cc, crystal=crystal, d_vdw=d_vdw)
    layer2 = layer1.copy()
    if stacking.upper() == 'AA':
        shift = np.array([0.0, 0.0, d_vdw])
    elif stacking.upper() == 'AB':
        pos = layer1.get_positions()
        shift_xy = pos[1] - pos[0]
        shift_xy[2] = 0.0
        shift = shift_xy + np.array([0.0, 0.0, d_vdw])
    else:
        raise ValueError("仅支持 'AA' 或 'AB' 堆垛方式")
    layer2.translate(shift)
    bilayer = layer1 + layer2
    bilayer.set_pbc([True, True, True])
    bilayer.wrap()

    return bilayer


def build_twisted_bilayer(m, n, stacking='AA', d_cc=1.42, d_vdw=3.35, crystal='graphene'):
    """
    构建扭转双层石墨烯/石墨结构

    参数:
        m, n: 扭转参数 (须互质)
        stacking: 初始堆垛方式 ('AA' 或 'AB')
        d_cc: 碳-碳键长 (Å)
        d_vdw: 层间距 (Å)
        crystal: 'graphene' 或 'graphite'

    石墨烯: c轴=30Å
    石墨: c轴=2*d_vdw
    """
    prim1 = build_graphene_primitive(d_cc=d_cc, crystal=crystal, d_vdw=d_vdw)
    prim2 = prim1.copy()
    if stacking.upper() == 'AB':
        pos = prim1.get_positions()
        shift_xy = pos[1] - pos[0]
        shift_xy[2] = 0.0
        prim2.translate(shift_xy)
    elif stacking.upper() != 'AA':
        raise ValueError("仅支持 'AA' 或 'AB' 堆垛方式")
    prim2.translate([0.0, 0.0, d_vdw])
    M1 = np.array([[m,  n, 0], [-n, m+n, 0], [0,  0, 1]])
    M2 = np.array([[n,  m, 0], [-m, n+m, 0], [0,  0, 1]])
    sc1 = make_supercell(prim1, M1)
    sc2 = make_supercell(prim2, M2)
    sc1 = align_to_x(sc1)
    sc2 = align_to_x(sc2)
    tbg = sc1 + sc2
    tbg.set_pbc([True, True, True])
    tbg.wrap()

    return tbg


def apply_inplane_strain(atoms, eps_xx=0.0, eps_yy=0.0):
    cell = atoms.get_cell()[:]
    deformation_matrix = np.array([
        [1.0 + eps_xx, 0.0,          0.0],
        [0.0,          1.0 + eps_yy, 0.0],
        [0.0,          0.0,          1.0]
    ])
    new_cell = np.dot(cell, deformation_matrix)
    atoms.set_cell(new_cell, scale_atoms=True)
    atoms.set_pbc([True, True, True])
    atoms.wrap()

    return atoms

def apply_wrinkles_auto_supercell(atoms, target_wavelength=50.0, amplitude=1.0, n_waves=(1, 0), wave_type='1D', phase=0.0):
    cell = atoms.get_cell()
    len_a = np.linalg.norm(cell[0])
    len_b = np.linalg.norm(cell[1])
    rep_x = max(1, math.ceil(target_wavelength / len_a))
    rep_y = max(1, math.ceil(target_wavelength / len_b))
    if rep_x > 1 or rep_y > 1:
        print(f"[*] 触发自适应扩胞：当前 a={len_a:.2f} Å, b={len_b:.2f} Å。")
        print(f"    自动扩胞为 {rep_x} x {rep_y} 结构，以容纳 {target_wavelength} Å 的物理波长。")
        working_atoms = atoms * (rep_x, rep_y, 1)
    else:
        print(f"[*] 晶胞尺寸充足：当前 a={len_a:.2f} Å, b={len_b:.2f} Å。无需额外扩胞。")
        working_atoms = atoms.copy()

    # 获取扩胞后的晶胞尺寸
    cell_expanded = working_atoms.get_cell()
    len_a_exp = np.linalg.norm(cell_expanded[0])
    len_b_exp = np.linalg.norm(cell_expanded[1])

    n1, n2 = int(n_waves[0]), int(n_waves[1])
    eps_xx = calc_iso_arclength_strain(len_a_exp, amplitude, n1) if n1 > 0 else 0.0
    eps_yy = calc_iso_arclength_strain(len_b_exp, amplitude, n2) if n2 > 0 else 0.0
    if eps_xx < 0 or eps_yy < 0:
        print(f"[*] 物理等弧长校正: 施加预压缩应变 eps_xx={eps_xx*100:.2f}%, eps_yy={eps_yy*100:.2f}%")
        working_atoms = apply_inplane_strain(working_atoms, eps_xx=eps_xx, eps_yy=eps_yy)
    scaled_pos = working_atoms.get_scaled_positions()
    u = scaled_pos[:, 0]
    v = scaled_pos[:, 1]
    if isinstance(phase, (int, float)):
        phase_u = phase_v = float(phase)
    elif isinstance(phase, (list, tuple)) and len(phase) == 2:
        phase_u, phase_v = phase
    else:
        raise ValueError("phase 必须是单一数值，或者长度为2的元组/列表 (例如: (0.0, 1.57))")  
    if wave_type.upper() == '1D':
        delta_z = amplitude * np.sin(2 * np.pi * (n1 * u + n2 * v) + phase_u)
    elif wave_type.upper() == '2D':
        if n1 == 0 or n2 == 0:
            raise ValueError("对于 2D 褶皱，n_waves 的两个分量都不能为 0。")
        delta_z = amplitude * np.sin(2 * np.pi * n1 * u + phase_u) * np.sin(2 * np.pi * n2 * v + phase_v)
    else:
        raise ValueError("wave_type 必须是 '1D' 或 '2D'")
    pos = working_atoms.get_positions()
    pos[:, 2] += delta_z
    working_atoms.set_positions(pos)
    working_atoms.set_pbc([True, True, True])
    working_atoms.wrap()
    
    return working_atoms


def export_structure_files(atoms, output_dir, m=None, n=None, stacking=None, d_cc=1.42, d_vdw=None, eps_xx=0.0, eps_yy=0.0, w_amp=0.0, w_lambda=50.0, w_type='1D', task=''):
    target_dir = Path(output_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    cell = atoms.get_cell()
    a = np.linalg.norm(cell[0])
    k_xy = max(1, int(np.round(45.0 / a)))
    k_grid = (k_xy, k_xy, 1)
    poscar_path = target_dir / 'POSCAR'
    write(poscar_path, atoms, format='vasp', vasp5=True, sort=True)
    json_path = target_dir / 'structure.json'
    structure_data = {
        "system_info": {
            "formula": atoms.get_chemical_formula(),
            "num_atoms": len(atoms),
            "lattice_constant_a": round(float(a), 6),
            "k_grid": k_grid
        },
        "metadata": {
            "task": task, 
            "m": m, 
            "n": n, 
            "stacking": stacking, 
            "d_cc": d_cc, 
            "d_vdw": d_vdw,
            "strain": {
                "eps_xx": eps_xx, "eps_yy": eps_yy
            },
            "wrinkle": {
                "amplitude": w_amp, "target_wavelength": w_lambda, "type": w_type
            }
        },
        "lattice_vectors": cell.tolist(),
        "atoms": []
    }
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    for idx, (symbol, pos) in enumerate(zip(symbols, positions)):
        structure_data["atoms"].append({
            "index": idx + 1,
            "symbol": symbol,
            "position": [round(p, 17) for p in pos]
        })
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(structure_data, f, indent=2)
    dat_path = target_dir / 'openmx.dat'
    with open(dat_path, 'w', encoding='utf-8') as f:
        f.write("System.Name                carbon\n")
        f.write("DATA.PATH                  /share/home/zeyungao/software/openmx3.9/DFT_DATA19\n")
        f.write("HS.fileout                 on\n")
        f.write("scf.maxIter                300\n")
        f.write("scf.energycutoff           300\n")
        f.write(f"scf.Kgrid                  {k_grid[0]} {k_grid[1]} {k_grid[2]}\n")
        f.write("scf.ElectronicTemperature  300.0\n")
        f.write("scf.Mixing.Type            rmm-diisk\n")
        f.write("scf.EigenvalueSolver       Band\n")
        f.write("scf.XcType                 GGA-PBE\n")
        f.write("scf.SpinPolarization       Off\n")
        f.write("scf.criterion              1.0e-6\n")
        f.write("\nMD.Type                    Opt\n")
        f.write("MD.maxIter                 1\n")
        f.write("MD.TimeStep                1\n")
        f.write("MD.Opt.criterion           1.0e-4\n")
        f.write("\nSpecies.Number             1\n")
        f.write("<Definition.of.Atomic.Species\n")
        f.write("  C   C6.0-s2p2d1   C_PBE19\n")
        f.write("Definition.of.Atomic.Species>\n\n")
        f.write("Atoms.UnitVectors.Unit  Ang\n")
        f.write("<Atoms.UnitVectors\n")
        for vec in cell:
            f.write(f"  {vec[0]:25.17f}  {vec[1]:25.17f}  {vec[2]:25.17f}\n")
        f.write("Atoms.UnitVectors>\n\n")
        f.write(f"Atoms.Number                {len(atoms)}\n")
        f.write("Atoms.SpeciesAndCoordinates.Unit  Ang\n")
        f.write("<Atoms.SpeciesAndCoordinates\n")
        for idx, (symbol, pos) in enumerate(zip(symbols, positions)):
            f.write(f"  {idx+1:<4}  {symbol:<2}  {pos[0]:25.17f}  {pos[1]:25.17f}  {pos[2]:25.17f}     2.0     2.0\n")
        f.write("Atoms.SpeciesAndCoordinates>\n")

def main():
    parser = argparse.ArgumentParser(
        description="石墨烯结构构建器",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-t", "--task", type=str, default="twisted", choices=["single", "bilayer", "twisted"], help="任务类型: single(单层), bilayer(双层), twisted(扭转双层)。")
    parser.add_argument("-o", "--outdir", type=str, default="./input", help="输出目录路径。")
    parser.add_argument("-c", "--crystal", type=str, default="graphene", help="晶体类型: graphene(石墨烯, c轴=30Å) 或 graphite(石墨, c轴=2*d_vdw)。")

    parser.add_argument("-s", "--stacking", type=str, default="AA", choices=["AA", "AB", "aa", "ab"], help="初始堆垛方式 ('AA' 或 'AB')。")
    parser.add_argument("--d_cc", type=float, default=1.42, help="碳-碳键长 (Å)。")
    parser.add_argument("--d_vdw", type=float, default=3.35, help="层间距 (Å)。")
    parser.add_argument("-m", type=int, default=None, help="扭转参数 m (须与 n 互质)。")
    parser.add_argument("-n", type=int, default=None, help="扭转参数 n (须与 m 互质)。")

    parser.add_argument("--eps_xx", type=float, default=0.0, help="x 方向面内应变 (例如 0.01 表示 1%% 拉应变)。")
    parser.add_argument("--eps_yy", type=float, default=0.0, help="y 方向面内应变。")

    parser.add_argument("--w_amp", type=float, default=0.0, help="褶皱振幅 (Å)，设为 > 0 以施加褶皱。")
    parser.add_argument("--w_lambda", type=float, default=50.0, help="褶皱目标波长 (Å)。")
    parser.add_argument("--w_type", type=str, default="1D", choices=["1D", "2D"], help="褶皱类型: '1D' (条纹) 或 '2D' (蛋盒状)。")
    args = parser.parse_args()

    print(f"[*] 开始构建任务 [{args.task}] ...")
    try:
        is_graphite = args.crystal.lower() == "graphite"
        crystal_name = "石墨" if is_graphite else "石墨烯"
        print(f"    晶体类型：{crystal_name}")

        if args.task == "single":
            print(f"    参数: 键长={args.d_cc} Å")
            atoms = build_graphene_primitive(d_cc=args.d_cc, crystal=args.crystal, d_vdw=args.d_vdw)
        elif args.task == "bilayer":
            print(f"    参数: 堆垛={args.stacking.upper()}, 键长={args.d_cc} Å, 层间距={args.d_vdw} Å")
            atoms = build_bilayer_graphene(stacking=args.stacking, d_cc=args.d_cc, d_vdw=args.d_vdw, crystal=args.crystal)
        elif args.task == "twisted":
            if args.m is None or args.n is None:
                parser.error("任务 'twisted' 需要同时指定 -m 和 -n 参数。")
            if math.gcd(args.m, args.n) != 1:
                print(f"[-] 错误: m={args.m} 和 n={args.n} 不互质")
                exit(1)
            print(f"    参数: m={args.m}, n={args.n}, 初始堆垛={args.stacking.upper()}")
            print(f"    键长={args.d_cc} Å, 层间距={args.d_vdw} Å")
            atoms = build_twisted_bilayer(m=args.m, n=args.n, stacking=args.stacking, d_cc=args.d_cc, d_vdw=args.d_vdw, crystal=args.crystal)
        if args.eps_xx != 0.0 or args.eps_yy != 0.0:
            print(f"[*] 施加面内应变: eps_xx={args.eps_xx}, eps_yy={args.eps_yy}")
            atoms = apply_inplane_strain(atoms, eps_xx=args.eps_xx, eps_yy=args.eps_yy)
        if args.w_amp > 0.0:
            print(f"[*] 施加面外褶皱: 振幅={args.w_amp} Å, 目标波长={args.w_lambda} Å, 类型={args.w_type}")
            n_waves = (1, 0) if args.w_type.upper() == '1D' else (1, 1)
            atoms = apply_wrinkles_auto_supercell(
                atoms, 
                target_wavelength=args.w_lambda, 
                amplitude=args.w_amp, 
                n_waves=n_waves, 
                wave_type=args.w_type
            )
        print(f"[+] 构建成功！总原子数: {len(atoms)}")

        m_val = args.m if args.task == "twisted" else None
        n_val = args.n if args.task == "twisted" else None
        stacking_val = args.stacking.upper() if args.task != "single" else None
        d_vdw_val = args.d_vdw if is_graphite or args.task != "single" else None
        export_structure_files(
            atoms=atoms,
            output_dir=args.outdir,
            m=m_val, n=n_val,
            stacking=stacking_val,
            d_cc=args.d_cc,
            d_vdw=d_vdw_val,
            eps_xx=args.eps_xx,
            eps_yy=args.eps_yy,
            w_amp=args.w_amp,
            w_lambda=args.w_lambda,
            w_type=args.w_type,
            task=args.task
        )
        print(f"[+] 文件已成功导出至目录: {os.path.abspath(args.outdir)}")

    except Exception as e:
        print(f"[-] 错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()