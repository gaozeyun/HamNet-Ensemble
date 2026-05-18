#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Universal Carbon Materials Structure Builder
================================================================================

All structures use lower-triangular lattice format:
    L1.x  L1.y  L1.z     a1   0   0
    L2.x  L2.y  L2.z  =  b1  b2   0
    L3.x  L3.y  L3.z     0   0   c3

For graphene/graphite:
    a = 2.46 A
    L1 = (a, 0, 0)
    L2 = (a/2, a*sqrt(3)/2, 0)
================================================================================
"""
import numpy as np
import os
import json
import argparse
import random
from ase import Atoms
from ase.build import nanotube, molecule
from ase.io import write

# ====================================================
# 工具函数
# ====================================================
def clean_small_values(arr, threshold=1e-5):
    """
    将数组中小于阈值的数值归零，避免浮点精度问题。
    """
    arr = np.array(arr, dtype=float)
    arr[np.abs(arr) < threshold] = 0.0
    return arr


def clean_atoms(atoms, threshold=1e-5):
    """
    清理atoms对象中的cell和positions，将小于阈值的数值归零。
    """
    cell = clean_small_values(atoms.get_cell()[:], threshold)
    positions = clean_small_values(atoms.get_positions(), threshold)
    atoms.set_cell(cell)
    atoms.set_positions(positions)
    return atoms


# ====================================================
# 物理常数
# ====================================================
A_GRAPHENE = 2.46  # 石墨烯晶格常数 (Å)
D_INTERLAYER = 3.35  # 层间距 (Å)

# ====================================================
# 辅助函数
# ====================================================
def calc_kgrid_2d(cell):
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    k_a = int(np.ceil(50.0 / a)) if a > 0 else 1
    k_b = int(np.ceil(50.0 / b)) if b > 0 else 1
    k = max(k_a, k_b, 1)
    return k

def calc_kgrid_cnt(cell):
    c = np.linalg.norm(cell[2])
    k = int(np.ceil(50.0 / c)) if c > 0 else 1
    return k

def write_structure_json(output_dir, structure_info):
    json_file = os.path.join(output_dir, "structure.json")
    with open(json_file, 'w') as f:
        json.dump(structure_info, f, indent=2)
    return json_file

def write_openmx(filename, atoms, init_file='init.dat', k_grid=None, task='graphite',
                  cell_2d=None, c_length=None, atoms_info=None):
    if atoms is None and atoms_info is not None and cell_2d is not None and c_length is not None:
        L1, L2 = cell_2d
        cell = np.array([
            [L1[0], L1[1], 0.0],
            [L2[0], L2[1], 0.0],
            [0.0, 0.0, c_length]
        ])
        symbols = [atom[1] for atom in atoms_info]
        positions = np.array([[atom[2], atom[3], atom[4]] for atom in atoms_info])
        n_atoms = len(atoms_info)
    else:
        cell = atoms.get_cell()[:]
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        n_atoms = len(atoms)

    # 清理微小数值
    cell = clean_small_values(cell)
    positions = clean_small_values(positions)

    if k_grid is None:
        if task == 'cnt':
            k = calc_kgrid_cnt(cell)
            k_grid = (1, 1, k)
        elif task in ('tbg_graphite', 'tbg_graphene'):
            a = np.linalg.norm(cell[0])
            b = np.linalg.norm(cell[1])
            c = np.linalg.norm(cell[2])
            k_a = int(np.ceil(50.0 / a)) if a > 0 else 1
            k_b = int(np.ceil(50.0 / b)) if b > 0 else 1
            if task == 'tbg_graphite':
                k_c = int(np.ceil(50.0 / c)) if c > 0 else 1
            else:
                k_c = 1
            k_grid = (k_a, k_b, k_c)
        else:
            k = calc_kgrid_2d(cell)
            k_grid = (k, k, 1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    init_path = os.path.join(script_dir, init_file)

    with open(filename, 'w') as f:
        if os.path.exists(init_path):
            with open(init_path, 'r') as init_f:
                for line in init_f:
                    if line.strip().startswith('scf.Kgrid'):
                        f.write(f"scf.Kgrid                   {k_grid[0]} {k_grid[1]} {k_grid[2]}\n")
                    else:
                        f.write(line)
        else:
            f.write("System.Name                 carbon\n")
            f.write("DATA.PATH                   /share/home/zeyungao/software/openmx3.9/DFT_DATA19\n")
            f.write("HS.fileout                  on\n")
            f.write("scf.maxIter                 300\n")
            f.write("scf.energycutoff            300\n")
            f.write(f"scf.Kgrid                   {k_grid[0]} {k_grid[1]} {k_grid[2]}\n")
            f.write("scf.ElectronicTemperature   300.0\n")
            f.write("scf.Mixing.Type             rmm-diisk\n")
            f.write("scf.EigenvalueSolver        Band\n")
            f.write("scf.XcType                  GGA-PBE\n")
            f.write("scf.SpinPolarization        Off\n")
            f.write("scf.criterion               1.0e-6\n")
            f.write("\nMD.Type                     Opt\n")
            f.write("MD.maxIter                  1\n")
            f.write("MD.TimeStep                 1\n")
            f.write("MD.Opt.criterion            1.0e-4\n")

        f.write("\nSpecies.Number              1\n")
        f.write("<Definition.of.Atomic.Species\n")
        f.write("  C   C6.0-s2p2d1   C_PBE19\n")
        f.write("Definition.of.Atomic.Species>\n\n")

        f.write("Atoms.UnitVectors.Unit  Ang\n")
        f.write("<Atoms.UnitVectors\n")
        for i in range(3):
            f.write(f"  {cell[i,0]:.17f}   {cell[i,1]:.17f}   {cell[i,2]:.17f}\n")
        f.write("Atoms.UnitVectors>\n\n")

        f.write(f"Atoms.Number                {n_atoms}\n")
        f.write("Atoms.SpeciesAndCoordinates.Unit  Ang\n")
        f.write("<Atoms.SpeciesAndCoordinates\n")
        for idx, (sym, p) in enumerate(zip(symbols, positions)):
            f.write(f" {idx+1:>3}    {sym}    {p[0]:.17f}    {p[1]:.17f}    {p[2]:.17f}      6.0      0.0\n")
        f.write("Atoms.SpeciesAndCoordinates>\n")

    return k_grid


# ====================================================
# 任务 1：碳纳米管 (CNT)
# ====================================================
def build_cnt(radius, length, walls=1, spacing=3.35, twist_angle=0.0):
    """
    构建碳纳米管结构，沿x轴方向，使用上三角晶格格式。

    晶格格式:
        L1 = (cnt_length, 0, 0)  - CNT轴向
        L2 = (0, vacuum_y, 0)    - 垂直于CNT
        L3 = (0, 0, vacuum_z)    - 垂直于CNT
    """
    a = A_GRAPHENE
    n_inner = max(1, int(round(radius * 2 * np.pi / (a * np.sqrt(3)))))
    actual_r_inner = a * np.sqrt(3) * n_inner / (2 * np.pi)

    repeats = max(1, int(round(length / a)))
    actual_length = repeats * a

    print(f"[CNT] Target radius: {radius} A -> chirality (n=m={n_inner})")
    print(f"     Actual radius: {actual_r_inner:.2f} A, length: {actual_length:.2f} A")

    # ASE nanotube 默认沿 z 轴
    tube = nanotube(n_inner, n_inner, length=repeats)
    actual_r_outer = actual_r_inner

    if walls == 2:
        target_r_outer = actual_r_inner + spacing
        n_outer = max(1, int(round(target_r_outer * 2 * np.pi / (a * np.sqrt(3)))))
        actual_r_outer = a * np.sqrt(3) * n_outer / (2 * np.pi)
        print(f"[DWCNT] Outer radius: {actual_r_outer:.2f} A, spacing: {actual_r_outer - actual_r_inner:.2f} A")

        outer_tube = nanotube(n_outer, n_outer, length=repeats)
        if twist_angle != 0.0:
            outer_tube.rotate(twist_angle, 'z', center=(0, 0, 0))
            print(f"       Twist angle: {twist_angle} deg")
        tube = tube + outer_tube

    # 旋转 CNT 使其沿 x 轴 (从 z 轴旋转到 x 轴)
    # 绕 y 轴旋转 90 度
    tube.rotate(90, 'y', center=(0, 0, 0))

    # 设置上三角晶格
    # L1 = (length, 0, 0) - CNT轴向
    # L2 = (0, vacuum_y, 0)
    # L3 = (0, 0, vacuum_z)
    vacuum = 15.0
    vacuum_y = 2 * actual_r_outer + 2 * vacuum
    vacuum_z = 2 * actual_r_outer + 2 * vacuum

    cell = np.array([
        [actual_length, 0.0, 0.0],
        [0.0, vacuum_y, 0.0],
        [0.0, 0.0, vacuum_z]
    ])

    tube.set_cell(cell)

    # 将 CNT 移到晶胞中心 (y 和 z 方向)
    tube.center(axis=(1, 2))

    print(f"     Lattice (upper triangular): L1=({actual_length:.2f}, 0, 0), L2=(0, {vacuum_y:.2f}, 0), L3=(0, 0, {vacuum_z:.2f})")

    return tube, {
        'task': 'cnt',
        'n_inner': n_inner,
        'n_outer': n_outer if walls == 2 else None,
        'radius_inner': actual_r_inner,
        'radius_outer': actual_r_outer,
        'length': actual_length,
        'twist_angle': twist_angle,
        'atom_count': len(tube)
    }


# ====================================================
# 任务 2：C60 Fullerene
# ====================================================
def build_c60_cluster(n_molecules):
    c60 = molecule('C60')
    c60_positions = c60.get_positions()
    c60_center = c60_positions.mean(axis=0)
    c60_positions_centered = c60_positions - c60_center

    c60_radius = 3.55
    min_distance = 2 * c60_radius + 1.0
    box_size = (n_molecules ** 0.5) * min_distance * 2.0

    molecule_centers = []
    max_attempts = 5000
    retry_count = 0
    max_retries = 5

    while len(molecule_centers) < n_molecules and retry_count < max_retries:
        attempts = 0
        molecule_centers = []

        while len(molecule_centers) < n_molecules and attempts < max_attempts:
            center = np.random.uniform(c60_radius + 2, box_size - c60_radius - 2, size=3)
            overlap = False
            for existing_center in molecule_centers:
                if np.linalg.norm(center - existing_center) < min_distance:
                    overlap = True
                    break
            if not overlap:
                molecule_centers.append(center)
            attempts += 1

        if len(molecule_centers) < n_molecules:
            retry_count += 1
            box_size *= 1.5

    if len(molecule_centers) < n_molecules:
        print(f"[Warning] Only placed {len(molecule_centers)}/{n_molecules} molecules")

    all_positions = []
    for center in molecule_centers:
        for pos in c60_positions_centered:
            all_positions.append(pos + center)

    all_positions = np.array(all_positions)
    vacuum = 15.0
    box = np.max(all_positions.max(axis=0) - all_positions.min(axis=0)) + 2 * vacuum

    cell = np.diag([box, box, box])
    atoms = Atoms('C'*len(all_positions), positions=all_positions)
    atoms.set_cell(cell)
    atoms.center()

    print(f"[C60 Cluster] Molecules: {len(molecule_centers)}, Atoms: {len(atoms)}")

    return atoms, {
        'task': 'fullerene',
        'n_molecules': len(molecule_centers),
        'atom_count': len(atoms)
    }


# ====================================================
# 任务 3：TBG (Twisted Bilayer Graphene)
# ====================================================

def rotate_to_lower_triangular(L1, L2, positions_2d):
    """
    旋转晶格使其满足下三角格式，且所有元素为正数:
        L1 = (a1, 0, 0)   with a1 > 0
        L2 = (b1, b2, 0)  with b1 > 0, b2 > 0

    通过旋转变换保持晶格形状不变，只改变坐标系。
    如果旋转后 L2.x < 0，通过 L2' = L2 + k*L1 来修正。

    Args:
        L1, L2: 原始晶格矢量 (2D)
        positions_2d: 原子位置 (N, 2)

    Returns:
        L1_new, L2_new: 旋转后的晶格矢量
        positions_new: 旋转后的原子位置
    """
    # 首先旋转使 L1 对齐 x 轴正方向 (L1.y = 0, L1.x > 0)
    angle = -np.arctan2(L1[1], L1[0])

    cos_a, sin_a = np.cos(angle), np.sin(angle)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])

    L1_new = R @ L1
    L2_new = R @ L2
    positions_new = (R @ positions_2d.T).T

    # 确保 L1.x > 0
    if L1_new[0] < 0:
        L1_new = -L1_new
        L2_new = -L2_new
        positions_new = -positions_new

    # 确保 L2.y > 0
    if L2_new[1] < 0:
        L2_new[1] = -L2_new[1]
        # 相应调整原子位置的 y 坐标
        positions_new[:, 1] = L2_new[1] - positions_new[:, 1]

    # 确保 L2.x >= 0: 如果 L2.x < 0，用 L2' = L2 + k*L1
    if L2_new[0] < 0:
        # 计算 k = ceil(-L2.x / L1.x)
        k = int(np.ceil(-L2_new[0] / L1_new[0]))
        # L2' = L2 + k * L1，原子位置不变（因为这只是晶格矢量的重新选择）
        L2_new = L2_new + k * L1_new

    return L1_new, L2_new, positions_new


def find_csl_hexagonal(theta_init, min_atoms, max_atoms, max_strain=1.0, R_grid=30):
    """
    使用六方CSL算法寻找扭转双层石墨烯的超胞参数。

    搜索策略：从请求的角度向两边搜索，优先返回角度最接近请求值的解。
    """
    a = A_GRAPHENE
    # 石墨烯原胞基矢 (非上三角形式，便于CSL计算)
    # a1 = (a, 0), a2 = (a/2, a*sqrt(3)/2)
    a1 = np.array([a, 0.0])
    a2 = np.array([a/2, a * np.sqrt(3) / 2])
    A = np.column_stack([a1, a2])

    # 存储所有找到的候选
    all_candidates = []

    # 向两边搜索，共搜索 ±15° 范围
    for delta in np.arange(0, 15.0, 0.1):
        for sign in [0, 1] if delta > 0 else [0]:  # delta=0时只搜索一次
            if sign == 0:
                current_theta = theta_init + delta
            else:
                current_theta = theta_init - delta
                if current_theta < 0.5:  # 避免角度太小
                    continue

            if current_theta > 30:
                continue

            theta_rad = np.radians(current_theta)
            R_mat = np.array([
                [np.cos(theta_rad), -np.sin(theta_rad)],
                [np.sin(theta_rad),  np.cos(theta_rad)]
            ])
            B = R_mat @ A  # 旋转后的晶格

            # 搜索 CSL 整数矩阵
            for i in range(1, R_grid):
                for j in range(-R_grid + 1, R_grid):
                    if i == 0 and j == 0:
                        continue

                    M1 = np.array([[i, -j], [j, i + j]])
                    Nc1 = i**2 + i*j + j**2

                    if Nc1 * 4 > max_atoms:
                        continue

                    # 超胞晶格 L = A @ M1
                    L_matrix = A @ M1

                    # 寻找第二层的匹配矩阵 M2
                    try:
                        M2_float = np.linalg.inv(B) @ L_matrix
                    except:
                        continue

                    M2 = np.round(M2_float).astype(int)
                    if np.max(np.abs(M2_float - M2)) > 0.1:
                        continue

                    Nc2 = int(round(np.linalg.det(M2)))
                    if Nc2 <= 0:
                        continue

                    total_atoms = Nc1 * 2 + Nc2 * 2

                    if min_atoms < total_atoms < max_atoms:
                        # 计算应变
                        B_prime = L_matrix @ np.linalg.inv(M2)
                        strain = np.linalg.norm(B_prime @ np.linalg.inv(B) - np.eye(2), ord=2) * 100

                        if strain < max_strain:
                            all_candidates.append({
                                'theta': current_theta,
                                'L1': L_matrix[:, 0],
                                'L2': L_matrix[:, 1],
                                'M1': M1,
                                'M2': M2,
                                'Nc1': Nc1,
                                'Nc2': Nc2,
                                'B': B,
                                'strain': strain,
                                'total_atoms': total_atoms,
                                'theta_diff': abs(current_theta - theta_init)
                            })

    if not all_candidates:
        raise ValueError(f"No CSL found for theta ~ {theta_init}° within {min_atoms}-{max_atoms} atoms")

    # 按角度偏差排序，返回最接近请求角度的解
    all_candidates.sort(key=lambda x: x['theta_diff'])
    return all_candidates[0]


def generate_layer_atoms(a1, a2, b1, b2, M_matrix, Nc):
    """
    在超胞内生成石墨烯层的原子位置。

    Args:
        a1, a2: 原胞基矢
        b1, b2: 原胞内两个原子的位置
        M_matrix: 从原胞到超胞的变换矩阵
        Nc: 原胞数 = det(M)

    Returns:
        positions: 原子位置数组
    """
    m11, m12 = M_matrix[0, 0], M_matrix[0, 1]
    m21, m22 = M_matrix[1, 0], M_matrix[1, 1]

    coords = []
    limit = abs(m11) + abs(m12) + abs(m21) + abs(m22) + 5

    for i in range(-limit, limit + 1):
        for j in range(-limit, limit + 1):
            F1 = m22 * i - m12 * j
            F2 = -m21 * i + m11 * j
            if 0 <= F1 < Nc and 0 <= F2 < Nc:
                origin = i * a1 + j * a2
                coords.append(origin + b1)
                coords.append(origin + b2)

    return np.array(coords)


def build_tbg(theta, min_atoms, max_atoms, mode='graphene', stacking='AB'):
    """
    构建扭转双层石墨烯，输出上三角晶格格式。
    """
    # 1. 使用CSL算法找到超胞参数
    csl = find_csl_hexagonal(theta, min_atoms, max_atoms)

    L1 = csl['L1']
    L2 = csl['L2']
    M1 = csl['M1']
    M2 = csl['M2']
    Nc1 = csl['Nc1']
    Nc2 = csl['Nc2']
    theta_deg = csl['theta']
    strain = csl['strain']
    B = csl['B']

    print(f"[TBG] Initial CSL found:")
    print(f"     Theta: {theta_deg:.2f} deg, Strain: {strain:.3f}%")
    print(f"     Atoms: {Nc1*2} + {Nc2*2} = {Nc1*2 + Nc2*2}")

    # 2. 生成第一层原子 (使用石墨烯原胞基矢)
    a = A_GRAPHENE
    a1 = np.array([a, 0.0])
    a2 = np.array([a/2, a * np.sqrt(3) / 2])

    b1_L1 = np.array([0.0, 0.0])
    b2_L1 = (a1 + a2) / 3.0
    layer1_pos = generate_layer_atoms(a1, a2, b1_L1, b2_L1, M1, Nc1)

    # 3. 生成第二层原子 (扭转后的晶格)
    a1_p, a2_p = B[:, 0], B[:, 1]
    b1_L2 = np.array([0.0, 0.0])
    b2_L2 = (a1_p + a2_p) / 3.0

    # AB堆垛偏移
    shift_L2 = (a1_p + a2_p) / 3.0 if stacking == 'AB' else np.array([0.0, 0.0])
    layer2_pos = generate_layer_atoms(a1_p, a2_p, b1_L2 + shift_L2, b2_L2 + shift_L2, M2, Nc2)

    # 4. 旋转整个系统到下三角格式
    L1_rot, L2_rot, layer1_pos_rot = rotate_to_lower_triangular(L1, L2, layer1_pos)
    _, _, layer2_pos_rot = rotate_to_lower_triangular(L1, L2, layer2_pos)

    # 验证下三角格式
    assert abs(L1_rot[1]) < 1e-10, f"L1.y must be 0, got {L1_rot[1]}"

    # 5. 设置 z 坐标
    if mode == 'graphene':
        z1 = 15.0
        c_length = 30.0
    else:
        c_length = 2 * D_INTERLAYER
        z1 = D_INTERLAYER / 2

    z2 = z1 + D_INTERLAYER

    # 6. 合并两层原子
    n_layer1 = len(layer1_pos_rot)
    n_layer2 = len(layer2_pos_rot)

    positions = np.zeros((n_layer1 + n_layer2, 3))
    positions[:n_layer1, :2] = layer1_pos_rot
    positions[:n_layer1, 2] = z1
    positions[n_layer1:, :2] = layer2_pos_rot
    positions[n_layer1:, 2] = z2

    # 7. 创建晶胞
    cell = np.array([
        [L1_rot[0], L1_rot[1], 0.0],
        [L2_rot[0], L2_rot[1], 0.0],
        [0.0, 0.0, c_length]
    ])

    print(f"[TBG] Lower-triangular lattice:")
    print(f"     L1 = ({L1_rot[0]:.6f}, 0, 0)")
    print(f"     L2 = ({L2_rot[0]:.6f}, {L2_rot[1]:.6f}, 0)")
    print(f"     L3 = (0, 0, {c_length:.6f})")

    atoms = Atoms(
        symbols=['C'] * (n_layer1 + n_layer2),
        positions=positions,
        cell=cell,
        pbc=[True, True, True]
    )

    metadata = {
        'task': 'tbg',
        'theta': theta_deg,
        'mode': mode,
        'stacking': stacking,
        'strain_percent': strain,
        'atom_count': n_layer1 + n_layer2,
        'layer_atoms': [n_layer1, n_layer2],
        'c_length': c_length,
        'lattice': {
            'L1': [L1_rot[0], L1_rot[1], 0.0],
            'L2': [0.0, L2_rot[1], 0.0],
            'L3': [0.0, 0.0, c_length]
        }
    }

    return atoms, metadata


def build_and_save_tbg(theta, min_atoms, max_atoms, mode='graphene', stacking='random', output_dir='./output'):
    """
    构建并保存TBG结构（有应变版本）。
    """
    os.makedirs(output_dir, exist_ok=True)

    if stacking == 'random':
        stacking_type = random.choice(['AA', 'AB'])
    else:
        stacking_type = stacking

    atoms, metadata = build_tbg(theta, min_atoms, max_atoms, mode, stacking_type)
    metadata['stacking'] = stacking_type

    # 清理微小数值
    atoms = clean_atoms(atoms)

    poscar_file = os.path.join(output_dir, "POSCAR.vasp")
    write(poscar_file, atoms, format='vasp', vasp5=True, direct=False)

    openmx_file = os.path.join(output_dir, "openmx.dat")
    task_type = 'tbg_graphite' if mode == 'graphite' else 'tbg_graphene'
    k_grid = write_openmx(openmx_file, atoms, task=task_type)

    metadata['k_grid'] = list(k_grid)

    print(f"[TBG] Stacking: {stacking_type}, K-grid: {k_grid}")

    return metadata


# ====================================================
# 任务 3.2：完美晶格 TBG (Perfect Commensurate TBG)
# ====================================================
def calc_perfect_theta(m, n):
    """
    计算完美公度扭转角度。

    公式：cos(theta) = (n² + 4mn + m²) / [2(n² + mn + m²)]

    Args:
        m, n: CSL整数参数

    Returns:
        theta: 扭转角度（度）
    """
    numerator = n**2 + 4*m*n + m**2
    denominator = 2 * (n**2 + m*n + m**2)
    cos_theta = numerator / denominator
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # 防止浮点误差
    theta_rad = np.arccos(cos_theta)
    return np.degrees(theta_rad)


def find_perfect_csl_params(m, n, min_atoms, max_atoms):
    """
    根据m, n参数计算完美晶格CSL参数。

    完美晶格不需要施加应变，扭转角度由m, n完全确定。

    Args:
        m, n: CSL整数参数
        min_atoms, max_atoms: 原子数范围

    Returns:
        dict: CSL参数
    """
    # 计算完美扭转角度
    theta_deg = calc_perfect_theta(m, n)
    theta_rad = np.radians(theta_deg)

    # 石墨烯晶格常数
    a = A_GRAPHENE
    a1 = np.array([a, 0.0])
    a2 = np.array([a/2, a * np.sqrt(3) / 2])
    A = np.column_stack([a1, a2])

    # 旋转矩阵
    R_mat = np.array([
        [np.cos(theta_rad), -np.sin(theta_rad)],
        [np.sin(theta_rad),  np.cos(theta_rad)]
    ])
    B = R_mat @ A  # 旋转后的晶格

    # CSL超胞矩阵: M = [[m, n], [-n, m+n]]
    # 对于完美晶格，两层使用相同的整数矩阵
    M1 = np.array([[m, n], [-n, m + n]])
    M2 = M1.copy()  # 完美晶格：两层使用相同的矩阵

    Nc = m**2 + m*n + n**2  # 原胞数 = det(M)
    total_atoms = Nc * 4  # 每层2原子 × 2层

    if not (min_atoms < total_atoms < max_atoms):
        raise ValueError(f"Total atoms {total_atoms} not in range ({min_atoms}, {max_atoms})")

    # 超胞晶格矢量
    L_matrix = A @ M1
    L1 = L_matrix[:, 0]
    L2 = L_matrix[:, 1]

    return {
        'm': m,
        'n': n,
        'theta': theta_deg,
        'L1': L1,
        'L2': L2,
        'M1': M1,
        'M2': M2,
        'Nc': Nc,
        'B': B,
        'strain': 0.0,  # 完美晶格无应变
        'total_atoms': total_atoms
    }


def build_tbg_perfect(m, n, min_atoms, max_atoms, mode='graphene', stacking='AB'):
    """
    构建完美公度扭转双层石墨烯。

    与build_tbg完全隔离的独立函数。
    不施加任何应变，扭转角度由m, n完全确定。

    Args:
        m, n: CSL整数参数
        min_atoms, max_atoms: 原子数范围
        mode: 'graphene' 或 'graphite'
        stacking: 'AA' 或 'AB'

    Returns:
        atoms: ASE Atoms对象
        metadata: 元数据字典
    """
    # 1. 获取完美晶格CSL参数
    csl = find_perfect_csl_params(m, n, min_atoms, max_atoms)

    L1 = csl['L1']
    L2 = csl['L2']
    M1 = csl['M1']
    M2 = csl['M2']
    Nc = csl['Nc']
    theta_deg = csl['theta']
    B = csl['B']

    print(f"[TBG-Perfect] CSL parameters: m={m}, n={n}")
    print(f"     Perfect theta: {theta_deg:.4f} deg (strain-free)")
    print(f"     Atoms: {Nc*2} + {Nc*2} = {Nc*4}")

    # 2. 生成第一层原子
    a = A_GRAPHENE
    a1 = np.array([a, 0.0])
    a2 = np.array([a/2, a * np.sqrt(3) / 2])

    b1_L1 = np.array([0.0, 0.0])
    b2_L1 = (a1 + a2) / 3.0
    layer1_pos = generate_layer_atoms(a1, a2, b1_L1, b2_L1, M1, Nc)

    # 3. 生成第二层原子（扭转后）
    a1_p, a2_p = B[:, 0], B[:, 1]
    b1_L2 = np.array([0.0, 0.0])
    b2_L2 = (a1_p + a2_p) / 3.0

    # AB堆垛偏移
    shift_L2 = (a1_p + a2_p) / 3.0 if stacking == 'AB' else np.array([0.0, 0.0])
    layer2_pos = generate_layer_atoms(a1_p, a2_p, b1_L2 + shift_L2, b2_L2 + shift_L2, M2, Nc)

    # 4. 旋转到下三角格式
    L1_rot, L2_rot, layer1_pos_rot = rotate_to_lower_triangular(L1, L2, layer1_pos)
    _, _, layer2_pos_rot = rotate_to_lower_triangular(L1, L2, layer2_pos)

    assert abs(L1_rot[1]) < 1e-10, f"L1.y must be 0, got {L1_rot[1]}"

    # 5. 设置z坐标
    if mode == 'graphene':
        z1 = 15.0
        c_length = 30.0
    else:
        c_length = 2 * D_INTERLAYER
        z1 = D_INTERLAYER / 2

    z2 = z1 + D_INTERLAYER

    # 6. 合并原子
    n_layer1 = len(layer1_pos_rot)
    n_layer2 = len(layer2_pos_rot)

    positions = np.zeros((n_layer1 + n_layer2, 3))
    positions[:n_layer1, :2] = layer1_pos_rot
    positions[:n_layer1, 2] = z1
    positions[n_layer1:, :2] = layer2_pos_rot
    positions[n_layer1:, 2] = z2

    # 7. 创建晶胞
    cell = np.array([
        [L1_rot[0], L1_rot[1], 0.0],
        [L2_rot[0], L2_rot[1], 0.0],
        [0.0, 0.0, c_length]
    ])

    print(f"[TBG-Perfect] Lower-triangular lattice:")
    print(f"     L1 = ({L1_rot[0]:.6f}, 0, 0)")
    print(f"     L2 = ({L2_rot[0]:.6f}, {L2_rot[1]:.6f}, 0)")
    print(f"     L3 = (0, 0, {c_length:.6f})")

    atoms = Atoms(
        symbols=['C'] * (n_layer1 + n_layer2),
        positions=positions,
        cell=cell,
        pbc=[True, True, True]
    )

    metadata = {
        'task': 'tbg_perfect',
        'm': m,
        'n': n,
        'theta': theta_deg,
        'mode': mode,
        'stacking': stacking,
        'strain_percent': 0.0,
        'atom_count': n_layer1 + n_layer2,
        'layer_atoms': [n_layer1, n_layer2],
        'c_length': c_length,
        'lattice': {
            'L1': [L1_rot[0], L1_rot[1], 0.0],
            'L2': [L2_rot[0], L2_rot[1], 0.0],
            'L3': [0.0, 0.0, c_length]
        }
    }

    return atoms, metadata


def build_and_save_tbg_perfect(m, n, min_atoms, max_atoms, mode='graphene', stacking='random', output_dir='./output'):
    """
    构建并保存完美晶格TBG结构。

    与build_and_save_tbg完全隔离的独立函数。
    """
    os.makedirs(output_dir, exist_ok=True)

    if stacking == 'random':
        stacking_type = random.choice(['AA', 'AB'])
    else:
        stacking_type = stacking

    atoms, metadata = build_tbg_perfect(m, n, min_atoms, max_atoms, mode, stacking_type)
    metadata['stacking'] = stacking_type

    # 清理微小数值
    atoms = clean_atoms(atoms)

    poscar_file = os.path.join(output_dir, "POSCAR.vasp")
    write(poscar_file, atoms, format='vasp', vasp5=True, direct=False)

    openmx_file = os.path.join(output_dir, "openmx.dat")
    task_type = 'tbg_graphite' if mode == 'graphite' else 'tbg_graphene'
    k_grid = write_openmx(openmx_file, atoms, task=task_type)

    metadata['k_grid'] = list(k_grid)

    print(f"[TBG-Perfect] Stacking: {stacking_type}, K-grid: {k_grid}")

    return metadata


# ====================================================
# 主程序
# ====================================================
def main():
    parser = argparse.ArgumentParser(description="Carbon materials builder")
    parser.add_argument('--task', type=str, required=True, choices=['cnt', 'fullerene', 'tbg', 'tbg_perfect'])
    parser.add_argument('--outdir', type=str, default="./output")
    parser.add_argument('--idx1', type=int, default=1)
    parser.add_argument('--idx2', type=int, default=1)

    # CNT parameters
    parser.add_argument('--radius', type=float, default=5.0)
    parser.add_argument('--length', type=float, default=20.0)
    parser.add_argument('--walls', type=int, default=1, choices=[1, 2])
    parser.add_argument('--spacing', type=float, default=3.35)
    parser.add_argument('--twist_angle', type=float, default=0.0)

    # Fullerene parameters
    parser.add_argument('--n_molecules', type=int, default=1)

    # TBG parameters
    parser.add_argument('--theta', type=float, default=10.0)
    parser.add_argument('--min_atoms', type=int, default=50)
    parser.add_argument('--max_atoms', type=int, default=500)
    parser.add_argument('--mode', type=str, default='graphene', choices=['graphene', 'graphite'])
    parser.add_argument('--stacking', type=str, default='random', choices=['AA', 'AB', 'random'])
    parser.add_argument('--outdir_suffix', type=str, default=None)

    # TBG-perfect parameters
    parser.add_argument('--m', type=int, default=1, help='CSL integer parameter m')
    parser.add_argument('--n', type=int, default=2, help='CSL integer parameter n')

    args = parser.parse_args()

    dir_name = args.outdir_suffix if args.outdir_suffix else args.task
    task_dir = os.path.join(args.outdir, dir_name)
    output_dir = os.path.join(task_dir, f"t-{args.idx1}-{args.idx2}")
    os.makedirs(output_dir, exist_ok=True)

    print("="*50)
    print(f"Task: {args.task.upper()}")
    print(f"Output: {dir_name}/t-{args.idx1}-{args.idx2}")
    print("="*50)

    if args.task == 'cnt':
        atoms, structure_info = build_cnt(args.radius, args.length, args.walls, args.spacing, args.twist_angle)
        atoms = clean_atoms(atoms)
        write(os.path.join(output_dir, "POSCAR.vasp"), atoms, format='vasp', vasp5=True, direct=False)
        k_grid = write_openmx(os.path.join(output_dir, "openmx.dat"), atoms, task='cnt')
        structure_info['k_grid'] = list(k_grid)
        write_structure_json(output_dir, structure_info)
        print(f"\n[Done] K-grid: {k_grid}")

    elif args.task == 'fullerene':
        atoms, structure_info = build_c60_cluster(args.n_molecules)
        atoms = clean_atoms(atoms)
        write(os.path.join(output_dir, "POSCAR.vasp"), atoms, format='vasp', vasp5=True, direct=False)
        k_grid = write_openmx(os.path.join(output_dir, "openmx.dat"), atoms, task='fullerene')
        structure_info['k_grid'] = list(k_grid)
        write_structure_json(output_dir, structure_info)
        print(f"\n[Done] K-grid: {k_grid}")

    elif args.task == 'tbg':
        structure_info = build_and_save_tbg(args.theta, args.min_atoms, args.max_atoms, args.mode, args.stacking, output_dir)
        write_structure_json(output_dir, structure_info)
        print(f"\n[Done] K-grid: {tuple(structure_info['k_grid'])}")

    elif args.task == 'tbg_perfect':
        structure_info = build_and_save_tbg_perfect(args.m, args.n, args.min_atoms, args.max_atoms, args.mode, args.stacking, output_dir)
        write_structure_json(output_dir, structure_info)
        print(f"\n[Done] K-grid: {tuple(structure_info['k_grid'])}")

    print("="*50)


if __name__ == "__main__":
    main()