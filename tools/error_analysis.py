#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MSE/MAE Analysis Script
Extract prediction error statistics from test_result.h5 file

Features:
1. Output MSE and MAE for each individual system
2. Output overall MSE and MAE for all systems
3. Support structure filtering (include/exclude)
4. Support exporting to CSV file
5. Parse hamiltonians_std.h5, calculate mean standard deviation of non-zero entries
6. Read structure.json, calculate d_cc, d_vdw change ratios and strain information
"""

import os
import re
import csv
import json
import argparse
import numpy as np
import h5py
from pathlib import Path
from datetime import datetime


def maskmse(input_arr, target, mask):
    """
    计算带掩码的均方误差 (MSE)

    Args:
        input_arr: 预测值数组
        target: 真实值数组
        mask: 有效数据掩码

    Returns:
        float: MSE 值
    """
    assert input_arr.shape == target.shape == mask.shape
    mse = np.power(np.abs(input_arr - target), 2)
    mse = mse[mask].mean()
    return mse


def maskmae(input_arr, target, mask):
    """
    计算带掩码的平均绝对误差 (MAE)

    Args:
        input_arr: 预测值数组
        target: 真实值数组
        mask: 有效数据掩码

    Returns:
        float: MAE 值
    """
    assert input_arr.shape == target.shape == mask.shape
    mae = np.abs(input_arr - target)
    mae = mae[mask].mean()
    return mae


def analyze_std_h5(std_file):
    """
    分析 hamiltonians_std.h5 文件，计算非零元素的平均标准差

    Args:
        std_file: hamiltonians_std.h5 文件路径

    Returns:
        dict: 包含统计信息的字典，若文件不存在返回 None
    """
    if not os.path.exists(std_file):
        return None

    try:
        total_sum = 0.0
        total_count = 0
        block_count = 0

        with h5py.File(std_file, 'r') as f:
            for key in f.keys():
                data = np.array(f[key])
                # 获取非零元素
                non_zero = data[data != 0]
                if len(non_zero) > 0:
                    total_sum += non_zero.sum()
                    total_count += len(non_zero)
                    block_count += 1

        if total_count == 0:
            return {"avg_std": 0.0, "total_nonzero": 0, "block_count": block_count}

        avg_std = total_sum / total_count
        return {
            "avg_std": float(avg_std),
            "total_nonzero": int(total_count),
            "block_count": int(block_count)
        }

    except Exception as e:
        print(f'[!] 解析 {std_file} 失败: {e}')
        return None


def read_structure_json(json_file, d_cc_0, d_vdw_0):
    """
    读取 structure.json 文件，计算 d_cc, d_vdw 变化率及应变信息

    Args:
        json_file: structure.json 文件路径
        d_cc_0: d_cc 的参考默认值
        d_vdw_0: d_vdw 的参考默认值

    Returns:
        dict: 包含结构变化信息的字典，若文件不存在返回 None
    """
    if not os.path.exists(json_file):
        return None

    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        result = {}

        # 获取 metadata 字段（可能在根级别或嵌套在 metadata 下）
        metadata = data.get('metadata', data)

        # 计算 d_cc 变化率
        if 'd_cc' in metadata and d_cc_0 is not None and d_cc_0 != 0:
            result['delta_cc'] = (metadata['d_cc'] - d_cc_0) / d_cc_0
        else:
            result['delta_cc'] = None

        # 计算 d_vdw 变化率
        if 'd_vdw' in metadata and d_vdw_0 is not None and d_vdw_0 != 0:
            result['delta_vdw'] = (metadata['d_vdw'] - d_vdw_0) / d_vdw_0
        else:
            result['delta_vdw'] = None

        # 读取应变信息
        if 'strain' in metadata:
            strain = metadata['strain']
            result['eps_xx'] = strain.get('eps_xx', None)
            result['eps_yy'] = strain.get('eps_yy', None)
            # 计算应变迹 (trace)，保留正负信息表示拉伸/压缩
            if result['eps_xx'] is not None and result['eps_yy'] is not None:
                result['delta_strain'] = result['eps_xx'] + result['eps_yy']
            else:
                result['delta_strain'] = None
        else:
            result['eps_xx'] = None
            result['eps_yy'] = None
            result['delta_strain'] = None

        return result

    except Exception as e:
        print(f'[!] 解析 {json_file} 失败: {e}')
        return None


def select_structures(all_structures, include, exclude):
    """
    根据正则表达式筛选结构

    Args:
        all_structures: 所有结构名称列表
        include: 包含的正则表达式列表
        exclude: 排除的正则表达式列表

    Returns:
        list: 筛选后的结构名称列表
    """
    assert not (len(include) > 0 and len(exclude) > 0), "不能同时使用 include 和 exclude"

    included = []
    if len(include) > 0:
        for s in all_structures:
            for inc in include:
                if re.match(inc, s):
                    included.append(s)
                    break
    else:
        included = list(all_structures)
        if len(exclude) > 0:
            for ix in range(len(included) - 1, -1, -1):
                for exc in exclude:
                    if re.match(exc, included[ix]):
                        included.pop(ix)
                        break

    return included


def analyze_errors(h5file, include=None, exclude=None, verbose=True, d_cc_0=None, d_vdw_0=None, struct_dir=None):
    """
    分析 test_result.h5 文件中的预测误差

    Args:
        h5file: test_result.h5 文件路径
        include: 包含的结构名称正则表达式列表
        exclude: 排除的结构名称正则表达式列表
        verbose: 是否打印详细信息
        d_cc_0: d_cc 的参考默认值，用于计算变化率
        d_vdw_0: d_vdw 的参考默认值，用于计算变化率
        struct_dir: structure.json 所在的基础目录，若为 None 则在 h5file 同目录查找

    Returns:
        dict: 包含各结构和总体误差统计的字典
    """
    include = include or []
    exclude = exclude or []

    results = {
        "structures": [],
        "summary": {
            "total_mse": 0.0,
            "total_mae": 0.0,
            "num_structures": 0,
            "avg_mse": 0.0,
            "avg_mae": 0.0
        }
    }

    # 获取 test_result.h5 所在目录
    h5_dir = os.path.dirname(h5file)

    with h5py.File(h5file, 'r') as f:
        all_structures = list(f.keys())
        included = select_structures(all_structures, include, exclude)

        if len(included) == 0:
            print('[!] 未找到符合条件的结构')
            return results

        if verbose:
            print(f'[*] 共发现 {len(all_structures)} 个结构，筛选后保留 {len(included)} 个')
            print('-' * 80)

        total_mse = 0.0
        total_mae = 0.0
        total_avg_std = 0.0
        std_count = 0

        for s in included:
            g = f[s]
            H_pred = np.array(g['H_pred'])
            label = np.array(g['label'])
            mask = np.array(g['mask'])

            mse = maskmse(H_pred, label, mask)
            mae = maskmae(H_pred, label, mask)

            total_mse += mse
            total_mae += mae

            structure_info = {
                "name": s,
                "mse": float(mse),
                "mae": float(mae),
                "avg_std": None,
                "delta_cc": None,
                "delta_vdw": None,
                "eps_xx": None,
                "eps_yy": None,
                "delta_strain": None
            }

            # 尝试解析对应的 hamiltonians_std.h5
            std_file = os.path.join(h5_dir, s, 'hamiltonians_std.h5')
            std_result = analyze_std_h5(std_file)
            if std_result is not None:
                structure_info["avg_std"] = std_result["avg_std"]
                structure_info["std_nonzero"] = std_result["total_nonzero"]
                structure_info["std_blocks"] = std_result["block_count"]
                total_avg_std += std_result["avg_std"]
                std_count += 1

            # 尝试解析对应的 structure.json
            if struct_dir is not None:
                json_file = os.path.join(struct_dir, s, 'structure.json')
            else:
                json_file = os.path.join(h5_dir, s, 'structure.json')
            struct_result = read_structure_json(json_file, d_cc_0, d_vdw_0)
            if struct_result is not None:
                structure_info["delta_cc"] = struct_result["delta_cc"]
                structure_info["delta_vdw"] = struct_result["delta_vdw"]
                structure_info["eps_xx"] = struct_result["eps_xx"]
                structure_info["eps_yy"] = struct_result["eps_yy"]
                structure_info["delta_strain"] = struct_result["delta_strain"]

            results["structures"].append(structure_info)

            if verbose:
                std_str = f"  STD: {structure_info['avg_std']:.6e}" if structure_info["avg_std"] is not None else ""
                delta_str = ""
                if structure_info["delta_cc"] is not None:
                    delta_str += f"  Δcc: {structure_info['delta_cc']:.4f}"
                if structure_info["delta_vdw"] is not None:
                    delta_str += f"  Δvdw: {structure_info['delta_vdw']:.4f}"
                if structure_info["delta_strain"] is not None:
                    delta_str += f"  Δstrain: {structure_info['delta_strain']:.4f}"
                print(f'  {s:40s}  MSE: {mse:.6e}  MAE: {mae:.6e}{std_str}{delta_str}')

        n = len(included)
        avg_mse = total_mse / n
        avg_mae = total_mae / n

        results["summary"]["total_mse"] = float(total_mse)
        results["summary"]["total_mae"] = float(total_mae)
        results["summary"]["num_structures"] = n
        results["summary"]["avg_mse"] = float(avg_mse)
        results["summary"]["avg_mae"] = float(avg_mae)

        # 添加 std 统计
        if std_count > 0:
            results["summary"]["avg_std"] = float(total_avg_std / std_count)
            results["summary"]["std_structures"] = std_count
        else:
            results["summary"]["avg_std"] = None
            results["summary"]["std_structures"] = 0

        if verbose:
            print('-' * 80)
            print(f'[*] 统计结果 ({n} 个结构):')
            print(f'    平均 MSE: {avg_mse:.6e}')
            print(f'    平均 MAE: {avg_mae:.6e}')
            if std_count > 0:
                print(f'    平均 STD: {results["summary"]["avg_std"]:.6e} ({std_count} 个结构有 std 数据)')

    return results


def save_to_csv(results, output_path):
    """
    将结果保存到 CSV 文件

    Args:
        results: analyze_errors 返回的结果字典
        output_path: 输出文件路径
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # 检查是否有 std 数据和结构信息
        has_std = any(s.get("avg_std") is not None for s in results["structures"])
        has_structure = any(s.get("delta_cc") is not None or s.get("delta_strain") is not None
                           for s in results["structures"])

        # 写入表头
        header = ['结构名称', 'MSE', 'MAE']
        if has_std:
            header.append('Avg_STD')
        if has_structure:
            header.extend(['Delta_cc', 'Delta_vdw', 'Delta_strain'])
        writer.writerow(header)

        # 写入各结构数据
        for s in results["structures"]:
            row = [s["name"], f'{s["mse"]:.10e}', f'{s["mae"]:.10e}']
            if has_std:
                if s.get("avg_std") is not None:
                    row.append(f'{s["avg_std"]:.10e}')
                else:
                    row.append('')
            if has_structure:
                # Delta_cc
                if s.get("delta_cc") is not None:
                    row.append(f'{s["delta_cc"]:.6f}')
                else:
                    row.append('')
                # Delta_vdw
                if s.get("delta_vdw") is not None:
                    row.append(f'{s["delta_vdw"]:.6f}')
                else:
                    row.append('')
                # Delta_strain
                if s.get("delta_strain") is not None:
                    row.append(f'{s["delta_strain"]:.6f}')
                else:
                    row.append('')
            writer.writerow(row)

        # 写入汇总行
        writer.writerow([])
        summary_header = ['统计项', 'MSE', 'MAE']
        if has_std:
            summary_header.append('Avg_STD')
        if has_structure:
            summary_header.extend(['Delta_cc', 'Delta_vdw', 'Delta_strain'])
        writer.writerow(summary_header)

        avg_row = ['平均值', f'{results["summary"]["avg_mse"]:.10e}', f'{results["summary"]["avg_mae"]:.10e}']
        if has_std:
            avg_std_str = f'{results["summary"]["avg_std"]:.10e}' if results["summary"]["avg_std"] is not None else ''
            avg_row.append(avg_std_str)
        if has_structure:
            avg_row.extend(['', '', ''])  # Delta 列不需要平均值统计
        writer.writerow(avg_row)

        total_row = ['总和', f'{results["summary"]["total_mse"]:.10e}', f'{results["summary"]["total_mae"]:.10e}']
        if has_std:
            total_row.append('')
        if has_structure:
            total_row.extend(['', '', ''])
        writer.writerow(total_row)

        count_row = ['结构数', results["summary"]["num_structures"], '', '']
        if has_std:
            count_row.append('')
        if has_structure:
            count_row.extend(['', '', ''])
        writer.writerow(count_row)

    print(f'[+] 结果已保存至: {output_path.absolute()}')


def save_to_json(results, output_path):
    """
    将结果保存到 JSON 文件

    Args:
        results: analyze_errors 返回的结果字典
        output_path: 输出文件路径
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f'[+] 结果已保存至: {output_path.absolute()}')


def main():
    parser = argparse.ArgumentParser(
        description='MSE/MAE Analysis Tool - Extract prediction error statistics from test_result.h5',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
example:
  %(prog)s -i ./test_result.h5
  %(prog)s -i ./test_result.h5 -o results.csv
  %(prog)s -i ./test_result.h5 --include "tbg_.*" --exclude ""
  %(prog)s -i ./test_result.h5 --json results.json
  %(prog)s -i ./test_result.h5 --d-cc 1.42 --d-vdw 3.35  # Calculate structure change ratio
        '''
    )

    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to test_result.h5 file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output CSV file path (default: not saved)')
    parser.add_argument('--json', type=str, default=None,
                        help='Output JSON file path (default: not saved)')
    parser.add_argument('--include', type=str, nargs='+', default=[],
                        help='Regex for structure names to include (e.g., "tbg_.*")')
    parser.add_argument('--exclude', type=str, nargs='+', default=[],
                        help='Regex for structure names to exclude')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Quiet mode, suppress detailed output')
    parser.add_argument('--d-cc', type=float, default=None,
                        help='Reference default value for d_cc, used to calculate Delta_cc change ratio')
    parser.add_argument('--d-vdw', type=float, default=None,
                        help='Reference default value for d_vdw, used to calculate Delta_vdw change ratio')
    parser.add_argument('--struct-dir', type=str, default=None,
                        help='Base directory containing structure.json (default: look in same directory as test_result.h5)')

    args = parser.parse_args()

    # 检查输入文件
    if not os.path.exists(args.input):
        print(f'[-] Error: File not found - {args.input}')
        return 1

    print(f'[*] Analyzing file: {args.input}')
    print(f'[*] Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

    # 显示参考值设置
    if args.d_cc is not None:
        print(f'[*] d_cc reference value: {args.d_cc}')
    if args.d_vdw is not None:
        print(f'[*] d_vdw reference value: {args.d_vdw}')
    if args.struct_dir is not None:
        print(f'[*] Structure info directory: {args.struct_dir}')

    # 执行分析
    results = analyze_errors(
        h5file=args.input,
        include=args.include,
        exclude=args.exclude,
        verbose=not args.quiet,
        d_cc_0=args.d_cc,
        d_vdw_0=args.d_vdw,
        struct_dir=args.struct_dir
    )

    # 保存结果
    if args.output:
        save_to_csv(results, args.output)

    if args.json:
        save_to_json(results, args.json)

    return 0


if __name__ == '__main__':
    exit(main())