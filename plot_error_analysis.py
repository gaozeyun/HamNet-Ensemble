#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
误差分析可视化脚本
读取 error_analyze.csv，绘制 MSE, MAE, Avg_STD 随 Delta 变化的图表
"""

import os
import csv
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def plot_error_analysis(csv_file, output_file, x_delta='Delta_cc', figsize=(24, 6), dpi=600):
    """
    绘制误差分析图表

    Args:
        csv_file: error_analyze.csv 文件路径
        output_file: 输出图片路径
        x_delta: 作为 x 轴的 Delta 列名 (Delta_cc, Delta_vdw, Delta_strain)
        figsize: 图片尺寸 (width, height)
        dpi: 图片分辨率
    """
    # 读取 CSV 文件
    data = []
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 过滤掉非数据行
            if row['结构名称'] in ['统计项', '平均值', '总和', '结构数', '']:
                continue
            try:
                data.append({
                    'name': row['结构名称'],
                    'MSE': float(row['MSE']),
                    'MAE': float(row['MAE']),
                    'Avg_STD': float(row['Avg_STD']) if row.get('Avg_STD') else 0.0,
                    'Delta_cc': float(row['Delta_cc']) if row.get('Delta_cc') else 0.0,
                    'Delta_vdw': float(row['Delta_vdw']) if row.get('Delta_vdw') else 0.0,
                    'Delta_strain': float(row['Delta_strain']) if row.get('Delta_strain') else 0.0
                })
            except (ValueError, KeyError):
                continue

    if len(data) == 0:
        print('[!] 错误: 没有有效数据')
        return

    # 按 x 轴排序
    data.sort(key=lambda d: d[x_delta])

    # 提取数据
    x = np.array([d[x_delta] for d in data]) * 100  # 转换为百分比 (%)
    mae = np.array([d['MAE'] for d in data]) * 1000  # 转换为 meV
    mse = np.array([d['MSE'] for d in data]) * 1e6  # 转换为 meV^2
    avg_std = np.array([d['Avg_STD'] for d in data]) * 1000  # 转换为 meV

    # 创建图表
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)

    # 字体大小
    fontsize = {
        'title': 28,
        'label': 24,
        'tick': 22,
        'legend': 22
    }

    # 线宽和标记大小
    line_width = 2.0
    marker_size = 8
    bar_width = 0.015 * (x.max() - x.min()) if len(x) > 1 else 0.01

    # x 轴位置（用于柱状图）
    x_bar = np.arange(len(x))

    # ========================================
    # 左轴: MAE 和 Avg_STD (单位: meV)
    # ========================================
    # 图层1: MAE 柱状图 (底层)
    bar_color = '#a8d5ba'  # 淡绿色
    bar_edge_color = '#2e7d32'  # 深绿色边框
    ax1.bar(x, mae, width=bar_width, color=bar_color, edgecolor=bar_edge_color,
            linewidth=1.0, alpha=0.7, label='MAE', zorder=1)

    # 图层2: Avg_STD 点线图 (顶层)
    std_color = '#d32f2f'  # 红色
    ax1.plot(x, avg_std, 's--', color=std_color, linewidth=line_width,
             markersize=marker_size, markerfacecolor='white',
             markeredgewidth=2.0, label=r'$\Delta\sigma$', zorder=3)

    # 左轴设置
    delta_label = x_delta.replace('Delta_', 'Δ')
    ax1.set_xlabel(f'$\\mathit{{{delta_label}}}$ (%)', fontsize=fontsize['label'])
    ax1.set_ylabel(r'MAE, $\Delta\sigma$ (meV)', fontsize=fontsize['label'])
    ax1.tick_params(axis='y', labelsize=fontsize['tick'],
                    direction='in', length=6, width=1.5)

    # ========================================
    # 右轴: MSE (单位: meV^2)
    # ========================================
    ax2 = ax1.twinx()
    mse_color = '#1976d2'  # 蓝色
    ax2.plot(x, mse, 'o-', color=mse_color, linewidth=line_width,
             markersize=marker_size, markerfacecolor='white',
             markeredgewidth=2.0, label='MSE', zorder=2)
    ax2.set_ylabel('MSE (meV$^2$)', fontsize=fontsize['label'])
    ax2.tick_params(axis='y', labelsize=fontsize['tick'],
                    direction='in', length=6, width=1.5)

    # ========================================
    # 共用坐标轴设置
    # ========================================
    # 设置 x 轴刻度
    ax1.tick_params(axis='x', which='major', labelsize=fontsize['tick'],
                    direction='in', length=6, width=1.5)

    # 启用所有边框
    for spine in ['top', 'bottom', 'left', 'right']:
        ax1.spines[spine].set_visible(True)
        ax1.spines[spine].set_linewidth(1.5)

    # 启用上边框刻度 (x轴)
    ax1.tick_params(axis='x', which='both',
                    top=True, bottom=True, direction='in', length=6, width=1.5)

    # y 轴科学计数法
    ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(-2, 3))
    ax2.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax2.ticklabel_format(axis='y', style='scientific', scilimits=(-2, 3))

    # ========================================
    # 图例 (横向放置于图片下方，无边框)
    # ========================================
    # 合并两个轴的图例
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc='upper center',
               bbox_to_anchor=(0.5, -0.12),
               ncol=3,
               fontsize=fontsize['legend'],
               frameon=False)

    # ========================================
    # 网格线
    # ========================================
    ax1.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=1.0, zorder=0)
    ax1.xaxis.grid(False)

    # ========================================
    # 保存图片
    # ========================================
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"[*] 图表已保存至: {output_file}")
    print(f"[*] 数据点数: {len(x)}")
    print(f"[*] x轴范围: [{x.min():.2f}%, {x.max():.2f}%]")


def main():
    parser = argparse.ArgumentParser(
        description='误差分析可视化工具 - 从 error_analyze.csv 绘制误差变化图表',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  %(prog)s -i ./error_analyze.csv -o error_plot.png
  %(prog)s -i ./error_analyze.csv -o error_plot.png --x-delta Delta_cc
  %(prog)s -i ./error_analyze.csv -o error_plot.png --x-delta Delta_strain --figsize 12,6
        '''
    )

    parser.add_argument('-i', '--input', type=str, required=True,
                        help='error_analyze.csv 文件路径')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='输出图片路径')
    parser.add_argument('--x-delta', type=str, default='Delta_cc',
                        choices=['Delta_cc', 'Delta_vdw', 'Delta_strain'],
                        help='作为 x 轴的 Delta 列名 (默认: Delta_cc)')
    parser.add_argument('--figsize', type=str, default='24,6',
                        help='图片尺寸，格式: width,height (默认: 24,6)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='图片分辨率 (默认: 150)')

    args = parser.parse_args()

    # 解析 figsize
    try:
        figsize = tuple(map(float, args.figsize.split(',')))
        if len(figsize) != 2:
            raise ValueError
    except ValueError:
        print(f'[!] 错误: figsize 格式无效，应为 "width,height"')
        return 1

    # 检查输入文件
    if not os.path.exists(args.input):
        print(f'[!] 错误: 文件不存在 - {args.input}')
        return 1

    # 绘图
    plot_error_analysis(
        csv_file=args.input,
        output_file=args.output,
        x_delta=args.x_delta,
        figsize=figsize,
        dpi=args.dpi
    )

    return 0


if __name__ == '__main__':
    exit(main())