#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structure Error Visualization Script
Read structure.csv and plot MSE, MAE, Avg_STD bar/line charts
"""

import os
import csv
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# 默认文件路径
DEFAULT_INPUT = r'D:\file\deep_ensemble\download\structure.csv'
DEFAULT_OUTPUT = r'D:\file\deep_ensemble\download\structure_plot.png'


def plot_structure_error(csv_file, output_file, figsize=(20, 6), dpi=300):
    """
    Plot structure error charts

    Args:
        csv_file: Path to structure.csv file
        output_file: Output image path
        figsize: Figure size (width, height)
        dpi: Image resolution
    """
    # 读取 CSV 文件
    data = []
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                data.append({
                    'name': row['Structure Name'],
                    'MSE': float(row['MSE']),
                    'MAE': float(row['MAE']),
                    'Avg_STD': float(row['Avg_STD']) if row.get('Avg_STD') else 0.0
                })
            except (ValueError, KeyError):
                continue

    if len(data) == 0:
        print('[!] Error: No valid data')
        return

    # 提取数据
    x_labels = [d['name'] for d in data]
    x = np.arange(len(x_labels))  # 柱状图位置
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
    bar_width = 0.2  # 柱状图宽度

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
             markeredgewidth=2.0, label=r'$\Delta \sigma$', zorder=3)

    # 左轴设置
    ax1.set_ylabel(r'MAE, $\Delta \sigma$ (meV)', fontsize=fontsize['label'])
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
    # 设置 x 轴刻度和标签
    ax1.set_xticks(x)
    ax1.set_xticklabels(x_labels, fontsize=fontsize['tick'])
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

    print(f"[*] Chart saved to: {output_file}")
    print(f"[*] Data points: {len(x)}")
    print(f"[*] Structure: {', '.join(x_labels)}")


def main():
    parser = argparse.ArgumentParser(
        description='Structure Error Visualization Tool - Plot error charts from structure.csv',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s -i ./structure.csv -o structure_plot.png
  %(prog)s -i ./structure.csv -o structure_plot.png --figsize 12,6
        '''
    )

    parser.add_argument('-i', '--input', type=str, default=DEFAULT_INPUT,
                        help=f'Path to structure.csv file (default: {DEFAULT_INPUT})')
    parser.add_argument('-o', '--output', type=str, default=DEFAULT_OUTPUT,
                        help=f'Output image path (default: {DEFAULT_OUTPUT})')
    parser.add_argument('--figsize', type=str, default='20,6',
                        help='Figure size, format: width,height (default: 20,6)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Image resolution (default: 150)')

    args = parser.parse_args()

    # Parse figsize
    try:
        figsize = tuple(map(float, args.figsize.split(',')))
        if len(figsize) != 2:
            raise ValueError
    except ValueError:
        print(f'[!] Error: Invalid figsize format, should be "width,height"')
        return 1

    # Check input file
    if not os.path.exists(args.input):
        print(f'[!] Error: File not found - {args.input}')
        return 1

    # Generate plot
    plot_structure_error(
        csv_file=args.input,
        output_file=args.output,
        figsize=figsize,
        dpi=args.dpi
    )

    return 0


if __name__ == '__main__':
    exit(main())