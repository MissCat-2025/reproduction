#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOOSE运行时间分析脚本
从所有案例的run.log文件中提取时间信息，汇总到CSV文件
"""

import os
import re
import csv
import glob
from collections import defaultdict
from datetime import datetime


target_times = [1000, 30000, 50000, 100000, 200000, 500000, 1000000, 1500000,99999999]





def extract_time_info_from_log(log_path, target_times):
    """
    从单个run.log文件中提取时间信息
    
    Args:
        log_path: log文件路径
        target_times: 目标时间列表 [1000, 30000, 50000, ...]
    
    Returns:
        dict: 包含时间信息的字典
    """
    time_info = {
        'case_name': '',
        'total_runtime': 0,
        'convergence_status': 'unknown',
        'return_code': -1,
        'time_steps': [],
        'wall_times': [],
        'errors': []
    }
    
    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            
        # 提取案例名称（从目录名）
        case_dir = os.path.dirname(log_path)
        time_info['case_name'] = os.path.basename(case_dir)
        
        # 提取总运行时间
        runtime_match = re.search(r'耗时:\s*([\d\.]+)s', content)
        if runtime_match:
            time_info['total_runtime'] = float(runtime_match.group(1))
        
        # 提取返回码
        rc_match = re.search(r'返回码:\s*(\d+)', content)
        if rc_match:
            time_info['return_code'] = int(rc_match.group(1))
            time_info['convergence_status'] = 'converged' if rc_match.group(1) == '0' else 'failed'
        
        # 提取所有时间步信息
        time_step_pattern = r'Time Step \d+, time = ([\d\.e+-]+)'
        time_steps = re.findall(time_step_pattern, content)
        time_info['time_steps'] = [float(t) for t in time_steps]
        
        # 提取wall time信息（如果存在）
        wall_time_pattern = r'Wall time:\s*([\d\.]+)s'
        wall_times = re.findall(wall_time_pattern, content)
        time_info['wall_times'] = [float(t) for t in wall_times]
        
        # 提取错误信息
        errors = []
        if 'Solve Did NOT Converge!' in content:
            errors.append('convergence_failed')
        if '*** ERROR ***' in content:
            errors.append('runtime_error')
        if 'dtmin' in content.lower():
            errors.append('dt_min_violation')
        if 'max steps' in content.lower():
            errors.append('max_steps_exceeded')
        time_info['errors'] = errors
        
        # 计算每个目标时间点的用时
        time_info['target_time_usage'] = {}
        for target_time in target_times:
            # 找到最接近目标时间的时间步
            if time_info['time_steps']:
                closest_time = min(time_info['time_steps'], 
                                 key=lambda x: abs(x - target_time))
                # 计算到该时间点的累计用时（线性插值）
                if len(time_info['wall_times']) > 0:
                    # 如果有wall time信息，使用wall time
                    if len(time_info['wall_times']) >= len(time_info['time_steps']):
                        # 找到对应时间步的wall time
                        step_index = time_info['time_steps'].index(closest_time)
                        if step_index < len(time_info['wall_times']):
                            time_info['target_time_usage'][target_time] = time_info['wall_times'][step_index]
                        else:
                            # 线性插值
                            time_info['target_time_usage'][target_time] = time_info['total_runtime'] * (closest_time / max(time_info['time_steps']))
                    else:
                        # 线性插值
                        time_info['target_time_usage'][target_time] = time_info['total_runtime'] * (closest_time / max(time_info['time_steps']))
                else:
                    # 没有wall time信息，使用总运行时间按比例估算
                    time_info['target_time_usage'][target_time] = time_info['total_runtime'] * (closest_time / max(time_info['time_steps']))
            else:
                time_info['target_time_usage'][target_time] = 0
        
    except Exception as e:
        print(f"警告：无法读取日志文件 {log_path}: {str(e)}")
        time_info['errors'].append(f'read_error: {str(e)}')
    
    return time_info

def analyze_all_cases(studies_dir, target_times):
    """
    分析所有案例的日志文件
    
    Args:
        studies_dir: 参数研究目录
        target_times: 目标时间列表
    
    Returns:
        list: 所有案例的时间信息列表
    """
    all_cases = []
    
    # 查找所有case目录
    case_dirs = glob.glob(os.path.join(studies_dir, "case_*"))
    
    print(f"找到 {len(case_dirs)} 个案例目录")
    
    for case_dir in case_dirs:
        # 查找run.log文件
        log_path = os.path.join(case_dir, 'run.log')
        
        if os.path.exists(log_path):
            print(f"分析案例: {os.path.basename(case_dir)}")
            time_info = extract_time_info_from_log(log_path, target_times)
            all_cases.append(time_info)
        else:
            print(f"警告：案例 {os.path.basename(case_dir)} 没有run.log文件")
    
    # 按案例名称排序（case_001, case_002, case_003...）
    def natural_sort_key(case):
        # 提取案例编号进行排序
        case_name = case['case_name']
        match = re.search(r'case_(\d+)', case_name)
        if match:
            return int(match.group(1))
        return float('inf')  # 如果没有编号，放到最后
    
    all_cases.sort(key=natural_sort_key)
    
    return all_cases

def generate_csv_report(all_cases, target_times, output_path):
    """
    生成CSV报告
    
    Args:
        all_cases: 所有案例的时间信息
        target_times: 目标时间列表
        output_path: 输出CSV文件路径
    """
    # 定义CSV列
    headers = ['case_name', 'convergence_status', 'total_runtime', 'return_code', 'errors']
    
    # 添加目标时间列
    for target_time in target_times:
        headers.append(f'time_{target_time}s')
    
    # 添加其他有用的时间信息列
    headers.extend(['max_time_reached', 'time_steps_count', 'wall_times_count'])
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        
        for case in all_cases:
            row = {
                'case_name': case['case_name'],
                'convergence_status': case['convergence_status'],
                'total_runtime': f"{case['total_runtime']:.2f}",
                'return_code': case['return_code'],
                'errors': ';'.join(case['errors']) if case['errors'] else 'none',
                'max_time_reached': f"{max(case['time_steps']):.2e}" if case['time_steps'] else '0',
                'time_steps_count': len(case['time_steps']),
                'wall_times_count': len(case['wall_times'])
            }
            
            # 添加每个目标时间的用时
            for target_time in target_times:
                usage = case['target_time_usage'].get(target_time, 0)
                row[f'time_{target_time}s'] = f"{usage:.2f}"
            
            writer.writerow(row)
    
    print(f"CSV报告已生成: {output_path}")

def generate_summary_report(all_cases, target_times, output_dir):
    """
    生成汇总报告
    
    Args:
        all_cases: 所有案例的时间信息
        target_times: 目标时间列表
        output_dir: 输出目录
    """
    # 统计信息
    total_cases = len(all_cases)
    converged_cases = sum(1 for case in all_cases if case['convergence_status'] == 'converged')
    failed_cases = total_cases - converged_cases
    
    # 计算平均运行时间
    runtimes = [case['total_runtime'] for case in all_cases if case['total_runtime'] > 0]
    avg_runtime = sum(runtimes) / len(runtimes) if runtimes else 0
    
    # 生成汇总报告
    summary_path = os.path.join(output_dir, 'time_analysis_summary.txt')
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write("MOOSE运行时间分析汇总报告\n")
        f.write("=" * 50 + "\n")
        f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"总案例数: {total_cases}\n")
        f.write(f"收敛案例数: {converged_cases}\n")
        f.write(f"失败案例数: {failed_cases}\n")
        f.write(f"成功率: {converged_cases/total_cases*100:.1f}%\n")
        f.write(f"平均运行时间: {avg_runtime:.2f}秒\n")
        f.write(f"目标时间点: {target_times}\n\n")
        
        # 显示案例列表（按顺序）
        f.write("案例列表（按编号排序）:\n")
        for case in all_cases:
            status_icon = "✓" if case['convergence_status'] == 'converged' else "✗"
            f.write(f"  {status_icon} {case['case_name']} - {case['convergence_status']} ({case['total_runtime']:.1f}s)\n")
        f.write("\n")
        
        # 每个目标时间点的统计
        for target_time in target_times:
            usages = [case['target_time_usage'].get(target_time, 0) for case in all_cases]
            valid_usages = [u for u in usages if u > 0]
            if valid_usages:
                avg_usage = sum(valid_usages) / len(valid_usages)
                f.write(f"时间 {target_time}s 平均用时: {avg_usage:.2f}秒\n")
        
        f.write("\n错误统计:\n")
        error_counts = defaultdict(int)
        for case in all_cases:
            for error in case['errors']:
                error_counts[error] += 1
        
        for error, count in sorted(error_counts.items()):
            f.write(f"  {error}: {count}次\n")
    
    print(f"汇总报告已生成: {summary_path}")

def main():
    """主函数"""
    # 配置参数
    script_dir = os.path.dirname(os.path.abspath(__file__))
    studies_dir = os.path.join(script_dir, 'parameter_studies_series')
    
    # 目标时间列表（秒）
    # target_times = [1000, 30000, 50000, 100000, 200000, 500000, 1000000, 1500000]
    
    # 检查目录是否存在
    if not os.path.exists(studies_dir):
        print(f"错误：参数研究目录不存在: {studies_dir}")
        return
    
    print(f"开始分析MOOSE运行时间...")
    print(f"目标时间点: {target_times}")
    print(f"研究目录: {studies_dir}")
    
    # 分析所有案例
    all_cases = analyze_all_cases(studies_dir, target_times)
    
    if not all_cases:
        print("未找到任何案例数据")
        return
    
    # 创建输出目录
    output_dir = os.path.join(studies_dir, 'time_analysis')
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成CSV报告
    csv_path = os.path.join(output_dir, 'runtime_analysis.csv')
    generate_csv_report(all_cases, target_times, csv_path)
    
    # 生成汇总报告
    generate_summary_report(all_cases, target_times, output_dir)
    
    print(f"\n分析完成！")
    print(f"共分析了 {len(all_cases)} 个案例")
    print(f"输出文件:")
    print(f"  - CSV报告: {csv_path}")
    print(f"  - 汇总报告: {os.path.join(output_dir, 'time_analysis_summary.txt')}")

if __name__ == '__main__':
    main()
