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
import sys
from collections import defaultdict
from datetime import datetime


BASE_DIR = os.environ.get("PROJECT_BASE_DIR", os.path.dirname(os.path.abspath(__file__)))
TARGET_STUDIES_DIR = os.environ.get("TARGET_STUDIES_DIR", os.path.join(BASE_DIR, "2Good-Gc2.25-3.5"))






TITLE_FONT_SIZE = int(os.environ.get("TIME_TITLE_FONT_SIZE", "72"))
LABEL_FONT_SIZE = int(os.environ.get("TIME_LABEL_FONT_SIZE", "28"))
TIME_FONT_SIZE = int(os.environ.get("TIME_AXIS_FONT_SIZE", "24"))
LABEL_WIDTH = int(os.environ.get("TIME_LABEL_WIDTH", "220"))
TITLE_HEIGHT = int(os.environ.get("TIME_TITLE_HEIGHT", "80"))
TIME_LABEL_HEIGHT = int(os.environ.get("TIME_AXIS_HEIGHT", "100"))
ROW_TIME_LABEL_HEIGHT = int(os.environ.get("TIME_ROW_AXIS_HEIGHT", "100"))

_env_time_picture_times = os.environ.get("TIME_PICTURE_TARGET_TIMES") or os.environ.get("TARGET_TIMES")
if _env_time_picture_times:
    try:
        target_times = [
            float(x) for x in _env_time_picture_times.split(",") if x.strip()
        ]
    except ValueError:
        target_times = [80000, 125000, 300000]
else:
    target_times = [80000, 125000, 300000]

def parse_image_time_from_name(filename):
    match = re.search(r'_(\d+p\d+)s(?=\.png$)', filename)
    if match:
        return float(match.group(1).replace('p', '.'))
    return None

def natural_case_sort_key(case_path):
    case_name = os.path.basename(case_path)
    match = re.search(r'case_(\d+)', case_name)
    if match:
        return int(match.group(1))
    return 10**9

def stitch_images_horizontally(image_paths):
    try:
        from PIL import Image
    except ImportError as exc:
        raise RuntimeError("需要安装 Pillow 才能拼接图片") from exc
    images = []
    widths = []
    heights = []
    times = []
    centers = []
    try:
        for path in image_paths:
            img = Image.open(path).convert('RGB')
            images.append(img)
            widths.append(img.width)
            heights.append(img.height)
            t = parse_image_time_from_name(os.path.basename(path))
            times.append(t)
        if not images:
            return None, [], []
        total_width = sum(widths)
        max_height = max(heights)
        canvas = Image.new('RGB', (total_width, max_height), (255, 255, 255))
        x_offset = 0
        for img, w in zip(images, widths):
            canvas.paste(img, (x_offset, 0))
            centers.append(x_offset + w // 2)
            x_offset += w
        return canvas, times, centers
    finally:
        for img in images:
            try:
                img.close()
            except Exception:
                pass

def stitch_images_vertically(images):
    try:
        from PIL import Image
    except ImportError as exc:
        raise RuntimeError("需要安装 Pillow 才能拼接图片") from exc
    if not images:
        return None
    max_width = max(img.width for img in images)
    total_height = sum(img.height for img in images)
    canvas = Image.new('RGB', (max_width, total_height), (255, 255, 255))
    y_offset = 0
    for img in images:
        canvas.paste(img, (0, y_offset))
        y_offset += img.height
    return canvas

def build_row_image_for_case(global_times, time_to_path, base_size):
    try:
        from PIL import Image
    except ImportError as exc:
        raise RuntimeError("需要安装 Pillow 才能拼接图片") from exc
    w, h = base_size
    row_width = w * len(global_times)
    row_image = Image.new("RGB", (row_width, h), (255, 255, 255))
    x_offset = 0
    for t in global_times:
        path = time_to_path.get(t)
        if path is not None:
            img = Image.open(path).convert("RGB")
            row_image.paste(img, (x_offset, 0))
            img.close()
        x_offset += w
    return row_image
def create_labeled_grid(row_entries, project_name, label_width=LABEL_WIDTH, title_height=TITLE_HEIGHT, time_label_height=TIME_LABEL_HEIGHT, row_time_label_height=ROW_TIME_LABEL_HEIGHT, time_centers=None, time_labels=None):
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise RuntimeError("需要安装 Pillow 才能拼接图片") from exc
    if not row_entries:
        return None
    def build_font(size):
        try:
            return ImageFont.truetype("DejaVuSans.ttf", size)
        except Exception:
            return ImageFont.load_default()
    font_title = build_font(TITLE_FONT_SIZE)
    font_label = build_font(LABEL_FONT_SIZE)
    font_time = build_font(TIME_FONT_SIZE)
    max_row_width = max(img.width for _, img in row_entries)
    if time_centers and time_labels:
        total_rows_height = sum(img.height + row_time_label_height for _, img in row_entries)
    else:
        total_rows_height = sum(img.height for _, img in row_entries)
    width = label_width + max_row_width
    extra_time_band = time_label_height if time_centers and time_labels else 0
    height = title_height + extra_time_band + total_rows_height
    grid = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(grid)
    def measure(text, font):
        try:
            return draw.textbbox((0, 0), text, font=font)[2:]
        except Exception:
            return draw.textsize(text, font=font)
    title_text = project_name
    title_size = measure(title_text, font_title)
    title_x = (width - title_size[0]) // 2
    title_y = (title_height - title_size[1]) // 2
    draw.text((title_x, title_y), title_text, fill=(0, 0, 0), font=font_title)
    y_offset = title_height
    if time_centers and time_labels:
        time_y_center = y_offset + time_label_height // 2
        for x, label in zip(time_centers, time_labels):
            if label is None or label == "":
                continue
            text = label
            size = measure(text, font_time)
            text_x = label_width + x - size[0] // 2
            text_y = time_y_center - size[1] // 2
            draw.text((text_x, text_y), text, fill=(0, 0, 0), font=font_time)
        y_offset += time_label_height
    for case_name, row_image in row_entries:
        grid.paste(row_image, (label_width, y_offset))
        case_text = case_name
        case_size = measure(case_text, font_label)
        text_x = max(0, label_width - 10 - case_size[0])
        text_y = y_offset + (row_image.height - case_size[1]) // 2
        draw.text((text_x, text_y), case_text, fill=(0, 0, 0), font=font_label)
        if time_centers and time_labels:
            row_time_y_center = y_offset + row_image.height + row_time_label_height // 2
            for x, label in zip(time_centers, time_labels):
                if label is None or label == "":
                    continue
                text = label
                size = measure(text, font_time)
                text_x = label_width + x - size[0] // 2
                text_y = row_time_y_center - size[1] // 2
                draw.text((text_x, text_y), text, fill=(0, 0, 0), font=font_time)
            y_offset += row_time_label_height
        y_offset += row_image.height
    return grid

def build_case_time_grid(studies_dir, output_prefix="combined"):
    if not os.path.isdir(studies_dir):
        print(f"错误：目录不存在 {studies_dir}")
        return []
    case_dirs = [
        os.path.join(studies_dir, name)
        for name in os.listdir(studies_dir)
        if name.startswith('case_') and os.path.isdir(os.path.join(studies_dir, name))
    ]
    case_dirs.sort(key=natural_case_sort_key)
    if not case_dirs:
        print("未找到case目录")
        return []
    project_names = set()
    for case_dir in case_dirs:
        post_dir = os.path.join(case_dir, 'post_results')
        if not os.path.isdir(post_dir):
            continue
        for name in os.listdir(post_dir):
            path = os.path.join(post_dir, name)
            if os.path.isdir(path):
                project_names.add(name)
    project_names = sorted(project_names)
    if not project_names:
        print("未找到post_results子目录")
        return []
    outputs = []
    for project in project_names:
        per_case_time_paths = {}
        global_times_set = set()
        for case_dir in case_dirs:
            project_dir = os.path.join(case_dir, "post_results", project)
            if not os.path.isdir(project_dir):
                continue
            image_paths = [
                os.path.join(project_dir, name)
                for name in os.listdir(project_dir)
                if name.lower().endswith(".png")
            ]
            if not image_paths:
                continue
            time_to_path = {}
            for path in image_paths:
                t = parse_image_time_from_name(os.path.basename(path))
                if t is None:
                    continue
                time_to_path[t] = path
                global_times_set.add(t)
            if not time_to_path:
                continue
            case_name = os.path.basename(case_dir)
            per_case_time_paths[case_name] = time_to_path
        if not per_case_time_paths or not global_times_set:
            continue
        available_times = sorted(global_times_set)
        if target_times:
            selected = []
            used = set()
            for t in target_times:
                if not available_times:
                    break
                closest = min(available_times, key=lambda x: abs(x - t))
                if closest not in used:
                    used.add(closest)
                    selected.append(closest)
            global_times = selected
        else:
            global_times = available_times
        if not global_times:
            continue
        base_size = None
        for time_paths in per_case_time_paths.values():
            if time_paths:
                sample_path = next(iter(time_paths.values()))
                try:
                    from PIL import Image
                except ImportError as exc:
                    raise RuntimeError("需要安装 Pillow 才能拼接图片") from exc
                with Image.open(sample_path) as img:
                    base_size = (img.width, img.height)
                break
        if base_size is None:
            continue
        row_entries = []
        for case_dir in case_dirs:
            case_name = os.path.basename(case_dir)
            time_to_path = per_case_time_paths.get(case_name)
            if not time_to_path:
                continue
            row_image = build_row_image_for_case(global_times, time_to_path, base_size)
            row_entries.append((case_name, row_image))
        if not row_entries:
            continue
        col_width = base_size[0]
        time_centers = [int(col_width * (i + 0.5)) for i in range(len(global_times))]
        time_labels = [f"{t:.0f}" for t in global_times]
        grid_image = create_labeled_grid(row_entries, project, time_centers=time_centers, time_labels=time_labels)
        if not grid_image:
            for _, img in row_entries:
                try:
                    img.close()
                except Exception:
                    pass
            continue
        output_name = f"{output_prefix}_{project}.png"
        output_path = os.path.join(studies_dir, output_name)
        grid_image.save(output_path)
        outputs.append(output_path)
        for _, img in row_entries:
            try:
                img.close()
            except Exception:
                pass
        try:
            grid_image.close()
        except Exception:
            pass
        print(f"已生成组合图: {output_path}")
    return outputs





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
    if len(sys.argv) > 1:
        studies_dir = sys.argv[1]
    else:
        studies_dir = TARGET_STUDIES_DIR
    
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
    print("开始生成组合图...")
    build_case_time_grid(studies_dir)

if __name__ == '__main__':
    main()
