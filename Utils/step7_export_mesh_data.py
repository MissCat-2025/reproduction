#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step7：导出目标时刻下的全网格点数据。
- 扫描 case_*，优先依据 file_base 定位 .e 文件
- 读取指定字段在目标时刻的所有网格点数据
- 输出为 CSV，每行表示一个网格点，列包含时刻对应的值
"""

import os
import sys
import re
import csv
import glob

def _parse_target_times():
    v = os.environ.get("DATA_TARGET_TIMES")
    if v:
        return [float(x.strip()) for x in v.split(",") if x.strip()]
    return []

def _parse_fields(env_value):
    if not env_value:
        return []
    fields = []
    for part in env_value.split(","):
        part = part.strip()
        if not part: continue
        parts = part.split(":")
        if len(parts) == 2:
            fields.append((parts[0].strip(), parts[1].strip()))
        else:
            fields.append((parts[0].strip(), parts[0].strip()))
    return fields

def _find_first_input(case_dir):
    for name in sorted(os.listdir(case_dir)):
        if name.endswith(".i") and not name.startswith("sub_"):
            return os.path.join(case_dir, name)
    return None

def _parse_file_base(input_path):
    if not input_path or not os.path.isfile(input_path):
        return None
    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("file_base"):
                parts = line.split("=", 1)
                if len(parts) == 2:
                    val = parts[1].strip()
                    val = val.strip("'\"")
                    return val
    return None

def _strip_placeholders(value):
    return re.sub(r"\$\{[^}]+\}", "", value)

def _find_exodus_by_file_base(case_dir, file_base):
    if not file_base:
        return None
    base_dir = os.path.dirname(file_base)
    base_name = os.path.basename(file_base)
    prefix = _strip_placeholders(base_name)
    output_dir = os.path.join(case_dir, base_dir)
    if not os.path.isdir(output_dir):
        return None
    matches = []
    for root, _, files in os.walk(output_dir):
        for f in files:
            if f.startswith(prefix) and f.endswith(".e"):
                matches.append(os.path.join(root, f))
    if not matches:
        return None
    return max(matches, key=os.path.getmtime)

def _find_case_exodus(case_dir):
    input_path = _find_first_input(case_dir)
    file_base = _parse_file_base(input_path)
    exodus_path = _find_exodus_by_file_base(case_dir, file_base)
    
    if exodus_path and os.path.isfile(exodus_path):
        return exodus_path
        
    for root, _, files in os.walk(case_dir):
        e_files = [os.path.join(root, f) for f in files if f.endswith(".e") and not f.endswith("_cp.e")]
        if e_files:
            return max(e_files, key=os.path.getmtime)
    return None

def _iter_leaf_datasets(obj):
    if obj is None:
        return
    if hasattr(obj, "GetNumberOfPartitionedDataSets") and hasattr(obj, "GetPartitionedDataSet"):
        for i in range(obj.GetNumberOfPartitionedDataSets()):
            yield from _iter_leaf_datasets(obj.GetPartitionedDataSet(i))
        return
    if hasattr(obj, "GetNumberOfPartitions") and hasattr(obj, "GetPartition"):
        for i in range(obj.GetNumberOfPartitions()):
            yield from _iter_leaf_datasets(obj.GetPartition(i))
        return
    if hasattr(obj, "GetNumberOfBlocks") and hasattr(obj, "GetBlock"):
        for i in range(obj.GetNumberOfBlocks()):
            yield from _iter_leaf_datasets(obj.GetBlock(i))
        return
    if hasattr(obj, "GetPoints") and hasattr(obj, "GetPointData"):
        yield obj

def _get_all_point_data(obj, field_name):
    """递归获取所有block的数据：节点坐标与对应的标量值"""
    coords = []
    values = []
    for ds in _iter_leaf_datasets(obj):
        pts = ds.GetPoints()
        pd = ds.GetPointData()
        if not pts or not pd:
            continue
            
        arr = pd.GetArray(field_name)
        if not arr:
            continue
            
        num_pts = pts.GetNumberOfPoints()
        for i in range(num_pts):
            coords.append(pts.GetPoint(i))
            values.append(arr.GetComponent(i, 0))
    return coords, values

def _closest_time(target, available_times):
    if not available_times:
        return None
    return min(available_times, key=lambda t: abs(t - target))

def export_mesh_data(exodus_path, case_dir, config):
    try:
        from paraview.simple import OpenDataFile, UpdatePipeline, CellDatatoPointData
        from paraview import servermanager
    except ImportError as e:
        print(f"错误: 导入ParaView模块失败。{e}")
        return

    output_dir_name = os.environ.get("PV_OUTPUT_DIR_SINGLE", "post_results")
    out_dir = os.path.join(case_dir, output_dir_name)
    os.makedirs(out_dir, exist_ok=True)
    
    reader = OpenDataFile(exodus_path)
    if not reader:
        print(f"[{case_dir}] 无法打开: {exodus_path}")
        return
        
    available_times = []
    if hasattr(reader, "TimestepValues") and reader.TimestepValues:
        available_times = [float(t) for t in reader.TimestepValues]
    elif hasattr(reader, "TimestepValue"):
        available_times = [float(reader.TimestepValue)]
        
    if not available_times:
        print(f"[{case_dir}] 没有时间步数据")
        return

    # 添加 CellDatatoPointData 过滤器，保证网格单元数据（CellData）转换为节点数据（PointData）
    c2p = CellDatatoPointData(Input=reader)
    c2p.ProcessAllArrays = 1

    target_times = config['target_times']
    fields = config['fields']
    
    for field_en, field_zh in fields:
        out_csv = os.path.join(out_dir, f"{field_en}_mesh_data.csv")
        
        all_coords = []
        all_time_vals = []
        headers = ["x", "y", "z"]
        
        for tgt_t in target_times:
            real_t = _closest_time(tgt_t, available_times)
            if real_t is None:
                continue
            
            headers.append(f"{tgt_t}")
            
            UpdatePipeline(time=real_t, proxy=c2p)
            data = servermanager.Fetch(c2p)
            
            coords, vals = _get_all_point_data(data, field_en)
            
            if not all_coords and coords:
                all_coords = coords
            
            # 如果点的数量不一致（可能自适应网格），此脚本暂时仅支持固定网格拓扑
            if len(vals) < len(all_coords):
                # 填充NaN
                vals.extend([None] * (len(all_coords) - len(vals)))
            all_time_vals.append(vals[:len(all_coords)])
        
        if not all_coords:
            print(f"[{case_dir}] 无法获取网格坐标和字段 '{field_en}' 数据")
            continue
            
        print(f"[{case_dir}] 导出点数据: {field_en} -> {out_csv}")
        with open(out_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            for i, coord in enumerate(all_coords):
                row = [coord[0], coord[1], coord[2]]
                for t_idx in range(len(target_times)):
                    if i < len(all_time_vals[t_idx]):
                        row.append(all_time_vals[t_idx][i])
                    else:
                        row.append(None)
                writer.writerow(row)
        return all_coords, all_time_vals, headers[3:]
    return None, None, None

def main():
    studies_subdir = os.environ.get("STUDIES_SUBDIR", "parameter_studies")
    script_dir = os.environ.get("PROJECT_BASE_DIR", os.path.dirname(os.path.abspath(__file__)))
    studies_dir = os.path.join(script_dir, studies_subdir)
    
    if not os.path.isdir(studies_dir):
        print(f"目录不存在: {studies_dir}")
        return

    config = {
        'target_times': _parse_target_times(),
        'fields': _parse_fields(os.environ.get("DATA_FIELDS", "sigma0:断裂强度"))
    }
    
    if not config['target_times']:
        print("未指定 DATA_TARGET_TIMES")
        return
    if not config['fields']:
        print("未指定 DATA_FIELDS")
        return

    case_dirs = [os.path.join(studies_dir, d) for d in os.listdir(studies_dir) 
                 if d.startswith("case_") and os.path.isdir(os.path.join(studies_dir, d))]
    
    # 按照 case 的自然顺序排序，保证合并时的顺序稳定
    case_dirs.sort(key=lambda x: [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', os.path.basename(x))])

    # 用于汇总的数据结构：{ field_en: [(case_name, time_labels, all_time_vals), ...] }
    merged_data = {f[0]: [] for f in config['fields']}
    # 记录第一个 case 的点数，以防止网格拓扑不同导致形状无法对齐，暂时以它作为坐标行的最大上限
    max_pts = 0
                 
    for case_dir in case_dirs:
        exo_path = _find_case_exodus(case_dir)
        if not exo_path:
            print(f"[{os.path.basename(case_dir)}] 找不到 .e 结果文件")
            continue
        coords, time_vals, t_labels = export_mesh_data(exo_path, case_dir, config)
        if coords and time_vals:
            max_pts = max(max_pts, len(coords))
            case_name = os.path.basename(case_dir)
            for field_en, field_zh in config['fields']:
                merged_data[field_en].append((case_name, t_labels, time_vals))

    # 生成各个 field 的汇总 CSV
    output_dir_name = os.environ.get("PV_OUTPUT_DIR_SINGLE", "post_results")
    data_n = int(os.environ.get("DATA_N", "200"))
    
    if max_pts > 0:
        try:
            import numpy as np
            HAS_NUMPY = True
        except ImportError:
            HAS_NUMPY = False
            
        for field_en in merged_data:
            case_records = merged_data[field_en]
            if not case_records:
                continue
                
            merged_csv_path = os.path.join(studies_dir, output_dir_name, f"step7_merged_{field_en}.csv")
            os.makedirs(os.path.dirname(merged_csv_path), exist_ok=True)
            
            headers = []
            for case_name, t_labels, _ in case_records:
                for t in t_labels:
                    headers.append(f"{case_name}_{t}")
            
            with open(merged_csv_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                
                for i in range(max_pts):
                    row = []
                    for case_name, t_labels, time_vals in case_records:
                        for t_idx in range(len(t_labels)):
                            # 如果该 case 网格点偏少，填充空字符串补齐
                            if i < len(time_vals[t_idx]):
                                row.append(time_vals[t_idx][i])
                            else:
                                row.append("")
                    writer.writerow(row)
            print(f"\n✅ 成功生成汇总 CSV: {merged_csv_path}")

            # ======= 扩展：计算全局统一区间的 PDF 和 CDF =======
            if HAS_NUMPY:
                all_vals = []
                for case_name, t_labels, time_vals in case_records:
                    for t_idx in range(len(t_labels)):
                        all_vals.extend([v for v in time_vals[t_idx] if v is not None])
                
                if not all_vals:
                    continue
                    
                v_min, v_max = min(all_vals), max(all_vals)
                if v_max <= v_min:
                    v_max = v_min + 1.0
                    
                # 设定全局范围的横坐标（区间和中心点）
                x_edges = np.linspace(v_min, v_max, data_n + 1)
                x_centers = (x_edges[:-1] + x_edges[1:]) * 0.5
                dx = (v_max - v_min) / data_n
                
                dist_headers = ["x"]
                for case_name, t_labels, _ in case_records:
                    for t in t_labels:
                        dist_headers.append(f"{case_name}_{t}_pdf")
                        dist_headers.append(f"{case_name}_{t}_cdf")
                        
                dist_rows = [[x_centers[i]] for i in range(data_n)]
                
                for case_name, t_labels, time_vals in case_records:
                    for t_idx in range(len(t_labels)):
                        vals = [v for v in time_vals[t_idx] if v is not None]
                        if vals:
                            pdf, _ = np.histogram(vals, bins=x_edges, density=True)
                            cdf = np.cumsum(pdf * dx)
                        else:
                            pdf = np.zeros(data_n)
                            cdf = np.zeros(data_n)
                            
                        # 组装数据
                        for i in range(data_n):
                            dist_rows[i].append(pdf[i])
                            dist_rows[i].append(cdf[i])
                            
                dist_csv_path = os.path.join(studies_dir, output_dir_name, f"step7_merged_{field_en}_dist.csv")
                with open(dist_csv_path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(dist_headers)
                    writer.writerows(dist_rows)
                print(f"✅ 成功生成分布组合 PDF/CDF CSV: {dist_csv_path}")

if __name__ == "__main__":
    main()
