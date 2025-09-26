#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
更稳健高性能的 ParaView 时间序列导出脚本（V2）
- 自动识别字段关联（Point/Cell/Field）并路由
- CellData 批量投到点上再按 ids 行号对齐
- FieldData 作为时间单值列导出（不按ID）
- 向量字段自动展开分量（如 disp:0/1/2）
- 每时间步批量 SaveData：Point/Cell→Point/Field 最多各一次
- 自动探测 ID 列（ids/GlobalNodeId/...）
"""

import os
import csv
from collections import namedtuple

# 用户可改：每隔多少个时间步导出一次（1 为每步）
TIME_STRIDE = 20

# 默认文件与配置
DEFAULT_E_FILE = "gap_conductance1/2D.e"
ExportConfig = namedtuple('ExportConfig', [
    'element_ids_point',   # 用于 PointData（以及 c2p 路径）
    'element_ids_cell',    # 用于直接按 Cell Data 提取（可选）
    'field_list',
    'output_filename'
])

# 字段列表可混合三类（Point/Cell/Field）。中文别名用于列名美化
DEFAULT_CONFIG = ExportConfig(
    element_ids_point=[702],
    element_ids_cell=[487],
    field_list=[
        ('burnup', '燃耗'),                # 常见 CellData
        ('burnup_avg', '燃耗'),
        # ('hoop_stress', '环向应力'),
        # ('effective_creep_strain', '有效蠕变应变'),
        # ('swelling_hoop_strain', '膨胀环向应变'),
        # ('densification_hoop_strain', '密实化环向应变'),
        # ('total_power', '总功率'),
        # ('radial_power_shape', '径向功率分布'),
        # ('qpoint_penetration', '穿透'),
        # ('paired_T', '对偶温度'),
        ('contact_pressure', '接触压力'),
        # ('disp', '位移'),
        # ('ids', 'ID'),
        # ('nodal_area', '节点面积'),
        # ('penetration', '穿透'),
        # ('strain_zz', '轴向应变'),
        # ('T', '温度'),
        # ('x', 'x坐标'),
        # ('h_eq_inner', '内壁热流'),
        # ('h_eq_outer', '外壁热流'),
        # ('h_gap_eq_inner', '内壁间隙热流'),
        # ('h_gap_eq_outer', '外壁间隙热流'),
        # ('inclad_outer_T_avg', '外包壳温度'),
        # ('pellet_inner_heat_rate', '内壁热流率'),
        # ('pellet_inner_perimeter', '内壁周长'),
        # ('pellet_inner_T_avg', '内壁温度'),
        # ('pellet_outer_heat_rate', '外壁热流率'),
        # ('pellet_outer_perimeter', '外壁周长'),
        # ('pellet_outer_T_avg', '外壁温度'),
        # ('pellet_T_max', '最大温度'),
        # ('QA Records', 'QA记录')
    ],
    output_filename="export_data.csv"
)

def detect_id_column(header):
    preferred = ['ids', 'GlobalNodeId', 'GlobalElementId', 'object_id', 'element_id', 'node_id']
    for name in preferred:
        if name in header:
            return name
    for name in header:
        lower = name.lower()
        if lower.endswith('id') or lower.endswith('ids'):
            return name
    return None


def build_row_index_map(csv_file, id_column=None):
    mapping = {}
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            name = id_column if id_column in header else detect_id_column(header)
            if not name or name not in header:
                return {}
            idx = header.index(name)
            for row_idx, row in enumerate(reader):
                if len(row) > idx:
                    try:
                        mapping[int(float(row[idx]))] = row_idx
                    except Exception:
                        continue
        return mapping
    except Exception:
        return {}


def read_scalar_from_field_csv(csv_file, key):
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            # 宽松匹配列
            target_idx = -1
            for i, name in enumerate(header):
                if name == key or name.lower().endswith(key.lower()):
                    target_idx = i
                    break
            if target_idx < 0:
                return float('nan')
            last_val = float('nan')
            for row in reader:
                if len(row) > target_idx:
                    try:
                        last_val = float(row[target_idx])
                    except Exception:
                        continue
            return last_val
    except Exception:
        return float('nan')


def extract_value_by_row(csv_file, col_name, row_index):
    try:
        if row_index < 0:
            return float('nan')
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            # 支持向量分量后缀匹配
            idx = -1
            for i, name in enumerate(header):
                if name == col_name or name.lower().endswith(col_name.lower()):
                    idx = i
                    break
            if idx < 0:
                return float('nan')
            for i, row in enumerate(reader):
                if i == row_index and len(row) > idx:
                    try:
                        return float(row[idx])
                    except Exception:
                        return float('nan')
        return float('nan')
    except Exception:
        return float('nan')


def get_all_fields(reader):
    fields = []
    for attr in ['PointData', 'CellData', 'FieldData']:
        if hasattr(reader, attr):
            data = getattr(reader, attr)
            if data:
                for arr in data:
                    name = arr.GetName()
                    if name:
                        fields.append((name, attr))
    return fields


def expand_vector_components(available_cols, base_name):
    # 返回存在的分量列名列表，如 ['disp:0','disp:1','disp:2'] 或空
    comps = []
    cand = [f"{base_name}:0", f"{base_name}:1", f"{base_name}:2"]
    for c in cand:
        if c in available_cols:
            comps.append(c)
    return comps


def main():
    try:
        from paraview.simple import (
            OpenDataFile, GetAnimationScene, UpdatePipeline, SaveData,
            CellDatatoPointData, GetTimeKeeper
        )
    except Exception:
        print("错误：无法导入 ParaView 模块。请在含 paraview.simple 的环境运行。")
        return

    epath = DEFAULT_E_FILE
    if not os.path.exists(epath):
        print(f"未找到文件: {epath}")
        return

    cfg = DEFAULT_CONFIG
    print(f"使用文件: {epath}")
    print(f"Point IDs: {cfg.element_ids_point}")
    print(f"Cell IDs: {cfg.element_ids_cell}")
    print(f"字段列表: {cfg.field_list}")

    reader = OpenDataFile(epath)
    if not reader:
        print("无法打开文件")
        return

    animation = GetAnimationScene()
    tk = GetTimeKeeper()
    timesteps = []
    try:
        if hasattr(tk, 'TimestepValues') and tk.TimestepValues:
            timesteps = list(tk.TimestepValues)
        elif hasattr(reader, 'TimestepValues') and reader.TimestepValues:
            timesteps = list(reader.TimestepValues)
    except Exception:
        timesteps = []

    num_steps = len(timesteps)
    if num_steps == 0:
        print("未发现时间步")
        return

    # 列出可用字段与可选 IDs 样例
    all_fields = get_all_fields(reader)
    try:
        tmp_ids_csv = "_v2_ids_preview.csv"
        SaveData(tmp_ids_csv, reader, FieldAssociation='Point Data', Precision=6)
        id_map = build_row_index_map(tmp_ids_csv)
        if os.path.exists(tmp_ids_csv):
            os.remove(tmp_ids_csv)
        preview = sorted(list(id_map.keys()))[:10]
        print(f"可选element_ids(点数据) 数量: {len(id_map)}，示例: {preview}")
    except Exception:
        pass
    print("可用字段(前50):")
    for name, atype in all_fields[:50]:
        print(f"  {atype}: {name}")

    # 将请求字段按关联分类
    names_point = set(n for n, t in all_fields if t == 'PointData')
    names_cell = set(n for n, t in all_fields if t == 'CellData')
    names_field = set(n for n, t in all_fields if t == 'FieldData')

    req_point = []
    req_cell = []
    req_field = []
    for fname, alias in cfg.field_list:
        if fname in names_point:
            req_point.append((fname, alias))
        elif fname in names_cell:
            req_cell.append((fname, alias))
        elif fname in names_field:
            req_field.append((fname, alias))
        else:
            # 宽松匹配后缀
            matched = False
            for pool, collector in ((names_point, req_point), (names_cell, req_cell), (names_field, req_field)):
                hit = next((n for n in pool if n.lower().endswith(fname.lower())), None)
                if hit:
                    collector.append((hit, alias))
                    matched = True
                    break
            if not matched:
                print(f"警告: 字段 '{fname}' 不存在，已跳过")

    # 输出 CSV
    out_csv = os.path.join(os.path.dirname(epath), cfg.output_filename)
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # 表头：时间 + FieldData单值列 + 按ID的Point/Cell→Point列
        header = ['时间 (s)']
        # FieldData 列（不带 _ID）
        for fname, alias in req_field:
            header.append(alias)
        # Point/Cell 列
        # 注意向量分量：列名用 "别名[:i]_ID"
        # 先不写，按首个时间步确认实际分量列后再扩展
        writer.writerow(header)

    # 预构建所需：CellData→PointData 过滤器
    c2p = None
    if req_cell:
        c2p = CellDatatoPointData(Input=reader)
        try:
            c2p.PassCellDataArrays = 1
        except Exception:
            pass

    # 时间循环（按步长）
    indices = list(range(0, num_steps, max(1, TIME_STRIDE)))
    total = len(indices)

    # 缓存首次时间步的分量列名（用于表头扩展）
    expanded_cols_cache = None

    for k, t_idx in enumerate(indices):
        is_last = (k == total - 1)
        if is_last:
            print(f"处理时间步 {k + 1}/{total} …")

        # 设置时间并刷新
        animation.AnimationTime = float(timesteps[t_idx])
        UpdatePipeline(time=float(timesteps[t_idx]), proxy=reader)
        if c2p is not None:
            try:
                UpdatePipeline(time=float(timesteps[t_idx]), proxy=c2p)
            except Exception:
                pass

        current_time = float(timesteps[t_idx])

        # 批量导出 CSV：PointData / Cell->Point / FieldData
        temp_point_csv = None
        temp_cellp_csv = None
        temp_cell_csv = None
        temp_field_csv = None

        # PointData（一次导出全部点列，后续按列名与 id 匹配）
        if req_point:
            temp_point_csv = f"_v2_point_{t_idx}.csv"
            # 仅写所需点数组，强制包含 ids
            point_arrays = list({fname for fname, _ in req_point})
            if 'ids' not in point_arrays:
                point_arrays.append('ids')
            try:
                SaveData(temp_point_csv, reader,
                         FieldAssociation='Point Data',
                         Precision=6,
                         PointDataArrays=point_arrays)
            except TypeError:
                # 某些版本用 ChooseArraysToWrite + PointDataArrays
                SaveData(temp_point_csv, reader,
                         FieldAssociation='Point Data',
                         Precision=6,
                         ChooseArraysToWrite=1,
                         PointDataArrays=point_arrays)

        # CellData：两种路径
        if req_cell:
            if cfg.element_ids_cell:  # 直接按 Cell IDs 提取
                temp_cell_csv = f"_v2_cell_{t_idx}.csv"
                cell_arrays = list({fname for fname, _ in req_cell})
                try:
                    SaveData(temp_cell_csv, reader,
                             FieldAssociation='Cell Data',
                             Precision=6,
                             CellDataArrays=cell_arrays)
                except TypeError:
                    SaveData(temp_cell_csv, reader,
                             FieldAssociation='Cell Data',
                             Precision=6,
                             ChooseArraysToWrite=1,
                             CellDataArrays=cell_arrays)
            else:  # 投到点上用 point ids
                if c2p is not None:
                    temp_cellp_csv = f"_v2_cellp_{t_idx}.csv"
                    cellp_arrays = list({fname for fname, _ in req_cell})
                    try:
                        SaveData(temp_cellp_csv, c2p,
                                 FieldAssociation='Point Data',
                                 Precision=6,
                                 PointDataArrays=cellp_arrays)
                    except TypeError:
                        SaveData(temp_cellp_csv, c2p,
                                 FieldAssociation='Point Data',
                                 Precision=6,
                                 ChooseArraysToWrite=1,
                                 PointDataArrays=cellp_arrays)

        # FieldData（单值，含时间列）
        if req_field:
            temp_field_csv = f"_v2_field_{t_idx}.csv"
            field_arrays = list({fname for fname, _ in req_field})
            try:
                SaveData(temp_field_csv, reader,
                         FieldAssociation='Field Data',
                         Precision=6,
                         AddTime=1,
                         FieldDataArrays=field_arrays)
            except TypeError:
                SaveData(temp_field_csv, reader,
                         FieldAssociation='Field Data',
                         Precision=6,
                         AddTime=1,
                         ChooseArraysToWrite=1,
                         FieldDataArrays=field_arrays)

        # 构造该时间步的数据行：时间 + field 单值 + point/cellp 按ID
        row = [current_time]

        # 1) FieldData 单值
        if req_field and temp_field_csv:
            for fname, alias in req_field:
                row.append(read_scalar_from_field_csv(temp_field_csv, fname))

        # 2) PointData 按 ID
        row_map_point = {}
        if temp_point_csv:
            row_map_point = build_row_index_map(temp_point_csv)

        # 3) CellData -> 直接按 Cell IDs 或投点按 Point IDs
        row_map_ids_point = None
        if temp_cellp_csv:
            if not temp_point_csv:
                tmp_ids = f"_v2_ids_{t_idx}.csv"
                SaveData(tmp_ids, reader, FieldAssociation='Point Data', Precision=6)
                row_map_ids_point = build_row_index_map(tmp_ids)
                if os.path.exists(tmp_ids):
                    os.remove(tmp_ids)
            else:
                row_map_ids_point = row_map_point

        # 在首个时间步确定向量分量展开，并把列头补全
        if expanded_cols_cache is None:
            expanded_cols_cache = []  # 列描述: (source, col_name, alias, element_id or None)

            # FieldData 已在表头写过

            # PointData：读取一次表头，判断向量分量
            point_header = []
            if temp_point_csv and os.path.exists(temp_point_csv):
                with open(temp_point_csv, 'r', encoding='utf-8') as fp:
                    point_header = next(csv.reader(fp), [])

            for fname, alias in req_point:
                comps = expand_vector_components(point_header, fname)
                target_cols = comps if comps else [fname]
                for col in target_cols:
                    for eid in cfg.element_ids_point:
                        expanded_cols_cache.append(('point', col, alias + (f":{col.split(':')[-1]}" if ':' in col else ''), eid))

            # Cell->Point：同理
            cellp_header = []
            if temp_cellp_csv and os.path.exists(temp_cellp_csv):
                with open(temp_cellp_csv, 'r', encoding='utf-8') as fp:
                    cellp_header = next(csv.reader(fp), [])

            for fname, alias in req_cell:
                comps = expand_vector_components(cellp_header, fname)
                if cfg.element_ids_cell:  # 直接 cell 路径
                    target_cols = comps if comps else [fname]
                    for col in target_cols:
                        for eid in cfg.element_ids_cell:
                            expanded_cols_cache.append(('cell', col, alias + (f":{col.split(':')[-1]}" if ':' in col else ''), eid))
                else:
                    target_cols = comps if comps else [fname]
                    for col in target_cols:
                        for eid in cfg.element_ids_point:
                            expanded_cols_cache.append(('cellp', col, alias + (f":{col.split(':')[-1]}" if ':' in col else ''), eid))

            # 把完整表头写回
            with open(out_csv, 'r+', newline='', encoding='utf-8') as f:
                reader_w = csv.reader(f)
                header = next(reader_w)
                f.seek(0)
                writer = csv.writer(f)
                writer.writerow(header + [f"{alias}_{eid}" if eid is not None else alias for _, _, alias, eid in expanded_cols_cache])
                # 重写已有数据行（无）
                for _ in reader_w:
                    pass

        # 取值写行
        # 顺序需与 expanded_cols_cache 对齐
        for source, col, alias, eid in expanded_cols_cache:
            if source == 'point' and temp_point_csv and eid in row_map_point:
                row_idx = row_map_point.get(eid, -1)
                row.append(extract_value_by_row(temp_point_csv, col, row_idx))
            elif source == 'cellp' and temp_cellp_csv and row_map_ids_point:
                row_idx = row_map_ids_point.get(eid, -1)
                row.append(extract_value_by_row(temp_cellp_csv, col, row_idx))
            elif source == 'cell' and temp_cell_csv:
                # 直接在单元数据里按单元ID匹配
                # 构建一次 cell id -> 行号 映射
                # 为避免多次开销，这里简单每次查，数据量通常较小
                with open(temp_cell_csv, 'r', encoding='utf-8') as fc:
                    reader_c = csv.reader(fc)
                    header_c = next(reader_c, [])
                    id_name = detect_id_column(header_c)
                    if id_name and id_name in header_c and col in header_c:
                        id_idx = header_c.index(id_name)
                        col_idx = header_c.index(col)
                        val = float('nan')
                        for row_c in reader_c:
                            try:
                                if int(float(row_c[id_idx])) == eid:
                                    val = float(row_c[col_idx]) if row_c[col_idx] != '' else float('nan')
                                    break
                            except Exception:
                                continue
                        row.append(val)
                    else:
                        row.append(float('nan'))
            else:
                row.append(float('nan'))

        # 追加到 CSV
        with open(out_csv, 'a', newline='', encoding='utf-8') as f:
            csv.writer(f).writerow(row)

        # 清理临时文件
        for path in (temp_point_csv, temp_cellp_csv, temp_field_csv, temp_cell_csv):
            if path and os.path.exists(path):
                os.remove(path)

    print("\n数据导出完成！")
    print(f"输出文件: {out_csv}")
    print(f"导出时间步数: {len(indices)} (步长={TIME_STRIDE})")


if __name__ == '__main__':
    main()


