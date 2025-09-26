#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用ParaView官方SaveData函数导出数据到CSV文件
参考：https://docs.paraview.org/en/latest/UsersGuide/savingResults.html
"""

import os
import sys
import csv
from collections import namedtuple

# 采样步长：每隔多少个时间步导出一次（1 表示每步都导出）
TIME_STRIDE = 50

# 配置参数
ExportConfig = namedtuple('ExportConfig', [
    'element_ids',       # 要导出的ID列表（自动匹配ID列）
    'field_list',        # 要导出的字段列表 [('burnup','燃耗'),('hoop_stress','环向应力')]
    'output_filename',   # 输出CSV文件名
    'include_burnup'     # 是否包含燃耗数据
])

# 默认.e文件路径
DEFAULT_E_FILE = "gap_conductance1/2D.e"

# 默认配置
DEFAULT_CONFIG = ExportConfig(
    element_ids=[497],
    field_list=[
        ('burnup', '燃耗'),
        ('burnup_avg', '燃耗'),
        ('contact_pressure', '接触压力'),
        ('disp', '位移'),
        ('ids', 'ID'),
        ('nodal_area', '节点面积'),
        ('penetration', '穿透'),
        ('strain_zz', '轴向应变'),
        ('T', '温度'),
        ('x', 'x坐标'),
        ('h_eq_inner', '内壁热流'),
        ('h_eq_outer', '外壁热流'),
        ('h_gap_eq_inner', '内壁间隙热流'),
        ('h_gap_eq_outer', '外壁间隙热流'),
        ('inclad_outer_T_avg', '外包壳温度'),
        ('pellet_inner_heat_rate', '内壁热流率'),
        ('pellet_inner_perimeter', '内壁周长'),
        ('pellet_inner_T_avg', '内壁温度'),
        ('pellet_outer_heat_rate', '外壁热流率'),
        ('pellet_outer_perimeter', '外壁周长'),
        ('pellet_outer_T_avg', '外壁温度'),
        ('pellet_T_max', '最大温度'),
        ('QA Records', 'QA记录')
    ],
    output_filename="export_data.csv",
    include_burnup=False
)

def get_all_fields(reader):
    """获取.e文件中所有可用的字段"""
    fields = []
    # 检查所有可能的数据属性
    for attr in ['PointData', 'CellData', 'FieldData']:
        if hasattr(reader, attr):
            data = getattr(reader, attr)
            if data:
                for arr in data:
                    field_name = arr.GetName()
                    if field_name:
                        fields.append((field_name, attr))
    return fields

def export_data_using_savedata(e_file_path, config):
    """使用ParaView的SaveData函数导出数据"""
    try:
        from paraview.simple import (
            OpenDataFile,
            GetAnimationScene,
            UpdatePipeline,
            SaveData,
            CellDatatoPointData,
            GetTimeKeeper,
        )
        
        # 打开.e文件
        reader = OpenDataFile(e_file_path)
        if not reader:
            print(f"无法打开文件: {e_file_path}")
            return None
        
        # 获取动画场景与TimeKeeper（更稳健的时间控制）
        animation_scene = GetAnimationScene()
        if not animation_scene:
            print("无法获取动画场景")
            return None
        time_keeper = GetTimeKeeper()
        
        # 获取时间步数组
        timesteps = []
        try:
            if hasattr(time_keeper, 'TimestepValues') and time_keeper.TimestepValues:
                timesteps = list(time_keeper.TimestepValues)
            elif hasattr(reader, 'TimestepValues') and reader.TimestepValues:
                timesteps = list(reader.TimestepValues)
        except:
            timesteps = []

        num_timesteps = len(timesteps) if timesteps else 0
        if num_timesteps == 0:
            # 兜底：维持原有逻辑
            try:
                if hasattr(animation_scene, 'GetNumberOfDataTimeSteps'):
                    num_timesteps = animation_scene.GetNumberOfDataTimeSteps()
                elif hasattr(animation_scene, 'GetNumberOfTimeSteps'):
                    num_timesteps = animation_scene.GetNumberOfTimeSteps()
                elif hasattr(reader, 'TimestepValues'):
                    num_timesteps = len(reader.TimestepValues)
                else:
                    num_timesteps = 0
            except:
                num_timesteps = 0
        
        # 静默：不打印每一步详细信息
        
        # 获取可用字段
        available_fields = get_all_fields(reader)
        # 启动时打印一次：可用ID与字段概览
        try:
            # 预探测点数据ID集合
            preview_ids_csv = "temp_preview_ids.csv"
            SaveData(preview_ids_csv, reader, FieldAssociation='Point Data', Precision=6)
            id_preview_map = build_row_index_map(preview_ids_csv)
            if os.path.exists(preview_ids_csv):
                os.remove(preview_ids_csv)
            preview_id_list = sorted(list(id_preview_map.keys()))
            print(f"可选element_ids(点数据) 数量: {len(preview_id_list)}，示例: {preview_id_list[:10]}")
        except Exception:
            pass
        print("可用字段:")
        for fname, ftype in available_fields[:50]:
            print(f"  {ftype}: {fname}")

        # 打印一次 burnup_avg 的类型信息，便于判断提取方式
        try:
            fdi = reader.GetFieldDataInformation()
            if fdi:
                arr = None
                for i in range(fdi.GetNumberOfArrays()):
                    ai = fdi.GetArrayInformation(i)
                    if ai and ai.GetName() == 'burnup':
                        arr = ai
                        break
                if arr:
                    dt = arr.GetDataTypeAsString() if hasattr(arr, 'GetDataTypeAsString') else 'Unknown'
                    num_comp = arr.GetNumberOfComponents() if hasattr(arr, 'GetNumberOfComponents') else -1
                    num_tuples = arr.GetNumberOfTuples() if hasattr(arr, 'GetNumberOfTuples') else -1
                    print(f"burnup 类型: {dt}, 分量: {num_comp}, 元组数: {num_tuples}")
        except Exception:
            pass
        
        # 过滤有效的字段
        valid_fields = []
        for field_name, field_desc in config.field_list:
            field_found = False
            for available_field, attr_type in available_fields:
                if available_field == field_name:
                    valid_fields.append((field_name, field_desc, attr_type))
                    # 静默
                    field_found = True
                    break
            if not field_found:
                # 静默
                pass
        
        if not valid_fields:
            print("没有找到任何有效字段")
            return None
        
        # 遍历时间步索引集合（支持步长抽样），不再限制时间步数量
        timestep_indices = list(range(0, num_timesteps, max(1, TIME_STRIDE))) if num_timesteps > 0 else []
        
        # 准备数据
        data_rows = []
        
        # 将CellData转换为PointData，便于按节点ID取值
        c2p = CellDatatoPointData(Input=reader)
        # 通常默认即可，这里确保传递所有Cell数组
        try:
            c2p.PassCellDataArrays = 1
        except Exception:
            pass

        # 处理每个时间步（按步长）
        total_steps = len(timestep_indices)
        for step_idx, timestep in enumerate(timestep_indices):
            is_last = (step_idx == total_steps - 1)
            if is_last:
                print(f"处理时间步 {step_idx + 1}/{total_steps}...")
            
            try:
                # 设置当前时间（优先使用TimeKeeper的时间数组）
                try:
                    if timesteps:
                        animation_scene.AnimationTime = float(timesteps[timestep])
                    else:
                        # 兜底：使用索引推进
                        if hasattr(animation_scene, 'SetCurrentTimeStep'):
                            animation_scene.SetCurrentTimeStep(timestep)
                except Exception:
                    pass
                
                # 用当前时间刷新 reader 与 c2p，确保SaveData导出对应时间步数据
                if timesteps:
                    UpdatePipeline(time=float(timesteps[timestep]), proxy=reader)
                    try:
                        UpdatePipeline(time=float(timesteps[timestep]), proxy=c2p)
                    except Exception:
                        pass
                else:
                    UpdatePipeline()
                
                # 获取当前时间
                current_time = 0.0
                try:
                    if timesteps:
                        current_time = float(timesteps[timestep])
                    elif hasattr(animation_scene, 'GetTimeValue'):
                        current_time = float(animation_scene.GetTimeValue())
                except:
                    pass
                
                # 构建行数据
                row_data = [current_time]
                
                # 获取燃耗值（如果启用）
                if config.include_burnup:
                    burnup_value = float('nan')
                    # 从Field Data导出并读取burnup_avg
                    try:
                        temp_field_csv = f"temp_fielddata_{timestep}.csv"
                        # 显式指定Field Data，并附加时间列，便于取对应时间步
                        SaveData(temp_field_csv, reader, FieldAssociation='Field Data', Precision=6, AddTime=1)
                        burnup_value = extract_fielddata_value(temp_field_csv, 'burnup_avg')
                        if os.path.exists(temp_field_csv):
                            os.remove(temp_field_csv)
                    except Exception:
                        pass
                    row_data.append(burnup_value)
                
                # 优化：对每个字段只导出一次，然后提取所有ID的数据
                for field_name, field_desc, attr_type in valid_fields:
                    if is_last:
                        print(f"  处理字段: {field_name}")
                    
                    # 创建临时CSV文件名
                    temp_csv = f"temp_{field_name}_{timestep}.csv"
                    
                    try:
                        # 使用SaveData导出数据
                        # 对CellData字段，改为使用c2p后的Point Data导出，便于按节点ID对齐
                        if attr_type == 'PointData':
                            # 直接按点数据导出
                            SaveData(temp_csv, reader, FieldAssociation='Point Data', Precision=6)
                            # 自动探测ID列并提取
                            for element_id in config.element_ids:
                                value = extract_value_from_csv(temp_csv, field_name, element_id, id_column=None)
                                row_data.append(value)
                        elif attr_type == 'CellData':
                            # 对CellData：
                            # 1) 导出reader点数据，构建 id -> 行号 映射
                            temp_ids_csv = f"temp_ids_{timestep}.csv"
                            SaveData(temp_ids_csv, reader, FieldAssociation='Point Data', Precision=6)
                            # 2) 导出c2p点数据（单元插值到点）
                            SaveData(temp_csv, c2p, FieldAssociation='Point Data', Precision=6)
                            # 3) 用行号在c2p文件中取该字段值
                            row_map = build_row_index_map(temp_ids_csv)
                            for element_id in config.element_ids:
                                row_index = row_map.get(element_id, -1)
                                value = extract_value_from_csv_by_row(temp_csv, field_name, row_index)
                                row_data.append(value)
                            if os.path.exists(temp_ids_csv):
                                os.remove(temp_ids_csv)
                        else:  # FieldData等
                            SaveData(temp_csv, reader, Precision=6)
                            for element_id in config.element_ids:
                                row_data.append(float('nan'))
                        
                        # 删除临时文件
                        if os.path.exists(temp_csv):
                            os.remove(temp_csv)
                            
                    except Exception as e:
                        if is_last:
                            print(f"    导出字段 {field_name} 时出错: {str(e)}")
                        # 如果出错，为所有ID添加nan值
                        for element_id in config.element_ids:
                            row_data.append(float('nan'))
                
                data_rows.append(row_data)
                if is_last:
                    print(f"  时间步 {timestep + 1} 完成")
                
            except Exception as e:
                print(f"处理时间步 {timestep} 时出错: {str(e)}")
                continue
        
        # 静默
        return data_rows, valid_fields
        
    except Exception as e:
        print(f"提取数据时出错: {str(e)}")
        return None

def detect_id_column(header):
    """在表头中自动探测ID列名称。返回列名或None。"""
    preferred = ['ids', 'GlobalNodeId', 'GlobalElementId', 'object_id', 'element_id', 'node_id']
    for name in preferred:
        if name in header:
            return name
    # 宽松：大小写不敏感匹配末尾
    for name in header:
        lower = name.lower()
        if lower.endswith('id') or lower.endswith('ids'):
            return name
    return None


def extract_value_from_csv(csv_file, field_name, element_id, id_column=None):
    """从临时CSV文件中提取特定字段和ID的值"""
    try:
        if not os.path.exists(csv_file):
            print(f"    临时文件不存在: {csv_file}")
            return float('nan')
        
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader)  # 读取表头
            
            # 查找字段列和ID列
            field_col = None
            id_col = None
            
            for i, col_name in enumerate(header):
                # 查找字段列 - 精确匹配
                if col_name == field_name:
                    field_col = i
                    # 静默
                
                # 查找ID列
                if id_column is not None and col_name == id_column:
                    id_col = i
                    # 静默
            
            # 如果精确匹配失败，尝试模糊匹配
            if field_col is None:
                for i, col_name in enumerate(header):
                    if (field_name.lower() in col_name.lower() or 
                        col_name.lower().endswith(field_name.lower())):
                        field_col = i
                        # 静默
                        break
            
            # 自动探测ID列
            if id_col is None:
                auto_name = detect_id_column(header)
                if auto_name is not None and auto_name in header:
                    id_col = header.index(auto_name)
            
            if field_col is None:
                # 静默
                return float('nan')
            
            if id_col is None:
                # 静默
                return float('nan')
            
            # 读取数据行，找到对应的ID
            row_count = 0
            for row in reader:
                row_count += 1
                if len(row) > max(field_col, id_col):
                    try:
                        current_id = int(row[id_col])
                        if current_id == element_id:
                            value = float(row[field_col])
                            # 静默
                            return value
                    except (ValueError, IndexError) as e:
                        # 静默
                        continue
            
            # 静默
            return float('nan')
        
    except Exception as e:
        print(f"    读取CSV文件时出错: {str(e)}")
        return float('nan')


def extract_fielddata_value(csv_file, field_key):
    """从导出的CSV的FieldData区域解析指定字段（如 burnup_avg）。
    由于不同版本的SaveData对FieldData写法不同，这里做宽松匹配：
    - 在表头中直接存在field_key列
    - 或在行中以"field_key,"形式出现且紧随其后的数值
    """
    try:
        if not os.path.exists(csv_file):
            return float('nan')
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            # 1) 直接在表头找列（或宽松匹配列名后缀）
            target_idx = -1
            for i, name in enumerate(header):
                if name == field_key or name.lower().endswith(field_key.lower()):
                    target_idx = i
                    break
            if target_idx >= 0:
                last_value = float('nan')
                for row in reader:
                    if len(row) > target_idx:
                        try:
                            last_value = float(row[target_idx])
                        except Exception:
                            continue
                return last_value
        # 2) 逐行扫描纯文本匹配（保底）
        with open(csv_file, 'r', encoding='utf-8') as f:
            for line in f:
                if field_key in line:
                    parts = [p.strip() for p in line.strip().split(',') if p.strip()]
                    # 尝试在该行中找第一个可解析为数值的字段
                    for token in parts:
                        try:
                            return float(token)
                        except Exception:
                            continue
        return float('nan')
    except Exception:
        return float('nan')

def find_row_index_by_id(csv_file, element_id):
    """在含ids列的点数据CSV中，查找给定ID的行号（从0开始，不含表头）。找不到返回-1。"""
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            if 'ids' not in header:
                return -1
            id_col = header.index('ids')
            for idx, row in enumerate(reader):
                if len(row) > id_col:
                    try:
                        if int(row[id_col]) == element_id:
                            return idx
                    except Exception:
                        continue
        return -1
    except Exception:
        return -1

def extract_value_from_csv_by_row(csv_file, field_name, row_index):
    """按行号读取某列值。row_index为数据行索引（不含表头）。失败返回nan。"""
    try:
        if row_index < 0:
            return float('nan')
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, [])
            if field_name not in header:
                # 宽松匹配
                for i, name in enumerate(header):
                    if name == field_name or name.lower().endswith(field_name.lower()):
                        field_idx = i
                        break
                else:
                    return float('nan')
            else:
                field_idx = header.index(field_name)
            for idx, row in enumerate(reader):
                if idx == row_index and len(row) > field_idx:
                    try:
                        return float(row[field_idx])
                    except Exception:
                        return float('nan')
        return float('nan')
    except Exception:
        return float('nan')


def build_row_index_map(csv_file, id_column=None):
    """从点数据CSV构建 id -> 行号 映射。失败返回空dict。"""
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
                        mapping[int(row[idx])] = row_idx
                    except Exception:
                        continue
        return mapping
    except Exception:
        return {}

def main():
    """主函数"""
    # 使用指定的.e文件路径
    e_file_path = DEFAULT_E_FILE
    if not os.path.exists(e_file_path):
        print(f"未找到文件: {e_file_path}")
        return
    
    print(f"使用文件: {e_file_path}")
    print(f"元素ID: {DEFAULT_CONFIG.element_ids}")
    print(f"字段列表: {DEFAULT_CONFIG.field_list}")
    
    # 使用SaveData方法提取数据
    result = export_data_using_savedata(e_file_path, DEFAULT_CONFIG)
    if result is None:
        print("数据提取失败")
        return
    
    data_rows, valid_fields = result
    
    # 生成输出文件名
    output_filename = DEFAULT_CONFIG.output_filename
    output_path = os.path.join(os.path.dirname(e_file_path), output_filename)
    
    # 写入CSV文件
    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        # 写入表头
        header = ['时间 (s)']
        if DEFAULT_CONFIG.include_burnup:
            header.append('燃耗 (MWd/tU)')
        
        for field_name, field_desc in DEFAULT_CONFIG.field_list:
            for element_id in DEFAULT_CONFIG.element_ids:
                header.append(f"{field_desc}_{element_id}")
        
        writer = csv.writer(csvfile)
        writer.writerow(header)
        
        # 写入数据行
        for row_data in data_rows:
            writer.writerow(row_data)
    
    print(f"\n数据导出完成！")
    print(f"输出文件: {output_path}")
    print(f"数据行数: {len(data_rows)}")
    print(f"有效字段: {[f[0] for f in valid_fields]}")

if __name__ == "__main__":
    main()
