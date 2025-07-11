#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ParaView/VTK数据提取脚本
从.e文件中提取指定节点/单元的数据并导出到Excel
"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd

# 设置matplotlib中文字体
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

# ================== 配置参数 ==================
# 请根据需要修改以下参数

# 数据文件路径
DATA_FILE = "elastoplasticity_out_1.5.e"

# 要提取的节点ID列表
NODE_IDS = [1, 10, 50, 100, 200, 500]

# 要提取的单元ID列表  
CELL_IDS = [1, 10, 50, 100, 200, 500]

# 时间范围 [开始时间, 结束时间]
TIME_RANGE = [0.0, 1000.0]  # 增大时间范围，获取后期的塑性应变数据

# 要提取的字段列表 [字段名, 数据类型]
# 数据类型: 'point' 或 'cell'
FIELD_LIST = [
    ('d', 'point'),                      # 相场变量
    ('effective_plastic_strain', 'cell'), # 有效塑性应变
]

# 输出Excel文件名
OUTPUT_FILE = "extracted_data.xlsx"

# ================== 数据读取模块 ==================

def extract_data_from_exodus(file_path):
    """从Exodus文件中提取数据"""
    try:
        from paraview.simple import (
            OpenDataFile, UpdatePipeline, GetAnimationScene, MergeBlocks
        )
        import paraview.simple as pv
        import paraview.servermanager as sm
        pv._DisableFirstRenderCameraReset()
        
        # 导入vtk用于数据处理
        import vtk
        from vtk.util import numpy_support
        
    except ImportError as e:
        print(f"❌ 无法导入ParaView模块: {e}")
        print("请确保ParaView环境已正确配置")
        return None

    try:
        # 1. 检查文件是否存在
        if not os.path.exists(file_path):
            print(f"❌ 数据文件不存在: {file_path}")
            return None
            
        # 2. 读取文件
        print(f"📖 正在读取文件: {os.path.basename(file_path)}")
        reader = pv.OpenDataFile(file_path)
        if not reader:
            print("❌ 无法读取文件")
            return None

        # 3. 获取时间步信息
        time_steps = []
        if hasattr(reader, 'TimestepValues'):
            time_steps = [float(t) for t in reader.TimestepValues]
        elif hasattr(reader, 'TimestepValue'):
            time_steps = [float(reader.TimestepValue)]
        
        if not time_steps:
            print("⚠️  未找到时间步信息")
            return None

        print(f"⏱️  找到 {len(time_steps)} 个时间步")
        print(f"⏱️  时间范围: {min(time_steps):.2e} - {max(time_steps):.2e}")

        # 4. 确定要提取的时间点
        start_time, end_time = TIME_RANGE
        target_times = [t for t in time_steps if start_time <= t <= end_time]
        
        if not target_times:
            print("⚠️  没有找到符合条件的时间点")
            return None

        print(f"📊 将提取 {len(target_times)} 个时间点的数据")

        # 5. 初始化数据存储
        extracted_data = {}
        
        # 6. 处理多块数据集 - 使用MergeBlocks
        print("🔄 合并多块数据集...")
        merged_data = MergeBlocks(Input=reader)
        pv.UpdatePipeline(proxy=merged_data)
        
        # 使用正确的API获取数据
        data_object = sm.Fetch(merged_data)
        
        if not data_object:
            print("❌ 无法获取数据对象")
            return None

        # 7. 检查可用字段
        point_arrays = []
        cell_arrays = []
        
        point_data = data_object.GetPointData()
        for i in range(point_data.GetNumberOfArrays()):
            array_name = point_data.GetArrayName(i)
            if array_name:
                point_arrays.append(array_name)
        
        cell_data = data_object.GetCellData()
        for i in range(cell_data.GetNumberOfArrays()):
            array_name = cell_data.GetArrayName(i)
            if array_name:
                cell_arrays.append(array_name)

        print(f"📋 可用点数据字段: {point_arrays}")
        print(f"📋 可用单元数据字段: {cell_arrays}")

        # 8. 处理每个时间步
        animationScene = pv.GetAnimationScene()
        
        for time_idx, target_time in enumerate(target_times):
            print(f"\n🔄 处理时间步 {time_idx+1}/{len(target_times)}: {target_time:.2e}")
            
            # 设置时间步
            animationScene.AnimationTime = float(target_time)
            pv.UpdatePipeline(time=target_time, proxy=merged_data)
            
            # 获取当前时间步的数据
            current_data = sm.Fetch(merged_data)
            
            # 9. 提取每个字段的数据
            for field_name, data_type in FIELD_LIST:
                if field_name not in extracted_data:
                    extracted_data[field_name] = {
                        'times': [],
                        'data': {},
                        'type': data_type
                    }
                
                # 记录时间
                if target_time not in extracted_data[field_name]['times']:
                    extracted_data[field_name]['times'].append(target_time)
                
                # 获取字段数据
                if data_type == 'point':
                    if field_name not in point_arrays:
                        print(f"  ⚠️  点数据字段 '{field_name}' 不存在")
                        continue
                    
                    data_array = current_data.GetPointData().GetArray(field_name)
                    if not data_array:
                        print(f"  ⚠️  无法获取点数据 '{field_name}'")
                        continue
                    
                    # 转换为numpy数组
                    np_array = numpy_support.vtk_to_numpy(data_array)
                    
                    # 提取指定ID的数据
                    for node_id in NODE_IDS:
                        if node_id < len(np_array):
                            if node_id not in extracted_data[field_name]['data']:
                                extracted_data[field_name]['data'][node_id] = []
                            extracted_data[field_name]['data'][node_id].append(np_array[node_id])
                        else:
                            print(f"  ⚠️  节点ID {node_id} 超出范围")
                
                elif data_type == 'cell':
                    if field_name not in cell_arrays:
                        print(f"  ⚠️  单元数据字段 '{field_name}' 不存在")
                        continue
                    
                    data_array = current_data.GetCellData().GetArray(field_name)
                    if not data_array:
                        print(f"  ⚠️  无法获取单元数据 '{field_name}'")
                        continue
                    
                    # 转换为numpy数组
                    np_array = numpy_support.vtk_to_numpy(data_array)
                    
                    # 添加调试信息，特别是对effective_plastic_strain
                    if field_name == 'effective_plastic_strain':
                        non_zero_count = np.count_nonzero(np_array)
                        print(f"  📊 {field_name}: 非零值数量 = {non_zero_count}/{len(np_array)}, max = {np_array.max():.8f}")
                        if non_zero_count > 0:
                            # 找出有非零值的单元ID
                            non_zero_indices = np_array.nonzero()[0]
                            print(f"    非零值单元ID: {non_zero_indices[:5]}...")
                    
                    # 提取指定ID的数据
                    for cell_id in CELL_IDS:
                        if cell_id < len(np_array):
                            if cell_id not in extracted_data[field_name]['data']:
                                extracted_data[field_name]['data'][cell_id] = []
                            value = np_array[cell_id]
                            extracted_data[field_name]['data'][cell_id].append(value)
                            
                            # 添加调试信息
                            if field_name == 'effective_plastic_strain' and value != 0:
                                print(f"    单元ID {cell_id}: {value:.8f}")
                        else:
                            print(f"  ⚠️  单元ID {cell_id} 超出范围")

        return extracted_data

    except Exception as e:
        print(f"❌ 数据提取失败: {str(e)}")
        return None

def export_to_excel(data, output_file):
    """导出数据到Excel文件"""
    try:
        print(f"\n📊 正在导出数据到Excel...")
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            
            for field_name, field_data in data.items():
                print(f"  📝 导出字段: {field_name}")
                
                # 获取时间列表
                times = sorted(field_data['times'])
                
                # 创建DataFrame
                df_data = {'Time': times}
                
                # 添加每个ID的数据
                for id_val, values in field_data['data'].items():
                    if len(values) == len(times):
                        df_data[f'ID_{id_val}'] = values
                    else:
                        print(f"    ⚠️  ID {id_val} 的数据长度不匹配")
                
                # 创建DataFrame并写入Excel
                df = pd.DataFrame(df_data)
                sheet_name = field_name[:31]  # Excel工作表名长度限制
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                print(f"    ✅ 工作表 '{sheet_name}' 已创建，数据形状: {df.shape}")
        
        print(f"\n✅ 数据已成功导出到: {output_file}")
        return True
        
    except Exception as e:
        print(f"❌ Excel导出失败: {str(e)}")
        return False

def main():
    """主函数"""
    print("=" * 60)
    print("🚀 ParaView数据提取工具")
    print("=" * 60)
    
    print("🔍 开始数据提取...")
    
    print(f"📁 数据文件: {DATA_FILE}")
    print(f"🔢 节点ID: {NODE_IDS}")
    print(f"🔢 单元ID: {CELL_IDS}")
    print(f"⏱️  时间范围: {TIME_RANGE}")
    print(f"📋 字段列表: {FIELD_LIST}")
    print(f"📊 输出文件: {OUTPUT_FILE}")
    
    # 获取脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_file_path = os.path.join(script_dir, DATA_FILE)
    output_file_path = os.path.join(script_dir, OUTPUT_FILE)
    
    # 检查数据文件是否存在
    if not os.path.exists(data_file_path):
        print(f"❌ 数据文件不存在: {data_file_path}")
        return
    
    # 提取数据
    extracted_data = extract_data_from_exodus(data_file_path)
    
    if extracted_data:
        # 导出到Excel
        if export_to_excel(extracted_data, output_file_path):
            print("\n🎉 数据提取和导出完成!")
            print(f"📂 输出文件位置: {output_file_path}")
        else:
            print("\n❌ 数据导出失败")
    else:
        print("\n❌ 数据提取失败")

if __name__ == "__main__":
    main() 