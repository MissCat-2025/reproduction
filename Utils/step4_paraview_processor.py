#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step4：ParaView 单例后处理脚本。
- 扫描 case_*，优先依据 file_base 定位 .e 文件
- 为每个字段导出目标时间截图
- 输出目录固定在 case 根目录下的 post_results
"""

import os
import re
import sys
import subprocess
from collections import namedtuple

# 配置参数
ExportConfig = namedtuple('ExportConfig', [
    'target_times',    # 目标时间列表 [4.0, 5.0, 6.0]
    'field_list',      # 要导出的字段配置 [('T', '温度'), ('stress', '应力')]
    'image_size',      # 图像分辨率 [1642, 1083]
    'output_dir_name'  # 输出目录名称后缀
])

def _parse_target_times():
    v = os.environ.get("TARGET_TIMES")
    if v:
        try:
            return [float(x) for x in v.split(",") if x.strip()]
        except ValueError:
            return [60000, 125000, 175000, 300000]
    return [60000, 125000, 175000, 300000]

def _parse_fields(env_value, default_fields):
    if not env_value:
        return default_fields
    fields = []
    for part in env_value.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            name, label = part.split(":", 1)
            fields.append((name.strip(), label.strip()))
        else:
            name = part.strip()
            fields.append((name, name))
    return fields or default_fields

def _parse_image_size():
    v = os.environ.get("PV_IMAGE_SIZE")
    if not v:
        return (1642, 1083)
    sep = "," if "," in v else "x"
    try:
        w, h = v.split(sep)
        return (int(w), int(h))
    except Exception:
        return (1642, 1083)

def _find_first_input(case_dir):
    for name in sorted(os.listdir(case_dir)):
        if name.endswith(".i"):
            return os.path.join(case_dir, name)
    return None

def _parse_file_base(input_path):
    if not input_path or not os.path.isfile(input_path):
        return None
    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "#" in line:
                line = line.split("#", 1)[0]
            if "file_base" not in line:
                continue
            m = re.search(r"file_base\s*=\s*['\"]([^'\"]+)['\"]", line)
            if m:
                return m.group(1).strip()
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
        for file in files:
            if not file.endswith(".e"):
                continue
            if prefix and not file.startswith(prefix):
                continue
            matches.append(os.path.join(root, file))
    if not matches:
        return None
    return max(matches, key=os.path.getmtime)

def _find_case_exodus(case_dir, output_dir_name):
    # 优先使用输入文件中的 file_base 进行精确定位
    input_path = _find_first_input(case_dir)
    file_base = _parse_file_base(input_path)
    exodus_path = _find_exodus_by_file_base(case_dir, file_base)
    if exodus_path:
        return exodus_path, file_base
    candidates = []
    for root, _, files in os.walk(case_dir):
        if output_dir_name and os.path.sep + output_dir_name in root:
            continue
        for file in files:
            if file.endswith(".e"):
                candidates.append(os.path.join(root, file))
    if not candidates:
        return None, file_base
    return max(candidates, key=os.path.getmtime), file_base

_default_fields_single = [
    ('d', '相场变量'),
    ('hoop_stress', '环向应力'),
]
_env_fields_single = os.environ.get("PV_FIELDS")
_fields_single = _parse_fields(_env_fields_single, _default_fields_single)
_image_size_single = _parse_image_size()
_output_dir_single = os.environ.get("PV_OUTPUT_DIR", "post_results")

DEFAULT_CONFIG = ExportConfig(
    target_times=_parse_target_times(),
    field_list=_fields_single,
    image_size=_image_size_single,
    output_dir_name=_output_dir_single
)

def export_field_data(file_path, config=DEFAULT_CONFIG, case_dir=None):
    """导出多字段数据"""
    try:
        from paraview.simple import (
            OpenDataFile, GetActiveViewOrCreate, GetDisplayProperties,
            ColorBy, GetColorTransferFunction, GetOpacityTransferFunction,
            GetAnimationScene, GetLayout, SaveScreenshot, UpdatePipeline,
            RenderAllViews,
            GetSettingsProxy, HideAll
        )
        import paraview.simple
        paraview.simple._DisableFirstRenderCameraReset()
    except ImportError:
        print("错误：无法导入ParaView模块")
        return False

    try:
        # 1. 读取文件
        print(f"\n处理文件: {os.path.basename(file_path)}")
        reader = OpenDataFile(file_path)
        if not reader:
            print("无法读取文件")
            return False

    # 2. 获取时间步列表
        time_steps = []
        if hasattr(reader, 'TimestepValues'):
            time_steps = [float(t) for t in reader.TimestepValues]
        elif hasattr(reader, 'TimestepValue'):
            time_steps = [float(reader.TimestepValue)]
        
        if not time_steps:
            print("警告: 未找到时间步信息")
            return False

        # 3. 创建输出目录（固定在 case 根目录）
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        if case_dir:
            output_root = os.path.join(case_dir, config.output_dir_name)
        else:
            output_root = os.path.join(os.path.dirname(file_path), config.output_dir_name)
        os.makedirs(output_root, exist_ok=True)

        # 为每个字段创建子目录
        field_dirs = {}
        for field, display_name in config.field_list:
            dir_path = os.path.join(output_root, f"{field}_images")
            os.makedirs(dir_path, exist_ok=True)
            field_dirs[field] = dir_path

        # 4. 设置视图参数（背景、相机、布局）
        renderView1 = GetActiveViewOrCreate('RenderView')
        layout1 = GetLayout()
        layout1.SetSize(*config.image_size)
        
        # 设置背景和文字颜色
        colorPalette = GetSettingsProxy('ColorPalette')
        colorPalette.Background = [1.0, 1.0, 1.0]  # 白色背景
        colorPalette.Text = [0.0, 0.0, 0.0]  # 黑色文字
        
        renderView1.CameraPosition = [0.0, 0.0, 0.03999187526517147]
        renderView1.CameraFocalPoint = [0.0, 0.0, 3e-05]
        renderView1.CameraParallelScale = 0.010342894396637723
        renderView1.OrientationAxesVisibility = 0
        # 5. 处理每个字段：设置颜色条并输出截图
        for field, display_name in config.field_list:
            print(f"\n正在处理字段: {display_name} ({field})")
            
            try:
                # 隐藏所有已有显示
                HideAll()
                # paraview.simple.HideAll(reader, renderView1)
                
                # 创建新的显示对象
                display = paraview.simple.Show(reader, renderView1)
                display.Representation = 'Surface'

                # 获取所有可用字段（包括点、单元和场数据）
                def get_all_fields(reader):
                    fields = []
                    # 检查所有可能的数据属性
                    for attr in ['PointData', 'CellData', 'FieldData']:
                        if hasattr(reader, attr):
                            data = getattr(reader, attr)
                            fields.extend([(arr.GetName(), attr) for arr in data])
                    return fields

                # 获取所有字段及其存储位置
                all_fields = get_all_fields(reader)
                available_fields = [f[0] for f in all_fields]
                print(f"可用字段列表: {available_fields}")

                # 查找字段实际存储名称（精确匹配）
                actual_field = next(
                    (f for f in all_fields if f[0] == field),  # 改为精确匹配
                    None
                )

                if not actual_field:
                    print(f"  警告: 字段 {field} 不存在，可用字段：{available_fields}")
                    continue

                # 获取数据数组
                data_attr = getattr(reader, actual_field[1])  # PointData/CellData
                data_array = data_attr[actual_field[0]]
                
                # 检查数据有效性
                if not data_array:
                    print(f"  警告: 字段 {field} 数据为空")
                    continue

                # 设置字段显示（使用新的display对象）
                association = 'POINTS' if actual_field[1] == 'PointData' else 'CELLS'
                ColorBy(display, (association, actual_field[0]))
                renderView1.Update()

                # 删除旧的颜色条（关键修改）
                display.SetScalarBarVisibility(renderView1, False)  # 先隐藏旧颜色条
                display.SetScalarBarVisibility(renderView1, True)   # 显示新颜色条

                # 安全获取数据范围
                try:
                    data_range = data_array.GetRange()
                    print(f"  数据范围: {data_range[0]:.2e} 至 {data_range[1]:.2e}")
                except Exception as e:
                    print(f"  无法获取数据范围: {str(e)}")
                    data_range = [0, 1]  # 使用默认范围

                # 设置颜色映射
                lut = GetColorTransferFunction(actual_field[0])
                pov = GetOpacityTransferFunction(actual_field[0])
                
                # 应用预设的颜色方案
                lut.ApplyPreset('Blue to Red Rainbow', True)
                lut.NumberOfTableValues = 20  # 将颜色表的值数量减少到20
                
                # 显示颜色条与字体样式
                display.SetScalarBarVisibility(renderView1, True)
                scalar_bar = paraview.simple.GetScalarBar(lut, renderView1)
                scalar_bar.Orientation = 'Horizontal'
                scalar_bar.WindowLocation = 'Any Location'
                scalar_bar.ScalarBarThickness = 16
                scalar_bar.ScalarBarLength = 0.146019
                scalar_bar.Title = actual_field[0]
                scalar_bar.TitleFontFamily = 'Times'
                scalar_bar.LabelFontFamily = 'Times'
                scalar_bar.TitleFontSize = 40
                scalar_bar.LabelFontSize = 40
                scalar_bar.TitleBold = 1
                scalar_bar.LabelBold = 1
                scalar_bar.TitleOpacity = 1.0
                scalar_bar.LabelOpacity = 1.0
                scalar_bar.Position = [0.4124640542898518, -0.1338517840805124]

                # 强制更新视图
                renderView1.ResetCamera()
                renderView1.Update()
                RenderAllViews()

            except Exception as e:
                print(f"  字段设置失败: {str(e)}")
                continue

            # 6. 处理每个时间点并写图像
            expected_outputs = []
            closest_times = []
            for target_time in config.target_times:
                closest_time = min(time_steps, key=lambda x: abs(x - target_time))
                closest_times.append(closest_time)
                time_str = f"{closest_time:.2f}s".replace('.', 'p')
                expected_outputs.append(
                    os.path.join(field_dirs[field], f'{base_name}_{field}_{time_str}.png')
                )
            if expected_outputs and all(os.path.isfile(p) for p in expected_outputs):
                print("  已存在目标时间截图，跳过该字段")
                continue

            for target_time, closest_time in zip(config.target_times, closest_times):
                print(f"  时间 {target_time}s → 使用 {closest_time:.2f}s")

                # 设置时间步
                animationScene1 = GetAnimationScene()
                animationScene1.AnimationTime = float(closest_time)

                # 强制更新数据
                UpdatePipeline(time=closest_time, proxy=reader)
                
                # 获取当前显示属性
                display = GetDisplayProperties(reader, view=renderView1)
                
                # 重新缩放传输函数（关键修改点）
                display.RescaleTransferFunctionToDataRange(False, True)  # 参数说明：(useVisibleRange=False, clamp=False)
                
                # 强制更新视图
                renderView1.Update()
                renderView1.ResetCamera()
                scalar_bar = paraview.simple.GetScalarBar(lut, renderView1)
                scalar_bar.WindowLocation = 'Any Location'
                scalar_bar.Position = [0.4124640542898518, 0.00338517840805124]
                RenderAllViews()

                # 生成输出路径
                time_str = f"{closest_time:.2f}s".replace('.', 'p')
                output_file = os.path.join(field_dirs[field], f'{base_name}_{field}_{time_str}.png')
                if os.path.isfile(output_file):
                    print(f"    已存在: {os.path.basename(output_file)}")
                    continue

                # 保存截图
                SaveScreenshot(
                    filename=output_file,
                    viewOrLayout=renderView1,
                    ImageResolution=config.image_size
                )
                print(f"    已保存: {os.path.basename(output_file)}")

        return True

    except Exception as e:
        print(f"\n错误: {str(e)}")
        return False

def main():
    """主函数"""
    script_dir = os.environ.get("PROJECT_BASE_DIR", os.path.dirname(os.path.abspath(__file__)))
    studies_subdir = os.environ.get("STUDIES_SUBDIR", "2Good-Gc2.25-3.5")
    studies_dir = os.path.join(script_dir, studies_subdir)
    
    if not os.path.exists(studies_dir):
        print("数据目录不存在")
        return
        
    output_dir_name = os.environ.get("PV_OUTPUT_DIR", "post_results")
    case_dirs = [
        os.path.join(studies_dir, name)
        for name in os.listdir(studies_dir)
        if name.startswith("case_") and os.path.isdir(os.path.join(studies_dir, name))
    ]
    case_dirs.sort()
    if not case_dirs:
        print("未找到case目录，改为全局扫描.e文件")
        e_files = []
        for root, _, files in os.walk(studies_dir):
            for file in files:
                if file.endswith(".e"):
                    e_files.append(os.path.join(root, file))
        if not e_files:
            print("未找到.e文件")
            return
        success_count = 0
        for idx, e_file in enumerate(e_files):
            print(f"\n处理文件 {idx+1}/{len(e_files)}: {os.path.basename(e_file)}")
            if export_field_data(e_file):
                success_count += 1
        print(f"\n处理完成：成功{success_count}个，失败{len(e_files)-success_count}个")
        return

    print(f"找到{len(case_dirs)}个case目录")
    success_count = 0
    for idx, case_dir in enumerate(case_dirs, start=1):
        case_name = os.path.basename(case_dir)
        e_file, file_base = _find_case_exodus(case_dir, output_dir_name)
        if not e_file:
            print(f"\n[{idx}/{len(case_dirs)}] {case_name}: 未找到.e文件")
            continue
        print(f"\n[{idx}/{len(case_dirs)}] {case_name}")
        if file_base:
            print(f"  file_base: {file_base}")
        print(f"  exodus: {os.path.relpath(e_file, studies_dir)}")
        if export_field_data(e_file, case_dir=case_dir):
            success_count += 1

    print(f"\n处理完成：成功{success_count}个，失败{len(case_dirs)-success_count}个")

if __name__ == "__main__":
    main()
