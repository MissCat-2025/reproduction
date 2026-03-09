#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
导出 Exodus 中的 pellet_total_strain_energy 时间序列。
输出格式：每行 "time value"（空格分隔）。
"""

import os
import sys


def _iter_vtk_blocks(obj):
    # 遍历 MultiBlock 结构，兼容嵌套块
    try:
        import vtk
    except Exception:
        return

    if obj is None:
        return

    if isinstance(obj, vtk.vtkMultiBlockDataSet):
        n = obj.GetNumberOfBlocks()
        for i in range(n):
            b = obj.GetBlock(i)
            if b is None:
                continue
            yield b
            yield from _iter_vtk_blocks(b)


def _get_field_data_scalar(obj, array_name: str):
    # 从 FieldData 读取标量（Postprocessor 写入通常在 FieldData）
    if obj is None:
        return None

    fd = obj.GetFieldData() if hasattr(obj, "GetFieldData") else None
    if fd is None:
        return None

    arr = fd.GetArray(array_name)
    if arr is None:
        return None

    if arr.GetNumberOfTuples() < 1 or arr.GetNumberOfComponents() < 1:
        return None

    return float(arr.GetComponent(0, 0))


def export_pellet_total_strain_energy_from_exodus(exodus_path: str, output_path: str):
    # 通过 ParaView 读取 Exodus 并逐时间步取值
    from paraview.simple import OpenDataFile, UpdatePipeline
    from paraview import servermanager

    reader = OpenDataFile(exodus_path)
    if reader is None:
        raise RuntimeError(f"无法读取文件: {exodus_path}")

    time_steps = []
    if hasattr(reader, "TimestepValues") and reader.TimestepValues:
        time_steps = [float(t) for t in reader.TimestepValues]
    elif hasattr(reader, "TimestepValue"):
        time_steps = [float(reader.TimestepValue)]

    if not time_steps:
        raise RuntimeError("未找到时间步信息")

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)

    # 逐时间步写出
    with open(output_path, "w", encoding="utf-8") as f:
        for t in time_steps:
            UpdatePipeline(time=t, proxy=reader)
            data = servermanager.Fetch(reader)

            value = _get_field_data_scalar(data, "pellet_total_strain_energy")
            if value is None:
                for b in _iter_vtk_blocks(data):
                    value = _get_field_data_scalar(b, "pellet_total_strain_energy")
                    if value is not None:
                        break

            if value is None:
                raise RuntimeError(
                    "在Exodus的FieldData中未找到 'pellet_total_strain_energy'，"
                    "请确认该后处理量确实写入了exodus输出。"
                )

            f.write(f"{t} {value}\n")


def _default_output_path(exodus_path: str):
    # 默认输出文件名与输入 .e 同目录
    base = os.path.splitext(os.path.basename(exodus_path))[0]
    return os.path.join(os.path.dirname(os.path.abspath(exodus_path)), f"{base}_pellet_total_strain_energy.txt")


def main(argv=None):
    # CLI 入口：支持可选输出文件名
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        raise SystemExit(
            "用法: python extract_pellet_total_strain_energy.py <main_out.e> [output.txt]"
        )

    exodus_path = argv[0]
    out_path = argv[1] if len(argv) > 1 else _default_output_path(exodus_path)
    export_pellet_total_strain_energy_from_exodus(exodus_path, out_path)
    print(out_path)


if __name__ == "__main__":
    main()
