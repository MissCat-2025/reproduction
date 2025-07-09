#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成Compact Tension (CT) 试样的Gmsh网格
基于ASTM E399标准几何

使用方法:
python generate_CT_mesh.py

需要安装: pip install gmsh
"""

import gmsh
import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

def create_CT_geometry():
    """创建CT试样几何"""
    
    # 几何参数 (单位: mm)
    W = 24.0      # 总宽度
    L = 25.0      # 总长度
    B = 7.1       # 厚度
    
    # 孔参数
    hole_diameter = 6.0
    hole_radius = hole_diameter / 2.0
    hole_spacing = 12.0  # 两孔间距
    
    # 孔的位置
    hole_x = 4.0
    hole_y1 = W/2 - hole_spacing/2  # 下孔
    hole_y2 = W/2 + hole_spacing/2  # 上孔
    
    # 缺口参数
    notch_length = 8.5
    notch_width = 1.0
    
    # 初始裂纹长度
    a0 = 2.0
    
    # 网格尺寸
    lc_fine = 0.2    # 裂纹前端细网格
    lc_coarse = 1.0  # 远场粗网格
    
    gmsh.initialize()
    gmsh.model.add("CT_specimen")
    
    # 创建外轮廓点
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc_coarse)
    p2 = gmsh.model.geo.addPoint(L, 0, 0, lc_coarse)
    p3 = gmsh.model.geo.addPoint(L, W, 0, lc_coarse)
    p4 = gmsh.model.geo.addPoint(0, W, 0, lc_coarse)
    
    # 创建外轮廓线
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    
    # 创建外轮廓环
    outer_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    
    # 创建上部圆孔
    hole_center1 = gmsh.model.geo.addPoint(hole_x, hole_y1, 0, lc_fine)
    hole_arc1_start = gmsh.model.geo.addPoint(hole_x + hole_radius, hole_y1, 0, lc_fine)
    hole_arc1_end = gmsh.model.geo.addPoint(hole_x - hole_radius, hole_y1, 0, lc_fine)
    
    hole_arc1_1 = gmsh.model.geo.addCircleArc(hole_arc1_start, hole_center1, hole_arc1_end)
    hole_arc1_2 = gmsh.model.geo.addCircleArc(hole_arc1_end, hole_center1, hole_arc1_start)
    
    hole_loop1 = gmsh.model.geo.addCurveLoop([hole_arc1_1, hole_arc1_2])
    
    # 创建下部圆孔
    hole_center2 = gmsh.model.geo.addPoint(hole_x, hole_y2, 0, lc_fine)
    hole_arc2_start = gmsh.model.geo.addPoint(hole_x + hole_radius, hole_y2, 0, lc_fine)
    hole_arc2_end = gmsh.model.geo.addPoint(hole_x - hole_radius, hole_y2, 0, lc_fine)
    
    hole_arc2_1 = gmsh.model.geo.addCircleArc(hole_arc2_start, hole_center2, hole_arc2_end)
    hole_arc2_2 = gmsh.model.geo.addCircleArc(hole_arc2_end, hole_center2, hole_arc2_start)
    
    hole_loop2 = gmsh.model.geo.addCurveLoop([hole_arc2_1, hole_arc2_2])
    
    # 创建缺口
    notch_p1 = gmsh.model.geo.addPoint(0, W/2 - notch_width/2, 0, lc_fine)
    notch_p2 = gmsh.model.geo.addPoint(notch_length, W/2 - notch_width/2, 0, lc_fine)
    notch_p3 = gmsh.model.geo.addPoint(notch_length, W/2 + notch_width/2, 0, lc_fine)
    notch_p4 = gmsh.model.geo.addPoint(0, W/2 + notch_width/2, 0, lc_fine)
    
    notch_l1 = gmsh.model.geo.addLine(notch_p1, notch_p2)
    notch_l2 = gmsh.model.geo.addLine(notch_p2, notch_p3)
    notch_l3 = gmsh.model.geo.addLine(notch_p3, notch_p4)
    notch_l4 = gmsh.model.geo.addLine(notch_p4, notch_p1)
    
    notch_loop = gmsh.model.geo.addCurveLoop([notch_l1, notch_l2, notch_l3, notch_l4])
    
    # 创建初始裂纹
    crack_p1 = gmsh.model.geo.addPoint(notch_length, W/2 - 0.05, 0, lc_fine)
    crack_p2 = gmsh.model.geo.addPoint(notch_length + a0, W/2 - 0.05, 0, lc_fine)
    crack_p3 = gmsh.model.geo.addPoint(notch_length + a0, W/2 + 0.05, 0, lc_fine)
    crack_p4 = gmsh.model.geo.addPoint(notch_length, W/2 + 0.05, 0, lc_fine)
    
    crack_l1 = gmsh.model.geo.addLine(crack_p1, crack_p2)
    crack_l2 = gmsh.model.geo.addLine(crack_p2, crack_p3)
    crack_l3 = gmsh.model.geo.addLine(crack_p3, crack_p4)
    crack_l4 = gmsh.model.geo.addLine(crack_p4, crack_p1)
    
    crack_loop = gmsh.model.geo.addCurveLoop([crack_l1, crack_l2, crack_l3, crack_l4])
    
    # 创建面（减去孔、缺口和裂纹）
    surface = gmsh.model.geo.addPlaneSurface([outer_loop, hole_loop1, hole_loop2, notch_loop, crack_loop])
    
    # 同步几何
    gmsh.model.geo.synchronize()
    
    # 添加物理组
    gmsh.model.addPhysicalGroup(2, [surface], 1, "specimen")
    gmsh.model.addPhysicalGroup(1, [hole_arc1_1, hole_arc1_2], 2, "upper_hole")
    gmsh.model.addPhysicalGroup(1, [hole_arc2_1, hole_arc2_2], 3, "lower_hole")
    gmsh.model.addPhysicalGroup(1, [crack_l1, crack_l3], 4, "crack_faces")
    gmsh.model.addPhysicalGroup(1, [crack_l2], 5, "crack_tip")
    
    # 在裂纹尖端添加细化点
    crack_tip_point = gmsh.model.geo.addPoint(notch_length + a0, W/2, 0, lc_fine/5)
    gmsh.model.geo.synchronize()
    
    # 生成网格
    gmsh.model.mesh.generate(2)
    
    # 细化裂纹尖端
    gmsh.model.mesh.refine()
    
    # 写入网格文件
    gmsh.write("CT_specimen.msh")
    
    print("✓ CT试样网格已生成: CT_specimen.msh")
    print(f"  - 试样尺寸: {L} × {W} mm")
    print(f"  - 孔径: {hole_diameter} mm")
    print(f"  - 缺口长度: {notch_length} mm")
    print(f"  - 初始裂纹长度: {a0} mm")
    
    # 显示网格
    if '-nopopup' not in gmsh.sys.argv:
        gmsh.fltk.run()
    
    gmsh.finalize()

def visualize_CT_geometry():
    """可视化CT试样几何"""
    
    # 几何参数
    W = 24.0
    L = 25.0
    hole_diameter = 6.0
    hole_radius = hole_diameter / 2.0
    hole_x = 4.0
    hole_y1 = W/2 - 6.0
    hole_y2 = W/2 + 6.0
    notch_length = 8.5
    notch_width = 1.0
    a0 = 2.0
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # 绘制外轮廓
    rect = plt.Rectangle((0, 0), L, W, linewidth=2, edgecolor='black', facecolor='lightgray', alpha=0.7)
    ax.add_patch(rect)
    
    # 绘制圆孔
    hole1 = plt.Circle((hole_x, hole_y1), hole_radius, facecolor='white', edgecolor='blue', linewidth=2)
    hole2 = plt.Circle((hole_x, hole_y2), hole_radius, facecolor='white', edgecolor='blue', linewidth=2)
    ax.add_patch(hole1)
    ax.add_patch(hole2)
    
    # 绘制缺口
    notch = plt.Rectangle((0, W/2 - notch_width/2), notch_length, notch_width, 
                         facecolor='white', edgecolor='red', linewidth=2)
    ax.add_patch(notch)
    
    # 绘制初始裂纹
    crack = plt.Rectangle((notch_length, W/2 - 0.05), a0, 0.1, 
                         facecolor='red', edgecolor='red', linewidth=2)
    ax.add_patch(crack)
    
    # 添加载荷箭头
    ax.annotate('', xy=(hole_x, hole_y1-hole_radius-2), xytext=(hole_x, hole_y1-hole_radius-1),
                arrowprops=dict(arrowstyle='->', color='green', lw=3))
    ax.annotate('', xy=(hole_x, hole_y2+hole_radius+2), xytext=(hole_x, hole_y2+hole_radius+1),
                arrowprops=dict(arrowstyle='->', color='green', lw=3))
    
    # 添加标注
    ax.text(hole_x-1, hole_y1-hole_radius-3, 'F', fontsize=14, ha='center', color='green', weight='bold')
    ax.text(hole_x-1, hole_y2+hole_radius+3, 'F', fontsize=14, ha='center', color='green', weight='bold')
    
    # 添加尺寸标注
    ax.text(L/2, -2, f'L = {L} mm', fontsize=12, ha='center')
    ax.text(-2, W/2, f'W = {W} mm', fontsize=12, ha='center', rotation=90)
    ax.text(notch_length/2, W/2-3, f'缺口长度 = {notch_length} mm', fontsize=10, ha='center')
    ax.text(notch_length + a0/2, W/2+1, f'a₀ = {a0} mm', fontsize=10, ha='center', color='red')
    
    ax.set_xlim(-3, L+3)
    ax.set_ylim(-5, W+5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_title('Compact Tension (CT) 试样几何', fontsize=14, weight='bold')
    
    plt.tight_layout()
    plt.savefig('CT_geometry.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✓ 几何图已保存: CT_geometry.png")

if __name__ == "__main__":
    print("正在生成CT试样几何...")
    
    # 可视化几何
    visualize_CT_geometry()
    
    # 生成网格
    create_CT_geometry()
    
    print("\n使用说明:")
    print("1. 在MOOSE输入文件中使用 FileMeshGenerator 读取 CT_specimen.msh")
    print("2. 使用物理组施加边界条件:")
    print("   - 物理组2 (upper_hole): 上部加载孔")
    print("   - 物理组3 (lower_hole): 下部加载孔")
    print("   - 物理组4 (crack_faces): 裂纹面")
    print("   - 物理组5 (crack_tip): 裂纹尖端") 