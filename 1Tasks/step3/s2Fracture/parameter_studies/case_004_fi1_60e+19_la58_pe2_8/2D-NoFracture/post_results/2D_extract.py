#!/usr/bin/env python3
import os
import paraview.simple as pv

# 设置环境变量
os.environ['DISPLAY'] = ':99'
os.environ['MESA_GL_VERSION_OVERRIDE'] = '3.3'
os.environ['MESA_GLSL_VERSION_OVERRIDE'] = '330'
os.environ['LIBGL_ALWAYS_SOFTWARE'] = '1'

try:
    # 读取文件
    reader = pv.OpenDataFile('/home/yp/projects/reproduction/1Tasks/step3/s2Fracture/parameter_studies/case_004_fi1_60e+19_la58_pe2_8/2D-NoFracture/2D.e')
    
    # 获取数据信息
    data_info = reader.GetDataInformation()
    num_points = data_info.GetNumberOfPoints()
    num_cells = data_info.GetNumberOfCells()
    
    print(f"数据点数: {num_points}")
    print(f"数据单元数: {num_cells}")
    
    # 处理字段
    fields = ['d', 'hoop_stress']
    for field in fields:
        try:
            # 获取字段数据
            point_data = reader.PointData
            if hasattr(point_data, field):
                data_array = point_data[field]
                data_range = data_array.GetRange()
                print(f"字段 {field} 范围: {data_range[0]:.2e} 至 {data_range[1]:.2e}")
                
                # 保存数据信息
                info_file = os.path.join('/home/yp/projects/reproduction/1Tasks/step3/s2Fracture/parameter_studies/case_004_fi1_60e+19_la58_pe2_8/2D-NoFracture/post_results', f'2D_{field}_info.txt')
                with open(info_file, 'w') as f:
                    f.write(f"字段: {field}\n")
                    f.write(f"时间: 当前时间步\n")
                    f.write(f"数据范围: {data_range[0]:.2e} 至 {data_range[1]:.2e}\n")
                    f.write(f"数据点数: {num_points}\n")
                    f.write(f"数据单元数: {num_cells}\n")
                print(f"已保存数据信息: {info_file}")
            else:
                print(f"字段 {field} 不存在")
                
        except Exception as e:
            print(f"处理字段 {field} 失败: {e}")
    
    # 清理
    pv.Delete(reader)
    print("数据提取完成")
    
except Exception as e:
    print(f"数据提取失败: {e}")
    sys.exit(1)
