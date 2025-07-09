# CT模型 - 自适应细化版本 (模拟混合网格效果)
# 通过多级细化在裂纹尖端创建极细的四边形网格
# 配合自适应细化可以达到类似三角形网格的效果
# 
# 运行命令：
# conda activate moose && mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i CT_adaptive_refined.i --mesh-only CT_adaptive_refined.e

# 几何参数 (单位: m)
L = 20.0e-3      # 总长度
W = 24.0e-3      # 总宽度  
notch_length = 8.5e-3    # 缺口长度
notch_width = 1.0e-3     # 缺口宽度
crack_length = 2.0e-3    # 初始裂纹长度

# 细化参数
crack_tip_x = '${fparse notch_length + crack_length}'  # 裂纹尖端X坐标
crack_tip_y = '${fparse W/2}'                          # 裂纹尖端Y坐标

[Mesh]
  # 步骤1：创建基础四边形网格
  [base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = ${W}
    nx = 40
    ny = 48
    elem_type = QUAD4
  []
  
  # 步骤2：创建缺口区域
  [notch_subdomain]
    type = SubdomainBoundingBoxGenerator
    input = base_mesh
    block_id = 1
    bottom_left = '0 ${fparse W/2 - notch_width/2} 0'
    top_right = '${notch_length} ${fparse W/2 + notch_width/2} 0'
    block_name = 'notch_region'
  []
  
  # 步骤3：删除缺口区域
  [remove_notch]
    type = BlockDeletionGenerator
    input = notch_subdomain
    block = '1'
  []
  
  # 步骤4：定义裂纹尖端大区域（第一级细化）
  [crack_region_large]
    type = SubdomainBoundingBoxGenerator
    input = remove_notch
    block_id = 2
    bottom_left = '${fparse crack_tip_x - 3.0e-3} ${fparse crack_tip_y - 3.0e-3} 0'
    top_right = '${fparse crack_tip_x + 1.0e-3} ${fparse crack_tip_y + 3.0e-3} 0'
    block_name = 'crack_region_large'
  []
  
  # 步骤5：第一级细化（大区域）
  [refine_level_1]
    type = RefineBlockGenerator
    input = crack_region_large
    block = '2'
    refinement_level = 1
    enable_neighbor_refinement = true
  []
  
  # 步骤6：定义裂纹尖端中区域（第二级细化）
  [crack_region_medium]
    type = SubdomainBoundingBoxGenerator
    input = refine_level_1
    block_id = 3
    bottom_left = '${fparse crack_tip_x - 1.5e-3} ${fparse crack_tip_y - 1.5e-3} 0'
    top_right = '${fparse crack_tip_x + 0.5e-3} ${fparse crack_tip_y + 1.5e-3} 0'
    block_name = 'crack_region_medium'
  []
  
  # 步骤7：第二级细化（中区域）
  [refine_level_2]
    type = RefineBlockGenerator
    input = crack_region_medium
    block = '3'
    refinement_level = 2
    enable_neighbor_refinement = true
  []
  
  # 步骤8：定义裂纹尖端小区域（第三级细化）
  [crack_region_small]
    type = SubdomainBoundingBoxGenerator
    input = refine_level_2
    block_id = 4
    bottom_left = '${fparse crack_tip_x - 0.5e-3} ${fparse crack_tip_y - 0.5e-3} 0'
    top_right = '${fparse crack_tip_x + 0.2e-3} ${fparse crack_tip_y + 0.5e-3} 0'
    block_name = 'crack_region_small'
  []
  
  # 步骤9：第三级细化（小区域）
  [refine_level_3]
    type = RefineBlockGenerator
    input = crack_region_small
    block = '4'
    refinement_level = 3
    enable_neighbor_refinement = true
  []
  
  # 步骤10：标记裂纹尖端点
  [mark_crack_tip]
    type = BoundingBoxNodeSetGenerator
    input = refine_level_3
    new_boundary = 'crack_tip'
    bottom_left = '${fparse crack_tip_x - 0.02e-3} ${fparse crack_tip_y - 0.02e-3} 0'
    top_right = '${fparse crack_tip_x + 0.02e-3} ${fparse crack_tip_y + 0.02e-3} 0'
  []
  
  # 步骤11：标记裂纹面
  [mark_crack_faces]
    type = BoundingBoxNodeSetGenerator
    input = mark_crack_tip
    new_boundary = 'crack_faces'
    bottom_left = '${notch_length} ${fparse crack_tip_y - 0.1e-3} 0'
    top_right = '${crack_tip_x} ${fparse crack_tip_y + 0.1e-3} 0'
  []
  
  # 步骤12：创建载荷点
  [upper_load_point]
    type = BoundingBoxNodeSetGenerator
    input = mark_crack_faces
    new_boundary = 'upper_load'
    bottom_left = '${fparse 4.0e-3 - 0.3e-3} ${fparse W/2 + 6.0e-3 - 0.3e-3} 0'
    top_right = '${fparse 4.0e-3 + 0.3e-3} ${fparse W/2 + 6.0e-3 + 0.3e-3} 0'
  []
  
  [lower_load_point]
    type = BoundingBoxNodeSetGenerator
    input = upper_load_point
    new_boundary = 'lower_load'
    bottom_left = '${fparse 4.0e-3 - 0.3e-3} ${fparse W/2 - 6.0e-3 - 0.3e-3} 0'
    top_right = '${fparse 4.0e-3 + 0.3e-3} ${fparse W/2 - 6.0e-3 + 0.3e-3} 0'
  []
  
  # 步骤13：重新组织块ID
  [final_blocks]
    type = RenameBlockGenerator
    input = lower_load_point
    old_blocks = '0 2 3 4'
    new_blocks = '1 2 3 4'
    old_block_names = 'any any any any'
    new_block_names = 'far_field coarse_refined medium_refined fine_refined'
  []
  
  # 步骤14：添加自适应细化（可选）
  [adaptive_refine]
    type = RefineBlockGenerator
    input = final_blocks
    block = '4'  # 只对最细的块进行额外细化
    refinement_level = 1
    enable_neighbor_refinement = false
  []
[]

# 只生成网格，不求解
[Problem]
  type = MooseProblem
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
  [mesh_info]
    type = Console
    execute_on = 'INITIAL'
    output_system_information = true
  []
[]

[Postprocessors]
  # 输出网格统计信息
  [num_nodes]
    type = NumNodes
    execute_on = 'INITIAL'
  []
  [num_elems]
    type = NumElems
    execute_on = 'INITIAL'
  []
  [mesh_volume]
    type = VolumePostprocessor
    execute_on = 'INITIAL'
  []
[] 