# CT模型 - 实用混合网格版本
# 使用XYDelaunayGenerator在裂纹尖端创建三角形网格
# 远场保持四边形网格
# 
# 运行命令：
# conda activate moose && mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i CT_hybrid_practical.i --mesh-only CT_hybrid_practical.e

# 几何参数 (单位: m)
L = 20.0e-3      # 总长度
W = 24.0e-3      # 总宽度  
notch_length = 8.5e-3    # 缺口长度
notch_width = 1.0e-3     # 缺口宽度
crack_length = 2.0e-3    # 初始裂纹长度

# 网格参数
crack_tip_x = '${fparse notch_length + crack_length}'  # 裂纹尖端X坐标
crack_tip_y = '${fparse W/2}'                          # 裂纹尖端Y坐标
refine_radius = 2.0e-3                                 # 细化半径

[Mesh]
  # 步骤1：创建远场四边形网格基础
  [base_quad]
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
    input = base_quad
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
  
  # 步骤4：定义裂纹尖端区域（将要三角化）
  [crack_tip_region]
    type = SubdomainBoundingBoxGenerator
    input = remove_notch
    block_id = 2
    bottom_left = '${fparse crack_tip_x - refine_radius} ${fparse crack_tip_y - refine_radius} 0'
    top_right = '${fparse crack_tip_x + refine_radius} ${fparse crack_tip_y + refine_radius} 0'
    block_name = 'crack_tip_region'
  []
  
  # 步骤5：先细化裂纹尖端区域
  [refine_crack]
    type = RefineBlockGenerator
    input = crack_tip_region
    block = '2'
    refinement_level = 2
    enable_neighbor_refinement = true
  []
  
  # 步骤6：创建裂纹尖端的三角形网格
  # 首先提取裂纹尖端区域的节点
  [extract_crack_nodes]
    type = BoundingBoxNodeSetGenerator
    input = refine_crack
    new_boundary = 'crack_nodes'
    bottom_left = '${fparse crack_tip_x - refine_radius} ${fparse crack_tip_y - refine_radius} 0'
    top_right = '${fparse crack_tip_x + refine_radius} ${fparse crack_tip_y + refine_radius} 0'
  []
  
  # 步骤7：在裂纹尖端附近进行更精细的三角化
  [triangulate_region]
    type = XYDelaunayGenerator
    input = extract_crack_nodes
    boundary = 'crack_nodes'
    stitch_holes = true
    desired_area = 5.0e-9  # 控制三角形大小
    smooth_triangulation = true
  []
  
  # 步骤8：标记重要边界
  [mark_crack_tip]
    type = BoundingBoxNodeSetGenerator
    input = triangulate_region
    new_boundary = 'crack_tip'
    bottom_left = '${fparse crack_tip_x - 0.05e-3} ${fparse crack_tip_y - 0.05e-3} 0'
    top_right = '${fparse crack_tip_x + 0.05e-3} ${fparse crack_tip_y + 0.05e-3} 0'
  []
  
  # 步骤9：创建载荷点
  [upper_load_point]
    type = BoundingBoxNodeSetGenerator
    input = mark_crack_tip
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
  
  # 步骤10：最终块重命名
  [rename_blocks]
    type = RenameBlockGenerator
    input = lower_load_point
    old_blocks = '0 2'
    new_blocks = '1 2'
    old_block_names = 'any any'
    new_block_names = 'far_field crack_tip'
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