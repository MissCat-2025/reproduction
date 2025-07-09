# CT模型 - 混合网格版本 (三角形+四边形)
# 裂纹尖端区域：三角形网格 (TRI3)
# 远场区域：四边形网格 (QUAD4)
# 
# 运行命令：
# conda activate moose && mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i CT_hybrid_mesh.i --mesh-only CT_hybrid_mesh.e

# 几何参数 (单位: m)
L = 20.0e-3      # 总长度
W = 24.0e-3      # 总宽度  
notch_length = 8.5e-3    # 缺口长度
notch_width = 1.0e-3     # 缺口宽度
crack_length = 2.0e-3    # 初始裂纹长度

# 网格细化参数
crack_tip_radius = 3.0e-3  # 裂纹尖端细化区域半径
fine_size = 0.1e-3         # 细网格尺寸
coarse_size = 0.5e-3       # 粗网格尺寸

[Mesh]
  # 步骤1：创建远场四边形网格
  [far_field]
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
    input = far_field
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
  
  # 步骤4：创建裂纹尖端细化区域 (准备转换为三角形)
  [crack_tip_region]
    type = SubdomainBoundingBoxGenerator
    input = remove_notch
    block_id = 2
    bottom_left = '${fparse notch_length - crack_tip_radius/2} ${fparse W/2 - crack_tip_radius/2} 0'
    top_right = '${fparse notch_length + crack_length + crack_tip_radius/2} ${fparse W/2 + crack_tip_radius/2} 0'
    block_name = 'crack_tip_region'
  []
  
  # 步骤5：在裂纹尖端区域细化网格
  [refine_crack_tip]
    type = RefineBlockGenerator
    input = crack_tip_region
    block = '2'
    refinement_level = 2
    enable_neighbor_refinement = true
  []
  
  # 步骤6：将裂纹尖端区域转换为三角形网格
  [triangulate_crack_tip]
    type = TriangleMeshGenerator
    input = refine_crack_tip
    triangulate_block = '2'
    desired_area = 1.0e-8  # 控制三角形大小
  []
  
  # 步骤7：标记裂纹尖端点（用于后续应力集中分析）
  [mark_crack_tip]
    type = BoundingBoxNodeSetGenerator
    input = triangulate_crack_tip
    new_boundary = 'crack_tip'
    bottom_left = '${fparse notch_length + crack_length - 0.1e-3} ${fparse W/2 - 0.1e-3} 0'
    top_right = '${fparse notch_length + crack_length + 0.1e-3} ${fparse W/2 + 0.1e-3} 0'
  []
  
  # 步骤8：创建载荷点（模拟圆孔位置）
  [upper_load_point]
    type = BoundingBoxNodeSetGenerator
    input = mark_crack_tip
    new_boundary = 'upper_load'
    bottom_left = '${fparse 4.0e-3 - 0.5e-3} ${fparse W/2 + 6.0e-3 - 0.5e-3} 0'
    top_right = '${fparse 4.0e-3 + 0.5e-3} ${fparse W/2 + 6.0e-3 + 0.5e-3} 0'
  []
  
  [lower_load_point]
    type = BoundingBoxNodeSetGenerator
    input = upper_load_point
    new_boundary = 'lower_load'
    bottom_left = '${fparse 4.0e-3 - 0.5e-3} ${fparse W/2 - 6.0e-3 - 0.5e-3} 0'
    top_right = '${fparse 4.0e-3 + 0.5e-3} ${fparse W/2 - 6.0e-3 + 0.5e-3} 0'
  []
  
  # 步骤9：重新标记块ID
  [rename_blocks]
    type = RenameBlockGenerator
    input = lower_load_point
    old_blocks = '0 2'
    new_blocks = '1 2'  # 1=四边形远场, 2=三角形裂纹尖端
    old_block_names = 'any any'
    new_block_names = 'far_field crack_tip'
  []
  
  # 步骤10：最终优化网格
  [final_mesh]
    type = MeshExtruderGenerator
    input = rename_blocks
    extrusion_vector = '0 0 0'  # 不挤出，只是为了触发网格优化
    num_layers = 0
  []
[]

# 可选：输出网格信息
[Problem]
  type = MooseProblem
  solve = false  # 只生成网格，不求解
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
  []
[] 