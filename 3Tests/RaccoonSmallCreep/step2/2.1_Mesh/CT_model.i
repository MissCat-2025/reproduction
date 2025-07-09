# 简化CT模型 - 实用版本
# 主要简化：
# 1. 将圆孔简化为在相应位置施加载荷
# 2. 缺口简化为矩形
# 3. 专注于裂纹扩展行为
# conda activate moose && mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i CT_model.i --mesh-only CT_model.e
# 几何参数 (单位: mm)
L = 20.0e-3      # 总长度
W = 24.0e-3      # 总宽度  
# B = 1.0       # 厚度 (简化为2D问题)
notch_length = 8.5e-3    # 缺口长度
notch_width = 1.0e-3     # 缺口宽度
crack_length = 2.0    # 初始裂纹长度

# 网格细化参数
crack_tip_radius = 3.0e-3  # 裂纹尖端细化区域半径
fine_size = 0.1e-3         # 细网格尺寸
coarse_size = 0.5e-3       # 粗网格尺寸



[Mesh]
  # 创建基本矩形网格
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = ${W}
    nx = 100
    ny = 100
    elem_type = QUAD4
  []
  
  # 创建缺口 - 使用SubdomainBoundingBox定义缺口区域
  [notch_subdomain]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    bottom_left = '0 ${fparse W/2 - notch_width/2} 0'
    top_right = '${notch_length} ${fparse W/2 + notch_width/2} 0'
    block_name = 'notch_region'
  []
  
  # 删除缺口区域
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
  # # 创建加载点的nodeset
  # [upper_load_point]
  #   type = BoundingBoxNodeSetGenerator
  #   input = remove_notch
  #   new_boundary = 'upper_load'
  #   bottom_left = '${fparse load_x - 0.5} ${fparse load_upper_y - 0.5} 0'
  #   top_right = '${fparse load_x + 0.5} ${fparse load_upper_y + 0.5} 0'
  # []
  
  # [lower_load_point]  
  #   type = BoundingBoxNodeSetGenerator
  #   input = upper_load_point
  #   new_boundary = 'lower_load'
  #   bottom_left = '${fparse load_x - 0.5} ${fparse load_lower_y - 0.5} 0'
  #   top_right = '${fparse load_x + 0.5} ${fparse load_lower_y + 0.5} 0'
  # []
  
  # # 为裂纹前端细化网格
  # [refine_crack]
  #   type = RefineBlockGenerator
  #   input = lower_load_point
  #   block = '0'
  #   # refinement_level = 2
  #   enable_neighbor_refinement = true
  #   refinement = 2
  # []
[]
