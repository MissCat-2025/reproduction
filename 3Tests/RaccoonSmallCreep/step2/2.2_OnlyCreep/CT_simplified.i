# 简化CT模型 - 实用版本
# 主要简化：
# 1. 将圆孔简化为在相应位置施加载荷
# 2. 缺口简化为矩形
# 3. 专注于裂纹扩展行为

# 几何参数 (单位: mm)
L = 25.0      # 总长度
W = 24.0      # 总宽度  
B = 1.0       # 厚度 (简化为2D问题)
notch_length = 8.5    # 缺口长度
notch_width = 2.0     # 缺口宽度
crack_length = 2.0    # 初始裂纹长度

# 加载点位置 (圆孔中心位置)
load_upper_y = 18.0   # 上部加载点Y坐标
load_lower_y = 6.0    # 下部加载点Y坐标  
load_x = 4.0          # 加载点X坐标

# 材料参数
E = 140e9     # 弹性模量 (Pa)
nu = 0.3      # 泊松比
Gc = 160e3    # 断裂韧性 (J/m^2)
l = 0.1       # 相场长度尺度 (mm)

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

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
    ny = 96
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
  
  # 创建加载点的nodeset
  [upper_load_point]
    type = BoundingBoxNodeSetGenerator
    input = remove_notch
    new_boundary = 'upper_load'
    bottom_left = '${fparse load_x - 0.5} ${fparse load_upper_y - 0.5} 0'
    top_right = '${fparse load_x + 0.5} ${fparse load_upper_y + 0.5} 0'
  []
  
  [lower_load_point]  
    type = BoundingBoxNodeSetGenerator
    input = upper_load_point
    new_boundary = 'lower_load'
    bottom_left = '${fparse load_x - 0.5} ${fparse load_lower_y - 0.5} 0'
    top_right = '${fparse load_x + 0.5} ${fparse load_lower_y + 0.5} 0'
  []
  
  # 为裂纹前端细化网格
  [refine_crack]
    type = RefineBlockGenerator
    input = lower_load_point
    block = '0'
    refinement_level = 2
    enable_neighbor_refinement = true
  []
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE  
  []
  [d]  # 相场损伤变量
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = FunctionIC
      # 在缺口尖端设置初始裂纹
      function = 'if(x>${notch_length} & x<=${fparse notch_length + crack_length} & abs(y-${fparse W/2})<0.2, 1.0, 0.0)'
    []
  []
[]

[AuxVariables]
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [vonmises]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_xx]
    type = ADRankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  []
  [stress_yy]
    type = ADRankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  [vonmises]
    type = ADRankTwoScalarAux
    variable = vonmises
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  # 力学平衡方程
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = false
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
  []
  
  # 相场方程
  [pff_diff]
    type = ADPFFDiffusion
    variable = d
  []
  [pff_source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[Materials]
  # 基本材料属性
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'E nu Gc l'
    prop_values = '${E} ${nu} ${Gc} ${l}'
  []
  
  # 降解函数
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = '(1-d)^2/((1-d)^2+d)'
    phase_field = d
  []
  
  # 弹性模型 (考虑相场降解)
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    youngs_modulus = E
    poissons_ratio = nu
    phase_field = d
    degradation_function = g
  []
  
  # 应变计算
  [strain]
    type = ADComputeSmallStrain
  []
  
  # 应力计算
  [stress]
    type = ComputeSmallDeformationStress  
    elasticity_model = elasticity
  []
  
  # 弹性能量密度
  [elastic_energy]
    type = ElasticEnergyDensity
    property_name = psie_active
  []
  
  # 相场能量
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'g*psie_active + Gc/4/l*(d^2 + l^2*(grad_d_x^2 + grad_d_y^2))'
    coupled_variables = 'd'
    material_property_names = 'g psie_active Gc l'
    derivative_order = 1
  []
[]

[BCs]
  # 对称约束 - 假设左边界对称
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  
  # 下边界固定Y方向
  [fix_y_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'  
    value = 0
  []
[]

[Functions]
  # 载荷函数
  [load_function]
    type = ParsedFunction
    expression = '1000*t'  # 线性加载，单位N
  []
[]

# CT试样的载荷施加
[DiracKernels]
  # 上部加载点 - 向上拉力
  [upper_load]
    type = ConstantPointSource
    variable = disp_y
    point = '${load_x} ${load_upper_y} 0'
    value = 1000  # 向上的力 (N)
  []
  
  # 下部加载点 - 向下拉力  
  [lower_load]
    type = ConstantPointSource
    variable = disp_y
    point = '${load_x} ${load_lower_y} 0'
    value = -1000  # 向下的力 (N)
  []
[]

[Executioner]
  type = Transient
  
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu superlu_dist'
  
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08
  nl_max_its = 15
  
  dt = 0.1
  end_time = 2.0
  
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
  
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'final'
    level = 2
  []
[]

[Postprocessors]
  # 最大应力
  [max_stress]
    type = ElementExtremeValue
    variable = vonmises
    value_type = max
  []
  
  # 最大损伤
  [max_damage]
    type = ElementExtremeValue
    variable = d
    value_type = max
  []
  
  # 载荷点位移
  [upper_disp]
    type = PointValue
    variable = disp_y
    point = '${load_x} ${load_upper_y} 0'
  []
  
  [lower_disp]
    type = PointValue
    variable = disp_y  
    point = '${load_x} ${load_lower_y} 0'
  []
  
  # 裂纹长度（基于损伤场）
  [crack_length]
    type = ElementIntegralVariablePostprocessor
    variable = d
    threshold = 0.5
    execute_on = 'TIMESTEP_END'
  []
[] 