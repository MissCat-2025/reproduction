# === 参数研究案例 ===
# endTime: 100050000
# 生成时间: 2025-08-26 00:06:03

# === 参数研究案例 ===
# end_time = 8.30e+6
# mesh_size: 8.00e-5
# 生成时间: 2025-08-15 21:21:25

# === 参数研究案例 ===
# end_time = 8.30e+6
# length_scale_paramete: 4.50e-5
# 生成时间: 2025-08-14 12:09:07

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]
# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)

# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
pellet_inner_diameter = 10.291         # 芯块内直径mm
pellet_outer_diameter = 14.627         # 芯块外直径mm
# length_scale_paramete = 4.50e-5
endTime = 100050000
endTime__50000 = '${fparse endTime-5000}'
endTime__100000 = '${fparse endTime-100000}'
mesh_size = 8.00e-5 #网格尺寸即可
# length_scale_paramete=${fparse mesh_size}
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
# 计算半径参数 (转换为米)
pellet_inner_radius = '${fparse pellet_inner_diameter/2*1e-3}'
pellet_outer_radius = '${fparse pellet_outer_diameter/2*1e-3}'
[Mesh]
[pellet1]
  type = AnnularMeshGenerator
  nr = ${n_radial_pellet}
  nt = ${n_azimuthal}
  rmin = ${pellet_inner_radius}
  rmax = ${pellet_outer_radius}
  growth_r = 1.006
  boundary_id_offset = 10
  boundary_name_prefix = 'pellet'
[]
[pellet]
  type = SubdomainIDGenerator
  input = pellet1
  subdomain_id = 1
[]
[rename1]
  type = RenameBoundaryGenerator
  input = pellet
  old_boundary = 'pellet_rmin pellet_rmax'
  new_boundary = 'pellet_inner pellet_outer'
[]
[cut_x]
  type = PlaneDeletionGenerator
  input = rename1
  point = '0 0 0'
  normal = '-1 0 0'  # 切割x>0区域
  new_boundary = 'y_axis'
[]
[cut_y]
  type = PlaneDeletionGenerator
  input = cut_x
  point = '0 0 0'
  normal = '0 -1 0'  # 切割y>0区域
  new_boundary = 'x_axis'
[]

[rename2]
  type = RenameBlockGenerator
  input = cut_y
  old_block = '1'
  new_block = 'pellet'
[]
[]

[Variables]
  [d]
    block = pellet
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [T]
    order = CONSTANT
    family = MONOMIAL
  []
  [a1]
    family = MONOMIAL
    order = CONSTANT
  []
  [Gc]
    family = MONOMIAL
    order = CONSTANT
  []
  [hoop_stress]
    family = MONOMIAL
    order = CONSTANT
  []
  [sigma0]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Bounds]
  [irreversibility]
    type = VariableOldValueBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
  []
[]

[BCs]
[]
[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
    block = pellet
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
    block = pellet
  []
[]

[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l'
    prop_values = '${l}'
    block = pellet
  []
[a11]
  type = ADParsedMaterial
  property_name = a1
  coupled_variables = 'a1'
  expression = 'a1'
  block = pellet
[]
[sigma0]
  type = ADParsedMaterial
  property_name = sigma0
  coupled_variables = 'sigma0'
  expression = 'sigma0'
  block = pellet
[]
[Gc1]
  type = ADParsedMaterial
  property_name = Gc
  coupled_variables = 'Gc'
  expression = 'Gc'
  block = pellet
[]
  #断裂力学-CZM模型
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d))*(1-eta)+eta
    phase_field = d
    material_property_names = 'a1'
    parameter_names = 'p a2 eta'
    parameter_values = '2 2 1e-6'
    block = pellet
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
    block = pellet
  [] 
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*(psie_active)'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    # output_properties = 'psi'
    # outputs = exodus
    block = pellet
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
  # solve_type = 'NEWTON'
  # solve_type = 'NEWTON'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type  -snes_type'
  petsc_options_value = '201                hypre    boomeramg  vinewtonrsls'  
  # petsc_options_iname = '-pc_type -ksp_type' 
  # petsc_options_value = 'lu gmres' 
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true  # 每个时间步都重新计算缩放
  # reuse_preconditioner = true
  # reuse_preconditioner_max_linear_its = 20
  nl_max_its = 20
  nl_rel_tol = 5e-9 # 非线性求解的相对容差
  nl_abs_tol = 5e-9 # 非线性求解的绝对容差
  l_tol = 5e-9  # 线性求解的容差
  l_abs_tol = 5e-9 # 线性求解的绝对容差
  l_max_its = 100 # 线性求解的最大迭代次数
  abort_on_solve_fail = true
  dtmin = 500
  dtmax = 50000
  end_time = ${endTime} # 总时间24h

  fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]
[Functions]
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 16000, 2000,
                   if(t < 105000, 750,
                   if(t < ${endTime__100000},50000,
                   if(t < (${endTime__50000}+10000), 750,10000))))'
  []
[]
[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = 2
  [Markers]
    [marker]
      type = PhasePiledFractureHSMarker
      von_mises_variable = hoop_stress
      sigma0 = sigma0
      x1 = 0.0001 #d变量小于x1时，标记为粗网格
      x2 = 0.05 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.1 #d变量大于xmax时，一定是细网格
      y1 = 0.4 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.5 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.02
      stress_change_threshold = 1e6
    []
  []
[]
[Outputs]
  # exodus = true
  print_linear_residuals = false
[]
