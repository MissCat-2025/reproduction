# === 参数研究案例 ===
# pellet_critical_energy: 9
# length_scale_paramete: 5.00e-5
# dt: 2000
# 生成时间: 2026-01-16 00:53:57

# === 参数研究案例 ===
# fission_rate: 1.60e+19
# largestPoreSize: 58
# pellet_critical_energy: 2.8
# 生成时间: 2025-10-20 20:26:05

# === 参数研究案例（对齐配对） ===
# LinearPower: 90
# initial_T_in: 570.7
# initial_T_out: 582.8
# 生成时间: 2025-09-23 18:08:10

# === 参数研究案例 ===
# end_time = 8.30e+6
# mesh_size: 8.00e-5
# 生成时间: 2025-08-15 21:21:25

# === 参数研究案例 ===
# end_time = 8.30e+6
# length_scale_paramete: 4.50e-5
# 生成时间: 2025-08-14 12:09:07
endTime = ${endTime}
endTime__50000 = '${fparse endTime-50000}'
endTime__100000 = '${fparse endTime-100000}'
[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

pellet_outer_radius = '${pellet_outer_radius}'         # 芯块外直径mm
mesh_size = '${mesh_size}' #网格尺寸即可
w = ${w} #裂纹尖端时，l是h的2**w倍
n_elems_azimuthal = '${fparse 2*ceil((3.1415*pellet_outer_radius/mesh_size)/2^w)}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int((pellet_outer_radius/mesh_size)/2^w)}'          # 芯块径向网格数（直接取整）

#相场断裂参数：
m = ${m}
a2 = ${a2}
a3 = ${a3}



[Mesh]
  [pellet_clad_gap]
    type = ConcentricCircleMeshGenerator
    num_sectors = '${n_elems_azimuthal}'  # 周向网格数
    radii = '${pellet_outer_radius}'
    rings = '${n_elems_radial_pellet}'
    has_outer_square = false
    preserve_volumes = true
    portion = top_right # 生成四分之一计算域
    smoothing_max_it=300 # 平滑迭代次数
  []
  [rename]
    type = RenameBoundaryGenerator
    input = pellet_clad_gap
    old_boundary = 'bottom left outer'
    new_boundary = 'yplane xplane pellet_outer' # 将边界命名为yplane xplane clad_outer

  []
  [rename2]
    type = RenameBlockGenerator
    input = rename
    old_block  = '1'
    new_block  = 'pellet' # 将block1和block3分别命名为pellet和clad
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
  [a1]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_I]
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
    prop_names = 'Gc sigma0 l'
    prop_values = '${Gc} ${sigma0} ${l}'
    block = pellet
  []
[a11]
  type = ADParsedMaterial
  property_name = a1
  coupled_variables = 'a1'
  expression = 'a1'
  block = pellet
[]
  #断裂力学-CZM模型
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'ksi*d+(1-ksi)*d*d'
    parameter_names = 'ksi'
    parameter_values = '${ksi}'
    phase_field = d
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))
    phase_field = d
    material_property_names = 'a1'
    parameter_names = 'p a2 a3'
    parameter_values = '${m} ${a2} ${a3}'
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
  # solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -ksp_type  -snes_type'
  # petsc_options_value = 'lu gmres  vinewtonrsls'  
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type  -snes_type'
  petsc_options_value = '201                hypre    boomeramg  vinewtonrsls'  
  # petsc_options_iname = '-pc_type -ksp_type' 
  # petsc_options_value = 'lu gmres' 
  # automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  # compute_scaling_once = true  # 每个时间步都重新计算缩放
  # reuse_preconditioner = true
  # reuse_preconditioner_max_linear_its = 20
  nl_max_its = 200
  nl_rel_tol = 5e-9 # 非线性求解的相对容差
  nl_abs_tol = 5e-10 # 非线性求解的绝对容差
  l_tol = 5e-9  # 线性求解的容差
  l_abs_tol = 5e-10 # 线性求解的绝对容差
  l_max_its = 100 # 线性求解的最大迭代次数
  abort_on_solve_fail = true
  dtmin = ${dtmin}
  dtmax = ${dtMax}
  end_time = ${endTime} # 总时间24h

  fixed_point_rel_tol =1e-8 # 固定点迭代的相对容差
  fixed_point_abs_tol = 1e-10 # 固定点迭代的绝对容差
  accept_on_max_fixed_point_iteration = true
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]
[Functions]
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 1000, 500,
                  if(t < 3000, 50,
                  if(t < 110000, ${dt},
                  if(t < (${endTime__100000}-5000),${dtMax},
                  if(t < (${endTime__50000}+10000), ${dt},10000)))))'
  []
[]
[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = ${w}
  [Markers]
    [marker]
      type = PhasePiledFractureHSMarker
      von_mises_variable = stress_I
      sigma0 = sigma0
      x1 = 0.000001 #d变量小于x1时，标记为粗网格
      x2 = 0.005 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.08 #d变量大于xmax时，一定是细网格
      y1 = 0.45 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.5 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.01
      stress_change_threshold = 1e6
    []
  []
[]
[Outputs]
  # exodus = true
  print_linear_residuals = false
[]
