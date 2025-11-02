# === 参数研究案例 ===
# end_time = 0.05
# dt: 1.00e-4
# 生成时间: 2025-11-01 20:48:52

# 陶瓷片热冲击实验 - 相场断裂部分
l = 0.10e-3                # 相场正则化长度 (m)
nh = 2
nx = '${fparse int(25e-3/(nh*nh*l/3))}'
ny = '${fparse int(5e-3/(nh*nh*l/3))}'
[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}            # 25mm / 0.05mm = 500
    ny = ${ny}            # 5mm / 0.05mm = 100
    xmax = 25e-3
    ymax = 5e-3
  []
[]

[Variables]
  [d]  # 相场变量
  []
[]

[AuxVariables]
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [bounds_dummy]
  []
  [MaxPrincipal]
    family = MONOMIAL
    order = CONSTANT
  []
  [sigma0]
    family = MONOMIAL
    order = CONSTANT
  []
  [a1]
    family = MONOMIAL
    order = CONSTANT
  []
[]
[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
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
[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'Gc l'
    prop_values = '${Gc} ${l}'
  []
  [sigma0]
  type = ADParsedMaterial
  property_name = sigma0
  coupled_variables = 'sigma0'
  expression = 'sigma0'
[]
  [a11]
    type = ADParsedMaterial
    property_name = a1
    coupled_variables = 'a1'
    expression = 'a1'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = '2*d-d*d'
    phase_field = d
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))
    phase_field = d
    material_property_names = 'a1'
    parameter_names = 'p a2 a3'
    parameter_values = '2 -0.5 0'
  []
  
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*psie_active'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
[]

[Adaptivity]
  initial_marker = boundary
  initial_steps = ${nh}
  marker = marker
  max_h_level = ${nh}
  [Markers]
    [marker]
      type = PhasePiledFractureHSMarker
      von_mises_variable = MaxPrincipal
      sigma0 = sigma0
      x1 = 1e-6 #d变量小于x1时，标记为粗网格
      x2 = 0.05 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.1 #d变量大于xmax时，一定是细网格
      y1 = 0.3 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.6 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.02
      stress_change_threshold = 1e6
    []
    [boundary]
      type = BoundaryMarker
      next_to = 'top left'
      distance = 1e-4
      mark = refine
    []
  []
[]
[Executioner]
  type = Transient
  
  solve_type = NEWTON
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -snes_type'
  petsc_options_value = '201                hypre    boomeramg vinewtonrsls'
  automatic_scaling = true
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  dt = 1.00e-4
  end_time = 50e-3
[]

[Outputs]
  exodus = false
  print_linear_residuals = false
[]