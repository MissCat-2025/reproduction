# === 平面应力子应用 (二维截面) ===
# 作为主应用的 MultiApp 子应用，提供断裂相场 d 等字段

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]
# 几何与网格
pellet_inner_diameter = 9.90
pellet_outer_diameter = 14.10
endTime = 2e7
mesh_size = 7.00e-5
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}'
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
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
    normal = '-1 0 0'
    new_boundary = 'y_axis'
  []
  [cut_y]
    type = PlaneDeletionGenerator
    input = cut_x
    point = '0 0 0'
    normal = '0 -1 0'
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
    block = pellet
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type  -snes_type'
  petsc_options_value = '201                hypre    boomeramg  vinewtonrsls'
  automatic_scaling = true
  compute_scaling_once = true
  nl_max_its = 20
  nl_rel_tol = 5e-9
  nl_abs_tol = 5e-9
  l_tol = 5e-9
  l_abs_tol = 5e-9
  l_max_its = 100
  abort_on_solve_fail = true
  dtmin = 500
  dtmax = 25000
  end_time = ${endTime}
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]

[Functions]
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 12000, 2000, if(t < 110000, 500, if(t < (${endTime}-500),100000, 500)))'
  []
[]

[Outputs]
  print_linear_residuals = false
[]