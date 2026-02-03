# 裂纹子算例，与第一阶段的 sub_dt1.i 结构保持一致
# 本阶段不做 solution restart，由 MultiApp 驱动重新计算 d 场

[Problem]
  kernel_coverage_check   = false
  material_coverage_check = false
[]

pellet_outer_radius = '${pellet_outer_radius}'
mesh_size           = '${mesh_size}'
w = ${w}
n_elems_azimuthal     = '${fparse 2*ceil((3.1415*pellet_outer_radius/mesh_size)/2^w)}'
n_elems_radial_pellet = '${fparse int((pellet_outer_radius/mesh_size)/2^w)}'

m  = ${m}
a2 = ${a2}
a3 = ${a3}

[Mesh]
  [pellet_clad_gap]
    type          = ConcentricCircleMeshGenerator
    num_sectors   = '${n_elems_azimuthal}'
    radii         = '${pellet_outer_radius}'
    rings         = '${n_elems_radial_pellet}'
    has_outer_square = false
    preserve_volumes = true
    portion       = top_right
    smoothing_max_it = 300
  []
  [rename]
    type  = RenameBoundaryGenerator
    input = pellet_clad_gap
    old_boundary = 'bottom left outer'
    new_boundary = 'yplane xplane pellet_outer'
  []
  [rename2]
    type  = RenameBlockGenerator
    input = rename
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
    order  = CONSTANT
    family = MONOMIAL
  []
  [a1]
    family = MONOMIAL
    order  = CONSTANT
  []
  [stress_I]
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[Bounds]
  [irreversibility]
    type            = VariableOldValueBounds
    variable        = bounds_dummy
    bounded_variable= d
    bound_type      = lower
  []
  [upper]
    type            = ConstantBounds
    variable        = bounds_dummy
    bounded_variable= d
    bound_type      = upper
    bound_value     = 1
  []
[]

[BCs]
[]

[Kernels]
  [diff]
    type                 = ADPFFDiffusion
    variable             = d
    fracture_toughness   = Gc
    regularization_length= l
    normalization_constant = c0
    block                = pellet
  []
  [source]
    type        = ADPFFSource
    variable    = d
    free_energy = psi
    block       = pellet
  []
[]

[Materials]
  [fracture_properties]
    type        = ADGenericConstantMaterial
    prop_names  = 'Gc sigma0 l'
    prop_values = '${Gc} ${sigma0} ${l}'
    block       = pellet
  []
  [a11]
    type              = ADParsedMaterial
    property_name     = a1
    coupled_variables = 'a1'
    expression        = 'a1'
    block             = pellet
  []
  [crack_geometric]
    type           = CrackGeometricFunction
    property_name  = alpha
    expression     = 'ksi*d+(1-ksi)*d*d'
    parameter_names  = 'ksi'
    parameter_values = '${ksi}'
    phase_field    = d
  []
  [degradation]
    type              = RationalDegradationFunction
    property_name     = g
    expression        = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))*(1-eta)+eta
    phase_field       = d
    material_property_names = 'a1'
    parameter_names   = 'p a2 a3 eta'
    parameter_values  = '${m} ${a2} ${a3} 1e-6'
  []
  [psi]
    type                  = ADDerivativeParsedMaterial
    property_name         = psi
    expression            = 'alpha*Gc/c0/l+g*(psie_active)'
    coupled_variables     = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    block                 = pellet
  []
[]

[Executioner]
  type       = Transient
  solve_type = 'PJFNK'

  petsc_options_iname  = '-ksp_gmres_restart -pc_type -pc_hypre_type  -snes_type'
  petsc_options_value  = '201                hypre    boomeramg       vinewtonrsls'

  nl_max_its = 200
  nl_rel_tol = 5e-9
  nl_abs_tol = 5e-10
  l_tol      = 5e-9
  l_abs_tol  = 5e-10
  l_max_its  = 100

  abort_on_solve_fail = true
  dtmin   = ${dtmin}
  end_time= ${endTime}

  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10
  accept_on_max_fixed_point_iteration = true

  [TimeStepper]
    type     = FunctionDT
    function = dt_limit_func
  []
[]

[Functions]
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < ${xTime}, 100,
                  if(t < (${xTime}+3), 1.5,
                  if(t < ${PowerTimeTotal},${dt1},
                  if(t < ${endTime},${dt1},${dt1}))))'
  []
[]

[Adaptivity]
  marker      = marker
  max_h_level = ${w}
  [Markers]
    [marker]
      type                 = PhasePiledFractureHSMarker
      von_mises_variable   = stress_I
      sigma0               = sigma0
      x1                   = 0.000001
      x2                   = 0.005
      xmax                 = 0.08
      y1                   = 0.45
      y2                   = 0.5
      variable             = d
      timeD                = 3
      timeStress           = 5
      d_change_threshold   = 0.01
      stress_change_threshold = 1e6
    []
  []
[]

[Outputs]
  print_linear_residuals = false
[]