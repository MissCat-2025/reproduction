#  mpirun -n 8 /home/yp/projects/pff/pff-opt -i ElasticityFracture.i
E = 2.1e5
nu = 0.3
Gc = 2.7
l = 0.02

[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
  planar_formulation = PLANE_STRESS
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 30
    ny = 15
    ymax = 0.5
  []
  [noncrack]
    type = BoundingBoxNodeSetGenerator
    input = gen
    new_boundary = noncrack
    bottom_left = '0.5 0 0'
    top_right = '1 0 0'
  []
  construct_side_list_from_node_list = true
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 2
  stop_time = 0
  max_h_level = 2
  [Markers]
    [marker]
      type = BoxMarker
      bottom_left = '0.4 0 0'
      top_right = '1 0.05 0'
      outside = DO_NOTHING
      inside = REFINE
    []
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [d]
  []
[]

[AuxVariables]
  [fy]
  []
  [bounds_dummy]
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
[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    save_in = fy
  []
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [pff_source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[BCs]
  [ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 't'
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = noncrack
    value = 0
  []
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = top
    value = 0
  []
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'Gc l E nu'
    prop_values = '${Gc} ${l} ${E} ${nu}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-10'
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd^2'
    phase_field = d
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticityThreshold
    youngs_modulus = E
    poissons_ratio = nu 
    phase_field = d
    degradation_function = g
    decomposition = spectral
    use_threshold = true
    tensile_strength = 0.0000001
    output_properties = 'psie_active psie'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+psie'
    coupled_variables = 'd'
    material_property_names = 'alpha(d) g(d) Gc c0 l psie(d)'
    derivative_order = 1
  []
[]

[Postprocessors]
  [Fy]
    type = NodalSum
    variable = fy
    boundary = top
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'lu'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '
  # petsc_options_iname = '-snes_type        -snes_qn_type   -snes_qn_scale_type' 
  # petsc_options_value = 'qn               lbfgs           DIAGONAL'
  #   petsc_options_iname = '-snes_type' 
  # petsc_options_value = 'qn'
  # petsc_options_iname = '-pc_type   -snes_type        -snes_qn_type   -snes_qn_scale_type' 
  # petsc_options_value = 'lu         qn               lbfgs           jacobian'  
  # petsc_options_iname = '-snes_type -snes_qn_type -snes_linesearch_type'
  # petsc_options_value = 'qn lbfgs bt'
  # petsc_options_iname = '-snes_qn_type' 
  # petsc_options_value = 'lbfgs'
#   petsc_options_iname = '-pc_type -snes_type_qn'
# petsc_options_value = 'lu lbfgs'
  petsc_options_iname = '-pc_type -snes_type   -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       vinewtonrsls NONZERO               1e-10'
  abort_on_solve_fail = true
  line_search = none
  automatic_scaling = true
  # line_search = bt
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8

  dt = 2e-5
  end_time = 3.5e-3
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
