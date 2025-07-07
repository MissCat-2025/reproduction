#  mpirun -n 8 /home/yp/projects/pff/pff-opt -i ElasticityFracture.i
E = 2.1e5
nu = 0.3
# K = '${fparse E/3/(1-2*nu)}'
# G = '${fparse E/2/(1+nu)}'

Gc = 2.7
l = 0.02

[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
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
  max_h_level = 6
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
  [strain_zz]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [fy]
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
  [./solid_z]
    type = ADWeakPlaneStress
    variable = strain_zz
  [../]
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
    prop_names = 'E nu Gc l'
    prop_values = '${E} ${nu} ${Gc} ${l}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [strain]
    type = ADComputePlaneSmallStrain
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd^2'
    phase_field = d
  []
  [elasticity]
    type = SmallDeformationHBasedElasticity
    youngs_modulus = E
    poissons_ratio = nu
    tensile_strength = 10
    fracture_energy = Gc
    # bulk_modulus = K
    # shear_modulus = G
    phase_field = d
    degradation_function = g
    # decomposition = SPECTRAL
    output_properties = 'psie_active'
    outputs = exodus
  []
  # [elasticity]
  #   type = SmallDeformationIsotropicElasticityThreshold
  #   youngs_modulus = E
  #   poissons_ratio = nu 
  #   phase_field = d
  #   degradation_function = g
  #   decomposition = spectral
  #   use_threshold = true
  #   tensile_strength = 0.0000001
  #   output_properties = 'psie_active psie'
  #   outputs = exodus
  # []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*psie_active'
    coupled_variables = 'd'
    material_property_names = 'alpha(d) g(d) Gc c0 l psie_active(d)'
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
  petsc_options_iname = '-pc_type -snes_type   -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       vinewtonrsls NONZERO               1e-10'
  automatic_scaling = true
  line_search = bt

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 2e-5
  end_time = 3.5e-3

  fixed_point_max_its = 20
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
