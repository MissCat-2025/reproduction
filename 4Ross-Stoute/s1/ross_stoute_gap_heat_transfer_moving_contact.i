## Units in the input file: m-Pa-s-K
#
# Two square blocks start separated by a gap. The left block is hotter and is
# translated slowly toward the fixed right block. This exercises the new
# Ross-Stoute-style gas, radiation, and contact gap flux models under a moving
# interface.

[GlobalParams]
  displacements = 'disp_x disp_y'
  block = '1 2'
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [left_rectangle]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 10
    xmin = 0.0
    xmax = 0.9
    ymin = 0.0
    ymax = 0.5
    boundary_name_prefix = moving_block
  []
  [left_block]
    type = SubdomainIDGenerator
    input = left_rectangle
    subdomain_id = 1
  []
  [right_rectangle]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 10
    xmin = 1.1
    xmax = 2.0
    ymin = 0.0
    ymax = 0.5
    boundary_name_prefix = fixed_block
    boundary_id_offset = 4
  []
  [right_block]
    type = SubdomainIDGenerator
    input = right_rectangle
    subdomain_id = 2
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'left_block right_block'
  []
  [mechanical_contact_primary_subdomain]
    type = LowerDBlockFromSidesetGenerator
    input = combined
    sidesets = 'moving_block_right'
    new_block_id = 3
    new_block_name = mechanical_contact_primary_subdomain
  []
  [mechanical_contact_secondary_subdomain]
    type = LowerDBlockFromSidesetGenerator
    input = mechanical_contact_primary_subdomain
    sidesets = 'fixed_block_left'
    new_block_id = 4
    new_block_name = mechanical_contact_secondary_subdomain
  []
[]

[Variables]
  [temperature]
    initial_condition = 525.0
  []
  [temperature_interface_lm]
    block = mechanical_contact_secondary_subdomain
  []
[]

[Functions]
  [move_left]
    type = PiecewiseLinear
    x = '0 1.5 3.0'
    y = '0 0.12 0.25'
  []

  [helium_k]
    type = ParsedFunction
    expression = '0.04679 + 3.81e-4*t - 6.786e-8*t^2'
  []

  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 2.3, 0.1,0.001)'
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [solid]
    add_variables = true
    strain = FINITE
    generate_output = 'vonmises_stress'
  []
[]

[Kernels]
  [heat_left]
    type = ADHeatConduction
    variable = temperature
    thermal_conductivity = left_thermal_conductivity
    block = 1
  []
  [heat_right]
    type = ADHeatConduction
    variable = temperature
    thermal_conductivity = right_thermal_conductivity
    block = 2
  []
[]

[BCs]
  [temperature_left]
    type = ADDirichletBC
    variable = temperature
    value = 800
    boundary = 'moving_block_left'
  []
  [temperature_right]
    type = ADDirichletBC
    variable = temperature
    value = 250
    boundary = 'fixed_block_right'
  []

  [move_left_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'moving_block_left moving_block_right moving_block_top moving_block_bottom'
    function = move_left
  []
  [move_left_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'moving_block_left moving_block_right moving_block_top moving_block_bottom'
    value = 0
  []

  [fix_right_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'fixed_block_left fixed_block_right fixed_block_top fixed_block_bottom'
    value = 0
  []
  [fix_right_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'fixed_block_left fixed_block_right fixed_block_top fixed_block_bottom'
    value = 0
  []
[]

[Materials]
  [left_elasticity]
    type = ComputeIsotropicElasticityTensor
    block = 1
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  []
  [left_stress]
    type = ComputeFiniteStrainElasticStress
    block = 1
  []

  [right_elasticity]
    type = ComputeIsotropicElasticityTensor
    block = 2
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  []
  [right_stress]
    type = ComputeFiniteStrainElasticStress
    block = 2
  []

  [left_thermal]
    type = ADHeatConductionMaterial
    block = 1
    thermal_conductivity = 16.2
    specific_heat = 1.0
  []
  [left_hardness]
    type = ADGenericConstantMaterial
    block = 1
    prop_names = 'left_hardness'
    prop_values = '129'
  []
  [right_thermal]
    type = ADHeatConductionMaterial
    block = 2
    thermal_conductivity = 210.0
    specific_heat = 1.0
  []
  [right_hardness]
    type = ADGenericConstantMaterial
    # block = 2
    prop_names = 'right_hardness left_thermal_conductivity right_thermal_conductivity'
    prop_values = '15 210 210'
  []
[]

[Contact]
  [mechanical_contact]
    primary = moving_block_right
    secondary = fixed_block_left
    model = frictionless
    formulation = penalty
    penalty = 1e12
    normalize_penalty = true
    tangential_tolerance = 1e-6
    normal_smoothing_distance = 1e-6
    generate_mortar_mesh = false
    # correct_edge_dropping = true
  []
[]

[UserObjects]
  [gas]
    type = GapFluxModelRossStouteGasConduction
    boundary = moving_block_right
    temperature = temperature
    gas_thermal_conductivity = helium_k
    fuel_roughness = 2.0e-6
    clad_roughness = 1.0e-6
    fuel_jump_distance = 0.0
    clad_jump_distance = 0.0
  []

  [radiation]
    type = GapFluxModelRadiation
    boundary = moving_block_right
    temperature = temperature
    primary_emissivity = 0.8
    secondary_emissivity = 0.6
  []

  [contact]
    type = GapFluxModelRossStouteContactConduction
    boundary = moving_block_right
    temperature = temperature
    contact_pressure = contact_pressure
    primary_conductivity = left_thermal_conductivity
    secondary_conductivity = right_thermal_conductivity
    primary_hardness = left_hardness
    secondary_hardness = right_hardness
    fuel_roughness = 2.0e-6
    clad_roughness = 1.0e-6
  []
[]

[Constraints]
  [thermal_contact]
    type = ModularGapConductanceConstraint
    variable = temperature_interface_lm
    secondary_variable = temperature
    use_displaced_mesh = true
    primary_boundary = moving_block_right
    primary_subdomain = mechanical_contact_primary_subdomain
    secondary_boundary = fixed_block_left
    secondary_subdomain = mechanical_contact_secondary_subdomain
    gap_flux_models = 'gas radiation contact'
  []
[]

[Postprocessors]
  [left_temperature]
    type = SideAverageValue
    boundary = moving_block_left
    variable = temperature
  []
  [right_temperature]
    type = SideAverageValue
    boundary = fixed_block_right
    variable = temperature
  []
  [contact_pressure_avg]
    type = ElementAverageValue
    variable = contact_pressure
    block = mechanical_contact_secondary_subdomain
  []
  [interface_heat_flux_left]
    type = ADSideDiffusiveFluxAverage
    variable = temperature
    boundary = moving_block_right
    diffusivity = left_thermal_conductivity
  []
  [interface_heat_flux_right]
    type = ADSideDiffusiveFluxAverage
    variable = temperature
    boundary = fixed_block_left
    diffusivity = right_thermal_conductivity
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu nonzero'
  dt = 0.01
  end_time = 3.0
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 20
  automatic_scaling = false
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]