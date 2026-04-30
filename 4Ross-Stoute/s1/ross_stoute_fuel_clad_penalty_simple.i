## Units in the input file: m-Pa-s-K
#
# Minimal Ross-Stoute fuel-clad interaction smoke test.
# - one fuel disk
# - one cladding ring
# - frictionless penalty contact for mechanics
# - Ross-Stoute gas + radiation + contact conduction for the thermal gap
# dos2unix ross_stoute_fuel_clad_penalty_simple.i &&mpirun -n 1 /home/yp/projects/reproduction/reproduction-opt -i ross_stoute_fuel_clad_penalty_simple.i
initial_T = 583.15
coolant_T = 550.0

pellet_density = 1.0431e4
clad_density = 6.59e3

fuel_youngs_modulus = 2.01e11
clad_youngs_modulus = 1.09e11
fuel_poissons_ratio = 0.345
clad_poissons_ratio = 0.334
fuel_thermal_expansion_coef = 1.0e-5
clad_thermal_expansion_coef = 5.0e-6
fuel_hardness = 15.0
clad_hardness = 129.0

fuel_heat_source = 4.0e8
coolant_conductance = 5.0e3
contact_penalty = 1.0e13

fuel_outer_diameter = 14.100e-3
clad_inner_diameter = 14.224e-3
clad_outer_diameter = 15.367e-3

fuel_outer_radius = '${fparse fuel_outer_diameter/2}'
clad_inner_radius = '${fparse clad_inner_diameter/2}'
clad_outer_radius = '${fparse clad_outer_diameter/2}'

n_azimuthal = 256
n_radial_fuel = 24
n_radial_clad = 4
growth_factor = 1.01

[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y'
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [fuel_disk]
    type = AnnularMeshGenerator
    nr = ${n_radial_fuel}
    nt = ${n_azimuthal}
    rmin = 3.100e-3
    rmax = ${fuel_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 10
    boundary_name_prefix = fuel
  []
  [fuel_block]
    type = SubdomainIDGenerator
    input = fuel_disk
    subdomain_id = 1
  []
  [clad_ring]
    type = AnnularMeshGenerator
    nr = ${n_radial_clad}
    nt = ${n_azimuthal}
    rmin = ${clad_inner_radius}
    rmax = ${clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 20
    boundary_name_prefix = clad
  []
  [clad_block]
    type = SubdomainIDGenerator
    input = clad_ring
    subdomain_id = 2
  []
  [combine]
    type = MeshCollectionGenerator
    inputs = 'fuel_block clad_block'
  []
  [rename_blocks]
    type = RenameBlockGenerator
    input = combine
    old_block = '1 2'
    new_block = 'fuel clad'
  []
  [rename_boundaries]
    type = RenameBoundaryGenerator
    input = rename_blocks
    old_boundary = 'fuel_rmax clad_rmin clad_rmax'
    new_boundary = 'fuel_outer clad_inner clad_outer'
  []
  [cut_x]
    type = PlaneDeletionGenerator
    input = rename_boundaries
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
  [fuel_interface]
    type = LowerDBlockFromSidesetGenerator
    input = cut_y
    sidesets = 'fuel_outer'
    new_block_id = 3
    new_block_name = fuel_interface
  []
  [clad_interface]
    type = LowerDBlockFromSidesetGenerator
    input = fuel_interface
    sidesets = 'clad_inner'
    new_block_id = 4
    new_block_name = clad_interface
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temperature]
    initial_condition = ${initial_T}
  []
  [temperature_interface_lm]
    block = clad_interface
  []
[]

[Functions]
  [helium_k]
    type = ParsedFunction
    expression = '0.04679 + 3.81e-4*t - 6.786e-8*t^2'
  []

  [gap_gas_pressure]
    type = ParsedFunction
    expression = '1.0e6 + 2.0e5*(1.0 - exp(-t/10.0))'
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    block = 'fuel clad'
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    block = 'fuel clad'
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = temperature
    thermal_conductivity = thermal_conductivity
    block = 'fuel clad'
  []
  [heat_time]
    type = ADHeatConductionTimeDerivative
    variable = temperature
    block = 'fuel clad'
  []
  [heat_source]
    type = ADMatHeatSource
    variable = temperature
    material_property = total_power
    block = fuel
  []
[]

[BCs]
  [sym_x]
    type = DirichletBC
    variable = disp_x
    boundary = y_axis
    value = 0
  []
  [sym_y]
    type = DirichletBC
    variable = disp_y
    boundary = x_axis
    value = 0
  []
  [coolant_bc]
    type = ConvectiveFluxFunction
    variable = temperature
    boundary = clad_outer
    T_infinity = ${coolant_T}
    coefficient = ${coolant_conductance}
  []
  [coolant_pressure_x]
    type = Pressure
    variable = disp_x
    boundary = clad_outer
    factor = 15.5e6
    use_displaced_mesh = true
  []
  [coolant_pressure_y]
    type = Pressure
    variable = disp_y
    boundary = clad_outer
    factor = 15.5e6
    use_displaced_mesh = true
  []
[]

[Materials]
  [fuel_properties]
    type = ADGenericConstantMaterial
    block = fuel
    prop_names = 'density specific_heat thermal_conductivity youngs_modulus poissons_ratio hardness total_power'
    prop_values = '${pellet_density} 300.0 3.0 ${fuel_youngs_modulus} ${fuel_poissons_ratio} ${fuel_hardness} ${fuel_heat_source}'
  []
  [clad_properties]
    type = ADGenericConstantMaterial
    block = clad
    prop_names = 'density specific_heat thermal_conductivity youngs_modulus poissons_ratio hardness'
    prop_values = '${clad_density} 300.0 16.2 ${clad_youngs_modulus} ${clad_poissons_ratio} ${clad_hardness}'
  []

  [fuel_thermal_expansion]
    type = ADComputeThermalExpansionEigenstrain
    eigenstrain_name = thermal_eigenstrain
    stress_free_temperature = ${initial_T}
    thermal_expansion_coeff = ${fuel_thermal_expansion_coef}
    temperature = temperature
    block = fuel
  []
  [fuel_strain]
    type = ADComputePlaneSmallStrain
    eigenstrain_names = thermal_eigenstrain
    block = fuel
  []
  [fuel_elasticity]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = fuel
  []

  [clad_thermal_expansion]
    type = ADComputeThermalExpansionEigenstrain
    eigenstrain_name = thermal_eigenstrain
    stress_free_temperature = ${initial_T}
    thermal_expansion_coeff = ${clad_thermal_expansion_coef}
    temperature = temperature
    block = clad
  []
  [clad_strain]
    type = ADComputePlaneSmallStrain
    eigenstrain_names = thermal_eigenstrain
    block = clad
  []
  [clad_elasticity]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = clad
  []

  [stress]
    type = ADComputeLinearElasticStress
    block = 'fuel clad'
  []
[]

[Contact]
  [mechanical_contact]
    primary = fuel_outer
    secondary = clad_inner
    model = frictionless
    formulation = penalty
    penalty = ${contact_penalty}
    normalize_penalty = true
    tangential_tolerance = 1e-6
    normal_smoothing_distance = 1e-6
    generate_mortar_mesh = false
  []
[]

[UserObjects]
  [gas]
    type = GapFluxModelRossStouteGasConduction
    boundary = fuel_outer
    temperature = temperature
    gas_thermal_conductivity = helium_k
    gas_pressure = gap_gas_pressure
    gas_mixture_mode = HELIUM_ONLY
    fuel_roughness = 2.0e-6
    clad_roughness = 1.0e-6
    fuel_jump_distance = 0.0
    clad_jump_distance = 0.0
  []
  [radiation]
    type = GapFluxModelRadiation
    boundary = fuel_outer
    temperature = temperature
    primary_emissivity = 0.8
    secondary_emissivity = 0.6
  []
  [contact]
    type = GapFluxModelRossStouteContactConduction
    boundary = fuel_outer
    temperature = temperature
    contact_pressure = contact_pressure
    primary_conductivity = thermal_conductivity
    secondary_conductivity = thermal_conductivity
    primary_hardness = hardness
    secondary_hardness = hardness
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
    displacements = 'disp_x disp_y'
    primary_boundary = fuel_outer
    primary_subdomain = fuel_interface
    secondary_boundary = clad_inner
    secondary_subdomain = clad_interface
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0 0 0'
    cylinder_axis_point_2 = '0 0 1'
    gap_flux_models = 'gas radiation contact'
  []
[]

[Postprocessors]
  [fuel_outer_T_avg]
    type = AverageNodalVariableValue
    variable = temperature
    block = fuel_interface
  []
  [clad_inner_T_avg]
    type = AverageNodalVariableValue
    variable = temperature
    block = clad_interface
  []
  [fuel_outer_heat_rate]
    type = ADSideDiffusiveFluxIntegral
    variable = temperature
    boundary = fuel_outer
    diffusivity = thermal_conductivity
    execute_on = 'TIMESTEP_END'
  []
  [clad_inner_heat_rate]
    type = ADSideDiffusiveFluxIntegral
    variable = temperature
    boundary = clad_inner
    diffusivity = thermal_conductivity
    execute_on = 'TIMESTEP_END'
  []
  [fuel_outer_perimeter]
    type = AreaPostprocessor
    boundary = fuel_outer
    execute_on = 'initial timestep_end'
  []
  [clad_inner_perimeter]
    type = AreaPostprocessor
    boundary = clad_inner
    execute_on = 'initial timestep_end'
  []
  [h_eq_gap_fuel]
    type = ParsedPostprocessor
    expression = 'abs(fuel_outer_heat_rate) / (fuel_outer_perimeter*max(fuel_outer_T_avg - clad_inner_T_avg, 1e-6))'
    pp_names = 'fuel_outer_heat_rate fuel_outer_perimeter fuel_outer_T_avg clad_inner_T_avg'
  []
  [h_eq_gap_clad]
    type = ParsedPostprocessor
    expression = 'abs(clad_inner_heat_rate) / (clad_inner_perimeter*max(fuel_outer_T_avg - clad_inner_T_avg, 1e-6))'
    pp_names = 'clad_inner_heat_rate clad_inner_perimeter fuel_outer_T_avg clad_inner_T_avg'
  []
  [contact_pressure_avg]
    type = AverageNodalVariableValue
    variable = contact_pressure
    block = clad_interface
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu preonly'
  # solve_type = 'NEWTON'
  line_search = contact
  automatic_scaling = true
  compute_scaling_once = true
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 20
  dt = 0.001
  dtmin = 1e-7
  dtmax = 10
  end_time = 20

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 8
    iteration_window = 4
  []
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
