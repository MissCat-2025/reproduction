# Ross-Stoute gap conductance template
#
# Drop the blocks below into your mortar thermal contact input.
# Assumptions:
# - primary side = fuel pellet
# - secondary side = cladding
# - contact_pressure is supplied by your mechanics/contact solve
# - replace the boundary and subdomain names with your own model names

[Functions]
  [helium_k]
    type = ParsedFunction
    expression = '0.04679 + 3.81e-4*t - 6.786e-8*t^2'
  []
[]

[Materials]
  [fuel_props]
    type = ADGenericConstantMaterial
    block = fuel
    prop_names = 'fuel_thermal_conductivity fuel_hardness'
    prop_values = '3.0 15.0'
  []

  [clad_props]
    type = ADGenericConstantMaterial
    block = clad
    prop_names = 'clad_thermal_conductivity clad_hardness'
    prop_values = '16.2 129.0'
  []
[]

[UserObjects]
  [gas]
    type = GapFluxModelRossStouteGasConduction
    boundary = fuel_outer
    temperature = temperature
    gas_thermal_conductivity = helium_k
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
    contact_pressure = interface_normal_lm
    primary_conductivity = fuel_thermal_conductivity
    secondary_conductivity = clad_thermal_conductivity
    primary_hardness = fuel_hardness
    secondary_hardness = clad_hardness
    fuel_roughness = 2.0e-6
    clad_roughness = 1.0e-6
  []
[]

[Constraints]
  [thermal_contact]
    type = ModularGapConductanceConstraint
    variable = temperature_interface_lm
    secondary_variable = temperature
    primary_boundary = fuel_outer
    primary_subdomain = fuel_interface
    secondary_boundary = clad_inner
    secondary_subdomain = clad_interface
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0 0 0'
    cylinder_axis_point_2 = '0 0 1'
    gap_flux_models = 'gas radiation contact'
    # Add displacements = 'ux uy uz' if your gap geometry is moving.
  []
[]

# The remaining mesh, variables, kernels, and boundary conditions are the same as in your
# existing mortar thermal contact input file.