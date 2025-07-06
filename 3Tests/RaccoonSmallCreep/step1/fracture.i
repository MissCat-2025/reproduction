a = 10e-3 #m 正方形边长
h = 0.05e-3 #m 网格尺寸 文中为0.02mm
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = ${a}
    ymax = ${a}
    nx = '${fparse ceil(a/h)}'
    ny = '${fparse ceil(a/h)}'
  []
[]

[Variables]
  [d]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [psip_active]
    order = CONSTANT
    family = MONOMIAL
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
  [pff_diff]
    type = ADPFFDiffusion
    variable = d
  []
  [pff_source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc'
    prop_values = '${l} ${Gc}'
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d*d))*(1-eta)+eta
    phase_field = d
    # material_property_names = 'a1'
    parameter_names = 'p a1 a2 a3 eta'
    parameter_values = '2 2 -0.5 0 1e-6'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd*d'
    phase_field = d
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*(psie_active+psip_active)'
    coupled_variables = 'd psie_active psip_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -snes_type'
  petsc_options_value = 'lu vinewtonrsls'

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  nl_max_its = 50
  dt = 3600
  end_time = 86400
  automatic_scaling = true

  abort_on_solve_fail = true
[]

[Outputs]
  print_linear_residuals = false
[]
