[Problem]
    kernel_coverage_check = false
    material_coverage_check = false
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 1000
  ymax = 1000
  elem_type = QUAD
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
    [sigma0_field]
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
    []
    [source]
      type = ADPFFSource
      variable = d
      free_energy = psi
    []
  []

  [Materials]
    [fracture_properties]
      type = ADGenericConstantMaterial
      prop_names = 'l Gc E0'
      prop_values = '${l} ${Gc} ${E0}'
    []
    [sigma0_mat]
      type = ADParsedMaterial
      property_name = sigma0
      coupled_variables = 'sigma0_field'
      expression = 'sigma0_field'
    []
    [degradation]
        type = RationalDegradationFunction
        property_name = g
        expression = (1-d)^p/((1-d)^p+(1.5*E0*Gc/sigma0^2)/l*d*(1+a2*d))*(1-eta)+eta
        phase_field = d
        material_property_names = 'Gc sigma0 l E0'
        parameter_names = 'p a2 eta'
        parameter_values = '2 2 1e-6'
      []
      [crack_geometric]
        type = CrackGeometricFunction
        property_name = alpha
        expression = 'd'
        phase_field = d
      []  
      [psi]
        type = ADDerivativeParsedMaterial
        property_name = psi
        expression = 'alpha*Gc/c0/l+g*psie_active'
        coupled_variables = 'd psie_active'
        material_property_names = 'alpha(d) g(d) Gc c0 l'
      []
[]

[Executioner]
    type = Transient
    solve_type = PJFNK
    petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type  -snes_type'
    petsc_options_value = '201                hypre    boomeramg  vinewtonrsls'  
    accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
    nl_max_its = 20
    nl_rel_tol = 1e-8   # 非线性求解的相对容差
    nl_abs_tol = 1e-10  # 非线性求解的绝对容差
    l_tol = 1e-8        # 线性求解的容差
    l_abs_tol = 1e-10   # 线性求解的绝对容差
    start_time = 0.0
    num_steps = 1000

    dt = 0.5  
[]

[Outputs]
  print_linear_residuals = false
[]