[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[Kernels]
  [TensorMechanics]
    displacements = 'disp_x disp_y'
  []
[]

[BCs]
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 't * 1e-4'  # 1e-4 strain rate
  []
[]

[Materials]
  [elastic]
    type = ComputeSmallDeformationStressMod
    displacements = 'disp_x disp_y'
    plasticity_model = creep_plasticity
    base_name = mech
  []
  
  [elasticity_model]
    type = ComputeSmallDeformationElasticityModelMod
    youngs_modulus = 2e11  # 200 GPa
    poissons_ratio = 0.3
    base_name = mech
  []
  
  [hardening_model]
    type = IsotropicHardeningModel
    yield_stress = 300e6  # 300 MPa
    hardening_modulus = 1e9  # 1 GPa
    base_name = mech
  []
  
  [creep_plasticity]
    type = SmallDeformationCreepPlasticitySeparatedMod
    hardening_model = hardening_model
    coefficient = 1e-20     # 蠕变系数 A
    exponent = 5.0          # 蠕变指数 n
    abs_tolerance = 1e-10   # 牛顿迭代容差
    max_iterations = 20     # 最大迭代次数
    base_name = mech
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       mumps'
  
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_max_its = 10
  
  start_time = 0.0
  end_time = 1000.0  # 1000 seconds
  dt = 100.0
  
  [TimeStepper]
    type = ConstantDT
    dt = 100.0
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    file_base = creep_plasticity_separated_out
  []
[]

[Postprocessors]
  [stress_xx]
    type = ElementAverageValue
    variable = mech_stress_xx
  []
  [stress_yy]
    type = ElementAverageValue
    variable = mech_stress_yy
  []
  [plastic_strain]
    type = ElementAverageValue
    variable = mech_effective_plastic_strain
  []
  [creep_strain]
    type = ElementAverageValue
    variable = mech_effective_creep_strain
  []
  [total_strain_yy]
    type = ElementAverageValue
    variable = mech_total_strain_yy
  []
[] 