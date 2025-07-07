# Power Law Creep Example with Inconel 625 Superalloy
# Creep model: dot_epsilon_c_eq = A * sigma_eq^n
# Material constants for Inconel 625: A = 9.541e-35, n = 11.56

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmax = 1.0
  ymax = 1.0
  zmax = 1.0
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_z]
    type = StressDivergenceTensors
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Materials]
  [./elasticity_model]
    type = SmallDeformationIsotropicElasticity
    lambda = 115.38e9  # Pa (Inconel 625)
    mu = 84.62e9       # Pa (Inconel 625)
  [../]
  
  [./creep_model]
    type = PowerLawCreepRate
    elasticity_model = elasticity_model
    # Power law creep parameters for Inconel 625
    A = 9.541e-35    # 蠕变系数 A
    n = 11.56        # 蠕变指数 n
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-8
    max_iterations = 50
  [../]
  
  [./stress]
    type = SmallDeformationStress
    displacements = 'disp_x disp_y disp_z'
    creep_model = creep_model
  [../]
[]

[BCs]
  [./fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  
  [./load_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '1e-3*t'  # 施加位移载荷
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'
  
  line_search = none
  
  l_max_its = 50
  l_tol = 1e-8
  nl_max_its = 20
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  
  start_time = 0.0
  end_time = 1000.0  # 1000秒
  dt = 10.0
  
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10.0
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 10
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    file_base = powerlaw_creep_out
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./csv]
    type = CSV
    file_base = powerlaw_creep_out
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Postprocessors]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./creep_strain_xx]
    type = ElementAverageValue
    variable = creep_strain_xx
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./effective_creep_strain]
    type = ElementAverageValue
    variable = effective_creep_strain
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[] 