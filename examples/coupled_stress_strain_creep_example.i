# Coupled Stress-Strain Creep Example
# Creep model: dot_epsilon_c_eq = A1*sigma_eq^n1*epsilon_c_eq^m1 + A2*sigma_eq^n2*epsilon_c_eq^m2 + A3*sigma_eq^n3*epsilon_c_eq^m3
# This model captures the interaction between stress and accumulated creep strain

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
    type = CoupledStressStrainCreepRate
    elasticity_model = elasticity_model
    
    # 第一项参数: A1 * sigma_eq^n1 * epsilon_c_eq^m1
    # 主导高应力下的蠕变
    A1 = 1.0e-20       # 第一项系数
    n1 = 5.0           # 第一项应力指数
    m1 = 0.5           # 第一项应变指数
    
    # 第二项参数: A2 * sigma_eq^n2 * epsilon_c_eq^m2
    # 主导中等应力下的蠕变
    A2 = 5.0e-15       # 第二项系数
    n2 = 3.0           # 第二项应力指数
    m2 = 1.0           # 第二项应变指数
    
    # 第三项参数: A3 * sigma_eq^n3 * epsilon_c_eq^m3
    # 主导低应力下的蠕变
    A3 = 1.0e-10       # 第三项系数
    n3 = 1.0           # 第三项应力指数
    m3 = 2.0           # 第三项应变指数
    
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
    function = '2e-3*t'  # 施加位移载荷
  [../]
[]

[Functions]
  [./ramp_function]
    type = ParsedFunction
    value = 'if(t<100, t/100, 1)'
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
  end_time = 2000.0  # 2000秒
  dt = 20.0
  
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 20.0
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 10
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    file_base = coupled_creep_out
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./csv]
    type = CSV
    file_base = coupled_creep_out
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
  [./von_mises_stress]
    type = ElementAverageValue
    variable = von_mises_stress
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[] 