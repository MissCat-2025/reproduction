# 1x1x1 unit cube with uniform pressure on top face
# conda activate moose && dos2unix power_law_creep.i  && mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i power_law_creep.i
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
[]

[Variables]
  [temp]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1000.0
  []
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [disp_z]
    order = FIRST
    family = LAGRANGE
  []
[]
[AuxVariables]
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_creep]
    order = CONSTANT
    family = MONOMIAL
  []
  [d]
    initial_condition = 0.0000
  []
[]
[AuxKernels]
  [vonMisesStress]
    type = ADRankTwoScalarAux
    variable = vonMises
    rank_two_tensor = stress
    execute_on = 'TIMESTEP_END'
    scalar_type = VonMisesStress
    # 不需要 index_i 和 index_j，因为我们使用 VonMisesStress 标量类型
  []
  [./creep_aux]
    type = ADMaterialRealAux
    property = effective_creep_strain
    variable = effective_creep
  [../]
[]

[Functions]
  [top_pull]
    type = PiecewiseLinear
    x = '0 1'
    y = '1 1'
  []
[]

[Kernels]
  [heat]
    type = ADDiffusion
    variable = temp
  []
  [heat_ie]
    type = ADTimeDerivative
    variable = temp
  []
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
  []
[]

[BCs]
  [u_top_pull]
    type = Pressure
    variable = disp_y
    boundary = top
    factor = -10.0e6
    function = top_pull
  []
  [u_bottom_fix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [u_yz_fix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [u_xy_fix]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []
  [temp_fix]
    type = DirichletBC
    variable = temp
    boundary = 'bottom top'
    value = 1000.0
  []
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'E nu l'
    prop_values = '2e11 0.3 0.1'
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = "(1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d*d))*(1-eta)+eta"
    phase_field = d
    # material_property_names = 'a1'
    parameter_names = 'p a1 a2 a3 eta'
    parameter_values = '2 2 -0.5 0 1e-6'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [hencky]
    type = IsotropicElasticity
    youngs_modulus = E
    poissons_ratio = nu
    phase_field = d
    degradation_function = g
    use_threshold = false
    output_properties = 'psie_active'
    tensile_strength = sigma_la
    outputs = exodus
  []

  [strain]
    type = ADComputeSmallStrain
  []

  [creep]
    type = PowerLawCreepRate2
    phase_field = d
    degradation_function = g
    # 可选参数（如果需要调整）
    A = 1.0e-15
    n = 4
    Q = 3.0e5
    temperature = temp
    relative_tolerance = 1e-10 #蠕变的相对残差
    absolute_tolerance = 1e-15 #蠕变的绝对残差
    # max_iter = 3
    # tolerance = 1e-8
    output_properties = 'effective_creep_strain psic_active'
    outputs = exodus 
  []
  [stress]
    type = ComputeCreepPlasticityDeformationStress
    elasticity_model = hencky
    creep_model = creep
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'

  line_search = 'none'

  l_max_its = 20
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-5
  start_time = 0.0
  end_time = 10.0
  num_steps = 100
  dt = 0.05
[]

[Outputs]
  exodus = true
[]
