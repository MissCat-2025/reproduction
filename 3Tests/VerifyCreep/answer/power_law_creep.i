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
[]
[AuxKernels]
  [vonMisesStress]
    type = RankTwoScalarAux
    variable = vonMises
    rank_two_tensor = stress
    execute_on = 'TIMESTEP_END'
    scalar_type = VonMisesStress
    # 不需要 index_i 和 index_j，因为我们使用 VonMisesStress 标量类型
  []
  [./creep_aux]
    type = MaterialRealAux
    property = effective_creep_strain
    variable = effective_creep
  [../]
[]
# [Physics/SolidMechanics/QuasiStatic]
#   [all]
#     strain = FINITE
#     incremental = true
#     add_variables = true
#     generate_output = 'stress_yy creep_strain_xx creep_strain_yy creep_strain_zz elastic_strain_yy'
#   []
# []

[Functions]
  [top_pull]
    type = PiecewiseLinear
    x = '0 1'
    y = '1 1'
  []
[]

[Kernels]
  [heat]
    type = Diffusion
    variable = temp
  []
  [heat_ie]
    type = TimeDerivative
    variable = temp
  []
  [solid_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [solid_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [solid_z]
    type = StressDivergenceTensors
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
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'power_law_creep'
    tangent_operator = elastic
    # output_properties = 'strain_increment'
    # outputs = exodus
  []
  [power_law_creep]
    type = PowerLawCreepStressUpdate
    coefficient = 1.0e-15
    n_exponent = 4
    activation_energy = 3.0e5
    temperature = temp
  []
  [strain]
    type = ComputeFiniteStrain
  []
  # [stress]
  #   type = ComputeLinearElasticStress
  # []
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
