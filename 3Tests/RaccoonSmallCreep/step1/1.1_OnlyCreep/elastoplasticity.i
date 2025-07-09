# conda activate moose && dos2unix elastoplasticity.i && dos2unix fracture.i &&mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i elastoplasticity.i
#Phase-field fracture modeling for creep crack
#4.1. Uniaxial tension model
A1 = 4.33e-32
A2 = 3.6e-25
A3 = 1.5e-25
n1 = 10.08
n2 = 8.9
n3 = 10.0
HHH = 1 #时间换算，
X = 1e-6 #Mpa换算
A11=${fparse A1*X^n1/HHH}
A22=${fparse A2*X^n2/HHH}
A33=${fparse A3*X^n3/HHH}

a = 10e-3 #m 正方形边长
h = 0.1e-3 #m 网格尺寸 文中为0.02mm
l = 0.2e-3 #m 


Gc = 160e3 #J/m^2
E = 140e9 #Pa
nu = 0.3
sigma_la = 588e6 #Pa 抗拉强度


Pressure1 = 290e6 #Pa 290, 308, 349, 359 and 366 MPa

[GlobalParams]
  displacements = 'disp_x disp_y'
  # volumetric_locking_correction = false 
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'fracture.i'
    cli_args = 'a=${a};Gc=${Gc};l=${l} '
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_fracture]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_fracture]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active psic_active'
    source_variable = 'psie_active psic_active'
  []
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = ${a}
    ymax = ${a}
    nx = '${fparse ceil(a/h)}'
    ny = '${fparse ceil(a/h)}'
  []
  # [nodeset1]
  #   type = ParsedGenerateNodeset
  #   input = gmg
  #   expression = 'x > 0.004995 & x < 0.005005'
  #   epsilon =  1e-9
  #   new_nodeset_name = leftpoint
  # []

[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [d]
  []
    [vonMises]
      order = CONSTANT
      family = MONOMIAL
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
[]

[Kernels]
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
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l sigma_la E nu'
    prop_values = '${l} ${sigma_la} ${E} ${nu}'
  []
  # [Gc1]
  #   type = ADDerivativeParsedMaterial
  #   property_name = Gc
  #   coupled_variables = 'd'
  #   expression = '${Gc} * (1-d)*(1-d)' 
  # []
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
  [hencky]
    type = IsotropicElasticityThreshold
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
    type = CoupledStressStrainCreepRate
    phase_field = d
    degradation_function = g
    A1 = ${A11}
    A2 = ${A22}
    A3 = ${A33}
    # 可选参数（如果需要调整）
    max_iter = 3
    tolerance = 1e-8
    output_properties = 'effective_creep_strain psic_active'
    outputs = exodus 
  []
  [stress]
    type = ComputeCreepPlasticityDeformationStress
    elasticity_model = hencky
    creep_model = creep
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.000
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [gap_pressure_fuel_x]
    type = Pressure
    variable = disp_y
    boundary = 'top'
    factor = -${Pressure1}
  []
[]
[Postprocessors]
  # 计算整个域的有效蠕变应变积分
  [total_effective_creep_strain]
    type = ADElementIntegralMaterialProperty
    mat_prop = effective_creep_strain
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # [total_effective_plastic_strain]
  #   type = ADElementIntegralMaterialProperty
  #   mat_prop = effective_plastic_strain
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []
  # 计算塑性应变和蠕变应变的总和
  # [total_effective_inelastic_strain]
  #   type = LinearCombinationPostprocessor
  #   pp_names = 'total_effective_plastic_strain total_effective_creep_strain'
  #   pp_coefs = '1.0 1.0'
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  nl_max_its = 50
  dt = 10
  end_time = 70000000

  automatic_scaling = true
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]

