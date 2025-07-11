# conda activate moose && dos2unix elastoplasticity.i && dos2unix fracture.i &&mpirun -n 12 /home/yp/projects/reproduction/reproduction-opt -i elastoplasticity.i
#Phase-field fracture modeling for creep crack
#4.1. Uniaxial tension model
A = 9.541e-35
n = 11.5
HHH = 1 #时间换算，
X = 1e-6 #Mpa换算
A11=${fparse A*X^n/HHH}

Gc = 160e3 #J/m^2
E = 165e9 #Pa
nu = 0.3
sigma_quf = 356e6 #Pa 屈服强度
l = 0.08e-3 #m

# Pressure1 = 8000 #Pa 290, 308, 349, 359 and 366 MPa

[GlobalParams]
  displacements = 'disp_x disp_y'
  # volumetric_locking_correction = false 
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'fracture.i'
    cli_args = 'Gc=${Gc};l=${l} '
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
  type = FileMesh
  file = CT_model.e
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
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l sigma_quf E nu'
    prop_values = '${l} ${sigma_quf} ${E} ${nu}'
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
    type = PowerLawCreepRate
    phase_field = d
    degradation_function = g
    A = ${A11}
    n = ${n}
    # 可选参数（如果需要调整）
    max_iter = 5
    tolerance = 1e-9
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
  [upper_load]
    type = Pressure
    boundary = 'upper_hole'
    variable = disp_y
    factor = 2.83e8
  []
  [lower_load]
    type = Pressure
    boundary = 'lower_hole'
    variable = disp_y
    factor = 2.83e8
  []
  # 约束刚体运动 - 固定右边界的x方向位移
  [fixed_x]
      type = DirichletBC
      variable = disp_x
      boundary = 'right'
      value = 0.0
  []
  
  # # 固定一个点的y方向位移（防止刚体运动）
  [fixed_y]
      type = DirichletBC
      variable = disp_y
      boundary = 'right'
      value = 0.0
  []
[]
[Postprocessors]
  [total_effective_creep_strain]
    type = ADElementIntegralMaterialProperty
    mat_prop = effective_creep_strain
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  nl_max_its = 50
  dt = 2.5
  end_time = 70000000

  automatic_scaling = true
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]

