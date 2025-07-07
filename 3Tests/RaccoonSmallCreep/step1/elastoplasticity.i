#  mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i elastoplasticity.i
#Phase-field fracture modeling for creep crack
#4.1. Uniaxial tension model

a = 10e-3 #m 正方形边长
h = 0.05e-3 #m 网格尺寸 文中为0.02mm
l = 0.2e-3 #m 


Gc = 160e3 #J/m^2
E = 140e9 #Pa
nu = 0.3
sigma_0 = 170e6 #Pa
psic = '${fparse E*Gc/sigma_0/sigma_0}'
H = '${fparse E/100}'


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
    variable = 'psie_active psip_active'
    source_variable = 'psie_active psip_active'
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
  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./hoop_stress]
    type = ADRankTwoScalarAux
    variable = hoop_stress
    rank_two_tensor = stress
    scalar_type = HoopStress
    point1 = '0 0 0'        # 圆心坐标
    point2 = '0 0 -0.0178'        # 定义旋转轴方向（z轴）
    execute_on = 'TIMESTEP_END'
  [../]
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
    prop_names = 'l Gc psic sigma_0 E nu H'
    prop_values = '${l} ${Gc} ${psic} ${sigma_0} ${E} ${nu} ${H}'
  []
  # [a1]
  #   type = ADDerivativeParsedMaterial
  #   property_name = a1
  #   material_property_names = 'Gc E l sigma0'
  #   expression = '4*E*Gc/sigma0/sigma0/l/3.14'
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
    tensile_strength = sigma_0
    outputs = exodus
  []
  [linear_hardening]
    type = LinearHardening
    phase_field = d
    yield_stress = sigma_0
    hardening_modulus = H
    degradation_function = g
    output_properties = 'psip_active'
    outputs = exodus
  []
  [strain]
    type = ADComputeSmallStrain
  []
  # [plasticity]
  #   type = J2Plasticity_C
  #   hardening_model = linear_hardening
  #   output_properties = 'effective_plastic_strain'
  #   outputs = exodus
  # []
  [creep]
    type = CoupledStressStrainCreepRate
    hardening_model = linear_hardening
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
    boundary = 'bottom'
    value = 0
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
    factor = ${Pressure1}
    use_displaced_mesh = true
  []
[]

[Postprocessors]
  [d]
    type = ElementAverageValue
    variable = d
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # [ep]
  #   type = ADElementAverageMaterialProperty
  #   mat_prop = effective_plastic_strain
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  nl_max_its = 50
  dt = 1
  end_time = 86400

  automatic_scaling = true
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]
