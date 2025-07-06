#  mpirun -n 12 /home/yp/projects/reproduction/reproduction-opt -i elastoplasticity.i


a = 250

Gc = 160e3 #J/m^2
l = 0.06e-3 #m

E = 140e9 #Pa
nu = 0.3

sigma_0 = 170e6 #Pa
psic = '${fparse E*Gc/sigma_0/sigma_0}'
H = '${fparse E/100}'

n = 11.56
creep_coef = 9.541e-35

Pressure = 8e3 #Pa

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = false 
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'fracture.i'
    cli_args = 'a=${a};psic=${psic};Gc=${Gc};l=${l} '
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
    dim = 3
    xmax = ${a}
    ymax = ${a}
    zmax = ${a}
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [d]
  []
  [stress]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress]
    type = ADRankTwoAux
    variable = stress
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    execute_on = 'INITIAL TIMESTEP_END'
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
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
  []
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc psic sigma_0 E nu'
    prop_values = '${l} ${Gc} ${psic} ${sigma_0} ${E} ${nu}'
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+(Gc/psic*xi/c0/l)*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta
    phase_field = d
    material_property_names = 'Gc psic xi c0 l '
    parameter_names = 'p a2 a3 eta '
    parameter_values = '2 -0.5 0 1e-6'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [hencky]
    type = SmallDeformationIsotropicElasticityThresholdMod
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
  [J2]
    type = SmallDeformationJ2PowerLawCreepMod
    hardening_model = linear_hardening
    coefficient = ${creep_coef}
    exponent = ${n}
  []
  [stress]
    type = ComputeSmallDeformationStressMod
    elasticity_model = hencky
    plasticity_model = J2
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  [ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    function = '${rate}*t'
    preset = false
  []
[]

[Postprocessors]
  [stress]
    type = ElementAverageValue
    variable = stress
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d]
    type = ElementAverageValue
    variable = d
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ep]
    type = ADElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  nl_max_its = 50
  dt = 0.0005
  end_time = 0.1

  automatic_scaling = true

  fixed_point_max_its = 100
  fixed_point_rel_tol = 1e-08
  fixed_point_abs_tol = 1e-10
  accept_on_max_fixed_point_iteration = true
  abort_on_solve_fail = true
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]
