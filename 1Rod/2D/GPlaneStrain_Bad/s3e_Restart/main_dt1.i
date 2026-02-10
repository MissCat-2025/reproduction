# === 第二阶段：Solution Restart 示例 ===
# 从第一阶段 case_001_dt1 的 exodus 文件 2.1-2D-New2026/1.e
# 读取解场作为初始条件，然后继续做瞬态计算
# mpirun -n 11 /home/yp/projects/reproduction/reproduction-opt -i main_dt1.i
# ----------------- 顶层参数（保持和第一阶段一致，可以按需要改） -----------------
LinearPowerAll = 35   # 总线功率 kW/m

# 功率升高阶段的时间节点
x      = 5                            # 功率从 0% 提到 x%
xTime  = '${fparse x*120}'            # s，从 0% 到 x% 的时间
# LinearPower0 = '${fparse LinearPowerAll*x*0.01}'
PowerTime     = 60                    # s，从 x% 到 100% 的时间
PowerTimeTotal = '${fparse xTime+PowerTime}'  # s，从 0% 到 100% 总时间
WatchTime      = 60                   # s，满功率后的观察时间
endTime        = '${fparse PowerTimeTotal+WatchTime+100}'

dt1   = 1
dtmin = 1e-3

initial_T = 550
pellet_nu = 0.316
pellet_thermal_expansion_coef = 1e-5
pellet_critical_fracture_strength = 6.0e7
density_percent = 0.95
fission_rate = 1.2e19
grain_size   = 10
pellet_critical_energy = 5

pellet_density = '${fparse density_percent*10980}'

length_scale_paramete = 6e-5

w         = 1
mesh_size = '${fparse 2*5e-5}'
pellet_outer_radius = 4.2e-3
# n_elems_azimuthal    = '${fparse 2*ceil((3.1415*pellet_outer_radius/mesh_size)/2^w)}'
# n_elems_radial_pellet = '${fparse int((pellet_outer_radius/mesh_size)/2^w)}'

creep_relative_tolerance = 1e-6
creep_absolute_tolerance = '${fparse creep_relative_tolerance*0.1}'

power_factor = '${fparse 1000*1/3.1415926/pellet_outer_radius/pellet_outer_radius}'

# 相场断裂参数
m   = 4
a2  = 0.5396842
a3  = 0
ksi = 2

# ----------------- Mesh：从第一阶段 exodus 读取网格 -----------------
[Mesh]
  # 使用第一阶段主算例的 exodus 网格
  file = '/home/yp/projects/reproduction/1Tasks/step2/3_2026/s2PowerOn/s1Dt/case_001_dt1/2.1-2D-New2026/1.e'
  # distribution = serial
[]

# ----------------- Problem：不做 checkpoint 重启，只是一个新的算例 -----------------
[Problem]
  kernel_coverage_check   = false
  material_coverage_check = false
[]

# ----------------- MultiApps：继续使用裂纹子算例 -----------------
[MultiApps]
  [fracture]
    type        = TransientMultiApp
    input_files = 'sub_dt1.i'
    cli_args = 'l=${length_scale_paramete};mesh_size=${mesh_size};Gc=${pellet_critical_energy};sigma0=${pellet_critical_fracture_strength};m=${m};w=${w};a2=${a2};a3=${a3};ksi=${ksi};pellet_outer_radius=${pellet_outer_radius};xTime=${xTime};PowerTimeTotal=${PowerTimeTotal};endTime=${endTime};dtmin=${dtmin};dt1=${dt1}'
    execute_on  = 'TIMESTEP_END'
    clone_parent_mesh=true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app   = 'fracture'
    variable         = d
    source_variable  = d
  []
  [to_ALL]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app     = 'fracture'
    variable         = 'psie_active a1 stress_I'
    source_variable  = 'psie_active a1 stress_I'
  []
[]

[GlobalParams]
  displacements       = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

# ----------------- AuxVariables -----------------
[AuxVariables]
  [hoop_stress]
    order  = CONSTANT
    family = MONOMIAL
  []
  [stress_I]
    order  = CONSTANT
    family = MONOMIAL
  []
  [radial_stress]
    order  = CONSTANT
    family = MONOMIAL
  []
  [d]
    block = pellet
    initial_from_file_var = d
        initial_from_file_timestep  = LATEST
  []
  [sigma0_field]
    family = MONOMIAL
    order  = CONSTANT
    # 从 exodus 读取上一阶段的随机强度场
    initial_from_file_var = sigma0_field
  []
[]

# ----------------- AuxKernels -----------------
[AuxKernels]
  [hoop_stress]
    type          = ADRankTwoScalarAux
    variable      = hoop_stress
    rank_two_tensor = stress
    scalar_type   = HoopStress
    point1        = '0 0 0.01'
    point2        = '0 0 -0.01'
    execute_on    = 'TIMESTEP_END'
  []
  [stress_I]
    type          = ADRankTwoScalarAux
    scalar_type   = MaxPrincipal
    rank_two_tensor = stress
    variable      = stress_I
    selected_qp   = 0
  []
  [radial_stress]
    type          = ADRankTwoScalarAux
    rank_two_tensor = stress
    variable      = radial_stress
    scalar_type   = RadialStress
    point1        = '0 0 0'
    point2        = '0 0 1'
  []
  [copy_sigma0]
    type         = ADMaterialRealAux
    variable     = sigma0_field
    property     = sigma0
    execute_on   = 'initial'
    block        = pellet
  []
[]

# ----------------- 原始变量 + 从 exodus 读取初始解 -----------------
[Variables]
  [disp_x]
    initial_from_file_var = disp_x
    initial_from_file_timestep  = LATEST
  []
  [disp_y]
    initial_from_file_var = disp_y
        initial_from_file_timestep  = LATEST
  []
  [T]
    initial_from_file_var = T
        initial_from_file_timestep  = LATEST
  []
  [x]
    initial_from_file_var = x
        initial_from_file_timestep  = LATEST
  []
  [strain_zz]
    initial_from_file_var = strain_zz
        initial_from_file_timestep  = LATEST
  []
[]

# ----------------- Kernels -----------------
[Kernels]
  # 力平衡
  [solid_x]
    type      = ADStressDivergenceTensors
    variable  = disp_x
    component = 0
  []
  [solid_y]
    type      = ADStressDivergenceTensors
    variable  = disp_y
    component = 1
  []
  [solid_z]
    type     = ADWeakPlaneStress
    variable = strain_zz
  []

  # 热传导
  [heat_conduction]
    type     = ADHeatConduction
    variable = T
  []
  [hcond_time]
    type     = ADHeatConductionTimeDerivative
    variable = T
  []

  [Fheat_source]
    type              = ADMatHeatSource
    variable          = T
    material_property = total_power
    block             = pellet
  []

  # 化学平衡
  [time_derivative]
    type     = ADTimeDerivative
    variable = x
    block    = pellet
  []
  [complex_diffusion]
    type       = ADComplexDiffusionKernel
    variable   = x
    temperature = T
    block      = pellet
  []
[]

# ----------------- BCs -----------------
[BCs]
  [y_zero_on_y_plane]
    type     = DirichletBC
    variable = disp_y
    boundary = 'yplane'
    value    = 0
  []
  [x_zero_on_x_plane]
    type     = DirichletBC
    variable = disp_x
    boundary = 'xplane'
    value    = 0
  []
  [coolant_bc_out]
    type        = ConvectiveFluxFunction
    variable    = T
    boundary    = 'pellet_outer'
    T_infinity  = T_infinity_out
    coefficient = gap_conductance
  []
  [gap_pressure_fuel_x]
    type      = Pressure
    variable  = disp_x
    boundary  = 'pellet_outer'
    factor    = 1e6
    function  = gap_pressure
  []
  [gap_pressure_fuel_y]
    type      = Pressure
    variable  = disp_y
    boundary  = 'pellet_outer'
    factor    = 1e6
    function  = gap_pressure
  []
[]

# ----------------- Materials -----------------
[Materials]
  [pellet_properties]
    type        = ADGenericConstantMaterial
    prop_names  = 'l Gc density nu'
    prop_values = '${length_scale_paramete} ${pellet_critical_energy} ${pellet_density} ${pellet_nu}'
    block       = pellet
  []
  [sigma0_mat]
    type              = ADParsedMaterial
    property_name     = sigma0
    coupled_variables = 'sigma0_field'
    expression        = 'sigma0_field'
    block             = pellet
  []
  [pellet_thermal_conductivity]
    type              = ADParsedMaterial
    property_name     = thermal_conductivity
    coupled_variables = 'T'
    expression        = '(1)*(100/(7.5408 + 17.692*T/1000 + 3.6142*(T/1000)^2) + 6400/((T/1000)^2.5)*exp(-16.35/(T/1000)))'
    block             = pellet
  []
  [pellet_specific_heat]
    type              = ADParsedMaterial
    property_name     = specific_heat
    coupled_variables = 'T x'
    expression        = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (x+2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
    block             = pellet
  []
  [pellet_elastic_constants]
    type              = ADParsedMaterial
    property_name     = E
    coupled_variables = 'T'
    expression        = '2.334*10^11*(1-2.752*(1-density/10960))*(1-1.0915*10^(-4)*T)'
    constant_names        = 'density'
    constant_expressions  = '${pellet_density}'
    block             = pellet
  []
  [total_power]
    type                = ADRimEffertPowerBurnupRod
    power_history       = 'power_history'
    pellet_outer_radius = ${pellet_outer_radius}
    output_properties   = 'total_power burnup radial_power_shape'
    outputs             = exodus
  []
  [pellet_thermal_eigenstrain]
    type                    = ADComputeThermalExpansionEigenstrain
    eigenstrain_name        = thermal_eigenstrain
    stress_free_temperature = ${initial_T}
    thermal_expansion_coeff = ${pellet_thermal_expansion_coef}
    temperature             = T
    block                   = pellet
  []
  [D_fickian]
    type              = ADParsedMaterial
    property_name     = D_fickian
    coupled_variables = 'x T d'
    expression        = '(1-0.99*d)*pow(10, -9.386 - 4260/(T) + 0.0012*T*x + 0.00075*T*log10(1+2/(x)))'
    block             = pellet
  []
  [D_soret]
    type                  = ADDerivativeParsedMaterial
    property_name         = D_soret
    coupled_variables     = 'x T d'
    material_property_names = 'D_fickian(x,T,d)'
    expression            = 'D_fickian * x * (-1380.8 - 134435.5*exp(-x/0.0261)) / ((2.0 + x)/(2.0 * (1.0 - 3.0*x) * (1.0 - 2.0*x)) * 8.314 * T * T)'
    block                 = pellet
  []
  [creep]
    type                 = UO2CreepRateBaseJ2Creep
    phase_field          = d
    degradation_function = g
    temperature          = T
    oxygen_ratio         = x
    fission_rate         = ${fparse fission_rate}
    grain_size           = ${grain_size}
    use_transition_stress   = false
    use_transient_creep     = false
    use_three_shear_modulus = false
    relative_tolerance      = ${creep_relative_tolerance}
    absolute_tolerance      = ${creep_absolute_tolerance}
    output_properties       = 'effective_creep_strain'
    outputs                 = exodus
  []
  [pellet_strain]
    type              = ADComputePlaneSmallStrain
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [crack_geometric]
    type           = CrackGeometricFunction
    property_name  = alpha
    expression     = 'ksi*d+(1-ksi)*d*d'
    parameter_names  = 'ksi'
    parameter_values = '${ksi}'
    phase_field    = d
  []
  [a1]
    type                  = ADDerivativeParsedMaterial
    property_name         = a1
    material_property_names = 'Gc E l sigma0'
    expression            = '4*E*Gc/sigma0/sigma0/l/3.14159'
    output_properties     = 'a1'
    outputs               = exodus
  []
  [degradation]
    type              = RationalDegradationFunction
    property_name     = g
    expression        = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))*(1-eta)+eta
    phase_field       = d
    material_property_names = 'a1'
    parameter_names   = 'p a2 a3 eta'
    parameter_values  = '${m} ${a2} ${a3} 1e-6'
  []
  [Elasticity]
    type                = IsotropicElasticity
    youngs_modulus      = E
    poissons_ratio      = nu
    phase_field         = d
    degradation_function= g
    kinematic_assumption= PLANE_STRESS
    decomposition       = SPECTRAL
    use_threshold       = false
    use_history_max     = false
    output_properties   = 'psie_active'
    tensile_strength    = sigma0
    outputs             = exodus
  []
  [stress]
    type            = ComputeCreepPlasticityDeformationStress
    elasticity_model= Elasticity
    creep_model     = creep
  []
[]

# ----------------- Functions -----------------
[Functions]
  [gap_conductance]
    type         = PiecewiseLinear
    x            = '0 ${endTime}'
    y            = '3400 3400'
    scale_factor = 1
  []
  [power_history]
    type = PiecewiseLinear
    x    = '0.0 ${xTime} ${PowerTimeTotal} ${endTime}'
    y    = '${LinearPowerAll} ${LinearPowerAll} ${LinearPowerAll} ${LinearPowerAll}'
    scale_factor = ${power_factor}
  []
  [gap_pressure]
    type = PiecewiseLinear
    x    = '0 ${endTime}'
    y    = '2.5 2.5'
  []
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < ${xTime}, 1,
                  if(t < (${xTime}+3), 1.5,
                  if(t < ${PowerTimeTotal},${dt1},
                  if(t < ${endTime},${dt1},${dt1}))))'
  []
  [T_infinity_out]
    type = PiecewiseLinear
    x    = '0 ${endTime}'
    y    = '${initial_T} ${initial_T}'
    scale_factor = 1
  []
[]

# ----------------- Executioner -----------------
[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  petsc_options_iname  = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value  = 'asm      100                lu           NONZERO'

  line_search          = 'bt'
  automatic_scaling    = true
  compute_scaling_once = true

  nl_max_its = 50
  nl_rel_tol = 5e-6
  nl_abs_tol = 5e-7

  dtmin    = ${dtmin}
  end_time = ${endTime}

  fixed_point_rel_tol  = 1e-8
  fixed_point_abs_tol  = 1e-10
  accept_on_max_fixed_point_iteration = true

  [TimeStepper]
    type     = FunctionDT
    function = dt_limit_func
  []
[]

# ----------------- Adaptivity（仅时间步自适应，不做初始 adapt） -----------------
[Adaptivity]
  marker     = marker
  max_h_level = ${w}
  [Markers]
    [marker]
      type              = PhasePiledFractureHSMarker
      von_mises_variable = stress_I
      sigma0            = sigma0
      x1                = 0.000001
      x2                = 0.005
      xmax              = 0.05
      y1                = 0.45
      y2                = 0.5
      variable          = d
      timeD             = 3
      timeStress        = 5
      d_change_threshold    = 0.01
      stress_change_threshold = 1e6
    []
  []
[]

# ----------------- Outputs -----------------
[Outputs]
  exodus = true
  print_linear_residuals = false
  file_base = 'RodRestart_e/1'
[]