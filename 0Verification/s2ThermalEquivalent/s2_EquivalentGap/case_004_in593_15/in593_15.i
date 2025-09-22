# === 参数研究案例 ===
# initial_T: 563.15
# 生成时间: 2025-08-18 00:42:26
# in：552.6、558.2、570.7、582.7、589.5
# out：553.7、560.8、582.8、602.5、611.8
# conda activate moose && dos2unix in593_15.i &&mpirun -n 4 /home/yp/projects/reproduction/reproduction-opt -i in593_15.i
initial_T = 582.7
initial_T_out = 602.5
LinePower = 70
EndTime = 110000000
pellet_nu = 0.345
# pellet_thermal_expansion_coef=1e-5#K-1
density_percent = 0.95
density_percent100 = '${fparse density_percent*100}'
pellet_density='${fparse density_percent*10980}'#10431.0*0.85#kg⋅m-3理论密度为10.980
# dt = 50000
fission_rate=3.6e19
grain_size = 10
#conda activate moose && dos2unix Complete2DQuarter.i&&mpirun -n 10 /home/yp/projects/reproduction/reproduction-opt -i Complete2DQuarter.i --mesh-only KAERI_HANARO_UpperRod1.e
#《《下面数据取自[1]Mechanism study and theoretical simulation on heat split phenomenon in dual-cooled annular fuel element》》

# 双冷却环形燃料几何参数 (单位：mm)
pellet_inner_diameter = 9.900         # 芯块内直径
pellet_outer_diameter = 14.100         # 芯块外直径               # 轴向长度(m)

# 网格控制参数n_azimuthal = 512时网格尺寸为6.8e-5m
mesh_size = 8e-5 #网格尺寸即可
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
growth_factor = 1.006       # 径向增长因子
# 计算半径参数 (转换为米)
pellet_inner_radius = '${fparse pellet_inner_diameter/2*1e-3}'
pellet_outer_radius = '${fparse pellet_outer_diameter/2*1e-3}'
[Mesh]
  [pellet1]
    type = AnnularMeshGenerator
    nr = ${n_radial_pellet}
    nt = ${n_azimuthal}
    rmin = ${pellet_inner_radius}
    rmax = ${pellet_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 20
    boundary_name_prefix = 'pellet'
  []
  [pellet]
    type = SubdomainIDGenerator
    input = pellet1
    subdomain_id = 1
  []
  [rename1]
    type = RenameBoundaryGenerator
    input = pellet
    old_boundary = 'pellet_rmin pellet_rmax'
    new_boundary = 'pellet_inner pellet_outer'
  []
  [cut_x]
    type = PlaneDeletionGenerator
    input = rename1
    point = '0 0 0'
    normal = '-1 0 0'  # 切割x>0区域
    new_boundary = 'y_axis'
  []
  [cut_y]
    type = PlaneDeletionGenerator
    input = cut_x
    point = '0 0 0'
    normal = '0 -1 0'  # 切割y>0区域
    new_boundary = 'x_axis'
  []
  [rename2]
    type = RenameBlockGenerator
    input = cut_y
    old_block = '1'
    new_block = 'pellet'
  []
[]




[GlobalParams]
    displacements = 'disp_x disp_y'
    out_of_plane_strain = strain_zz
[]
[AuxVariables]
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
  []
  [hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [x]
    initial_condition = 0.01
  []
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
  [vonMisesStress]
    type = ADRankTwoScalarAux
    variable = vonMises
    rank_two_tensor = stress
    execute_on = 'TIMESTEP_END'
    scalar_type = VonMisesStress
    # 不需要 index_i 和 index_j，因为我们使用 VonMisesStress 标量类型
  []
[]

[Variables]
    [disp_x]
    []
    [disp_y]
    []
    [T]
      initial_condition = ${initial_T}
    []
    [strain_zz]
    []

[]
[Kernels]
  #力平衡方程
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
    [./solid_z]
      type = ADWeakPlaneStress
      variable = strain_zz
    [../]
    #热传导方程
    [heat_conduction]
      type = ADHeatConduction
      variable = T
    []
    [hcond_time]
      type = ADHeatConductionTimeDerivative
      variable = T
    []
    [Fheat_source]
      type = ADMatHeatSource
      variable = T
      material_property = total_power
      block = pellet
    []
[]
[BCs]
  #固定平面
  [y_zero_on_y_plane]
    type = DirichletBC
    variable = disp_y
    boundary = 'x_axis'
    value = 0
  []
  [x_zero_on_x_plane]
    type = DirichletBC
    variable = disp_x
    boundary = 'y_axis'
    value = 0
  []

  #芯块包壳间隙压力
  [gap_pressure_fuel_x]
    type = Pressure
    variable = disp_x
    boundary = 'pellet_inner'
    factor = 1.0e6
    function = gap_pressure_inner
    use_displaced_mesh = true
  []
  [gap_pressure_fuel_y]
    type = Pressure
    variable = disp_y
    boundary = 'pellet_inner'
    factor = 1.0e6
    function = gap_pressure_inner
    use_displaced_mesh = true
  []

    #芯块包壳间隙压力
    [gap_pressure_pellet_outerx]
      type = Pressure
      variable = disp_x
      boundary = 'pellet_outer'
      factor = 1.0e6
      function = gap_pressure_outer
      use_displaced_mesh = true
    []
    [gap_pressure_pellet_outery]
      type = Pressure
      variable = disp_y
      boundary = 'pellet_outer'
      factor = 1.0e6
      function = gap_pressure_outer
      use_displaced_mesh = true
    []

  [coolant_bc_in]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'pellet_inner'
    T_infinity = ${initial_T}
    coefficient = coolant_equivalent_conductance_in#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  [coolant_bc_out]#对流边界条件
  type = ConvectiveFluxFunction
  variable = T
  boundary = 'pellet_outer'
  T_infinity = ${initial_T_out}
  coefficient = coolant_equivalent_conductance_out#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
[]
[Materials]
    #定义芯块热导率、密度、比热等材料属性
    [pellet_properties2]
      type = ADGenericConstantMaterial
      prop_names = 'density nu'
      prop_values = '${pellet_density} ${pellet_nu}'
      block = pellet
    []

    [pellet_thermal_conductivity] #新加的！！！！！！！！！！！！！！！！！！！！！！
      type = ADParsedMaterial
      property_name = thermal_conductivity #参考某论文来的，不是Fink-Lukuta model（非常复杂）
      coupled_variables = 'T'
      expression = '(100/(7.5408 + 17.692*T/1000 + 3.6142*(T/1000)^2) + 6400/((T/1000)^2.5)*exp(-16.35/(T/1000)))'
      block = pellet
    []
    [pellet_specific_heat]
      type = ADParsedMaterial
      property_name = specific_heat #Fink model
      coupled_variables = 'T'  # 需要在AuxVariables中定义Y变量
      expression = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
      block = pellet
    []
    [pellet_elastic_constants]
      type = ADParsedMaterial
      property_name = E #Fink model
      coupled_variables = 'T'  # 需要在AuxVariables中定义Y变量
      expression = '2.334*10^11*(1-2.752*(1-D))*(1-1.0915*10^(-4)*T)'
      constant_names = 'D'
      constant_expressions = '${density_percent}'
      block = pellet
    []
    # 肿胀应变函数
    [total_power]
      type = ADRimEffertPowerBurnup
      # property_name = total_power
      power_history = 'power_history'  # 声明使用的函数
      pellet_inner_radius = ${pellet_inner_radius}  # 为函数指定符号名称
      pellet_outer_radius = ${pellet_outer_radius}  # 直接使用函数符号进行计算
      initial_density = ${pellet_density}
      block = pellet
      output_properties = 'total_power burnup radial_power_shape'
      outputs = exodus
    []

    [swelling_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = swelling_coef
      coupled_variables = 'T'
      material_property_names = 'burnup'
      expression = '(${pellet_density}*5.577e-5*burnup + 1.101e-29*pow(2800-T,11.73)*exp(-0.0162*(2800-T))*(1-exp(-0.0178*${pellet_density}*burnup)))/3'
      block = pellet
    []
    # 密实化温度因子函数
    [CD_factor]
      type = ADParsedMaterial
      property_name = CD_factor
      coupled_variables = 'T'
      expression = 'if(T < 1023.15, 7.2-0.0086*(T-298.15),1)'
      block = pellet
    []
    [densification_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = densification_coef
      coupled_variables = 'T'
      material_property_names = 'CD_factor(T) burnup'
      expression = '0.001 * (exp(-4.605 * (burnup*1) / (CD_factor * 0.006024)) - 1)/3'# 0.6024是5000MWd/tU的转换系数
      block = pellet
    []
    [thermal_eigenstrain_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = thermal_eigenstrain_coef
      coupled_variables = 'T'
      expression = '-4.972e-4+7.107e-6*T+2.581e-9*T^2+1.14e-13*T^3'# 0.6024是5000MWd/tU的转换系数
      block = pellet
    []
    # 肿胀应变计算
    [swelling_eigenstrain]
      type = ADComputeVariableFunctionEigenstrain
      eigen_base = '1 1 1 0 0 0'
      prefactor = swelling_coef
      eigenstrain_name = swelling_eigenstrain
      block = pellet
    [../]
    # 肿胀应变计算
    [thermal_eigenstrain]
      type = ADComputeVariableFunctionEigenstrain
      eigen_base = '1 1 1 0 0 0'
      prefactor = thermal_eigenstrain_coef
      eigenstrain_name = thermal_eigenstrain
      block = pellet
    [../]

    # 密实化应变计算
    [densification_eigenstrain]
      type = ADComputeVariableFunctionEigenstrain
      eigen_base = '1 1 1 0 0 0'
      prefactor = densification_coef
      eigenstrain_name = densification_eigenstrain
      
      block = pellet
    [../]
      [creep_rate]
        type = UO2CreepRateExplicit
        temperature = T
        oxygen_ratio = x
        fission_rate = ${fparse fission_rate}
        theoretical_density = ${fparse density_percent100}
        grain_size = ${grain_size}
        vonMisesStress = vonMises
        block = pellet
      []
      [creep_eigenstrain]
        type = UO2CreepEigenstrain
        eigenstrain_name = creep_eigenstrain
        output_properties = 'effective_creep_strain'
        outputs = exodus
        block = pellet
      []
    [pellet_strain]
      type = ADComputePlaneSmallStrain 
      eigenstrain_names = 'thermal_eigenstrain swelling_eigenstrain densification_eigenstrain creep_eigenstrain'
      block = 'pellet'
    []
    [pellet_elasticity_tensor]
      type = ADComputeVariableIsotropicElasticityTensor
      youngs_modulus = E
      poissons_ratio = nu
      block = pellet
    []
    [stress]
      type = ADComputeLinearElasticStress
    []
    
[]
# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/(pellet_outer_radius^2-pellet_inner_radius^2)}' #新加的！！！！！！！！！！！！！！！！！！！！！！
[Functions]
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = PiecewiseLinear
    x = '0   ${EndTime}'
    y = '${LinePower} ${LinePower}'
    scale_factor = ${power_factor}         # 保持原有的转换因子
    # 论文中只给了线密度，需要化为体积密度
  []
  [gap_pressure_inner] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
        x = '0   ${EndTime}'
        y = '1.5 1.5'
    scale_factor = 1
  []
  [gap_pressure_outer] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
        x = '0   ${EndTime}'
        y = '1.5 1.5'
    scale_factor = 1
  []

  # reproduction/0Verification/s2ThermalEquivalent/s2_EquivalentGap/case_001_in563_15/in563_15.i
  # reproduction/0Verification/s2ThermalEquivalent/s1_OutputGap/case_001_in563_15/gap_conductance1/0.0.csv
  [coolant_equivalent_conductance_in]
    type = PiecewiseLinear
    data_file = '../.././s1_OutputGap/case_004_in593_15/gap_conductance1/0.75.csv'
    format = columns  # CSV 数据格式为列式
    x_index_in_file = 0  # 0.5innerBu 列（第3列，索引从0开始）
    y_index_in_file = 2  # 0.5innerConductance 列（第4列，索引从0开始）
    xy_in_file_only = false
    scale_factor = 1
  []
  [coolant_equivalent_conductance_out]
    type = PiecewiseLinear
    data_file = '../.././s1_OutputGap/case_004_in593_15/gap_conductance1/0.75.csv'
    format = columns  # CSV 数据格式为列式
    x_index_in_file = 0  # 0.5outerBu 列（第3列，索引从0开始）
    y_index_in_file = 3  # 0.5outerConductance 列（第4列，索引从0开始）
    xy_in_file_only = false
    scale_factor = 1
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
  # solve_type = 'NEWTON'
  # solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type   -snes_type        -snes_qn_type   -snes_qn_scale_type -snes_linesearch_type' 
  # petsc_options_value = 'lu         qn               lbfgs           jacobian           bt'
  petsc_options_iname = '-pc_type -ksp_type' 
  petsc_options_value = 'lu gmres' 
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true  # 每个时间步都重新计算缩放
  # reuse_preconditioner = true
  # reuse_preconditioner_max_linear_its = 20
  nl_max_its = 150
  nl_rel_tol = 5e-6 # 非线性求解的相对容差
  nl_abs_tol = 5e-6 # 非线性求解的绝对容差
  l_tol = 5e-8  # 线性求解的容差
  l_abs_tol = 5e-9 # 线性求解的绝对容差
  l_max_its = 500 # 线性求解的最大迭代次数
  accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
  dtmin = 1e-6
  dtmax = 100000
  end_time = ${EndTime} # 总时间24h

  fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
    growth_factor = 10
    cutback_factor = 0.5
    optimal_iterations = 100
    iteration_window = 25
  [../]
[]

[Postprocessors]
  # 已有
  [burnup_avg]
    type = ElementAverageValue
    variable = burnup
    block = pellet
    execute_on = 'TIMESTEP_END'
  []
  [pellet_inner_heat_rate]
    type = ADSideDiffusiveFluxIntegral
    variable = T
    boundary = pellet_inner
    diffusivity = thermal_conductivity
    execute_on = 'TIMESTEP_END'
  []
  [pellet_outer_heat_rate]
    type = ADSideDiffusiveFluxIntegral
    variable = T
    boundary = pellet_outer
    diffusivity = thermal_conductivity
    execute_on = 'TIMESTEP_END'
  []
  [pellet_inner_T_avg]
    type = SideAverageValue
    variable = T
    boundary = pellet_inner
    execute_on = 'TIMESTEP_END'
  []
  [pellet_outer_T_avg]
    type = SideAverageValue
    variable = T
    boundary = pellet_outer
    execute_on = 'TIMESTEP_END'
  []
  [pellet_inner_perimeter]
    type = AreaPostprocessor
    boundary = pellet_inner
    execute_on = 'initial timestep_end'
  []
  
  [pellet_outer_perimeter]
    type = AreaPostprocessor
    boundary = pellet_outer
    execute_on = 'initial timestep_end'
  []
  [h_eq_inner]
    type = ParsedPostprocessor
    expression = 'abs(pellet_inner_heat_rate) / (pellet_inner_perimeter*max(pellet_inner_T_avg - ${initial_T}, 1e-6))'
    pp_names = 'pellet_inner_heat_rate pellet_inner_T_avg pellet_inner_perimeter'
  []
  [h_eq_outer]
    type = ParsedPostprocessor
    expression = 'abs(pellet_outer_heat_rate) / (pellet_outer_perimeter*max(pellet_outer_T_avg - ${initial_T_out}, 1e-6))'
    pp_names = 'pellet_outer_heat_rate pellet_outer_T_avg pellet_outer_perimeter'
  []
  [pellet_T_max]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = pellet
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [my_checkpoint]
    type = Checkpoint
    time_step_interval = 5    # 每5个时间步保存
    num_files = 2            # 保留最近4个检查点
    wall_time_interval = 600 # 每10分钟保存一次（秒）
  []
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = false
  file_base = 'gap_conductance1/0.75'
  # csv = true
  [./csv]
    type = CSV
    precision = 5  # 默认保留1位小数
    execute_on = 'TIMESTEP_END'
    show = 'burnup_avg h_eq_inner h_eq_outer pellet_inner_T_avg pellet_T_max pellet_outer_T_avg'
  [../]
  
[]
