# conda activate moose && dos2unix 2D.i &&mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i 2D.i
pellet_nu = 0.316
pellet_thermal_expansion_coef=1e-5#K-1
density_percent = 0.95
theoretical_density = '${fparse density_percent*100}'

clad_density=6.59e3#kg⋅m-3
# clad_elastic_constants=7.52e10#Pa
clad_nu = 0.33
# clad_specific_heat=264.5
# clad_thermal_conductivity = 16
clad_thermal_expansion_coef=5.0e-6#K-1
# dt = 50000
# pellet_critical_fracture_strength=6.0e7#Pa
# fission_rate=2.0e19
# (已验证)这第一步就是测试生成网格文件KAERI_HANARO_UpperRod1.e   注意：这文件生成的是完整的3D模型，而不是1/4模型3D模型
# 语法要求：仅仅为了生成网格文件：Run with --mesh-only:
#https://mooseframework.inl.gov/source/meshgenerators/ConcentricCircleMeshGenerator.html
#conda activate moose && dos2unix Complete2DQuarter.i&&mpirun -n 10 /home/yp/projects/reproduction/reproduction-opt -i Complete2DQuarter.i --mesh-only KAERI_HANARO_UpperRod1.e
#《《下面数据取自[1]Thermomechanical Analysis and Irradiation Test of Sintered Dual-Cooled Annular pellet》》

# 双冷却环形燃料几何参数 (单位：mm)
inclad_inner_diameter = 9.0      # 内包壳内直径
inclad_outer_diameter = 10.14    # 内包壳外直径
pellet_inner_diameter = 10.271         # 芯块内直径
pellet_outer_diameter = 14.685         # 芯块外直径
outclad_inner_diameter = 14.76    # 外包壳内直径
outclad_outer_diameter = 15.9     # 外包壳外直径                   # 轴向长度(m)

# 网格控制参数n_azimuthal = 512时网格尺寸为6.8e-5m
n_radial_inner_clad = 3    # 内包壳径向单元数
mesh_size = 8e-5 #网格尺寸即可
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
n_radial_outer_clad = 3    # 外包壳径向单元数
growth_factor = 1.006       # 径向增长因子
# 计算半径参数 (转换为米)
inner_clad_inner_radius = '${fparse inclad_inner_diameter/2*1e-3}'
inner_clad_outer_radius = '${fparse inclad_outer_diameter/2*1e-3}'
pellet_inner_radius = '${fparse pellet_inner_diameter/2*1e-3}'
pellet_outer_radius = '${fparse pellet_outer_diameter/2*1e-3}'
outer_clad_inner_radius = '${fparse outclad_inner_diameter/2*1e-3}'
outer_clad_outer_radius = '${fparse outclad_outer_diameter/2*1e-3}'
[Mesh]
  [inner_clad1]
    type = AnnularMeshGenerator
    nr = ${n_radial_inner_clad}
    nt = ${n_azimuthal}
    rmin = ${inner_clad_inner_radius}
    rmax = ${inner_clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 10
    boundary_name_prefix = 'inclad'
  []
  [inner_clad]
    type = SubdomainIDGenerator
    input = inner_clad1
    subdomain_id = 1
  []
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
    subdomain_id = 2
  []
  [outer_clad1]
    type = AnnularMeshGenerator
    nr = ${n_radial_outer_clad}
    nt = ${n_azimuthal}
    rmin = ${outer_clad_inner_radius}
    rmax = ${outer_clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 30
    boundary_name_prefix = 'outclad'
  []
  [outer_clad]
    type = SubdomainIDGenerator
    input = outer_clad1
    subdomain_id = 3
  []
  [combine]
    type = CombinerGenerator
    inputs = 'inner_clad pellet outer_clad'
  []
  [rename1]
    type = RenameBoundaryGenerator
    input = combine
    old_boundary = 'inclad_rmin inclad_rmax pellet_rmin pellet_rmax outclad_rmin outclad_rmax'
    new_boundary = 'inclad_inner inclad_outer pellet_inner pellet_outer outclad_inner outclad_outer'
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
    old_block = '1 2 3'
    new_block = 'inclad pellet outclad'
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
      initial_condition = 293.15
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
      type = HeatSource
      variable = T
      function = power_history
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
    boundary = 'pellet_inner inclad_outer'
    factor = 1.0e6
    function = gap_pressure_inner
    use_displaced_mesh = true
  []
  [gap_pressure_fuel_y]
    type = Pressure
    variable = disp_y
    boundary = 'pellet_inner inclad_outer'
    factor = 1.0e6
    function = gap_pressure_inner
    use_displaced_mesh = true
  []

    #芯块包壳间隙压力
    [gap_pressure_pellet_outerx]
      type = Pressure
      variable = disp_x
      boundary = 'pellet_outer outclad_inner'
      factor = 1.0e6
      function = gap_pressure_outer
      use_displaced_mesh = true
    []
    [gap_pressure_pellet_outery]
      type = Pressure
      variable = disp_y
      boundary = 'pellet_outer outclad_inner'
      factor = 1.0e6
      function = gap_pressure_outer
      use_displaced_mesh = true
    []

  [coolant_bc_in]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'inclad_inner'
    T_infinity = 313.15
    coefficient = coolant_conductance_in#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  [coolant_bc_out]#对流边界条件
  type = ConvectiveFluxFunction
  variable = T
  boundary = 'outclad_outer'
  T_infinity = 313.15
  coefficient = coolant_conductance_out#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  #冷却剂压力
  [colden_pressure_fuel_x]
    type = Pressure
    variable = disp_x
    boundary = 'inclad_inner outclad_outer'
    factor = 15.5e6
    use_displaced_mesh = true
  []
  [colden_pressure_fuel_y]
    type = Pressure
    variable = disp_y
    boundary = 'inclad_inner outclad_outer'
    factor = 15.5e6
    use_displaced_mesh = true
  []
[]
[Materials]
    #定义芯块热导率、密度、比热等材料属性
    [pellet_properties2]
      type = ADGenericConstantMaterial
      prop_names = 'density nu'
      prop_values = '${theoretical_density} ${pellet_nu}'
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
    [E_clad_elastic_constants]
      type = ADParsedMaterial
      property_name = E_clad #Fink model
      coupled_variables = 'T'  # 需要在AuxVariables中定义Y变量
      expression = '-46.67e6 * T + 1.09e11'
      block = 'inclad outclad' 
    []

    [clad_properties]
      type = ADGenericConstantMaterial
      prop_names = ' density nu'
      prop_values = '${clad_density} ${clad_nu}'
      block = 'inclad outclad' 
    []
    [clad_specific_heat_T]
      type = ADDerivativeParsedMaterial
      property_name = specific_heat
      coupled_variables = 'T'
      functor_names = 'Cp'
      functor_symbols = 'Cp'
      expression = 'Cp'
      block = 'inclad outclad'
    []
    [clad_thermal_conductivity]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = thermal_conductivity
      coupled_variables = 'T'
      expression = '7.51+2.09e-2*T-1.45e-5*T^2 + 7.67e-9*T^3'
      block = 'inclad outclad'
    []



    [clad_thermal_eigenstrain]
      type = ADComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = 293.15
      thermal_expansion_coeff = ${clad_thermal_expansion_coef}
      temperature = T
      block = 'inclad outclad'
    []

    [pellet_thermal_eigenstrain]
      type = ADComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = 293.15
      thermal_expansion_coeff = ${pellet_thermal_expansion_coef}
      temperature = T
      block = pellet
    []
    # 肿胀应变函数
    [total_power]
      type = ADDerivativeParsedMaterial
      property_name = total_power
      functor_names = 'power_history'  # 声明使用的函数
      functor_symbols = 'P'  # 为函数指定符号名称
      expression = 'P'  # 直接使用函数符号进行计算
      derivative_order = 1  # 需要计算导数时指定
      block = pellet
      output_properties = 'total_power'
      outputs = exodus
    []
    [burnup]
      type = ADBurnupMaterial
      total_power = total_power
      initial_density = ${theoretical_density}
      block = pellet
      output_properties = 'burnup'
      outputs = exodus
    []

    [swelling_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = swelling_coef
      coupled_variables = 'T'
      material_property_names = 'burnup'
      expression = '(${theoretical_density}*5.577e-5*burnup + 1.101e-29*pow(2800-T,11.73)*exp(-0.0162*(2800-T))*(1-exp(-0.0178*${theoretical_density}*burnup)))/3'
      block = pellet
    []
    # 密实化温度因子函数
    [densification_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = densification_coef
      coupled_variables = 'T'
      material_property_names = 'CD_factor burnup'
      expression = '0.005 * (exp(-4.605 * (burnup*1) / (CD_factor * 0.006024)) - 1)/3'# 0.6024是5000MWd/tU的转换系数
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

    # 密实化应变计算
    [densification_eigenstrain]
      type = ADComputeVariableFunctionEigenstrain
      eigen_base = '1 1 1 0 0 0'
      prefactor = densification_coef
      eigenstrain_name = densification_eigenstrain
      block = pellet
    [../]

    [pellet_strain]
      type = ADComputePlaneSmallStrain 
      eigenstrain_names = 'thermal_eigenstrain swelling_eigenstrain densification_eigenstrain'
      block = 'pellet'
    []
    [clad_strain]
      type = ADComputePlaneSmallStrain
      eigenstrain_names = thermal_eigenstrain
      block = 'inclad outclad'
    []
    [pellet_elasticity_tensor]
      type = ADComputeVariableIsotropicElasticityTensor
      youngs_modulus = E
      poissons_ratio = nu
      block = pellet
    []
    [clad_elasticity_tensor]
      type = ADComputeVariableIsotropicElasticityTensor
      youngs_modulus = E_clad
      poissons_ratio = nu
      block = 'inclad outclad'
  []
    [stress]
      type = ADComputeLinearElasticStress
    []
[]

[ThermalContact]
  [./thermal_contact1]
    type = GapHeatTransfer
    variable = T
    primary = inclad_outer
    secondary = pellet_inner
    emissivity_primary = 0.8
    emissivity_secondary = 0.8
    gap_conductivity = 6000
    quadrature = true
    gap_geometry_type = CYLINDER
    cylinder_axis_point_1 = '0 0 -0.0001'
    cylinder_axis_point_2 = '0 0 0.0001'
  [../]
    [./thermal_contact2]
      type = GapHeatTransfer
      variable = T
      primary = outclad_inner
      secondary = pellet_outer
      emissivity_primary = 0.8
      emissivity_secondary = 0.8
      gap_conductivity = 6000
      quadrature = true
      gap_geometry_type = CYLINDER
      cylinder_axis_point_1 = '0 0 -0.0001'
      cylinder_axis_point_2 = '0 0 0.0001'
    [../]
    []
# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/(pellet_outer_radius^2-pellet_inner_radius^2)}' #新加的！！！！！！！！！！！！！！！！！！！！！！
[Functions]
  [coolant_conductance_in]
    type = PiecewiseLinear
    x = '0 10000000'
    y = '3200 3200'
    scale_factor = 1         # 保持原有的转换因子
  []
  [coolant_conductance_out]
    type = PiecewiseLinear
    x = '0 10000000'
    y = '3200 3200'
    scale_factor = 1         # 保持原有的转换因子
  []
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = PiecewiseLinear
    data_file = power_history.csv    # 创建一个包含上述数据的CSV文件，数据为<s,w/m>
    format = columns                 # 指定数据格式为列式
    scale_factor = ${power_factor}         # 保持原有的转换因子
    # 论文中只给了线密度，需要化为体积密度
  []
  [Cp] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = PiecewiseLinear
    data_file = clad_Cp.csv    # 创建一个包含上述数据的CSV文件，数据为<s,w/m>
    format = columns                 # 指定数据格式为列式
    scale_factor = 1         # 保持原有的转换因子
    # 论文中只给了线密度，需要化为体积密度
  []
  [gap_pressure_inner] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
        x = '0   10000000'
        y = '1.1 5' #论文的1.1
    scale_factor = 1
  []
  [gap_pressure_outer] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
        x = '0   10000000'
        y = '1.1 5'
    scale_factor = 1
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  # solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
  # solve_type = 'NEWTON'
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type   -snes_type        -snes_qn_type   -snes_qn_scale_type -snes_linesearch_type' 
  petsc_options_value = 'lu         qn               lbfgs           jacobian           bt'
  # petsc_options_iname = '-pc_type -ksp_type' 
  # petsc_options_value = 'lu gmres' 
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
  dtmax = 10000
  end_time = 10000000 # 总时间24h

  fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
    growth_factor = 1.1
    cutback_factor = 0.5
    optimal_iterations = 100
    iteration_window = 25
  [../]
[]

[Postprocessors]
  # 已有
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
    expression = 'abs(pellet_inner_heat_rate) / (pellet_inner_perimeter*max(pellet_inner_T_avg - 313.15, 1e-6))'
    pp_names = 'pellet_inner_heat_rate pellet_inner_T_avg pellet_inner_perimeter'
  []
  [h_eq_outer]
    type = ParsedPostprocessor
    expression = 'abs(pellet_outer_heat_rate) / (pellet_outer_perimeter*max(pellet_outer_T_avg - 313.15, 1e-6))'
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
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = false
  file_base = '2D-NoFracture/2D'
  # csv = true
  [./csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
    show = 'h_eq_inner h_eq_outer pellet_inner_T_avg pellet_T_max pellet_outer_T_avg'
  [../]
[]
