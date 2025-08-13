# === 参数研究案例 ===
# end_time = 1.10e+6
# mesh_size: 4.00e-5
# Gc: 10
# 生成时间: 2025-03-28 11:16:02

# conda activate moose && dos2unix 1_main.i&& dos2unix 1_Sub.i &&mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i 1_main.i
pellet_density=10431.0#10431.0*0.85#kg⋅m-3
pellet_nu = 0.345
pellet_thermal_expansion_coef=1e-5#K-1
Gc = 5 #断裂能
density_percent = 0.95
theoretical_density = '${fparse density_percent*100}'
# dt = 50000
fission_rate=3.0e19
grain_size = 10
initial_T = 393.15
pellet_critical_fracture_strength='${fparse 1.7*10^8*(1-2.62*(1-density_percent))^0.5*exp(-1590/8.314/initial_T)}'#10431.0*0.85#kg⋅m-3理论密度为10.980
length_scale_paramete=4e-5
mesh_size = 8e-5 #网格尺寸即可
pellet_critical_energy=${fparse Gc} #J⋅m-2
#《《下面数据取自[1]Thermomechanical Analysis and Irradiation Test of Sintered Dual-Cooled Annular pellet》》

# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
pellet_inner_diameter = 10.291         # 芯块内直径mm
pellet_outer_diameter = 14.627         # 芯块外直径mm

# length_scale_paramete=${fparse mesh_size}
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
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
    growth_r = 1.006
    boundary_id_offset = 10
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

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = '1_Sub.i'
    cli_args = 'l=${length_scale_paramete}'
    execute_on = 'TIMESTEP_END'
    # 强制同步参数
    sub_cycling = false          # 禁止子循环
    catch_up = false             # 禁止追赶步
    max_failures = 0             # 严格同步模式
  []
[]

[Transfers]
  [from_d]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_ALL]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active a1 hoop_stress Gc'
    source_variable = 'psie_active a1 hoop_stress Gc'
  []
[]


[GlobalParams]
    displacements = 'disp_x disp_y'
    out_of_plane_strain = strain_zz
[]
[AuxVariables]
  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./creep_hoop_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
  []
  [d]
    block = pellet
  []

  [sigma0_field]
    family = MONOMIAL
    order = CONSTANT
    [InitialCondition]
      type = WeibullIC
      scale = ${pellet_critical_fracture_strength}
      shape = 50
      location = 0.0
      seed = 0
      block = pellet
    []
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
    [./creep_strain]
      type = ADRankTwoAux
      variable = creep_hoop_strain
      rank_two_tensor = creep_eigenstrain
      index_i = 2
      index_j = 2  # zz分量对应环向
      execute_on = 'TIMESTEP_END'
      block = pellet
    [../]
    [copy_sigma0]
      type = ADMaterialRealAux
      variable = sigma0_field
      property = sigma0
      execute_on = 'initial'
      block = pellet
    []
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
    [strain_zz]
    []
    [T]
      initial_condition = 393.15
    []
    [x]
      initial_condition = 0.01
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
    #化学平衡方程
    [time_derivative]
      type = ADTimeDerivative
      variable = x
      block = pellet
    []
    [complex_diffusion]
      type = ADComplexDiffusionKernel
      variable = x
      temperature = T
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

  # #芯块包壳间隙压力
  # [gap_pressure_fuel_x]
  #   type = Pressure
  #   variable = disp_x
  #   boundary = 'pellet_inner pellet_outer'
  #   factor = 2.0e6
  #   use_displaced_mesh = false
  # []
  # [gap_pressure_fuel_y]
  #   type = Pressure
  #   variable = disp_y
  #   boundary = 'pellet_inner pellet_outer'
  #   factor = 2.0e6
  #   use_displaced_mesh = false
  # []
  [coolant_bc]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'pellet_inner pellet_outer'
    T_infinity = 393.15
    coefficient = gap_conductance#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
[]
[Materials]
    #定义芯块热导率、密度、比热等材料属性
    [pellet_properties2]
      type = ADGenericConstantMaterial
      prop_names = 'l density nu'
      prop_values = '${length_scale_paramete} ${pellet_density} ${pellet_nu}'
      block = pellet
    []
    [pellet_Gc]
      type = ADDerivativeParsedMaterial
      property_name = Gc
      material_property_names = 'burnup'
      expression = 'Gc1*(1-0.9*burnup/0.0012)'
      constant_names = 'Gc1'
      constant_expressions = '${pellet_critical_energy}'
      output_properties = 'Gc'
      outputs = exodus
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
    [total_power]
      type = ADDerivativeParsedMaterial
      property_name = total_power
      functor_names = 'power_history'  # 声明使用的函数
      functor_symbols = 'P'  # 为函数指定符号名称
      expression = 'P'  # 直接使用函数符号进行计算
      derivative_order = 2  # 需要计算导数时指定
      block = pellet
      output_properties = 'total_power'
      outputs = exodus
    []
    [burnup]
      type = ADBurnupMaterial
      total_power = total_power
      initial_density = ${pellet_density}
      block = pellet
      output_properties = 'burnup'
      outputs = exodus
    []
    # 为临界断裂强度生成威布尔分布
    [sigma0_mat]
      type = ADDerivativeParsedMaterial
      property_name = sigma00
      coupled_variables = 'sigma0_field T'
      expression = 'if(T < 1000,sigma0_field/F*1.7*10^8*(1-2.62*(1-D))^0.5*exp(-1590/8.314/T),sigma0_field/F*(1-2.62*(1-D))^0.5*1.4040832037*10^8)'  # 直接使用辅助变量的值
      constant_names = 'F D'
      constant_expressions = '${pellet_critical_fracture_strength} ${density_percent}'
      block = pellet
    []
    [sigma0]
      type = ADDerivativeParsedMaterial
      property_name = sigma0
      material_property_names = 'sigma00 burnup'
      expression = 'sigma00*(1-0.7*burnup/0.001)'
      output_properties = 'sigma0'
      outputs = exodus
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
      coupled_variables = 'T x'  # 需要在AuxVariables中定义Y变量
      expression = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (x+2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
      block = pellet
    []
    [pellet_thermal_eigenstrain]
      type = ADComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = 393.15
      thermal_expansion_coeff = ${pellet_thermal_expansion_coef}
      temperature = T
      block = pellet
    []

    #化学相关
    [D_fickian]
      type = ADParsedMaterial
      property_name = D_fickian
      coupled_variables = 'x T d'
      expression = '(1-0.99*d)*pow(10, -9.386 - 4260/(T) + 0.0012*T*x + 0.00075*T*log10(1+2/(x)))'
      block = pellet
    []
    [D_soret]
      type = ADDerivativeParsedMaterial
      property_name = D_soret
      coupled_variables = 'x T d'
      material_property_names = 'D_fickian(x,T,d)'
      expression = 'D_fickian * x * (-1380.8 - 134435.5*exp(-x/0.0261)) / ((2.0 + x)/(2.0 * (1.0 - 3.0*x) * (1.0 - 2.0*x)) * 8.314 * T * T)'
      block = pellet
    []
    # # # # 蠕变相关
  [creep_rate]
    type = UO2CreepRateExplicit
    temperature = T
    oxygen_ratio = x
    fission_rate = ${fparse fission_rate}
    theoretical_density = ${fparse theoretical_density}
    grain_size = ${grain_size}
    vonMisesStress = vonMises
    block = pellet
  []
    # 蠕变特征应变
    [creep_eigenstrain]
      type = UO2CreepEigenstrain
      eigenstrain_name = creep_eigenstrain
      block = pellet
    []
    [pellet_strain]
      type = ADComputePlaneSmallStrain
      eigenstrain_names = 'thermal_eigenstrain creep_eigenstrain'
      block = pellet
    []
    [a1]
      type = ADDerivativeParsedMaterial
      property_name = a1
      material_property_names = 'Gc E l sigma0'
      expression = '4*E*Gc/sigma0/sigma0/3.14159/l'
      output_properties = 'a1'
      outputs = exodus
      block = pellet
    []
    [degradation]
      type = RationalDegradationFunction
      property_name = g
      expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d))*(1-eta)+eta
      phase_field = d
      material_property_names = 'a1'
      parameter_names = 'p a2 eta'
      parameter_values = '2 -0.5 1e-6'
      block = pellet
    []
    [crack_geometric]
      type = CrackGeometricFunction
      property_name = alpha
      expression = '2*d-d*d'
      phase_field = d
      block = pellet
    []
    [pellet_elasticity]
      type = SmallDeformationH
      youngs_modulus = E
      poissons_ratio = nu
      tensile_strength = sigma0
      # fracture_energy = Gc
      phase_field = d
      degradation_function = g
      output_properties = 'psie_active'
      outputs = exodus
      block = pellet
    []
    [pellet_stress]
      type = ComputeSmallDeformationStress
      elasticity_model = pellet_elasticity
      block = pellet
    []
[]
# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/(pellet_outer_radius^2-pellet_inner_radius^2)}' #新加的！！！！！！！！！！！！！！！！！！！！！！
[Functions]
  [gap_conductance]
    type = PiecewiseLinear
    x = '0 1100000'
    y = '3500 2000'
    scale_factor = 1         # 保持原有的转换因子
  []
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
  type = PiecewiseLinear #论文的功率历史
  x = '0 800000 950000 1000000   '
  y = '0 105 100 0'
  scale_factor = ${power_factor}         # 保持原有的转换因子
  # 论文中只给了线密度，需要化为体积密度
  []
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 800000, 5000,
                   if(t < 950000, 500, 1000))'
  []
[]


[Executioner]
    type = Transient
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
    petsc_options_value = 'lu superlu_dist gmres'
    accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
    automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
    compute_scaling_once = true  # 每个时间步都重新计算缩放
    nl_max_its = 30
    nl_rel_tol = 1e-6 # 非线性求解的相对容差
    nl_abs_tol = 1e-7 # 非线性求解的绝对容差
    l_tol = 1e-7  # 线性求解的容差
    l_abs_tol = 1e-8 # 线性求解的绝对容差
    l_max_its = 150 # 线性求解的最大迭代次数
    dtmin = 500
    # dt = ${dt}
    end_time = 1100000
    [TimeStepper]
      type = FunctionDT
      function = dt_limit_func
    []
[]
[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = 2
  [Markers]
    [marker]
      type = PhasePiledFractureMarker
      von_mises_variable = hoop_stress
      x1 = 1e-4 #d变量小于x1时，标记为粗网格
      x2 = 0.05 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.1 #d变量大于xmax时，一定是细网格
      y1 = 1e7 #vonMises应力小于y1时，标记为粗网格
      y2 = 1e7 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      enable_time_check = true
      time_steps = 3
      d_change_threshold = 1e-2
    []
  []
[]
[Outputs]
  exodus = true
  print_linear_residuals = false
  file_base = '2D/2D'
[]