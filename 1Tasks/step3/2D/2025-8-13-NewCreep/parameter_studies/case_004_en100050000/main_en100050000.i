# === 参数研究案例 ===
# endTime: 100050000
# 生成时间: 2025-08-26 11:52:26


# conda activate moose && dos2unix 2D_Main.i&& dos2unix 2D_Sub.i &&mpirun -n 2 /home/yp/projects/reproduction/reproduction-opt -i 2D_Main.i
initial_T = 293.15
initial_T_in = 570.7
initial_T_out = 582.8
LinearPower = 90
LinearPower0_2 = '${fparse LinearPower*0.2}'
endTime = 100050000
endTime__50000 = '${fparse endTime-5000}'
endTime__100000 = '${fparse endTime-150000}'
Pressure1 = 1.5
Pressure2 = 1.8
pellet_nu = 0.345
pellet_thermal_expansion_coef=1e-5#K-1
density_percent = 0.95
Gc = 5#断裂能

# pellet_critical_fracture_strength=9.0e7#Pa
# pellet_critical_fracture_strength='${fparse 1.7*10^8*(1-2.62*(1-density_percent))^0.5*exp(-1590/8.314/initial_T)}'#10431.0*0.85#kg⋅m-3理论密度为10.980
#(1-2.62*(1-0.98))^0.5 = 0.9734,exp(-1590/8.314/293.15)=0.521,pellet_critical_fracture_strength = 83.9MPa
#(1-2.62*(1-0.98))^0.5 = 0.9734,exp(-1590/8.314/293.15)=0.727,pellet_critical_fracture_strength = 120MPa
#(1-2.62*(1-0.95))^0.5 = 0.9322,exp(-1590/8.314/293.15)=0.521,pellet_critical_fracture_strength = 82.5MPa
fission_rate=7.20e19
grain_size =10
pellet_critical_energy=${fparse Gc} #J⋅m-2
pellet_density='${fparse density_percent*10980}'#10431.0*0.85#kg⋅m-3理论密度为10.980
GcX = 0.5
largestPoreSize = 30

#《《下面数据取自[1]Thermomechanical Analysis and Irradiation Test of Sintered Dual-Cooled Annular pellet》》

# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
pellet_inner_diameter = 10.291         # 芯块内直径mm
pellet_outer_diameter = 14.627         # 芯块外直径mm
length_scale_paramete = 4e-5
mesh_size = 8.00e-5 #网格尺寸即可
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
    input_files = 'sub_en100050000.i'
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
    variable = 'psie_active a1 hoop_stress Gc sigma0'
    source_variable = 'psie_active a1 hoop_stress Gc sigma0'
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
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
  []
  [d]
    block = pellet
    initial_condition = 0.0
  []
  [sigma0_F_Density]
    family = MONOMIAL
    order = CONSTANT
    [InitialCondition]
      type = WeibullICDensity
      scale = 6.0e7
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
    point1 = '0 0 0.01'        # 圆心坐标
    point2 = '0 0 -0.01'        # 定义旋转轴方向（z轴）
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
    [x]
      initial_condition = 0.01
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
  #芯块包壳间隙压力
  [gap_pressure_fuel_x]
    type = Pressure
    variable = disp_x
    boundary = 'pellet_inner'
    factor = 1.0e6
    function = gap_pressure_inner
    # use_displaced_mesh = true
  []
  [gap_pressure_fuel_y]
    type = Pressure
    variable = disp_y
    boundary = 'pellet_inner'
    factor = 1.0e6
    function = gap_pressure_inner
    # use_displaced_mesh = true
  []

    #芯块包壳间隙压力
    [gap_pressure_pellet_outerx]
      type = Pressure
      variable = disp_x
      boundary = 'pellet_outer'
      factor = 1.0e6
      function = gap_pressure_outer
      # use_displaced_mesh = true
    []
    [gap_pressure_pellet_outery]
      type = Pressure
      variable = disp_y
      boundary = 'pellet_outer'
      factor = 1.0e6
      function = gap_pressure_outer
      # use_displaced_mesh = true
    []

  [coolant_bc_in]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'pellet_inner'
    T_infinity = T_infinity_in
    coefficient = gap_conductance_in#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  [coolant_bc_out]#对流边界条件
  type = ConvectiveFluxFunction
  variable = T
  boundary = 'pellet_outer'
  T_infinity = T_infinity_out
  coefficient = gap_conductance_out#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
[]
[]
[Materials]
    #定义芯块热导率、密度、比热等材料属性
    [pellet_properties2]
      type = ADGenericConstantMaterial
      prop_names = 'l nu'
      prop_values = '${length_scale_paramete} ${pellet_nu}'
      block = pellet
    []
    [strain_adjusted_density]
      type = ADStrainAdjustedDensity
      displacements = 'disp_x disp_y'
      strain_free_density = ${pellet_density}
      block = pellet
      output_properties = 'density'
      outputs = exodus
    []
    [largestPoreSize]
      type = ADDerivativeParsedMaterial
      property_name = largestPoreSize
      material_property_names = 'burnup'
      expression = 'PS*(1+burnup/0.06)'
      constant_names = 'PS'
      constant_expressions = '${largestPoreSize}'
      block = pellet
    []
    [porosity]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = porosity
      coupled_variables = 'T'
      material_property_names = 'density'
      expression = '(1-(2.58-5.8e-4*T)*(1-density/10960))/(1-(2.58-5.8e-4*T)*(1-${density_percent}))'
      output_properties = 'porosity'
      outputs = exodus
      block = pellet
    []
    [sigma0]
      type = ADDerivativeParsedMaterial
      property_name = sigma0
      coupled_variables = 'T sigma0_F_Density'
      material_property_names = 'porosity(T) largestPoreSize'
      expression = 'sigma0_F_Density*626e6*(largestPoreSize+0.5*GS)^(-0.5)*exp(-0.057*porosity)'  # 直接使用辅助变量的值
      constant_names = 'GS'
      constant_expressions = '${grain_size}'
      output_properties = 'sigma0'
      outputs = exodus
      block = pellet
    []
    # [sigma0_mat]
    #   type = ADDerivativeParsedMaterial
    #   property_name = sigma00
    #   coupled_variables = 'sigma0_field T'
    #   material_property_names = 'density'
    #   expression = 'if(T < 1000,sigma0_field/F*1.7*10^8*(1-2.62*(1-density/10960))^0.5*exp(-1590/8.314/T),sigma0_field/F*(1-2.62*(1-density/10960))^0.5*1.4040832037*10^8)'  # 直接使用辅助变量的值
    #   constant_names = 'F'
    #   constant_expressions = '${pellet_critical_fracture_strength}'
    #   block = pellet
    # []


    [pellet_thermal_conductivity] #新加的！！！！！！！！！！！！！！！！！！！！！！
      type = ADParsedMaterial
      property_name = thermal_conductivity #参考某论文来的，不是Fink-Lukuta model（非常复杂）
      coupled_variables = 'T'
      expression = '(1)*(100/(7.5408 + 17.692*T/1000 + 3.6142*(T/1000)^2) + 6400/((T/1000)^2.5)*exp(-16.35/(T/1000)))'
      block = pellet
    []
    [pellet_specific_heat]
      type = ADParsedMaterial
      property_name = specific_heat #Fink model
      coupled_variables = 'T x'  # 需要在AuxVariables中定义Y变量
      expression = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (x+2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
      block = pellet
    []
    [pellet_elastic_constants]
      type = ADParsedMaterial
      property_name = E #Fink model
      coupled_variables = 'T'  # 需要在AuxVariables中定义Y变量
      material_property_names = 'density'
      expression = '2.334*10^11*(1-2.752*(1-density/10960))*(1-1.0915*10^(-4)*T)'
      block = pellet
    []
    # 肿胀应变函数
    [total_power]
      type = ADRimEffertPowerBurnup
      power_history = 'power_history'  # 声明使用的函数
      pellet_inner_radius = ${pellet_inner_radius}  # 为函数指定符号名称
      pellet_outer_radius = ${pellet_outer_radius}  # 直接使用函数符号进行计算
      block = pellet
      output_properties = 'total_power burnup radial_power_shape'
      outputs = exodus
    []
    [pellet_Gc]
      type = ADDerivativeParsedMaterial
      property_name = Gc
      material_property_names = 'burnup'
      expression = 'Gc1*(1-GcX*burnup/0.1)'
      constant_names = 'Gc1 GcX'
      constant_expressions = '${pellet_critical_energy} ${GcX}'
      output_properties = 'Gc'
      outputs = exodus
      block = pellet
    []

    [pellet_thermal_eigenstrain]
      type = ADComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = ${initial_T}
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
    [a1]
      type = ADDerivativeParsedMaterial
      property_name = a1
      material_property_names = 'Gc E l sigma0'
      expression = '1.5*E*Gc/sigma0/sigma0/l'
      output_properties = 'a1'
      outputs = exodus
      block = pellet
    []
    [creep]
      type = UO2CreepRateBaseJ2Creep
      phase_field = d
      degradation_function = g
      temperature = T
      oxygen_ratio = x
      fission_rate = ${fparse fission_rate}
      grain_size = ${grain_size}
      # 相场断裂相关参数
      use_transition_stress = false
      # use_transient_creep = true
      use_three_shear_modulus = false

      relative_tolerance = 1e-8 #蠕变的相对残差
      absolute_tolerance = 1e-10 #蠕变的绝对残差

      output_properties = 'effective_creep_strain psic_active'
      outputs = exodus 
    []
    [swelling_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = swelling_coef
      coupled_variables = 'T'
      material_property_names = 'burnup density'
      expression = '(density*5.577e-5*burnup + 1.101e-29*pow(2800-T,11.73)*exp(-0.0162*(2800-T))*(1-exp(-0.0178*density*burnup)))/3'
      block = pellet
    []

    [CD_factor]
      type = ADParsedMaterial
      property_name = CD_factor
      coupled_variables = 'T'
      expression = 'if(T < 1023.15, 7.2-0.0086*(T-298.15),1)'
      block = pellet
    []
    # 密实化温度因子函数
    [densification_coef]
      type = ADDerivativeParsedMaterial  # 改为ADParsedMaterial
      property_name = densification_coef
      coupled_variables = 'T'
      material_property_names = 'CD_factor(T) burnup'
      expression = '0.01 * (exp(-4.605 * (burnup*1) / (CD_factor * 0.006024)) - 1)/3'# 0.6024是5000MWd/tU的转换系数
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
      block = pellet
    []
    [crack_geometric]
      type = CrackGeometricFunction
      property_name = alpha
      expression = 'd'
      phase_field = d
      block = pellet
    []

    [degradation]
      type = RationalDegradationFunction
      property_name = g
      expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d))*(1-eta)+eta
      phase_field = d
      material_property_names = 'a1'
      parameter_names = 'p a2 eta'
      parameter_values = '2 2 1e-6'
      block = pellet
    []
    [Elasticity]
      type = IsotropicElasticity
      youngs_modulus = E
      poissons_ratio = nu
      phase_field = d
      degradation_function = g
      decomposition = SPECTRAL
      use_threshold = false
      use_history_max = false
      output_properties = 'psie_active'
      tensile_strength = sigma0
      outputs = exodus
    []
    [stress]
      type = ComputeCreepPlasticityDeformationStress
      elasticity_model = Elasticity
      creep_model = creep
    []

[]



# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/(pellet_outer_radius^2-pellet_inner_radius^2)}' #新加的！！！！！！！！！！！！！！！！！！！！！！
[Functions]
  [gap_conductance_in]
    type = PiecewiseLinear
    x = '0 ${endTime}'
    y = '5000 6000'
    scale_factor = 1         # 保持原有的转换因子
  []
  [gap_conductance_out]
    type = PiecewiseLinear
    x = '0 1000000 ${endTime}'
    y = '8063.5 6016.1 13642'
    scale_factor = 1         # 保持原有的转换因子
  []
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
  type = PiecewiseLinear
  # data_file = '../../../.././power_history2.csv'    # 创建一个包含上述数据的CSV文件，数据为<s,w/m>
  # format = columns                 # 指定数据格式为列式
  # scale_factor = ${power_factor}         # 保持原有的转换因子
  # 论文中只给了线密度，需要化为体积密度
    x = '0.0 2400.0 98400.0 ${endTime__100000} ${endTime__50000} ${endTime}'
    y = '0.0 ${LinearPower0_2} ${LinearPower} ${LinearPower} 0 0'
    scale_factor = ${power_factor}
  []
  [gap_pressure_inner] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
    x = '0   ${endTime}'
    y = '${Pressure1} ${Pressure2}'
    scale_factor = 1
  []
  [gap_pressure_outer] #新加的！！！！！！！！！！！！！！！！！！！！！！
  #间隙压力随时间的变化
  type = PiecewiseLinear
    x = '0 ${endTime}'
    y = '${Pressure1} ${Pressure2}'
  scale_factor = 1
[]
[T_infinity_in]
  #间隙压力随时间的变化
  type = PiecewiseLinear
  x = '0 98400 ${endTime__100000} ${endTime__50000} ${endTime}'
  y = '293.15 ${initial_T_in} ${initial_T_in} 293.15 293.15'
  scale_factor = 1
[]
[T_infinity_out]
  #间隙压力随时间的变化
  type = PiecewiseLinear
  x = '0 98400 ${endTime__100000} ${endTime__50000} ${endTime}'
  y = '293.15 ${initial_T_out} ${initial_T_out} 293.15 293.15'
  scale_factor = 1
[]
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 16000, 2000,
                   if(t < 105000, 750,
                   if(t < ${endTime__100000},50000,
                   if(t < (${endTime__50000}+20000), 750,10000))))'
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  # solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type -ksp_type' 
  # petsc_options_value = 'lu gmres' 
#=================================================
  # solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -ksp_type' 
  # petsc_options_value = 'lu gmres' 
#=================================================
  # solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  # petsc_options_value = 'lu superlu_dist gmres'
#=================================================
  # solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  # petsc_options_value = 'lu superlu_dist gmres'
  #=================================================
  solve_type = 'NEWTON'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  petsc_options_value = '201                hypre    boomeramg'  
#=================================================
  # solve_type = 'PJFNK'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
#=================================================
  line_search = 'bt'
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true  # 每个时间步都重新计算缩放
  nl_max_its = 150
  nl_rel_tol = 5e-6 # 非线性求解的相对容差
  nl_abs_tol = 5e-6 # 非线性求解的绝对容差
  l_tol = 5e-8  # 线性求解的容差
  l_abs_tol = 5e-9 # 线性求解的绝对容差
  l_max_its = 500 # 线性求解的最大迭代次数
  # abort_on_solve_fail = true
  dtmin = 500
  dtmax = 50000
  end_time = ${endTime} # 总时间24h

  fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
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
      type = PhasePiledFractureHSMarker
      von_mises_variable = hoop_stress
      sigma0 = sigma0
      x1 = 0.0001 #d变量小于x1时，标记为粗网格
      x2 = 0.05 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.1 #d变量大于xmax时，一定是细网格
      y1 = 0.6 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.7 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.02
      stress_change_threshold = 1e6
    []
  []
[]

[Postprocessors]
  [pellet_area]
    type = VolumePostprocessor
    block = pellet
    use_displaced_mesh = true
  []
  [burnup_avg]
    type = ElementAverageValue
    variable = burnup
    block = pellet
    execute_on = 'TIMESTEP_END'
  []
  [pellet_area0_full]
    type = ParsedPostprocessor
    constant_names = 'R r'
    constant_expressions = '${fparse pellet_outer_radius} ${fparse pellet_inner_radius}'
    expression = '3.141592653589793*(R*R - r*r)'
  []
  [pellet_area_change_full]
    type = ParsedPostprocessor
    pp_names = 'pellet_area pellet_area0_full'
    expression = '4*pellet_area - pellet_area0_full'
  []
  [pellet_area_ratio_full]
    type = ParsedPostprocessor
    pp_names = 'pellet_area pellet_area0_full'
    expression = '4*pellet_area/pellet_area0_full'
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
  file_base = '2D-NoFracture/2D'
  [csv_out]
    type = CSV
    execute_on = 'TIMESTEP_END'
    show = 'burnup_avg pellet_area_change_full pellet_area_ratio_full'
    file_base = 'Volume/data'
  []
[]