# === 参数研究案例 ===

# conda activate moose && dos2unix 3.1.1main.i && dos2unix 3.1.1_Sub.i && mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i 3.1.1main.i
# mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i 1.1.2main.i --mesh-only KAERI_HANARO_UpperRod1.e
pellet_E=201.3e9
pellet_density=10431.0#10431.0*0.85#kg⋅m-3
pellet_nu = 0.345
pellet_thermal_expansion_coef=1e-5#K-1
Gc = 2.896 #断裂能
grain_size = 10
pellet_critical_fracture_strength=6.0e7#Pa
fission_rate = 1.2e19

length_scale_paramete=3.0e-5
grid_sizes = 8e-5

#几何与网格参数
pellet_outer_radius = 4.1e-3#直径变半径，并且单位变mm
#将下列参数转化为整数
n_elems_azimuthal = '${fparse 2*ceil(3.1415*2*(pellet_outer_radius/(4*grid_sizes)/2))}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int(pellet_outer_radius/(4*grid_sizes))}'          # 芯块径向网格数（直接取整）
# pellet_critical_energy=${fparse (1+grid_sizes/2/length_scale_paramete)*Gc} #J⋅m-2
pellet_critical_energy = ${fparse Gc}

[Mesh]
  [pellet_clad_gap]
    type = ConcentricCircleMeshGenerator
    num_sectors = '${n_elems_azimuthal}'  # 周向网格数
    radii = '${pellet_outer_radius}'
    rings = '${n_elems_radial_pellet}'
    has_outer_square = false
    preserve_volumes = true
    portion = full # 生成四分之一计算域
    smoothing_max_it=666 # 平滑迭代次数
  []
  [rename]
    type = RenameBoundaryGenerator
    input = pellet_clad_gap
    old_boundary = 'outer'
    new_boundary = 'pellet_outer' # 将边界命名为yplane xplane clad_outer
  []
  [rename2]
    type = RenameBlockGenerator
    input = rename
    old_block  = '1'
    new_block  = 'pellet' # 将block1和block3分别命名为pellet和clad
  []
  [center_point]
    type = ExtraNodesetGenerator
    input = rename2
    coord = '0 0 0'
    new_boundary  = 'center_point'
  []
[]


[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = '3.1.1_Sub.i'
    cli_args = 'Gc=${pellet_critical_energy};l=${length_scale_paramete}'
    execute_on = 'TIMESTEP_END'
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
    variable = 'psie_active a1 vonMises'
    source_variable = 'psie_active a1 vonMises'
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
    [radial_stress]    # 径向应力，用于检查压力在径向的传递
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
  [vonMises]
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
    [radial_stress]
      type = ADRankTwoScalarAux
      variable = radial_stress
      rank_two_tensor = stress
      scalar_type = RadialStress
      point1 = '0 0 0'
      point2 = '0 0 -1'
      execute_on = 'TIMESTEP_END'
    []
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
    [T]
      initial_condition = 500
    []
    [x]
      initial_condition =0.01
      block = pellet
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
    boundary = 'center_point'
    value = 0
  []
  [x_zero_on_x_plane]
    type = DirichletBC
    variable = disp_x
    boundary = 'center_point'
    value = 0
  []
  [coolant_bc]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'pellet_outer'
    T_infinity = 500
    coefficient = gap_conductance#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  #氧浓度边界条件
  # [x_dirichlet_bc]
  #   type = DirichletBC
  #   variable = x
  #   boundary = pellet_outer
  #   value = 0.01
  # []

    #芯块包壳间隙压力边界条件
    [gap_pressure_fuel_x]
      type = Pressure
      variable = disp_x
      boundary = 'pellet_outer'
      factor = 1e6 # 间隙压力2.5MPa
      function = gap_pressure #新加的！！！！！！！！！！！！！！！！！！！！！！
      # use_displaced_mesh = true
    []
    [gap_pressure_fuel_y]
      type = Pressure
      variable = disp_y
      boundary = 'pellet_outer'
      factor = 1e6
      function = gap_pressure #新加的！！！！！！！！！！！！！！！！！！！！！！
      # use_displaced_mesh = true
    []

[]
[Materials]
    #定义芯块热导率、密度、比热等材料属性
    [pellet_properties]
      type = ADGenericConstantMaterial
      prop_names = 'l Gc density nu E'
      prop_values = '${length_scale_paramete} ${pellet_critical_energy} ${pellet_density} ${pellet_nu} ${pellet_E}'
      block = pellet
    []
    # 为临界断裂强度生成威布尔分布
    [sigma0_mat]
      type = ADParsedMaterial
      property_name = sigma0
      coupled_variables = 'sigma0_field'
      expression = 'sigma0_field'  # 直接使用辅助变量的值
      block = pellet
    []
    [pellet_thermal_conductivity] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = ADParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '(100/(7.5408 + 17.692*T/1000 + 3.6142*(T/1000)^2) + 6400/((T/1000)^2.5)*exp(-16.35/(T/1000)))'
    block = pellet    []
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
      stress_free_temperature = 500
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
    [creep]
      type = UO2CreepRate
      phase_field = d
      degradation_function = g
      temperature = T
      oxygen_ratio = x
      fission_rate = ${fparse fission_rate}
      theoretical_density = 95.0
      grain_size = '${fparse grain_size}'
      # 相场断裂相关参数
      use_transition_stress = false
      # use_transient_creep = true
      use_three_shear_modulus = true
      # full_three_shear_modulus_strategy = true

      relative_tolerance = 1e-8 #蠕变的相对残差
      absolute_tolerance = 1e-10 #蠕变的绝对残差

      output_properties = 'effective_creep_strain psic_active'
      outputs = exodus 
    []
    #力学属性
    #力学属性
    [strain]
      type = ADComputePlaneSmallStrain
      eigenstrain_names = thermal_eigenstrain
    []

    [a1]
      type = ADDerivativeParsedMaterial
      property_name = a1
      material_property_names = 'Gc E l sigma0'
      expression = '1.5*E*Gc/sigma0/sigma0/l'
      output_properties = 'a1'
      outputs = exodus
      block = pellet
    []
    # [a1]
    #   type = ADDerivativeParsedMaterial
    #   property_name = a1
    #   material_property_names = 'Gc E l sigma0'
    #   expression = '4*E*Gc/sigma0/sigma0/3.14159/l'
    #   output_properties = 'a1'
    #   outputs = exodus
    #   block = pellet
    # []
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
power_factor = '${fparse 1000*1/3.1415926/pellet_outer_radius/pellet_outer_radius}' #新加的！！！！！！！！！！！！！！！！！！！！！！
[Functions]
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = PiecewiseLinear
    x = '0.0 100000 125000 175000 300000 370000'
    y = '0.0 18.0 35.0 35.0 0.0 0'
    scale_factor = ${power_factor}
  []
  [gap_conductance]
    type = PiecewiseLinear
    x = '0 370000'
    y = '3400 3400'
    scale_factor = 1         # 保持原有的转换因子
  []
  [gap_pressure] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
    x = '0          125000   175000   350000'
    y = '2.5  6  10  14'
    scale_factor = 1
  []
  # [dt_limit_func]
  #   type = ParsedFunction
  #   expression = 'if(t < 25000, 5000,
  #                   if(t < 100000, 1000,
  #                   if(t < 125000, 1000,
  #                   if(t < 175000, 5000,1000))))'
  # []
[]

[Executioner]
  type = Transient # 瞬态求解器
  solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
  # solve_type = 'NEWTON'
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
  dtmax = 1000
  end_time = 3.7e5 # 总时间24h

  fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1000
    growth_factor = 2
    cutback_factor = 0.8
    optimal_iterations = 100
    iteration_window = 25
  [../]
[]
[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = 2
  [Markers]
    [marker]
      type = PhasePiledFractureMarker
      von_mises_variable = vonMises
      x1 = 0.001 #d变量小于x1时，标记为粗网格
      x2 = 0.2 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.95 #d变量大于xmax时，标记为粗网格
      y1 = 4.0e7 #vonMises应力小于y1时，标记为粗网格
      y2 = 4.0e7 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      enable_time_check = true
      time_steps = 3
      d_change_threshold = 0.01
      # invert = true
      # third_state = coarsen
    []
  []
[]
[Outputs]
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = false
  file_base = '2.1-2D-NewCreep1/1'
[]
