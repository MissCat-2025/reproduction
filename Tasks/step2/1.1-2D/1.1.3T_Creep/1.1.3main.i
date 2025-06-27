# === 参数研究案例 ===

# mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i 1.1.3main.i
# mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i 1.1.2main.i --mesh-only KAERI_HANARO_UpperRod1.e
pellet_E=201.3e9
pellet_density=10431.0#10431.0*0.85#kg⋅m-3
pellet_nu = 0.345
pellet_thermal_expansion_coef=1e-5#K-1

length_scale_paramete=6e-5
NNN=1.5
grid_sizes = '${fparse length_scale_paramete/NNN}'

#几何与网格参数
pellet_outer_radius = 4.1e-3#直径变半径，并且单位变mm
#将下列参数转化为整数
n_elems_azimuthal = '${fparse 2*ceil(3.1415*2*(pellet_outer_radius/(4*grid_sizes)/2))}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int(pellet_outer_radius/(4*grid_sizes))}'          # 芯块径向网格数（直接取整）

[Mesh]
  [pellet_clad_gap]
    type = ConcentricCircleMeshGenerator
    num_sectors = '${n_elems_azimuthal}'  # 周向网格数
    radii = '${pellet_outer_radius}'
    rings = '${n_elems_radial_pellet}'
    has_outer_square = false
    preserve_volumes = true
    portion = top_right # 生成四分之一计算域
    smoothing_max_it=666 # 平滑迭代次数
  []
  [rename]
    type = RenameBoundaryGenerator
    input = pellet_clad_gap
    old_boundary = 'bottom left outer'
    new_boundary = 'yplane xplane pellet_outer' # 将边界命名为yplane xplane clad_outer

  []
  [rename2]
    type = RenameBlockGenerator
    input = rename
    old_block  = '1'
    new_block  = 'pellet' # 将block1和block3分别命名为pellet和clad
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
    [effective_creep]
      order = CONSTANT
      family = MONOMIAL
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
    [effective_creep]
      type = ADMaterialRealAux
      variable = effective_creep
      property = creep_effective_strain_value
      execute_on = 'TIMESTEP_END'
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
    [strain_zz]
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
    boundary = 'yplane'
    value = 0
  []
  [x_zero_on_x_plane]
    type = DirichletBC
    variable = disp_x
    boundary = 'xplane'
    value = 0
  []
  [T_0BC]
    type = NeumannBC
    variable = T
    boundary = 'xplane yplane'
    value = 0
  []
  [coolant_bc]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'pellet_outer'
    T_infinity = 500
    coefficient = gap_conductance#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []

    #芯块包壳间隙压力边界条件
    [gap_pressure_fuel_x]
      type = Pressure
      variable = disp_x
      boundary = 'pellet_outer'
      factor = 1e6 # 间隙压力2.5MPa
      function = gap_pressure #新加的！！！！！！！！！！！！！！！！！！！！！！
      use_displaced_mesh = true
    []
    [gap_pressure_fuel_y]
      type = Pressure
      variable = disp_y
      boundary = 'pellet_outer'
      factor = 1e6
      function = gap_pressure #新加的！！！！！！！！！！！！！！！！！！！！！！
      use_displaced_mesh = true
    []

[]
[Materials]
    #定义芯块热导率、密度、比热等材料属性
    [pellet_properties]
      type = ADGenericConstantMaterial
      prop_names = 'density nu E'
      prop_values = '${pellet_density} ${pellet_nu} ${pellet_E}'
      block = pellet
    []

    [pellet_thermal_conductivity] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = ADParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '(100/(7.5408 + 17.692*T/1000 + 3.6142*(T/1000)^2) + 6400/((T/1000)^2.5)*exp(-16.35/(T/1000)))'
    block = pellet
    []
    [pellet_specific_heat]
      type = ADParsedMaterial
      property_name = specific_heat #Fink model
      coupled_variables = 'T'  # 需要在AuxVariables中定义Y变量
      expression = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (0+2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
      block = pellet
    []
        #化学相关
        [D_fickian]
          type = ADParsedMaterial
          property_name = D_fickian
          coupled_variables = 'x T'
          expression = 'pow(10, -9.386 - 4260/(T) + 0.0012*T*x + 0.00075*T*log10(1+2/(x)))'
          block = pellet
        []
        [D_soret]
          type = ADDerivativeParsedMaterial
          property_name = D_soret
          coupled_variables = 'x T'
          material_property_names = 'D_fickian(x,T)'
          expression = 'D_fickian * x * (-1380.8 - 134435.5*exp(-x/0.0261)) / ((2.0 + x)/(2.0 * (1.0 - 3.0*x) * (1.0 - 2.0*x)) * 8.314 * T * T)'
          block = pellet
        []
    [elasticity_tensor]
      type = ADComputeVariableIsotropicElasticityTensor
      youngs_modulus = E
      poissons_ratio = nu
    []
    [pellet_thermal_eigenstrain]
      type = ADComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = 500
      thermal_expansion_coeff = ${pellet_thermal_expansion_coef}
      temperature = T
      block = pellet
    []
    [g]
      type = ADParsedMaterial
      property_name = g
      coupled_variables = 'x'
      expression = '1-x'
      block = pellet
    []
    [creep]
      type = UO2DegradablePowerLawCreepStressUpdate
      temperature = T
      oxygen_ratio = x
      fission_rate = 1.2e19
      theoretical_density = 95.0
      grain_size = 10
      
      # 相场断裂相关参数
      degradation_function = g
      use_stress_degradation = false
      # 输出设置
    []
    [creep_effective_strain]
      type = ADRankTwoInvariant
      rank_two_tensor = creep_strain
      property_name = creep_effective_strain_value
      invariant = EffectiveStrain
      block = pellet
    []
    [strain]
      type = ADComputePlaneSmallStrain
      eigenstrain_names = 'thermal_eigenstrain creep_strain'
    []
    [stress]
      type = ADComputeLinearElasticStress
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
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 20000, 5000,
                   if(t < 100000, 2500,
                   if(t < 125000, 1000,
                   if(t < 175000, 1000,2500))))'
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  petsc_options_value = 'lu superlu_dist gmres'
  # petsc_options_iname = '-pc_type -ksp_type'
  # petsc_options_value = 'lu gmres'
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true  # 每个时间步都重新计算缩放
  reuse_preconditioner = true
  reuse_preconditioner_max_linear_its = 20
  line_search = BT
  nl_max_its = 100
  nl_rel_tol = 1e-8 # 非线性求解的相对容差
  nl_abs_tol = 5e-9 # 非线性求解的绝对容差
  l_tol = 1e-8  # 线性求解的容差
  l_abs_tol = 5e-9 # 线性求解的绝对容差
  l_max_its = 500 # 线性求解的最大迭代次数
  accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
  dtmin = 100
  end_time = 3.7e5 # 总时间24h

  fixed_point_rel_tol =1e-5 # 固定点迭代的相对容差
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]

[Outputs]
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = false
[]
