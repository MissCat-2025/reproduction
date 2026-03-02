# === 参数研究案例 ===
# end_time = 3.70e+5
# mesh_size: 5.00e-5
# Gc: 3
# length_scale_paramete: 2.00e-5
# 生成时间: 2025-06-24 09:31:06
#conda activate moose &&mpirun -n 15 /home/yp/projects/reproduction/reproduction-opt -i main_Pr3_00e+6.i --recover main_Pr3_00e+6_my_checkpoint_cp/0110

pellet_E=201.3e9
pellet_density=10431.0#10431.0*0.85#kg⋅m-3
pellet_nu = 0.345
pellet_thermal_expansion_coef=1e-5#K-1
pellet_critical_energy = 3.25 #断裂能
grain_size = 10.0
pellet_critical_fracture_strength=6.0e7#Pa
length_scale_paramete = 4e-5
mesh_size = '${fparse 3.5e-5}' #mm,最大网格尺寸
# length_scale_paramete = '${fparse 1.2*mesh_size}'
fission_rate = 1.2e19
power_factor_mod = 1.0
PressureFactor = 3.00e+6 # 间隙压力2.5MPa
# 各种参数都取自[1]Multiphysics phase-field modeling of quasi-static cracking in urania ceramic nuclear fuel
#几何与网格参数
pellet_outer_radius = 4.1e-3#直径变半径，并且单位变mm
#将下列参数转化为整数
n_elems_azimuthal = '${fparse 2*ceil(3.1415*2*(pellet_outer_radius/(4*mesh_size)/2))}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int(pellet_outer_radius/(4*mesh_size))}'          # 芯块径向网格数（直接取整）
#相场断裂参数：
degradation_factor = 1e-6

#功率突降的时间参数
LinearPowerAll=35   #总线功率kw/m
#注意提第一个时间节点一定要大于125000s，因为125000s功率才刚刚稳定，温度可能还没稳定

#第一个节点：功率突降的开始时的节点
power1 = 0.5  # 功率从1降低到power1
TimeDown1 = 60  #TimeDown1秒后功率从1降低到power1
Time1 = '${fparse 167000+TimeDown1}'  #s 功率从0%增加到x%的时间,假设正常工况下的温升，即2400s升到20%
LinearPower1 =  '${fparse LinearPowerAll*power1}'  
#第二个节点：停堆时的节点
TimeDown2 = 30  #TimeDown2秒后功率从power1降低到0
Time2 = '${fparse Time1+TimeDown2}'  #s 功率从0%增加到x%的时间,假设正常工况下的温升，即2400s升到20%
#第二个节点：观察结束时的节点
WatchTime = 40   #s 停堆后的观察时间(一般以裂纹稳定后为准)
endTime = '${fparse Time2+WatchTime}'  #s 观察结束时的时间

# timePlus1='${fparse 1}'
# dt = '${fparse 1}' #s 时间步长
# # initial_T_in: 570.7
dtmin = 1e-3
# length_scale_paramete = '${fparse 1.2*mesh_size}'

# 各种参数都取自[1]Multiphysics phase-field modeling of quasi-static cracking in urania ceramic nuclear fuel
#几何与网格参数#相场断裂参数：
# 双曲
m = 4
a2 = 0.5396842
a3 = 0
ksi = 2

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


[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'sub_Pr3_00e+6.i'
    cli_args = 'Gc=${pellet_critical_energy};l=${length_scale_paramete};pellet_outer_radius=${pellet_outer_radius};mesh_size=${mesh_size};m=${m};a2=${a2};a3=${a3};ksi=${ksi};degradation_factor=${degradation_factor};Time1=${Time1};Time2=${Time2};endTime=${endTime};dtmin=${dtmin}'
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
    variable = 'psie_active a1'
    source_variable = 'psie_active a1'
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
    [radial_stress]
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
  [effective_creep]
    order = CONSTANT
    family = MONOMIAL
  []
  [creep_strain_hoop]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [creep_strain_radial]
    order = CONSTANT
    family = MONOMIAL
  [../]
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
    rank_two_tensor = stress
    variable = radial_stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = 'TIMESTEP_END'
  []
  [copy_sigma0]
      type = ADMaterialRealAux
      variable = sigma0_field
      property = sigma0
      execute_on = 'initial'
      block = pellet
  []
  [effective_creep]
      type = ADMaterialRealAux
      variable = effective_creep
      property = creep_effective_strain_value
      execute_on = 'TIMESTEP_END'
      block = pellet
  []
      [./creep_strain_hoop]
    type = ADRankTwoScalarAux
    variable = creep_strain_hoop
    rank_two_tensor = creep_strain
    scalar_type = HoopStress
    point1 = '0 0 0.01'        # 圆心坐标
    point2 = '0 0 -0.01'        # 定义旋转轴方向（z轴）
    execute_on = 'TIMESTEP_END'
  [../]
  [./creep_strain_radial]
    type = ADRankTwoScalarAux
    rank_two_tensor = creep_strain
    variable = creep_strain_radial
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
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
    initial_from_file_var = disp_x
  []
  [disp_y]
    initial_from_file_var = disp_y
  []
  [T]
    initial_from_file_var = T
  []
  [x]
    initial_from_file_var = x
  []
  [strain_zz]
    initial_from_file_var = strain_zz
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
    boundary = 'yplane'
    value = 0
  []
  [x_zero_on_x_plane]
    type = DirichletBC
    variable = disp_x
    boundary = 'xplane'
    value = 0
  []
  [coolant_bc]#对流边界条件
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'pellet_outer'
    T_infinity = 500
    coefficient = 3400#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
    #芯块包壳间隙压力边界条件
    [gap_pressure_fuel_x]
      type = Pressure
      variable = disp_x
      boundary = 'pellet_outer'
      factor = ${PressureFactor} # 间隙压力2.5MPa
    []
    [gap_pressure_fuel_y]
      type = Pressure
      variable = disp_y
      boundary = 'pellet_outer'
      factor = ${PressureFactor}
    []
    [x_on_OUTER]
    type = DirichletBC
    variable = x
    boundary = 'pellet_outer'
    value = 0.01
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
    block = pellet
    []
    [pellet_specific_heat]
      type = ADParsedMaterial
      property_name = specific_heat #Fink model
      coupled_variables = 'T'  # 需要在AuxVariables中定义Y变量
      expression = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
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

    # # # # 蠕变相关
    [creep]
      type = UO2PowerLawCreepStressUpdate
      temperature = T
      oxygen_ratio = x
      fission_rate = ${fparse fission_rate}
      theoretical_density = 95.0
      grain_size = '${fparse grain_size}'
      degradation_function = g
      use_new_three_shear_modulus = false

    []
    [creep_effective_strain]
      type = ADRankTwoInvariant
      rank_two_tensor = creep_strain
      property_name = creep_effective_strain_value
      invariant = EffectiveStrain
      block = pellet
    []
    #力学属性
    #力学属性
    [pellet_strain]
      type = ADComputePlaneSmallStrain
      eigenstrain_names = 'thermal_eigenstrain'
    []
        # 肿胀应变函数
    [total_power]
      type = ADRimEffertPowerBurnupRod
      power_history = 'power_history'  # 声明使用的函数
      pellet_outer_radius = ${pellet_outer_radius}  # 直接使用函数符号进行计算
      output_properties = 'burnup'
      outputs = exodus
    []
    #断裂
    [crack_geometric]
      type = CrackGeometricFunction
      property_name = alpha
      expression = 'ksi*d+(1-ksi)*d*d'
      parameter_names = 'ksi'
      parameter_values = '${ksi}'
      phase_field = d
    []
    [a1]
      type = ADParsedMaterial
      property_name = a1
      material_property_names = 'Gc E l sigma0'
      expression = '4*E*Gc/sigma0/sigma0/l/3.14159'
      output_properties = 'a1'
      outputs = exodus
    []
    #   [a1]
    #   type = ADDerivativeParsedMaterial
    #   property_name = a1
    #   material_property_names = 'Gc E l sigma0'
    #   expression = '1.5*E*Gc/sigma0/sigma0/l'
    #   output_properties = 'a1'
    #   outputs = exodus
    #   block = pellet
    # []
    [degradation]
      type = RationalDegradationFunction
      property_name = g
      expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))*(1-eta)+eta
      phase_field = d
      material_property_names = 'a1'
      parameter_names = 'p a2 a3 eta'
      parameter_values = '${m} ${a2} ${a3} ${degradation_factor}'
      block = pellet
    []
    [elasticity]
      type = SmallDeformationIsotropicElasticityThresholdModMax
      youngs_modulus = E
      poissons_ratio = nu
      phase_field = d
      degradation_function = g
      decomposition = SPECTRAL
      use_history_max = true
      # 阈值参数
      use_threshold = true
      tensile_strength = sigma0
      
      output_properties = 'psie_active'
      outputs = exodus
      block = pellet
    []
    [stress]
      type = ComputeSmallDeformationStressAll
      inelastic_models = 'creep'
      elasticity_model = elasticity
    []
[]
# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/pellet_outer_radius/pellet_outer_radius/power_factor_mod}' #新加的！！！！！！！！！！！！！！！！！！！！！！
[Functions]
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
    type = PiecewiseLinear
    x = '0.0 167000 ${Time1} ${Time2} ${endTime}'
    y = '0.0 35.0 ${LinearPower1} 0 0'
    scale_factor = ${power_factor}
  []
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 167000, 1000,
                   if(t < ${Time1}, 1,
                   if(t < ${Time2}, 1,
                   if(t < ${endTime}, 1,1))))'
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  solve_type = 'NEWTON' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm      100                lu           NONZERO'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  # petsc_options_value = 'lu superlu_dist gmres'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = 'lu       superlu_dist                  51'
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true  # 每个时间步都重新计算缩放
  nl_max_its = 100
  nl_rel_tol = 1e-5 # 非线性求解的相对容差
  nl_abs_tol = 1e-6 # 非线性求解的绝对容差
  l_tol = 1e-7  # 线性求解的容差
  l_abs_tol = 1e-8 # 线性求解的绝对容差
  l_max_its = 150 # 线性求解的最大迭代次数
  accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
  dtmin = 100
  end_time =  3.50e+5# 总时间24h3.7e5

  # fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]

[Outputs]
  [my_checkpoint]
    type = Checkpoint
    time_step_interval = 10    # 每5个时间步保存
    num_files = 500            # 保留最近4个检查点
    wall_time_interval = 600 # 每10分钟保存一次（秒）
  []
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = false
[]
