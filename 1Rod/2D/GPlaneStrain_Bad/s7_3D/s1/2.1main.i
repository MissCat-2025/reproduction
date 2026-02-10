#conda activate moose && mpirun -n 11 /home/yp/projects/reproduction/reproduction-opt -i 2.1main.i --recover 2.1-2D-New2026/1_cp/0100
# conda activate moose &&mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i 2.1main.i --recover
# conda activate moose &&mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i 2.1main.i
LinearPower = 35
LinearPower0_2 = '${fparse LinearPower*0.2}'
endTime = 3e5
dtmin = 1
dt = 2000
dtMax = 5000
endTime__50000 = '${fparse endTime-5000}'
endTime__100000 = '${fparse endTime-10000}'
initial_T = 550
pellet_E=201.3e9
pellet_nu = 0.345   #RELAP5
pellet_thermal_expansion_coef=1e-5#K-1
pellet_critical_fracture_strength=6.0e7#Pa
density_percent = 0.95
# Gc = 6#断裂能
# fission_rate = 1.2e19
# grain_size =10
pellet_critical_energy = 8# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
pellet_density='${fparse density_percent*10980}'#10431.0*0.85#kg⋅m-3理论密度为10.980
#几何与网格参数
# density_percent100 = '${fparse density_percent*100}'
length_scale_paramete = 6e-5

w = 0 #裂纹尖端时，l是mesh_size的2**w倍
mesh_size = '${fparse 2*length_scale_paramete}' #网格尺寸即可
#将下列参数转化为整数
pellet_outer_radius = 4.2e-3#直径变半径，并且单位变mm
length = 0.05e-3 # 芯块长度17.78mm
n_elems_axial = 1 # 轴向网格数
n_elems_azimuthal = '${fparse 2*ceil((3.1415*pellet_outer_radius/mesh_size)/2^w)}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int((pellet_outer_radius/mesh_size)/2^w)}'          # 芯块径向网格数（直接取整）

# creep_relative_tolerance = 1e-7 #蠕变的相对残差
# creep_absolute_tolerance = '${fparse creep_relative_tolerance}' #蠕变的绝对残差
# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/pellet_outer_radius/pellet_outer_radius}' #新加的！！！！！！！！！！！！！！！！！！！！！！

#相场断裂参数：
# m = 2
# a2 = 2
# a3 = 0
# ksi = 2

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
    smoothing_max_it=300 # 平滑迭代次数
  []
  [rename]
    type = RenameBoundaryGenerator
    input = pellet_clad_gap
    old_boundary = 'bottom left outer'
    new_boundary = 'yplane xplane pellet_outer' # 将边界命名为yplane xplane clad_outer
  []
[extrude]
  type = AdvancedExtruderGenerator
  input = rename                   # 修改输入为切割后的网格
  heights = '${length}'
  num_layers = '${n_elems_axial}'
  direction = '0 0 1'
  bottom_boundary = '100'
  top_boundary = '101'
  subdomain_swaps = '1 1'
[]
[rename_extrude]
  type = RenameBoundaryGenerator
  input = rename_extrude
  old_boundary = '100 101'
  new_boundary = 'bottom top' # 最终边界命名
[]
  [rename2]
    type = RenameBlockGenerator
    input = rename_extrude
    old_block  = '1'
    new_block  = 'pellet' # 将block1和block3分别命名为pellet和clad
  []
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp    input_files = '2.1_Sub.i'
    cli_args = 'l=${length_scale_paramete};mesh_size=${mesh_size};Gc=${pellet_critical_energy};sigma0=${pellet_critical_fracture_strength};m=${m};w=${w};a2=${a2};a3=${a3};ksi=${ksi};endTime=${endTime};dtmin=${dtmin};dt=${dt};pellet_outer_radius=${pellet_outer_radius};dtMax=${dtMax}'
    execute_on = 'TIMESTEP_END'
    sub_cycling = false
    catch_up = false
    max_failures = 0
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
    variable = 'psie_active a1 stress_I'
    source_variable = 'psie_active a1 stress_I'
  []
[]

[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
[]
[AuxVariables]
  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [stress_I]
    order = CONSTANT
    family = MONOMIAL
  []
  [radial_stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
  []
  # [creep_strain_hoop]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [creep_strain_radial]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [d]
    block = pellet
    initial_condition = 0.0
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
  # [total_power_aux]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   block = pellet
  # []
[]

[AuxKernels]
  [./hoop_stress]
    type = RankTwoScalarAux
    variable = hoop_stress
    rank_two_tensor = stress
    scalar_type = HoopStress
    point1 = '0 0 0.01'        # 圆心坐标
    point2 = '0 0 -0.01'        # 定义旋转轴方向（z轴）
    execute_on = 'TIMESTEP_END'
  [../]
  [./stress_I]
    type = RankTwoScalarAux
    scalar_type = MaxPrincipal
    rank_two_tensor = stress
    variable = stress_I
    selected_qp = 0
    execute_on = 'TIMESTEP_END'
    block = pellet
  [../]
  [radial_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = radial_stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = 'TIMESTEP_END'
  []
    [vonMisesStress]
      type = RankTwoScalarAux
      variable = vonMises
      rank_two_tensor = stress
      execute_on = 'TIMESTEP_END'
      scalar_type = VonMisesStress
      # 不需要 index_i 和 index_j，因为我们使用 VonMisesStress 标量类型
    []
  #   [./creep_strain_hoop]
  #   type = RankTwoScalarAux
  #   variable = creep_strain_hoop
  #   rank_two_tensor = creep_strain
  #   scalar_type = HoopStress
  #   point1 = '0 0 0.01'        # 圆心坐标
  #   point2 = '0 0 -0.01'        # 定义旋转轴方向（z轴）
  #   execute_on = 'TIMESTEP_END'
  # [../]
  # [./creep_strain_radial]
  #   type = RankTwoScalarAux
  #   rank_two_tensor = creep_strain
  #   variable = creep_strain_radial
  #   scalar_type = RadialStress
  #   point1 = '0 0 0'
  #   point2 = '0 0 1'
  # [../]
    [copy_sigma0]
    type = MaterialRealAux
    variable = sigma0_field
    property = sigma0
    execute_on = 'initial'
    block = pellet
  []
  # [total_power_aux]
  #   type = ADMaterialRealAux
  #   variable = total_power_aux
  #   property = total_power
  #   execute_on = 'INITIAL TIMESTEP_BEGIN'
  #   block = pellet
  # []
[]

[Variables]
    [disp_x]
    []
    [disp_y]
    []
    [disp_z]
    []
    [T]
      initial_condition = ${initial_T}
    []
    [x]
      initial_condition = 0.01
    []
    [scalar_strain_zz]
      family = SCALAR
      order = FIRST
    []
[]
[Kernels]
    #热传导方程
    [heat_conduction]
      type = HeatConduction
      variable = T
    []
    [hcond_time]
      type = HeatConductionTimeDerivative
      variable = T
    []
    [Fheat_source]
      type = HeatSource
      variable = T
      function = power_history
    []
    #化学平衡方程
    [time_derivative]
      type = TimeDerivative
      variable = x
      block = pellet
    []
    [complex_diffusion]
      type = ComplexDiffusionKernel
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
  [x_on_OUTER]
    type = DirichletBC
    variable = x
    boundary = 'pellet_outer'
    value = 0.01
  []
  [coolant_bc_out]#对流边界条件
  type = ConvectiveFluxFunction
  variable = T
  boundary = 'pellet_outer'
  T_infinity = ${initial_T}
  coefficient = 3400#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  #芯块包壳间隙压力边界条件
  [gap_pressure_fuel_x]
    type = Pressure
    variable = disp_x
    boundary = 'pellet_outer'
    factor = 2e6 # 间隙压力2.5MPa
  []
  [gap_pressure_fuel_y]
    type = Pressure
    variable = disp_y
    boundary = 'pellet_outer'
    factor = 2e6
    # function = 2 #新加的！！！！！！！！！！！！！！！！！！！！！！
    # use_displaced_mesh = true
  []
[]



[Materials]
    [pellet_properties]
      type = GenericConstantMaterial
      prop_names = 'density E nu Gc l'
      prop_values = '${pellet_density} ${pellet_E} ${pellet_nu} ${pellet_critical_energy} ${length_scale_paramete}'
      block = pellet
    []
    # 为临界断裂强度生成威布尔分布
    [sigma0_mat]
      type = ParsedMaterial
      property_name = sigma0
      coupled_variables = 'sigma0_field'
      expression = 'sigma0_field'  # 直接使用辅助变量的值
      block = pellet
    []
    [pellet_thermal_conductivity] #新加的！！！！！！！！！！！！！！！！！！！！！！
      type = ParsedMaterial
      property_name = thermal_conductivity #参考某论文来的，不是Fink-Lukuta model（非常复杂）
      coupled_variables = 'T'
      expression = '(1)*(100/(7.5408 + 17.692*T/1000 + 3.6142*(T/1000)^2) + 6400/((T/1000)^2.5)*exp(-16.35/(T/1000)))'
      block = pellet
    []
    [pellet_specific_heat]
      type = ParsedMaterial
      property_name = specific_heat #Fink model
      coupled_variables = 'T x'  # 需要在AuxVariables中定义Y变量
      expression = '(296.7 * 535.285^2 * exp(535.285/T))/(T^2 * (exp(535.285/T) - 1)^2) + 2.43e-2 * T + (x+2) * 8.745e7 * 1.577e5 * exp(-1.577e5/(8.314*T))/(2 * 8.314 * T^2)'
      block = pellet
    []
    # 肿胀应变函数
    # [total_power]
    #   type = ADRimEffertPowerBurnupRod
    #   power_history = 'power_history'  # 声明使用的函数
    #   pellet_outer_radius = ${pellet_outer_radius}  # 直接使用函数符号进行计算
    #   output_properties = 'burnup'
    #   outputs = exodus
    # []
    [pellet_thermal_eigenstrain]
      type = ComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = ${initial_T}
      thermal_expansion_coeff = ${pellet_thermal_expansion_coef}
      temperature = T
    []
    #化学相关
        [D_fickian]
      type = ParsedMaterial
      property_name = D_fickian
      coupled_variables = 'x T d'
      expression = '(1-0.99*d)*pow(10, -9.386 - 4260/(T) + 0.0012*T*x + 0.00075*T*log10(1+2/(x)))'
      block = pellet
    []
        [D_soret]
      type = DerivativeParsedMaterial
      property_name = D_soret
      coupled_variables = 'x T d'
      material_property_names = 'D_fickian(x,T,d)'
      expression = 'D_fickian * x * (-1380.8 - 134435.5*exp(-x/0.0261)) / ((2.0 + x)/(2.0 * (1.0 - 3.0*x) * (1.0 - 2.0*x)) * 8.314 * T * T)'
      block = pellet
    []
    [pellet_elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus =  ${pellet_E}
      poissons_ratio = ${pellet_nu} 
      # args = 'T'
    []
    [stress]
      type = ComputeCreepPlasticityDeformationStressNew
      elasticity_model = Elasticity
      use_elasticity_model_for_stress = true
      # inelastic_models = 'creep_update'
      # debug = true
      # debug_qp = 0
      # debug_step_interval = 10
      # block = pellet
    []
    [strain_energy_density]
      type = StrainEnergyDensity
      incremental = true
      output_properties = 'strain_energy_density'
      outputs = exodus
    []
    # [creep_eigenstrain]
    #   type = UO2CreepEigenstrainNonAD
    #   # block = pellet
    #   eigenstrain_name = creep_eigenstrain
    #   output_properties = 'effective_creep_strain'
    #   outputs = exodus
    # []
    # # # # 蠕变相关

    # [creep]
    #   type = UO2CreepRateBaseJ2Creep
    #   phase_field = d
    #   degradation_function = g
    #   temperature = T
    #   oxygen_ratio = x
    #   fission_rate = ${fparse fission_rate}
    #   grain_size = ${grain_size}
    #   # 相场断裂相关参数
    #   use_transition_stress = false
    #   use_transient_creep = false
    #   use_three_shear_modulus = false

    #   relative_tolerance = ${creep_relative_tolerance} #蠕变的相对残差
    #   absolute_tolerance = ${creep_absolute_tolerance} #蠕变的绝对残差

    #   # output_properties = 'effective_creep_strain psic_active'
    #   output_properties = 'effective_creep_strain'
    #   outputs = exodus 
    # []
    # [pellet_strain]
    #   type = ComputePlaneSmallStrain
    #   scalar_out_of_plane_strain = scalar_strain_zz
    #   eigenstrain_names = 'thermal_eigenstrain creep_eigenstrain'
    # []

      # 相场断裂模型材料
    [crack_geometric]
    type = DerivativeParsedMaterial
    property_name = alpha
    coupled_variables = 'd'
    expression = 'ksi*d+(1-ksi)*d*d'
    constant_names = 'ksi'
    constant_expressions = '${ksi}'
    derivative_order = 2
  []
  [a1]
    type = ParsedMaterial
    property_name = a1
    material_property_names = 'Gc E l sigma0'
    expression = '4*E*Gc/sigma0/sigma0/l/3.14159'
    output_properties = 'a1'
    outputs = exodus 
  []
    #   [a1]
    #   type = ParsedMaterial
    #   property_name = a1
    #   material_property_names = 'Gc E l sigma0'
    #   expression = '1.5*E*Gc/sigma0/sigma0/l'
    #   # output_properties = 'a1'
    #   # outputs = exodus
    #   block = pellet
    # []
  [degradation]
    type = DerivativeParsedMaterial
    property_name = g
    coupled_variables = 'd'
    material_property_names = 'a1'
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))*(1-eta)+eta
    constant_names = 'p a2 a3 eta'
    constant_expressions = '${m} ${a2} ${a3} 1e-6'
    derivative_order = 2

  []
    [Elasticity]
      type = IsotropicElasticityNonAD
      youngs_modulus = ${pellet_E}
      poissons_ratio = ${pellet_nu}
      phase_field = d
      degradation_function = g
      kinematic_assumption = PLANE_STRAIN
      decomposition = SPECTRAL
      use_threshold = true
      use_history_max = true
      tensile_strength = sigma0
      
      output_properties = 'psie_active'
      outputs = exodus 
      # debug = true
      # debug_qp = 0
      # debug_step_interval = 10
    []
  #   [gps_E]
  #     type = ADParsedMaterial
  #     property_name = gps_E
  #     material_property_names = 'E g'
  #     expression = 'E*g'
  #     block = pellet
  #   []
  #   [gps_elasticity_tensor]
  #     type = ADComputeVariableIsotropicElasticityTensor
  #     base_name = gps
  #     youngs_modulus = gps_E
  #     poissons_ratio = nu
  #     block = pellet
  #   []
  #   [convert_gps_stress]
  #     type = RankTwoTensorMaterialADConverter
  #     ad_props_in = 'stress'
  #     reg_props_out = 'gps_stress'
  #     block = pellet
  #   []
  #   [convert_gps_jacobian_mult]
  #     type = RankFourTensorMaterialADConverter
  #     ad_props_in = 'gps_elasticity_tensor'
  #     reg_props_out = 'gps_Jacobian_mult'
  #     block = pellet
  #   []
  #   [convert_gps_creep]
  #     type = RankTwoTensorMaterialADConverter
  #     ad_props_in = 'creep_strain'
  #     reg_props_out = 'gps_creep_strain'
  #     block = pellet
  #   []
[]


[Functions]
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
  type = PiecewiseLinear
    x = '0.0 2400.0 98400.0 ${endTime__100000} ${endTime__50000} ${endTime}'
    y = '0.0 ${LinearPower0_2} ${LinearPower} ${LinearPower} 0 0'
    scale_factor = ${power_factor}
  []
  # [gap_pressure] #新加的！！！！！！！！！！！！！！！！！！！！！！
  #   #间隙压力随时间的变化
  #   type = PiecewiseLinear
  #   x = '0 ${endTime}'
  #   y = '2.5  14'
  # []
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < 1000, 500,
                  if(t < 3000, 100,
                  if(t < 150000, ${dt},
                  if(t < (${endTime__100000}),${dtMax},
                  if(t < (${endTime__50000}),${dt},10000)))))'
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
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu superlu_dist'52000
#=================================================
  solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  # petsc_options_value = 'lu superlu_dist gmres'
  #=================================================
  # solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type' 
  # petsc_options_value = 'lu'# 60000
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm      100                lu           NONZERO'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
#=================================================
  # solve_type = 'PJFNK'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  petsc_options_value = '201                hypre    boomeramg'  
#=================================================
#李伟老师的配置
  # solve_type = 'PJFNK'
  #开启 SNES–KSP 的 Eisenstat–Walker (EW) 非线性/线性容差耦合策略 ，用来自动调整线性求解器的容差，提高整体效率和鲁棒性。
  # petsc_options = '-snes_ksp_ew -snes_converged_reason -ksp_converged_reason -snes_monitor -ksp_monitor_true_residual'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  # petsc_options_value = 'lu       superlu_dist                  51'
#=================================================



  line_search = 'bt'
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true
  nl_max_its = 200
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  dtmin = ${dtmin}
  dt = ${dt}
  end_time = ${endTime}
  # fixed_point_max_its = 6
  # fixed_point_rel_tol =1e-5 # 固定点迭代的相对容差
  # fixed_point_abs_tol = 1e-6
  accept_on_max_fixed_point_iteration = true
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]
[Outputs]
  [my_checkpoint]
    type = Checkpoint
    time_step_interval = 5    # 每5个时间步保存
    num_files = 50            # 保留最近4个检查点
    wall_time_interval = 600 # 每10分钟保存一次（秒）
  []
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = true
  file_base = '2.1-2D-New2026/1'
[]
