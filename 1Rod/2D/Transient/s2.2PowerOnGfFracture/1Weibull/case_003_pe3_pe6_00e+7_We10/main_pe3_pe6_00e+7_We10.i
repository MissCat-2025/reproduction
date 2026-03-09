# === 参数研究案例 ===
# pellet_critical_energy: 3
# pellet_critical_fracture_strength: 6.00e+7
# WeibullShape: 10
# 生成时间: 2026-03-04 21:17:20







# # === 参数研究案例（对齐配对） ===
LinearPowerAll=35   #总线功率kw/m

#第一个节点：功率突升的开始时的节点
x = 5  # 功率从0%增加到x%
xTime = '${fparse x*120}'  #s 功率从0%增加到x%的时间,假设正常工况下的温升，即2400s升到20%
LinearPower0 =  '${fparse LinearPowerAll*x*0.01}'  
#第二个节点：功率突升的结束时的节点
PowerTime = 60  #s 功率从x%增加到100%的时间
PowerTimeTotal = '${fparse xTime+PowerTime}'  #s 功率从0%增加到100%的时间
#第二个节点：观察结束时的节点
WatchTime = 40   #s 功率上升后的观察时间(一般以裂纹稳定后为准)
endTime = '${fparse PowerTimeTotal+WatchTime}'  #s 观察结束时的时间

# timePlus1='${fparse 0.05*PowerTime}'
dt1 = '${fparse 0.25/60*PowerTime}' #s 时间步长
# # initial_T_in: 570.7
dtmin = 1e-5
# # initial_T_out: 582.8
# # 生成时间: 2025-09-23 18:08:10
# conda activate moose &&mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i Main.i
initial_T = 550
pellet_nu = 0.345   #RELAP5
pellet_thermal_expansion_coef=1e-5#K-1
pellet_critical_fracture_strength = 6.00e+7#Pa
density_percent = 0.95
# Gc = 6#断裂能
fission_rate = 1.2e19
grain_size =10
pellet_critical_energy = 3# 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
pellet_density='${fparse density_percent*10980}'#10431.0*0.85#kg⋅m-3理论密度为10.980
#几何与网格参数

length_scale_paramete = 4e-5
length_mesh = 2 #裂纹尺度/网格尺寸=(4/length_mesh)
w = 1 #裂纹尖端时，l是mesh_size的2**w倍
mesh_size = '${fparse length_mesh*length_scale_paramete}' #网格尺寸即可
#将下列参数转化为整数
pellet_outer_radius = 4.2e-3#直径变半径，并且单位变mm
n_elems_azimuthal = '${fparse 2*ceil((3.1415*pellet_outer_radius/mesh_size)/2^w)}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int((pellet_outer_radius/mesh_size)/2^w)}'          # 芯块径向网格数（直接取整）

creep_relative_tolerance = 1e-7 #蠕变的相对残差
creep_absolute_tolerance = '${fparse creep_relative_tolerance*0.1}' #蠕变的绝对残差
# 线密度转为体积密度的转换系数
power_factor = '${fparse 1000*1/3.1415926/pellet_outer_radius/pellet_outer_radius}' #新加的！！！！！！！！！！！！！！！！！！！！！！
WeibullSeed = 0
WeibullShape = 10
#相场断裂参数：
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
    input_files = 'sub_pe3_pe6_00e+7_We10.i'
    cli_args = 'l=${length_scale_paramete};mesh_size=${mesh_size};m=${m};w=${w};a2=${a2};a3=${a3};ksi=${ksi};pellet_outer_radius=${pellet_outer_radius};xTime=${xTime};endTime=${endTime};dtmin=${dtmin};dt1=${dt1};PowerTimeTotal=${PowerTimeTotal}'
    execute_on = 'TIMESTEP_END'
        # 强制同步参数
        sub_cycling = false          # 禁止子循环
        catch_up = false              # 子程序失败后尝试用更小dt追赶，失败则触发主程序cutback
        max_catch_up_steps =0      # 允许追赶的最大步数
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
    variable = 'psie_active a1 stress_I Gc sigma0'
    source_variable = 'psie_active a1 stress_I Gc sigma0'
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
  [stress_I]
    order = CONSTANT
    family = MONOMIAL
  []
  [radial_stress]
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
      scale = 1
      shape = ${WeibullShape}
      location = 0.0
      seed = ${WeibullSeed}
      block = pellet
    []
  []
    [grad_T_vec]
    order = CONSTANT
    family = MONOMIAL_VEC
  []
  [grad_T_mag_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [grad_T_mag_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [grad_T_mag]
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
    point1 = '0 0 0.01'        # 圆心坐标
    point2 = '0 0 -0.01'        # 定义旋转轴方向（z轴）
    execute_on = 'TIMESTEP_END'
  [../]
  [./stress_I]
    type = ADRankTwoScalarAux
    scalar_type = MaxPrincipal
    rank_two_tensor = stress
    variable = stress_I
    selected_qp = 0
  [../]
  [radial_stress]
    type = ADRankTwoScalarAux
    rank_two_tensor = stress
    variable = radial_stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
  []
  [grad_T_vec]
    type = ADFunctorElementalGradientAux
    variable = grad_T_vec
    functor = T
  []
  [grad_T_mag_x]
    type = VectorVariableComponentAux
    variable = grad_T_mag_x
    vector_variable = grad_T_vec
    component = x
  []
  [grad_T_mag_y]
    type = VectorVariableComponentAux
    variable = grad_T_mag_y
    vector_variable = grad_T_vec
    component = y
  []
  [grad_T_mag]
    type = VectorMagnitudeAux
    variable = grad_T_mag
    x = grad_T_mag_x
    y = grad_T_mag_y
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
    boundary = 'yplane'
    value = 0
  []
  [x_zero_on_x_plane]
    type = DirichletBC
    variable = disp_x
    boundary = 'xplane'
    value = 0
  []
  [coolant_bc_out]#对流边界条件
  type = ConvectiveFluxFunction
  variable = T
  boundary = 'pellet_outer'
  T_infinity = T_infinity_out
  coefficient = gap_conductance#3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
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
      prop_names = 'l density nu'
      prop_values = '${length_scale_paramete} ${pellet_density} ${pellet_nu}'
      block = pellet
    []
    # 为临界断裂强度生成威布尔分布
    [GcXXX]
      type = ADDerivativeParsedMaterial
      property_name = Gc
      coupled_variables = 'sigma0_F_Density'
      expression = 'sigma0_F_Density*sigma0_F_Density*Gc0'
      constant_names = 'Gc0'
      constant_expressions = '${pellet_critical_energy}'
      output_properties = 'Gc'
      outputs = exodus
      block = pellet
    []
    [sigma0]
      type = ADDerivativeParsedMaterial
      property_name = sigma0
      coupled_variables = 'sigma0_F_Density'
      expression = 'sigma0_F_Density*Strength0'
      constant_names = 'Strength0'
      constant_expressions = '${pellet_critical_fracture_strength}'
      output_properties = 'sigma0'
      outputs = exodus
      block = pellet
    []
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
      expression = '2.334*10^11*(1-2.752*(1-density/10960))*(1-1.0915*10^(-4)*T)'
      constant_names = 'density'
      constant_expressions = '${pellet_density}'
      block = pellet
    []
    # 肿胀应变函数
    [total_power]
      type = ADRimEffertPowerBurnupRod
      power_history = 'power_history'  # 声明使用的函数
      pellet_outer_radius = ${pellet_outer_radius}  # 直接使用函数符号进行计算
      output_properties = 'total_power burnup radial_power_shape'
      outputs = exodus
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
      use_transient_creep = false
      use_three_shear_modulus = false

      relative_tolerance = ${creep_relative_tolerance} #蠕变的相对残差
      absolute_tolerance = ${creep_absolute_tolerance} #蠕变的绝对残差

      # output_properties = 'effective_creep_strain psic_active'
      output_properties = 'effective_creep_strain'
      outputs = exodus 
    []
    [pellet_strain]
      type = ADComputePlaneSmallStrain
      eigenstrain_names = 'thermal_eigenstrain'
    []

      # 相场断裂模型材料
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'ksi*d+(1-ksi)*d*d'
    parameter_names = 'ksi'
    parameter_values = '${ksi}'
    phase_field = d
  []
  [a1]
    type = ADDerivativeParsedMaterial
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
    parameter_values = '${m} ${a2} ${a3} 1e-6'
  []
    [Elasticity]
      type = IsotropicElasticity
      youngs_modulus = E
      poissons_ratio = nu
      phase_field = d
      degradation_function = g
      kinematic_assumption = PLANE_STRESS
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
        [strain_energy_density]
      type = ADStrainEnergyDensity
      incremental = false
      output_properties = 'strain_energy_density'
      outputs = exodus
    []
[]


[Functions]
    [gap_conductance]
    type = PiecewiseLinear
    x = '0 ${endTime}'
    y = '3400 3400'
    scale_factor = 1         # 保持原有的转换因子
  []
  [power_history] #新加的！！！！！！！！！！！！！！！！！！！！！！
  type = PiecewiseLinear
    x = '0.0 ${xTime} ${PowerTimeTotal} ${endTime}'
    y = '0.0 ${LinearPower0} ${LinearPowerAll} ${LinearPowerAll}'
    scale_factor = ${power_factor}
  []
  [gap_pressure] #新加的！！！！！！！！！！！！！！！！！！！！！！
    #间隙压力随时间的变化
    type = PiecewiseLinear
    x = '0 ${endTime}'
    y = '2.5  2.5'
  []
  [dt_limit_func]
    type = ParsedFunction
    expression = 'if(t < ${xTime}, 100,
                  if(t < (${PowerTimeTotal}),${dt1},
                  if(t < ${endTime},0.025,${dt1})))'
  []
  [T_infinity_out]
      #间隙压力随时间的变化
    type = PiecewiseLinear
    x = '0 ${endTime}'
    y = '${initial_T} ${initial_T}'
    scale_factor = 1
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
  # petsc_options_value = 'lu gmres' 57000
#=================================================
  solve_type = 'NEWTON'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = 'lu       superlu_dist                  51'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu superlu_dist'52000
#=================================================
  # solve_type = 'PJFNK'
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
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  # petsc_options_value = '201                hypre    boomeramg'  
#=================================================
  line_search = 'bt'
  automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  compute_scaling_once = true  # 每个时间步都重新计算缩放
  nl_max_its = 50
  nl_rel_tol = 5e-6
  nl_abs_tol = 5e-7
  dtmin = ${dtmin}
  # dt = ${dt}
  end_time = ${endTime}
  fixed_point_rel_tol =1e-8 # 固定点迭代的相对容差
  fixed_point_abs_tol = 1e-10
  accept_on_max_fixed_point_iteration = true
  [TimeStepper]
    type = FunctionDT
    function = dt_limit_func
  []
[]
[Postprocessors]
  [pellet_total_strain_energy]
    type = ElementIntegralMaterialProperty
    mat_prop = strain_energy_density
    block = pellet
    execute_on = 'TIMESTEP_END'
  []
[]
[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = ${w}
  [Markers]
    [marker]
      type = PhasePiledFractureHSMarker
      von_mises_variable = stress_I
      sigma0 = sigma0
      x1 = 0.000001 #d变量小于x1时，标记为粗网格
      x2 = 0.005 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.05 #d变量大于xmax时，一定是细网格
      y1 = 0.45 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.5 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.01
      stress_change_threshold = 1e6
    []
  []
[]
[Outputs]
 
  exodus = true #表示输出exodus格式文件
  print_linear_residuals = false
  hide = 'grad_T_vec grad_T_mag_x grad_T_mag_y'
  file_base = 'Output/C_${PowerTime}'
[]
