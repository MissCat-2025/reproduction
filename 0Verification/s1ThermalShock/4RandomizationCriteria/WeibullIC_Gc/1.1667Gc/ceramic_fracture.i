# 陶瓷片热冲击实验 - 热弹性模拟部分
# conda activate moose &&mpirun -n 4 /home/yp/projects/reproduction/reproduction-opt -i ceramic_fracture.i --timing > timing.txt

[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

# 陶瓷材料参数
E_ceramic = 370e9       # 陶瓷杨氏模量 (Pa)
nu_ceramic = 0.3       # 陶瓷泊松比
K = '${fparse E_ceramic/(3*(1-2*nu_ceramic))}'
G = '${fparse E_ceramic/(2*(1+nu_ceramic))}'
alpha_ceramic = 7.5e-6    # 陶瓷热膨胀系数 (1/°C)
k_ceramic = 31          # 陶瓷导热系数 (W/m·K)
cp_ceramic = 880        # 比热容 (J/kg·K)
rho_ceramic = 3980      # 密度 (kg/m³)

T_initial_condition = 573.15
# 断裂参数
Gc0 = 42.47               # 断裂能 (J/m^2)
# Gc=${fparse (1+grid_sizes/2/length_scale_paramete)*Gc0}
l = 0.1e-3 #裂纹附近加密4倍)                # 相场正则化长度 (m)
nh = 2 # 加密次数
nx = '${fparse int(25e-3/(nh*nh*l/3))}'
ny = '${fparse int(5e-3/(nh*nh*l/3))}'
ft = 180e6                # 抗拉强度 (Pa)
Gc=${fparse (1+1/6)*Gc0}
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}            # 25mm / 0.05mm = 500
    ny = ${ny}            # 5mm / 0.05mm = 100
    xmax = 25e-3
    ymax = 5e-3
  []
[]
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'ceramic_fracture_Sub.i'
    cli_args = 'Gc=${Gc};l=${l}'
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
  [to_psie_active]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active MaxPrincipal a1 sigma0'
    source_variable = 'psie_active MaxPrincipal a1 sigma0'
  []
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temp]
    initial_condition = '${T_initial_condition}'  # 初始温度300°C
  []
  [strain_zz]
  []
[]

[AuxVariables]
  [d]                         # 相场变量
  []
  [MaxPrincipal]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  # 力学平衡
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
  
  # 热传导
  [heat_conduction]
    type = ADHeatConduction
    variable = temp
  []
  [heat_dt]
    type = ADHeatConductionTimeDerivative
    variable = temp
  []
[]

[AuxKernels]
  [MaxPrincipalStress]
    type = ADRankTwoScalarAux
    variable = MaxPrincipal
    rank_two_tensor = stress
    scalar_type = MaxPrincipal
    execute_on = 'TIMESTEP_END'
  []
[]
[AuxVariables]
  [sigma0_F_Density]
    family = MONOMIAL
    order = CONSTANT
    [InitialCondition]
      type = WeibullICDensity
      scale = 1
      shape = 30
      location = 0.0
      seed = 0
    []
  []
[]
[BCs]
  # 力学边界条件 - 右侧对称面
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  [symm_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  
  # 热边界条件
  [left_temp]
    type = DirichletBC
    variable = temp
    boundary = 'top left'
    value = 293.15  # 水淬温度20°C
  []
  # 右侧为绝热边界 - 不需要额外的边界条件
[]

[Materials]
  # 热物理属性
  [thermal]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat density K G'
    prop_values = '${k_ceramic} ${cp_ceramic} ${rho_ceramic} ${K} ${G}'
  []
  
  # 断裂属性
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'E nu l Gc'
    prop_values = '${E_ceramic} ${nu_ceramic} ${l} ${Gc}'
  []
  [sigma0]
    type = ADDerivativeParsedMaterial
    property_name = sigma0
    coupled_variables = 'sigma0_F_Density'
    expression = 'sigma0_F_Density*ft'  # 直接使用辅助变量的值
    constant_names = 'ft'
    constant_expressions = ' ${ft}'
    output_properties = 'sigma0'
    outputs = exodus
  []
  # 相场断裂模型材料
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = '2*d-d*d'
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
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))
    phase_field = d
    material_property_names = 'a1'
    parameter_names = 'p a2 a3'
    parameter_values = '2 -0.5 0'
  []
  
  [eigenstrain]
    type = ADComputeThermalExpansionEigenstrain
    eigenstrain_name = thermal_eigenstrain
    stress_free_temperature = 573.15  # 应力自由温度为初始温度
    thermal_expansion_coeff = ${alpha_ceramic}
    temperature = temp
  []

  [strain]
    type = ADComputePlaneSmallStrain
    eigenstrain_names = thermal_eigenstrain
  []
  [elasticity]
    type = SmallDeformationH
    youngs_modulus = E
    poissons_ratio = nu
    tensile_strength = sigma0
    phase_field = d
    degradation_function = g
    output_properties = 'psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
  []
[]

[Adaptivity]
  initial_marker = boundary
  initial_steps = ${nh}
  marker = marker
  max_h_level = ${nh}
  [Markers]
    [marker]
      type = PhasePiledFractureHSMarker
      von_mises_variable = MaxPrincipal
      sigma0 = sigma0
      x1 = 1e-6 #d变量小于x1时，标记为粗网格
      x2 = 0.05 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.1 #d变量大于xmax时，一定是细网格
      y1 = 0.3 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.6 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.02
      stress_change_threshold = 1e6
    []
    [boundary]
      type = BoundaryMarker
      next_to = 'top left'
      distance = 1e-4
      mark = refine
    []
  []
[]

[Executioner]
  type = Transient
  
  solve_type = 'NEWTON'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
  petsc_options_value = '201                hypre    boomeramg'
  automatic_scaling = true
  
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  
  # 时间步长设置
  dt = 0.1e-3  # 小的时间步长以捕捉快速的温度变化
  end_time = 50e-3  # 总模拟时间
  fixed_point_max_its = 4
  fixed_point_rel_tol = 1e-5
  fixed_point_abs_tol = 1e-6
  accept_on_max_fixed_point_iteration = true
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
  file_base = 'outputs/${l}'
[]