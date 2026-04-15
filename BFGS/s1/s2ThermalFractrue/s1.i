#This input uses PhaseField-Nonconserved Action to add phase field fracture bulk rate kernels
#弹性
E = 200e9
Nu = 0.33
#断裂
c0 = 3.1415
Ksi = 2
Gc = 5
Sigma = 6e7
b = 6e-5

m = 2

a2 = -0.5
a3 = 0

#几何与网格参数
w = 1 #裂纹尖端时，l是mesh_size的2**w倍
mesh_size = '${fparse 12e-5}' #网格尺寸即可
#将下列参数转化为整数
pellet_outer_radius = 4.2e-3#直径变半径，并且单位变mm
n_elems_azimuthal = '${fparse 2*ceil((3.1415*pellet_outer_radius/mesh_size)/2^w)}'  # 周向网格数（向上取整）
n_elems_radial_pellet = '${fparse int((pellet_outer_radius/mesh_size)/2^w)}'  

# 芯块径向网格数（直接取整）
LinearPower = 35
power_factor = '${fparse LinearPower*1000*1/3.1415926/pellet_outer_radius/pellet_outer_radius}' #新加的！！！！！！！！！！！！！！！！！！！！！！

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

[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[Modules]
  [./PhaseField]
    [./Nonconserved]
      [./d]
        free_energy = F
        kappa = kappa_op
        mobility = b
      [../]
    [../]
  [../]
[]
# [Physics]
#   [./SolidMechanics]
#     [./QuasiStatic]
#       [./mech]
#         add_variables = true
#         strain = SMALL
#         planar_formulation = WEAK_PLANE_STRESS
#         eigenstrain_names = thermal_eigenstrain
#         temperature = T
#       [../]
#     [../]
#   [../]
# []

[Kernels]
    #力平衡方程
    [solid_x2]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
    []
    [solid_y2]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
    []
    [./solid_z]
      type = WeakPlaneStress
      variable = strain_zz
    [../]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = d
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = d
  [../]

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
      value = ${power_factor}
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
  T_infinity = 293.15
  coefficient = 3500 #3500 W·m-2 K-1！！！！！！！！！！！！！！！！！！！！！！！！！！！
  []
  #芯块包壳间隙压力边界条件
  [gap_pressure_fuel_x]
    type = Pressure
    variable = disp_x
    boundary = 'pellet_outer'
    factor = 2e6 # 间隙压力2.5MPa
    # function = gap_pressure #新加的！！！！！！！！！！！！！！！！！！！！！！
    # use_displaced_mesh = true
  []
  [gap_pressure_fuel_y]
    type = Pressure
    variable = disp_y
    boundary = 'pellet_outer'
    factor = 2e6
    # function = gap_pressure #新加的！！！！！！！！！！！！！！！！！！！！！！
    # use_displaced_mesh = true
  []
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
    [sigma0_field]
    family = MONOMIAL
    order = CONSTANT
    [InitialCondition]
      type = WeibullIC
      scale = ${Sigma}
      shape = 50
      location = 0.0
      seed = 0
      block = pellet
    []
  []
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
  [../]
  [radial_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = radial_stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
  []
    [copy_sigma0]
    type = MaterialRealAux
    variable = sigma0_field
    property = sigma0
    execute_on = 'initial'
    block = pellet
  []
[]

[Materials]
    [pellet_thermal_eigenstrain]
      type = ComputeThermalExpansionEigenstrain
      eigenstrain_name = thermal_eigenstrain
      stress_free_temperature = 293.15
      thermal_expansion_coeff = 1e-5
      temperature = T
      block = pellet
    []



  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'Gc b visco c0 gc_prop l density thermal_conductivity specific_heat'
    prop_values = '${Gc} ${b} 1e-1 ${c0} ${Gc} ${b} 10431.0 5 280'
  [../]
        # 为临界断裂强度生成威布尔分布
    [sigma0_mat]
      type = ParsedMaterial
      property_name = sigma0
      coupled_variables = 'sigma0_field'
      expression = 'sigma0_field'  # 直接使用辅助变量的值
      block = pellet
    []
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'Gc visco'
    property_name = L
    expression = '1.0/(Gc * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'Gc b c0'
    property_name = kappa_op
    expression = 'Gc * b * 2 / c0'
  [../]
  [elast_tensor2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${Nu}
  []
      [pellet_strain]
      type = ComputePlaneSmallStrain
      eigenstrain_names = 'thermal_eigenstrain'
    []
  [./elastic]
    type = ComputeLinearElasticPFFractureStress
    c = d
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'fracture_energy'
    I_name = 'degradation'
    barrier_energy = 'barrier'
    use_current_history_variable = true
    decomposition_type = strain_spectral
    # outputs = exodus
  [../]
  [a1]
        # a1 = '${fparse 2*Ksi/c0/b*E/Sigma/Sigma}'
    type = ParsedMaterial
    property_name = a1
    material_property_names = 'Gc l sigma0'
    expression = '2*Ksi/l/c0*E*Gc/sigma0/sigma0'
    constant_names = 'E Ksi c0' 
    constant_expressions = '${E} ${Ksi} ${c0}'
  []
    [./degradation]
    type = DerivativeParsedMaterial
    property_name = degradation
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))*(1-eta)+eta
    coupled_variables = 'd'
    constant_names = 'p a2 a3 eta '
    constant_expressions = '${m} ${a2} ${a3} 1e-7'
    material_property_names = 'a1'
    derivative_order = 2
  [../]
  [./CrackGeometricFunction]
    type = DerivativeParsedMaterial
    property_name = alpha
    coupled_variables = 'd'
    expression = 'Ksi*d+d*d*(1-Ksi)'
    constant_names       = 'Ksi'
    constant_expressions = '${Ksi}'
    derivative_order = 2
  [../]
  [./fracture_energy]
    type = DerivativeParsedMaterial
    property_name = fracture_energy
    coupled_variables = 'd'
    material_property_names = 'Gc b alpha(d) c0'
    expression = 'Gc*alpha/(c0*b)'
    derivative_order = 2
    output_properties = 'fracture_energy'
    outputs = exodus
  [../]
  [./fracture_driving_energy]
    type = DerivativeParsedMaterial
    coupled_variables = d
        # expression = 'elastic_energy + fracture_energy + 1e1 * (max(0, -d)^2 + max(0, d-1)^2)'
    expression = 'elastic_energy + fracture_energy'
    material_property_names = 'elastic_energy(d) fracture_energy(d) '
    # expression = 'elastic_energy + fracture_energy + 1e10 * (max(0, -d)^2 + max(0, d-1)^2)'
    derivative_order = 2
    # output_properties = 'dFdd'
    outputs = exodus
    property_name = F
  [../]
  [./barrier_energy]
    type = ParsedMaterial
    property_name = barrier
    material_property_names = 'sigma0'
    expression = 'sigma0*sigma0/E'
    constant_names       = 'E'
    constant_expressions = '${E}'
    output_properties = 'barrier'
    outputs = exodus
  [../]
[]
[Executioner]
  type = Transient

  # solve_type = PJFNK
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm      31                  preonly       lu           1'
  solve_type = NEWTON
  #   petsc_options_iname = '-snes_type        -snes_qn_type   -snes_qn_scale_type' 
  # petsc_options_value = 'qn               lbfgs           DIAGONAL'
  #   petsc_options_iname = '-pc_type   -snes_type        -snes_qn_type ' 
  # petsc_options_value = 'lu         qn               lbfgs'  
  petsc_options_iname = '-pc_type   -snes_type        -snes_qn_type   -snes_qn_scale_type' 
  petsc_options_value = 'lu         qn               lbfgs           JACOBIAN'  
    # automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
  # compute_scaling_once = true  # 每个时间步都重新计算缩放
  nl_rel_tol = 1e-3
  l_max_its = 10
  nl_max_its = 20
  # line_search = 'bt'
  dt = 0.05
  dtmin = 1e-9
  num_steps = 500
[]
[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = ${w}
  [Markers]
    [marker]
      type = PhasePiledFractureHSMarkerNoAD
      von_mises_variable = stress_I
      sigma0 = sigma0
      x1 = 0.000001 #d变量小于x1时，标记为粗网格
      x2 = 0.005 #d变量在x1和x2之间时，标记为细网格
      xmax = 0.12 #d变量大于xmax时，一定是细网格
      y1 = 0.45 #vonMises应力小于y1时，标记为粗网格
      y2 = 0.6 #vonMises应力大于y2之间时，标记为细网格
      variable = d
      timeD = 3
      timeStress = 5
      d_change_threshold = 0.02
      stress_change_threshold = 1e6
    []
  []
[]
[Outputs]
  exodus = true
  file_base = 'Output/C_1'
[]
