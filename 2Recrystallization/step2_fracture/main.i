# mpirun -n 14 /home/yp/projects/reproduction/reproduction-opt -i main.i
# mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i main.i --mesh-only 1.e
# 材料参数
E_modulus = 2.1e11  # 弹性模量 (Pa)
poisson_ratio = 0.3  # 泊松比
K = '${fparse E_modulus/3/(1-2*poisson_ratio)}'
G = '${fparse E_modulus/2/(1+poisson_ratio)}'
Gf = 3 #断裂能
critical_fracture_strength=6.0e7#Pa
length_scale_paramete=1e-5
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'fracture_sub.i'
    cli_args = 'Gc=${Gf};l=${length_scale_paramete};E0=${E_modulus}'
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    # type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    # type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = psie_active
    source_variable = psie_active
  []
  [to_sigma0]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = 'fracture'
    variable = sigma0_field
    source_variable = sigma0_field
  []
[]



[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 1000
  ymax = 1000
  elem_type = QUAD
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [d]
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
[]

[BCs]
  # 左下边界固定（x=0处固定x方向，y=0处固定y方向）
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  # 右边界向上位移
  [top_displacement]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'top'
    function = displacement_function
  []
[]

[Functions]
  [displacement_function]
    type = ParsedFunction
    expression = 't'
  []
[]

[Materials]
  [properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc E0'
    prop_values = '${K} ${G} ${length_scale_paramete} ${Gf} ${E_modulus}'
  []
  # 直接读取CSV文件的AD材料
  [csv_bnds_material]
    type = ADCSVBndsMaterial2
    csv_file = '../step1_bnds/results/final_bnds.csv'
    data_column = 'bnds'
    x_column = 'x'
    y_column = 'y'
    property_name = 'bnds'
    # verbose = true  # 关闭调试输出
    outputs = exodus
  []
  # 使用函数读取晶界数据并计算弹性模量
  # 使用bnds材料属性计算弹性模量
  # 材料属性
  [sigma0_field]
    type = ADDerivativeParsedMaterial
    property_name = sigma0_field
    material_property_names = 'bnds'
    constant_names = 'critical_fracture_strength'
    constant_expressions = '${critical_fracture_strength}'
    expression = 'critical_fracture_strength * (bnds)'  # 直接使用函数符号进行计算
    outputs = exodus
  []
  # 小应变计算
  [strain]
    type = ADComputeSmallStrain
    displacements = 'disp_x disp_y'
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+(1.5*E0*Gc/sigma0_field^2)/l*d*(1+a2*d))*(1-eta)+eta
    phase_field = d
    material_property_names = 'Gc sigma0_field l E0'
    parameter_names = 'p a2 eta'
    parameter_values = '2 2 1e-6'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []  
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = None
    output_properties = 'psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
  []

[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  petsc_options_value = 'lu superlu_dist gmres'
  accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
 
  nl_max_its = 20
  nl_rel_tol = 1e-8   # 非线性求解的相对容差
  nl_abs_tol = 1e-10  # 非线性求解的绝对容差
  l_tol = 1e-8        # 线性求解的容差
  l_abs_tol = 1e-10   # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 1000

  dt = 0.1
[]

[Outputs]
  exodus = true
  file_base = 'results/1'
[]
