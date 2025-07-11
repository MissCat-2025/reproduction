E = 165e9
nu = 0.3

# conda activate moose && dos2unix 2.2.mechanical.i && mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i 2.2.mechanical.i
[Mesh]
    type = FileMesh
    file = CT_model.e
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
    [vonmises]
      order = CONSTANT
      family = MONOMIAL
    []
  []

  [AuxKernels]
    [vonmises]
        type = ADRankTwoScalarAux
        variable = vonmises
        rank_two_tensor = stress
        scalar_type = VonMisesStress
        execute_on = 'TIMESTEP_END'
    []
[]

[Kernels]
    [solid_x]
        type = ADStressDivergenceTensors
        variable = disp_x
        component = 0
        use_displaced_mesh = false
      []
      [solid_y]
        type = ADStressDivergenceTensors
        variable = disp_y
        component = 1
        use_displaced_mesh = false
    []
[]

[BCs]
    # CT试件标准边界条件：在两个孔上施加相反的压力载荷
    # [upper_load]
    #     type = ADFunctionDirichletBC
    #     variable = disp_y
    #     boundary = 'upper_hole'
    #     # factor = 1e7  # 向上的压力载荷 (Pa)，相当于1 MPa
    #     function = func_upper_load

    # []
    # [lower_load]
    #     type = ADFunctionDirichletBC
    #     variable = disp_y
    #     boundary = 'lower_hole'
    #     # factor = 1e7  # 向下的压力载荷 (Pa)，相当于-1 MPa
    #     function = func_lower_load
    # []
    [upper_load]
      type = Pressure
      boundary = 'upper_hole'
      variable = disp_y
      factor = 5e6
    []
    [lower_load]
      type = Pressure
      boundary = 'lower_hole'
      variable = disp_y
      factor = -5e6
    []
    # # 约束刚体运动 - 固定右边界的x方向位移
    [fixed_x]
        type = DirichletBC
        variable = disp_x
        boundary = 'right'
        value = 0.0
    []
    
    # # 固定一个点的y方向位移（防止刚体运动）
    [fixed_y]
        type = DirichletBC
        variable = disp_y
        boundary = 'right'
        value = 0.0
    []
[]

# [Functions]
#   [func_upper_load]
#     type = ParsedFunction
#     expression = '0.00001*t'
#   []
#   [func_lower_load]
#     type = ParsedFunction
#     expression = '-0.00001*t'
#   []
# []

[Materials]
      # 基本材料属性
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'E nu'
    prop_values = '${E} ${nu}'
  []
  [elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = E
    poissons_ratio = nu
  []
  # 应变计算
  [strain]
    type = ADComputeSmallStrain
  []
  # 应力计算
  [stress]
    type = ADComputeLinearElasticStress  
  []
[]

[Executioner]
  type = Transient # 瞬态求解器
  solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  petsc_options_value = 'lu superlu_dist gmres'
  nl_max_its = 100
  nl_rel_tol = 1e-8 # 非线性求解的相对容差
  nl_abs_tol = 5e-9 # 非线性求解的绝对容差
  l_tol = 1e-8  # 线性求解的容差
  l_abs_tol = 1e-9 # 线性求解的绝对容差
  l_max_its = 500 # 线性求解的最大迭代次数
  accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
  dt = 1
  end_time = 3.7e5 # 总时间24h

[]

[Outputs]
  exodus = true
[]