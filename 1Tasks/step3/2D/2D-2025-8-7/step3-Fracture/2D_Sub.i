[Problem]
    kernel_coverage_check = false
    material_coverage_check = false
  []
  # 双冷却环形燃料几何参数 (单位：mm)(无内外包壳)
# (已验证)这第一步就是测试生成网格文件KAERI_HANARO_UpperRod1.e   注意：这文件生成的是完整的3D模型，而不是1/4模型3D模型
# 语法要求：仅仅为了生成网格文件：Run with --mesh-only:
#https://mooseframework.inl.gov/source/meshgenerators/ConcentricCircleMeshGenerator.html
#conda activate moose && dos2unix Complete2DQuarter.i&&mpirun -n 10 /home/yp/projects/reproduction/reproduction-opt -i Complete2DQuarter.i --mesh-only KAERI_HANARO_UpperRod1.e
#《《下面数据取自[1]Thermomechanical Analysis and Irradiation Test of Sintered Dual-Cooled Annular pellet》》
# 双冷却环形燃料几何参数 (单位：mm)
inclad_inner_diameter = 9.0      # 内包壳内直径
inclad_outer_diameter = 10.14    # 内包壳外直径
pellet_inner_diameter = 10.271         # 芯块内直径
pellet_outer_diameter = 14.685         # 芯块外直径
outclad_inner_diameter = 14.76    # 外包壳内直径
outclad_outer_diameter = 15.9     # 外包壳外直径                   # 轴向长度(m)

# 网格控制参数n_azimuthal = 512时网格尺寸为6.8e-5m
n_radial_inner_clad = 1    # 内包壳径向单元数
mesh_size = 8e-5 #网格尺寸即可
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
n_azimuthal_clad = '${fparse int(3*3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_outer_clad = 1    # 外包壳径向单元数
growth_factor = 1.006       # 径向增长因子
# 计算半径参数 (转换为米)
inner_clad_inner_radius = '${fparse inclad_inner_diameter/2*1e-3}'
inner_clad_outer_radius = '${fparse inclad_outer_diameter/2*1e-3}'
pellet_inner_radius = '${fparse pellet_inner_diameter/2*1e-3}'
pellet_outer_radius = '${fparse pellet_outer_diameter/2*1e-3}'
outer_clad_inner_radius = '${fparse outclad_inner_diameter/2*1e-3}'
outer_clad_outer_radius = '${fparse outclad_outer_diameter/2*1e-3}'

[Mesh]
  [inner_clad1]
    type = AnnularMeshGenerator
    nr = ${n_radial_inner_clad}
    nt = ${n_azimuthal_clad}
    rmin = ${inner_clad_inner_radius}
    rmax = ${inner_clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 10
    boundary_name_prefix = 'inclad'
  []
  [inner_clad]
    type = SubdomainIDGenerator
    input = inner_clad1
    subdomain_id = 1
  []
  [pellet1]
    type = AnnularMeshGenerator
    nr = ${n_radial_pellet}
    nt = ${n_azimuthal}
    rmin = ${pellet_inner_radius}
    rmax = ${pellet_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 20
    boundary_name_prefix = 'pellet'
  []
  [pellet]
    type = SubdomainIDGenerator
    input = pellet1
    subdomain_id = 2
  []
  [outer_clad1]
    type = AnnularMeshGenerator
    nr = ${n_radial_outer_clad}
    nt = ${n_azimuthal_clad}
    rmin = ${outer_clad_inner_radius}
    rmax = ${outer_clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 30
    boundary_name_prefix = 'outclad'
  []
  [outer_clad]
    type = SubdomainIDGenerator
    input = outer_clad1
    subdomain_id = 3
  []
  [combine]
    type = CombinerGenerator
    inputs = 'inner_clad pellet outer_clad'
  []
  [rename1]
    type = RenameBoundaryGenerator
    input = combine
    old_boundary = 'inclad_rmin inclad_rmax pellet_rmin pellet_rmax outclad_rmin outclad_rmax'
    new_boundary = 'inclad_inner inclad_outer pellet_inner pellet_outer outclad_inner outclad_outer'
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
    old_block = '1 2 3'
    new_block = 'inclad pellet outclad'
  []
[]


  
  [Variables]
    [d]
      block = pellet
    []
  []
  
  [AuxVariables]
    [bounds_dummy]
      block = pellet
    []
    [psie_active]
      block = pellet
      order = CONSTANT
      family = MONOMIAL
    []
    [T]
      order = CONSTANT
      family = MONOMIAL
    []
    [a1]
      block = pellet
      family = MONOMIAL
      order = CONSTANT
    []
    [vonMises]
      family = MONOMIAL
      order = CONSTANT
    []
  []
  
  [Bounds]
    [irreversibility]
      type = VariableOldValueBounds
      variable = bounds_dummy
      bounded_variable = d
      bound_type = lower
      block = pellet
    []
    [upper]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = d
      bound_type = upper
      bound_value = 1
      block = pellet
    []
  []
  
  [BCs]
  []
  [Kernels]
    [diff]
      type = ADPFFDiffusion
      variable = d
      fracture_toughness = Gc
      regularization_length = l
      normalization_constant = c0
      block = pellet
    []
    [source]
      type = ADPFFSource
      variable = d
      free_energy = psi
      block = pellet
    []
  []
  
  [Materials]
    [fracture_properties]
      type = ADGenericConstantMaterial
      prop_names = 'l Gc'
      prop_values = '${l} ${Gc}'
      block = pellet
    []
  [a11]
    type = ADParsedMaterial
    property_name = a1
    coupled_variables = 'a1'
    expression = 'a1'
    block = pellet
  []
    #断裂力学-CZM模型
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
    [crack_geometric]
      type = CrackGeometricFunction
      property_name = alpha
      expression = 'd'
      phase_field = d
      block = pellet
    [] 
    [psi]
      type = ADDerivativeParsedMaterial
      property_name = psi
      expression = 'alpha*Gc/c0/l+g*(psie_active)'
      coupled_variables = 'd psie_active'
      material_property_names = 'alpha(d) g(d) Gc c0 l'
      # output_properties = 'psi'
      # outputs = exodus
      block = pellet
    []
  []
  
  [Executioner]
    type = Transient # 瞬态求解器
    solve_type = 'PJFNK' #求解器，PJFNK是预处理雅可比自由牛顿-克雷洛夫方法
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
    # petsc_options_value = '201                hypre    boomeramg'  
    # solve_type = 'NEWTON'
    # solve_type = 'NEWTON'
    petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type  -snes_type'
    petsc_options_value = '201                hypre    boomeramg  vinewtonrsls'  
    # petsc_options_iname = '-pc_type -ksp_type' 
    # petsc_options_value = 'lu gmres' 
    automatic_scaling = true # 启用自动缩放功能，有助于改善病态问题的收敛性
    compute_scaling_once = true  # 每个时间步都重新计算缩放
    # reuse_preconditioner = true
    # reuse_preconditioner_max_linear_its = 20
    nl_max_its = 150
    nl_rel_tol = 5e-6 # 非线性求解的相对容差
    nl_abs_tol = 5e-6 # 非线性求解的绝对容差
    l_tol = 5e-7  # 线性求解的容差
    l_abs_tol = 5e-8 # 线性求解的绝对容差
    l_max_its = 1000 # 线性求解的最大迭代次数
    accept_on_max_fixed_point_iteration = true # 达到最大迭代次数时接受解
    dtmin = 1e-6
    dtmax = 25000
    end_time = 10000000 # 总时间24h
  
    fixed_point_rel_tol =1e-4 # 固定点迭代的相对容差
    [./TimeStepper]
      type = IterationAdaptiveDT
      dt = 10
      growth_factor = 1.1
      cutback_factor = 0.5
      optimal_iterations = 20
      iteration_window = 10
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
        xmax = 0.90 #d变量大于xmax时，标记为粗网格
        y1 = 3.0e7 #vonMises应力小于y1时，标记为粗网格
        y2 = 3.0e7 #vonMises应力大于y2之间时，标记为细网格
        variable = d
        enable_time_check = true
        time_steps = 3
        d_change_threshold = 0.05
        # invert = true
        # third_state = coarsen
      []
    []
  []
  [Outputs]
    # exodus = true
    print_linear_residuals = false
  []
  