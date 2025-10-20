# mpirun -n 6 /home/yp/projects/reproduction/reproduction-opt -i cahn_hilliardAH.i
# 位错密度参数 (SI单位制)

f = 1.0e20
T = 373
A_SI = 1.11e7 #J/m^3
B_SI = 1.11e7 #J/m^3
C_SI = 1.665e7 #J/m^3
Sl = 1 


M_SI = 4.55e-30 #m^5/(Js)原为4.55e-31
kappa_c_SI = 2.74e-7 #J/m
kappa_eta_SI = 5.55e-8 #J/m
L_SI = 1.82e-10 #m^3/(Js)论文中为1.82e-14
W_SI = 0.1 #J/m^3
burgers_vector_SI = 3.42e-10    # 伯格斯矢量 (m)
shear_modulus_SI = 36.0e9          # 剪切模量 (Pa)=kg/(m·s²)

JtoeV = 6.24150974e18 # Joule to eV conversion
length_scale = 1e9  #模拟的尺度为纳米
time_scale = 1  #模拟的尺度为纳秒

shear_modulus = '${fparse shear_modulus_SI/(length_scale*time_scale^2)}'
# shear_modulus = '${fparse shear_modulus_SI}'
burgers_vector = '${fparse burgers_vector_SI*length_scale}'
M = ${fparse M_SI*length_scale^5/(JtoeV*time_scale)}
kappa_c = ${fparse kappa_c_SI/length_scale*JtoeV}
kappa_eta = ${fparse kappa_eta_SI/length_scale*JtoeV}
L = ${fparse L_SI*length_scale^3/(JtoeV*time_scale)}
W = ${fparse W_SI*JtoeV/length_scale^3}
A = ${fparse A_SI*JtoeV/length_scale^3}
B = ${fparse B_SI*JtoeV/length_scale^3}
C = ${fparse C_SI*JtoeV/length_scale^3}

[GlobalParams]
  op_num = 3
  var_name_base = gr
  grain_num = 3 # Number of grains
  derivative_order = 2
[]
[Variables]
  [./PolycrystalVariables]
  [../]
[]
[AuxVariables]
  # [./bnds]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
  # [./bnds_old]  # 专门存储上一时间步的bnds值
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
  [./bnds_new]  # 合并的晶界
    order = FIRST
    family = LAGRANGE
  [../]

  [./effective_gr0]
    order = FIRST
    family = LAGRANGE
  [../]
  [./effective_gr1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./effective_gr2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./effective_regr0]
    order = FIRST
    family = LAGRANGE
  [../]
  [./effective_regr0_old]
    order = FIRST
    family = LAGRANGE
  [../]
  [./effective_regr_combined]
    order = FIRST
    family = LAGRANGE
  [../]
[]
[AuxKernels]
  [effective_gr0_aux]
    type = ParsedAux
    variable = effective_gr0
    coupled_variables = 'gr0 effective_regr_combined'
    constant_names = 'center sharpness'
    constant_expressions = '0.2 20'  # center是转换中心，sharpness控制锐度
    expression = 'gr0*(1/(1+exp(sharpness*(effective_regr_combined-center))))'  # sigmoid抑制函数
  []
  [effective_gr1_aux]
    type = ParsedAux
    variable = effective_gr1
    coupled_variables = 'gr1 effective_regr_combined'
    constant_names = 'center sharpness'
    constant_expressions = '0.2 20'
    expression = 'gr1*(1/(1+exp(sharpness*(effective_regr_combined-center))))'
  []
  [effective_gr2_aux]
    type = ParsedAux
    variable = effective_gr2
    coupled_variables = 'gr2 effective_regr_combined'
    constant_names = 'center sharpness'
    constant_expressions = '0.2 20'
    expression = 'gr2*(1/(1+exp(sharpness*(effective_regr_combined-center))))'
  []
  [effective_regr0_aux]
    type = ParsedAux
    variable = effective_regr0
    coupled_variables = 'regr0'
    constant_names = 'threshold_high threshold_low'
    constant_expressions = '2e-8 1e-11'  # 高阈值和低阈值
    expression = 'if(abs(regr0)>=threshold_high, 0.83,max(0.5*(1+tanh(2*(abs(regr0)-0.5*(threshold_high+threshold_low))/(threshold_high-threshold_low))),0))' # 线性插值
  []
  [effective_regr0_aux_old]
    type = ParsedAux
    variable = effective_regr0_old
    coupled_variables = 'effective_regr_combined'
    expression = 'effective_regr_combined'
    execute_on = 'timestep_begin'
  []
  # 合并当前晶界和历史晶界
  [effective_regr_combined_aux]
    type = ParsedAux
    variable = effective_regr_combined
    coupled_variables = 'effective_regr0 effective_regr0_old'
    constant_names = 'threshold_high threshold_low'
    constant_expressions = '0.8 0.2'
    expression = '
    if(effective_regr0_old>threshold_high,effective_regr0_old,
    if(effective_regr0_old<threshold_low,effective_regr0,
    if(threshold_low<effective_regr0<threshold_high,effective_regr0,
    if(threshold_high<effective_regr0,effective_regr0_old,effective_regr0_old))))'  # 取并集
    execute_on = 'timestep_end'
  []

  # 首先计算当前时刻的晶界
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds_new
    v = 'effective_gr0 effective_gr1 effective_gr2 effective_regr_combined'
    execute_on = 'initial timestep_end'
  []

[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 10e3
  ymax = 10e3
  elem_type = QUAD
[]

[Modules]
  [./PhaseField]
    [./Conserved]
      [./c]
        free_energy = f_local0
        mobility = M
        kappa = kappa_c
        solve_type = REVERSE_SPLIT
        args = 'gr0 gr1 gr2 regr0'
      [../]
    [../]
    [./Nonconserved]
      [./regr0]
        args = 'c gr0 gr1 gr2'
        free_energy = f_local0
        kappa = kappa_eta
        mobility = L
      [../]
    [../]
  [../]
[]



[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]
[Kernels]
  # 添加以下NullKernel配置阻止晶粒变量演化
  [gr0_null]
    type = NullKernel
    variable = gr0
  []
  [gr1_null]
    type = NullKernel
    variable = gr1
  []
  [gr2_null]
    type = NullKernel
    variable = gr2
  []
[]
[Materials]
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names  = 'M kappa_c L kappa_eta'
    prop_values = '${M} ${kappa_c} ${L} ${kappa_eta}'
  [../]
    #ρ_d=6340744985690.971((f) ) ̇^(1/6)*exp⁡(-1044.406/T)
  [./mobility]
    type = DerivativeParsedMaterial
    property_name = rho
    constant_names       = 'f T length_scale'
    constant_expressions = '${f} ${T} ${length_scale}'
    expression = 6340744985690.971*f^(1/6)*exp(-1044.406/T)/length_scale^2
    [../]
      [./f_stored]
        type = DerivativeParsedMaterial
        property_name = f_stored
        coupled_variables = 'gr0 gr1 gr2 regr0'
        material_property_names = 'rho'
        constant_names = 'G b_v2'
        constant_expressions = '${shear_modulus} ${burgers_vector}'
        expression = 0.5*G*b_v2*rho*(gr0^2+gr1^2+gr2^2+regr0^2)
        derivative_order = 2
        output_properties = f_stored
        outputs = exodus
      [../]
      [./bulk_free_energy]
        type = DerivativeParsedMaterial
        property_name = f_bulk
        coupled_variables = 'c gr0 gr1 gr2 regr0'
        constant_names = 'W A B C Sl'
        constant_expressions = '${W} ${A} ${B} ${C} ${Sl}'
        expression = W*c^2*(1-c)^2+(-A/2*gr0^2+B/4*gr0^4)+(-A/2*gr1^2+B/4*gr1^4)+(-A/2*gr2^2+B/4*gr2^4)+(-A/2*regr0^2+B/4*regr0^4)+C*(gr0^2*gr1^2+gr0^2*gr2^2+gr1^2*gr0^2+gr1^2*gr2^2+gr2^2*gr0^2+gr2^2*gr1^2+gr0^2*regr0^2+regr0^2*gr0^2+gr1^2*regr0^2+regr0^2*gr1^2+gr2^2*regr0^2+regr0^2*gr2^2)+1/4+Sl*c*(1-(gr0^2+gr1^2+gr2^2+regr0^2))
        derivative_order = 2
        output_properties = f_bulk
        outputs = exodus
      [../]
      [./f_local]
        type = DerivativeParsedMaterial
        property_name = f_local
        coupled_variables = 'c gr0 gr1 gr2 regr0'
        expression = 'f_bulk + f_stored'
        material_property_names = 'f_bulk(c,gr0,gr1,gr2,regr0) f_stored(gr0,gr1,gr2,regr0)'
        # outputs = exodus
      [../]
    [./probability]
      type = ParsedMaterial
      property_name = P
      coupled_variables = 'bnds_new effective_regr_combined'  # 如果还要考虑c的值
      constant_names = 'min_regr0 base_prob'
      constant_expressions = '0.9 1e-3'
      expression = 'max(max((0.7-bnds_new)*base_prob,0) + 
                    if(0.0001<effective_regr_combined<min_regr0, 0.5, 0)*base_prob-if(0.0001>effective_regr_combined, 0.5, 0)*base_prob,0)'
      outputs = exodus
    [../]
  [./nucleation]
    # The nucleation material is configured to insert nuclei into the free energy
    # tht force the concentration to go to 0.95, and holds this enforcement for 500
    # time units.
    type = DiscreteNucleation
    property_name = Fn
    op_names  = c
    op_values = 0.95      # 成核瞬间直接达到0.95
    penalty = 1e6         # 足够大的penalty确保瞬间到达
    penalty_mode = MIN
    map = map
    outputs = exodus
  [../]
    [./nucleation_regr0]
      # The nucleation material is configured to insert nuclei into the free energy
      # tht force the concentration to go to 0.95, and holds this enforcement for 500
      # time units.
      type = DiscreteNucleation
      property_name = Fn_regr0
      op_names  = regr0
      op_values = 0.93      # 成核瞬间直接达到0.95
      penalty = 10#1.5e6         # 足够大的penalty确保瞬间到达
      penalty_mode = MIN
      map = map_regr0
      outputs = exodus
    [../]
  [./free_energy]
    # add the chemical and nucleation free energy contributions together
    type = DerivativeSumMaterial
    f_name = f_local0
    coupled_variables = 'c bnds_new gr0 gr1 gr2 regr0'
    sum_materials = 'f_local Fn Fn_regr0'
  [../]
[]

[UserObjects]
  [./inserter]
    # The inserter runs at the end of each time step to add nucleation events
    # that happend during the timestep (if it converged) to the list of nuclei
    type = DiscreteNucleationInserter
    hold_time = 0.0001#这个参数不错0.001
    probability = P
    radius = 0.2e3
  [../]
  [./map]
    # The map UO runs at the beginning of a timestep and generates a per-element/qp
    # map of nucleus locations. The map is only regenerated if the mesh changed or
    # the list of nuclei was modified.
    # The map converts the nucleation points into finite area objects with a given radius.
    type = DiscreteNucleationMap
    periodic = c
    inserter = inserter
  [../]
    [./inserter_regr0]
      # The inserter runs at the end of each time step to add nucleation events
      # that happend during the timestep (if it converged) to the list of nuclei
      type = DiscreteNucleationInserter
      hold_time = 0.0001#这个参数不错0.001
      probability = P
      radius = 0.2e3
    [../]
    [./map_regr0]
      # The map UO runs at the beginning of a timestep and generates a per-element/qp
      # map of nucleus locations. The map is only regenerated if the mesh changed or
      # the list of nuclei was modified.
      # The map converts the nucleation points into finite area objects with a given radius.
      type = DiscreteNucleationMap
      periodic = regr0
      inserter = inserter_regr0
    [../]
  [voronoi]
    type = PolycrystalVoronoi
    rand_seed = 2
    int_width = 0.2e3
  []
  [grain_tracker]
    type = GrainTracker
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'
  petsc_options = '-ksp_type=gmres -pc_type=lu'
  line_search = bt

  nl_max_its = 20

  nl_rel_tol = 1e-6 # 非线性求解的相对容差
  nl_abs_tol = 1e-5 # 非线性求解的绝对容差
  l_tol = 1e-6  # 线性求解的容差
  l_abs_tol = 1e-5 # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 100

  dt = 0.001
  # [./Adaptivity]
  #   max_h_level = 2
  #   initial_adaptivity = 1
  #   refine_fraction = 0.9
  #   coarsen_fraction = 0.1
  # [../]
[]
[VectorPostprocessors]
  [./bnds_new_export]
    type = NodalValueSampler
    variable = 'bnds_new'
    sort_by = 'id'
    execute_on = 'final'  # 只在最后执行
  [../]
[]
[Outputs]
  exodus = true
  file_base = 'results/2'
    # CSV输出最终的bnds_new数据
  [./csv_final]
    type = CSV
    file_base = 'results/final_bnds_data'
    execute_on = 'final'
    # vectorpostprocessors = 'bnds_new_export'
  [../]
[]
