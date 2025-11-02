# conda activate moose && mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i ss4_phi_eta_cv_cg.i
# ========== 晶粒序参量演化 ==========
# 本输入文件实现晶粒序参量演化:
#
# 晶粒序参量 (φi):
#    ∂φi/∂t = -Lφ(δF/δφi)
#    实现位置: [Modules]/Nonconserved/phi*


# 参数设置
# --- 基本常数 ---
JtoeV = 6.24150974e18           # J → eV 转换
length_scale = 1e9              # nm 为长度尺度
time_scale = 1.0               # ns 为时间尺度
T = 1000                       # K 温度
k = 8.617333262e-5             # eV/K 玻尔兹曼常数

# --- 自由能势垒 / 势阱参数 ---
omega_SI = 1.54e8              # J/m³ (Potential height)
kp_SI = 6.4e11                 # J/m³ (Prefactor)
Va_SI = 0.04092e-27            # m³/atom (0.04092 nm³)
# b_SI = 0.085e-27               # m³/atom (van der Waals constant)

# --- 梯度能系数 ---
kappa_v_SI = 3.38e-8            # J/m (vacancy gradient)
kappa_phi_SI = 1.67e-9          # J/m (order parameter gradient)
kappa_eta_SI = 1.67e-9          # J/m (bubble order parameter gradient)
# --- 迁移率与界面迁移率 ---
Mv_SI = 2.69e-26                # m⁵/(J·s) vacancy mobility
L_SI = 4.33e-6                 # m³/(J·s) interface mobility
# --- UN 物理参数 ---
# Ev_f_eV = 7.85                  # Vacancy formation energy (eV)
# cv0 = ${fparse exp(-Ev_f_eV/(k*T))}                    # 初始空位浓度 c0ν=exp(-Efv/kBT)
cv_eq = 1.0e-5
cg =  0
cg_eq = 1e-5 #the yield of gas atoms (Xe and Kr) in the fission process of one uranium atom is 0.251 (Olander, 1976),
# nQ = ${fparse 1.0e32/length_scale^3} #量子浓度???
P_bubble = 0.0005 #气泡成核概率

# Pcasc = 0.01  # 级联发生概率（单位体积单位时间的发生概率）

# VG = 1e-1                       # 级联空位浓度最大增加 (无量纲) 假设值 1e-8

# --- 单位换算到 MOOSE 内部（eV、nm尺度） ---
omega = ${fparse omega_SI * JtoeV / length_scale^3}  #  J/m³ (Potential height)
kp = ${fparse kp_SI * JtoeV / length_scale^3}  # J/m³ (Prefactor)
Va = ${fparse Va_SI * length_scale^3}  # m³/atom (0.04092 nm³)
# b = ${fparse b_SI * length_scale^3}  # m³/atom (van der Waals constant) 
kappa_v = ${fparse kappa_v_SI / length_scale * JtoeV}  # J/m (vacancy gradient)
kappa_phi = ${fparse kappa_phi_SI / length_scale * JtoeV} # J/m (order parameter gradient)
kappa_eta = ${fparse kappa_eta_SI / length_scale * JtoeV} # J/m (bubble order parameter gradient)
L = ${fparse L_SI * length_scale^3 / (JtoeV * time_scale)} # m³/(J·s) interface mobility
Mv = ${fparse Mv_SI * length_scale^5 / (JtoeV * time_scale)} # m⁵/(J·s) vacancy mobility


[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  xmax = 512
  ymax = 512
  elem_type = QUAD4
[]
  
[GlobalParams]
  # phi_i 为多晶序参量
  op_num = 3
  var_name_base = phi
  grain_num = 3
#   derivative_order = 2
[]

[Variables]
  # 多晶序参量
  [./PolycrystalVariables]
  [../]
  
  # 气泡序参量 eta（气泡内=1，基体=0）
  [./eta]
    family = LAGRANGE
    order = FIRST
  [../]

  # 守恒相场：空位浓度 cv、气体原子浓度 cg、固溶原子浓度 cs
  [./cv]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[AuxVariables]
  # 晶界掩膜
  [./bnds]
    family = LAGRANGE
    order = FIRST
  [../]
  
  # 可视化变量 Φ = ∑_(j=1)^p φ_j^2 - η^2
  [./Phi]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[AuxKernels]
  # 晶界掩膜计算
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    v = 'phi0 phi1 phi2'
    execute_on = 'initial timestep_end'
  [../]
  
  # 计算可视化变量 Φ = ∑_(j=1)^p φ_j^2 - η^2
  [./Phi_calc]
    type = ParsedAux
    variable = Phi
    coupled_variables = 'phi0 phi1 phi2 eta'
    expression = 'phi0^2 + phi1^2 + phi2^2 - eta^2'
    execute_on = 'initial timestep_end'
  [../]
[]



[Modules]
    [./PhaseField]
      [./Conserved]
        # 空位演化: ∂cν/∂t = ∇·(Mν∇(δF/δcν)) + Pν(基体级联生成)
        [./cv]
          free_energy = f_total
          mobility    = Mv
          kappa       = kappa_cv
          solve_type  = REVERSE_SPLIT
          args = 'eta phi0 phi1 phi2'
        [../]
      []
      [./Nonconserved]
        # 气泡序参量演化: ∂η/∂t = -Lη(δF/δη)
        [./eta]
          free_energy = f_total
          mobility    = L_eta
          kappa       = kappa_eta
          coupled_variables = 'cv phi0 phi1 phi2'
        [../]
      [../]
    [../]
[]

[ICs]
  # 多晶结构：Voronoi
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
  
  # 气泡序参量初始条件
  [./eta_init]
    type = ConstantIC
    variable = eta
    value = 1.0e-30         # 初始无气泡
  [../]

    # 守恒场
    # 初始条件实现动力学方程的初始态：cg=0, c0ν=exp(-Efv/kBT), cs为局部U原子浓度, η=0
    # [./cv_init]
    #   type = ConstantIC
    #   variable = cv
    #   value = ${cv0}   # 从材料属性读取：cv_eq = exp(-Ev/kT)
    # [../]
      [./cv_init]
        type = RandomIC
        variable = cv
        min = ${fparse 1e-10}
        max = ${fparse 1e-1}
        seed = 42
      [../]
[]

[Kernels]
  # 冻结多晶序参量（phi0/phi1/phi2 不参与演化）：∂φi/∂t = 0
  [./phi0_null]
    type = NullKernel
    variable = phi0
  [../]
  [./phi1_null]
    type = NullKernel
    variable = phi1
  [../]
  [./phi2_null]
    type = NullKernel
    variable = phi2
  [../]

  # [./Pv_kernel]
  #   type = VacancySourceKernel
  #   variable = cv
  #   eta = eta
  #   Pcasc = ${Pcasc}
  #   VG = ${VG}
  #   eta_threshold = 0.8
  #   sign = -1
  # []
  # === 空位方程源项：∂cν/∂t = ... + P_ν ===

[]
[Materials]
  # === PhaseField 动力学常数，供 Modules/PhaseField 读取 ===
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names  = 'L_phi kappa_phi L_eta kappa_eta kappa_cv Mv'
    prop_values = '${L} ${kappa_phi} ${L} ${kappa_eta} ${kappa_v} ${Mv}'
  [../]
  # === 基体自由能 f^m（简化版本，用于演示） ===
  [./f_m]
    type = DerivativeParsedMaterial
    property_name = f_m
    coupled_variables = 'eta cv'
    constant_names = 'k T Valpha cg cv_eq cg_eq'
    constant_expressions = '${k} ${T} ${Va} ${cg} ${cv_eq} ${cg_eq}'
    # kT*[ cv ln(cv) + cg ln(cg) + (1-cv-cg) ln(1-cv-cg) ] + (kp/2)*(1 - cv - cg)^2
    expression = '(1-eta)^2*((k*T/Valpha) * ( cv * ( log(cv) - log(cv_eq) ) + (1-cv-cg) * ( log(1-cv-cg) - log(1-cv_eq-cg_eq) ) ))'
    derivative_order = 2
  [../]
    # [./f_m]
    #   type = DerivativeParsedMaterial
    #   property_name = f_m
    #   coupled_variables = 'eta cv'
    #   constant_names = 'k T Valpha cg cv_eq cg0'
    #   constant_expressions = '${k} ${T} ${Va} ${cg} ${cv_eq} ${cg0}'
    #   # kT*[ cv ln(cv) + cg ln(cg) + (1-cv-cg) ln(1-cv-cg) ] + (kp/2)*(1 - cv - cg)^2
    #   expression = '(1-eta)^2*((k*T/Valpha) * ( cv * ( log(cv) - log(cv_eq) ) + cg * ( log(cg) - log(cg0) ) + (1-cv-cg) * ( log(1-cv-cg) - log(1-cv_eq-cg0) ) ))'
    #   derivative_order = 2
    # [../]
  # === 气泡自由能 f^b（简化版本，用于演示） ===
  [./f_b]
    type =DerivativeParsedMaterial
    property_name = f_b
    coupled_variables = 'eta cv'
    constant_names = 'kp cg'
    constant_expressions = '${kp} ${cg}'
    # 常见实现：cg * [ -ln( nQ * (Va/cg - b) ) - 1 ] ；与文献一致的结构
    expression = 'eta^2*(0.5*kp*(1 - cv - cg)^2)'
    derivative_order = 2
  [../]
    # # === 气泡自由能 f^b（简化版本，用于演示） ===
    # [./f_b]
    #   type =DerivativeParsedMaterial
    #   property_name = f_b
    #   coupled_variables = 'eta cv'
    #   constant_names = 'k T Vi b kp nQ cg'
    #   constant_expressions = '${k} ${T}  ${Va} ${b} ${kp} ${nQ} ${cg}'
    #   # 常见实现：cg * [ -ln( nQ * (Va/cg - b) ) - 1 ] ；与文献一致的结构
    #   expression = 'eta^2*((k*T/Vi) * cg * ( -log(nQ*(Vi/cg - b)) - 1.0 ) + 0.5*kp*(1 - cv - cg)^2)'
    #   derivative_order = 2
    # [../]
  # === 多相势 f^{poly}(eta,phi_i)（标准多项式势） ===
  [./f_poly]
    type = DerivativeParsedMaterial
    property_name = f_poly
    coupled_variables = 'eta phi0 phi1 phi2'
    constant_names = 'a_gb a_s'  # 晶界能，表面能
    constant_expressions = '1.5 1.8'
    derivative_order = 2
    # 采用标准多项式：每个序参量的双/四项式 + 取向间二次耦合 + 晶界/表面能耦合
    expression = '(phi0^4+phi1^4+phi2^4)/5 - 0.5*(phi0^2+phi1^2+phi2^2) + (eta^4/4 - 0.5*eta^2) + a_gb*(phi0^2*phi1^2 + phi0^2*phi2^2 + phi1^2*phi2^2) + a_s*eta^2*(phi0^2+phi1^2+phi2^2) + 0.25'
  [../]
  
  # === 化学耦合 f^{chem}（按需，可为0）===
  
  # === 组装总自由能 F = h(η)f^m + j(η)f^b + ωf^poly + f^chem ===
  [./f_0]
    type = DerivativeParsedMaterial
    property_name = f_0
    coupled_variables = 'eta phi0 phi1 phi2 cv'
    material_property_names = 'f_m(eta,cv) f_b(eta,cv) f_poly(eta,phi0,phi1,phi2)'
    constant_names = 'omega'
    constant_expressions = '${omega}'
    expression = 'f_m+ f_b+ omega*f_poly'
    derivative_order = 2
  [../]

  # === 成核概率场 P_bubble ===
  [./nucleation_probability]
    type = ParsedMaterial
    property_name = P_bubble
    coupled_variables = 'eta'
    constant_names = 'prob_base'
    constant_expressions = '${P_bubble}'
    # 简化的成核概率：在基体中均匀分布
    expression = 'if(eta<0.1, prob_base, 0)'
    outputs = exodus
  [../]

  # === 供 PhaseField 使用的最终自由能 ===
  [./f_sum]
    type = DerivativeSumMaterial
    property_name = f_total
    coupled_variables = 'eta phi0 phi1 phi2 cv'
    # 将基体/气泡/多相等基础能量 f0，与两类离散成核惩罚 Fn_* 一并汇总
    sum_materials = 'f_0 Fn_eta'
  [../]
  [./nucleation_eta]
    type = DiscreteNucleation
    property_name = Fn_eta
    op_names  = eta
    op_values = 1
    penalty   = 1e3
    penalty_mode = MIN
    map = map_eta
    outputs = exodus
  [../]
[]

[UserObjects]
  # 晶界 Voronoi 生成
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 2
    int_width = 6.0
  [../]
  
  [./grain_tracker]
    type = GrainTracker
  [../]

  # --- 气泡成核：对 eta（把气泡序参量拉到1） ---
  [./inserter_eta]
    type = DiscreteNucleationInserter
    hold_time  = 1
    probability = P_bubble     # 概率: 基于 ln(cv/cv_eq)
    radius     = 3.2
  [../]
  
  [./map_eta]
    type = DiscreteNucleationMap
    periodic = eta
    inserter = inserter_eta
  [../]
[]

[Adaptivity]
  initial_steps = 2
  initial_marker = combined_marker
  marker = combined_marker
  max_h_level = 3
  [Markers]
    [bnds_marker]
       type = ValueRangeMarker
       lower_bound = 0
       upper_bound = 0.7
       variable = bnds
    []
    [eta_marker]
       type = ValueRangeMarker
       lower_bound = 0.01
       upper_bound = 1
       variable = eta
    []
    [combined_marker]
       type = ComboMarker
       markers = 'bnds_marker eta_marker'
    []
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
  # solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
petsc_options_value = 'asm      31                  preonly       lu           2'
# line_search = bt
  solve_type = 'NEWTON'
  # petsc_options = '-ksp_type=preonly -pc_type=lu -pc_factor_mat_solver_type=mumps'
  # line_search = bt
  # scheme = bdf2

  nl_rel_tol = 1e-10 # 非线性求解的相对容差
  nl_abs_tol = 1e-10 # 非线性求解的绝对容差
  l_tol = 1e-10  # 线性求解的容差
  l_abs_tol = 1e-10 # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 5000

  dt = 0.01
[]

[Outputs]
  exodus = true
  file_base = 'ss4_phi_eta_cv_cg/ss4_phi_eta_cv_cg'
[]
