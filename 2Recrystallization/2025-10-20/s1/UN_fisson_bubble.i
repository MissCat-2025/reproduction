# conda activate moose && mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i UN_fisson_bubble.i
# 参数设置
# --- 基本常数 ---
JtoeV = 6.24150974e18           # J → eV 转换
length_scale = 1e6              # nm 为长度尺度
time_scale = 1                  # ns 为时间尺度
T = 1000                        # K 温度

# --- 自由能势垒 / 势阱参数 ---
omega_SI = 1.54e8               # J/m³ (Potential height)
kp_SI = 6.4e11                  # J/m³ (Prefactor)
Va_SI = 0.04092e-27             # m³/atom (0.04092 nm³)
b_SI = 0.085e-27                # m³/atom (van der Waals constant)

# —— 化学耦合项：S_l（先用 0.05*omega_SI，数值可扫参） —— （先用0.05，后面可以改为0.01或0.1）
S_l_SI = ${fparse 0.05*omega_SI}                     # J/m^3
S_l     = ${fparse S_l_SI * JtoeV / length_scale^3}  # eV/nm^3

# --- 梯度能系数 ---
kappa_v_SI = 3.38e-8            # J/m (vacancy gradient)
kappa_eta_SI = 1.67e-9          # J/m (order parameter gradient)

# --- 迁移率与界面迁移率 ---
Mv_SI = 2.69e-26                # m⁵/(J·s) vacancy mobility
Mg_SI = 1.33e-28                # m⁵/(J·s) gas atom mobility
L_SI = 4.33e-6                  # m³/(J·s) interface mobility

# --- UN 物理参数 ---
# Ev_f_eV = 7.85                  # Vacancy formation energy
# Ev_m_eV = 2.2                   # Vacancy migration energy
# Eg_m_eV = 1.79                  # Gas atom migration energy
# Dv_m_SI = 4.0e-7                # Vacancy diffusivity (m²/s)
# Dg_m_SI = 4.05e-7               # Gas diffusivity (m²/s)
# gamma_gb_SI = 0.70              # J/m² Grain boundary energy
# gamma_s_SI = 1.62               # J/m² Surface energy

# --- 单位换算到 MOOSE 内部（eV尺度） ---
omega = ${fparse omega_SI * JtoeV / length_scale^3}
kp = ${fparse kp_SI * JtoeV / length_scale^3}
Va = ${fparse Va_SI * length_scale^3}
b = ${fparse b_SI * length_scale^3}

kappa_v = ${fparse kappa_v_SI / length_scale * JtoeV}
kappa_eta = ${fparse kappa_eta_SI / length_scale * JtoeV}

Mv = ${fparse Mv_SI * length_scale^5 / (JtoeV * time_scale)}
Mg = ${fparse Mg_SI * length_scale^5 / (JtoeV * time_scale)}
L = ${fparse L_SI * length_scale^3 / (JtoeV * time_scale)}

# gamma_gb = ${fparse gamma_gb_SI * JtoeV / length_scale^2}
# gamma_s = ${fparse gamma_s_SI * JtoeV / length_scale^2}

# --- 注释 ---
# Ev_f_eV, Ev_m_eV, Eg_m_eV 保留在 eV 单位用于计算扩散激活能。
# 如果要耦合温度依赖扩散率，可按 D = D0 * exp(-Em/kT) 转换为有效迁移率。

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 300
  ny = 300
  xmax = 10e3
  ymax = 10e3
  elem_type = QUAD4
[]
  
[GlobalParams]
  # phi_i 为多晶序参量；cv、cg、cs 守恒；eta 为气泡序参量（非守恒）
  op_num = 3
  var_name_base = phi
  grain_num = 3
  derivative_order = 2
[]

[Variables]
  # 多晶序参量（冻结与否在 Kernels 里再加 NullKernel）
  [./PolycrystalVariables]
  [../]

  # --- 气泡模型中的场变量（论文 2.1、2.2 节）---
  # 守恒相场：空位浓度 cv、气体原子浓度 cg、固溶原子浓度 cs
  [./cv]
    family = LAGRANGE
    order = FIRST
  [../]
  [./cg]
    family = LAGRANGE
    order = FIRST
  [../]
  [./cs]
    family = LAGRANGE
    order = FIRST
  [../]

  # 非守恒相场：气泡序参量 eta（气泡内=1，基体=0）
  [./eta]
    family = LAGRANGE
    order = FIRST
  [../]
[]
[AuxVariables]
  # 晶界掩膜（用于后续成核概率/界面能等）：按你模板的 bnds 逻辑
  [./bnds]
    family = LAGRANGE
    order = FIRST
  [../]

  # 气泡掩膜：把连续的 eta 转成 0~1 的有效“是否为气泡”场（供成核与禁核用）
  [./bubble_mask]
    family = LAGRANGE
    order = FIRST
  [../]
  [./bubble_mask_old]
    family = LAGRANGE
    order = FIRST
  [../]
  [./bubble_mask_combined]
    family = LAGRANGE
    order = FIRST
  [../]
[]
[AuxKernels]
  # 1) 晶界掩膜（供成核概率 P(bnds) 使用）：对 phi_i 计算 bnds
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    v = 'phi0 phi1 phi2'
    execute_on = 'initial timestep_end'
  [../]

  # 2) 气泡掩膜（把连续的 eta → 0~1 掩膜；用于禁在已成核处重复成核）
  [./bubble_mask_aux]
    type = ParsedAux
    variable = bubble_mask
    coupled_variables = 'eta'
    constant_names = 'thigh tlow eps'
    constant_expressions = '0.8 0.2 1e-12'
    # eta ≥ thigh 视为稳定气泡区；中间用 tanh 过渡；避免除零加 eps
    expression = 'if(eta>=thigh,0.83,max(0.5*(1+tanh(2*(eta-0.5*(thigh+tlow))/(thigh-tlow+eps))),0))'
  [../]

  # 3) 历史掩膜（时间步开始时拷贝上一时刻合并后的气泡掩膜）
  [./bubble_mask_old_aux]
    type = ParsedAux
    variable = bubble_mask_old
    coupled_variables = 'bubble_mask_combined'
    expression = 'bubble_mask_combined'
    execute_on = 'timestep_begin'
  [../]

  # 4) 合并掩膜（时间步结束：保持已成核区域，新增则更新，其它延续）
  [./bubble_mask_combined_aux]
    type = ParsedAux
    variable = bubble_mask_combined
    coupled_variables = 'bubble_mask bubble_mask_old'
    constant_names = 'thigh tlow'
    constant_expressions = '0.8 0.2'
    expression = 'if(bubble_mask_old>thigh,bubble_mask_old,if(bubble_mask_old<tlow,bubble_mask,if(tlow<bubble_mask<thigh,bubble_mask,bubble_mask_old)))'
    execute_on = 'timestep_end'
  [../]
[]
[Modules]
  [./PhaseField]
    [./Conserved]
      [./cv]
        free_energy = f_bubble_total
        mobility    = Mv
        kappa       = kappa_cv
        solve_type  = REVERSE_SPLIT
        coupled_variables = 'cg eta phi0 phi1 phi2'
      [../]
      [./cg]
        free_energy = f_bubble_total
        mobility    = Mg
        kappa       = kappa_cg
        solve_type  = REVERSE_SPLIT
        coupled_variables = 'cv eta phi0 phi1 phi2'
      [../]
      # [./cs]
      #   free_energy = f_bubble_total
      #   mobility    = Ms
      #   kappa       = kappa_cs
      #   solve_type  = REVERSE_SPLIT
      #   coupled_variables = 'cv cg eta phi0 phi1 phi2'
      # [../]
    [../]

    [./Nonconserved]
      [./eta]
        free_energy = f_bubble_total
        mobility    = L_eta
        kappa       = kappa_eta
        coupled_variables = 'cv cg phi0 phi1 phi2'
      [../]
    [../]
  [../]
[]
[ICs]
    # 多晶结构：Voronoi 与你现有用法一致
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
  
    # 守恒场与气泡序参量的初始值（可按文献/表格替换）
  [./cv_init]
    type = ConstantIC
    variable = cv
    value = 1.0e-6      # cv0：空位初始浓度
  [../]
  [./cg_init]
    type = ConstantIC
    variable = cg
    value = 1.0e-8      # 初始极低游离气体原子
  [../]
  [./cs_init]
    type = ConstantIC
    variable = cs
    value = 1.0         # cs_eq：单相 UN 的基体占位
  [../]
  [./eta_init]
    type = ConstantIC
    variable = eta
    value = 0.0         # 初始无气泡
  [../]
  
    # （可选）给 cv/cg/eta 加极小随机扰动以避免完全对称：
    # [./cv_jitter] type=RandomIC variable=cv min=0.0 max=1e-8 [../]
[]
[Kernels]
    # 冻结多晶序参量（phi0/phi1/phi2 不参与演化）
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
  
    # 冻结 cs：不参与演化但保持为活动变量
  [./cs_null]
    type = NullKernel
    variable = cs
  [../]
[]
[Materials]
    # === 常量/参数（eV–nm–ns 体系）；请把数值替换成你表里的 ===
    # [./params]
    # #   type = GenericConstantMaterial
    #   # omega: 多相势系数；kp: 体积分数约束前因子；Va: 基质原子体积 (nm^3/atom)；b: vdW 常数 (nm^3/atom)
    #   # S_l: 化学耦合系数（如需）；kB_eVK: eV/K 中的 Boltzmann 常数
    # #   prop_names  = 'omega kp Va b S_l kB_eVK'
    # #   prop_values = '1.0   4000.0 0.04092 0.085 1.0  8.617333262e-5'
    #   constant_names       = 'Valpha cnu0 cg0 eps'
    #   constant_expressions = '${Va}  1.0e-6 1.0e-8 1.0e-12'   # cnu0/cg0 请换成你文档里的参考值
    #   constant_names       = 'Vi b kp nQ eps'
    #   constant_expressions = '${Va} ${b} ${kp} 1.0 1.0e-12'
    # [../]

    # 若你已在别处定义 T（温度，K），这里把 kT 提供给自由能表达式用
    [./kT_mat]
      type = ParsedMaterial
      property_name = kT
      constant_names = 'kB_eVK T'
      constant_expressions = '8.617333262e-5 ${T}'
      expression = 'kB_eVK*T'
    [../]
  
    # === PhaseField 动力学常数，供 Modules/PhaseField 读取 ===
    [./pfmobility]
      type = GenericConstantMaterial
    prop_names  = 'Mv Mg L_eta kappa_cv kappa_cg kappa_eta'
    prop_values = '${Mv} ${Mg} ${L} ${kappa_v} ${kappa_v} ${kappa_eta}'
      # 注：若暂不演化 cs，可将 Ms=0 且 kappa_cs=0；反之请给出实际值
    [../]
  
    # === 基体自由能 f^m(cv,cg,cs) ===
    # 采用理想混合熵 + 体积分数约束（kp）常见写法；与你 Word 文档思路一致
    [./f_m]
      type = DerivativeParsedMaterial
      property_name = f_m
      coupled_variables = 'cv cg cs'
      material_property_names = 'kT()'
      constant_names = 'kp Valpha cnu0 cg0 eps'
      constant_expressions = '${kp} ${Va} 1.0e-6 1.0e-8 1.0e-12'
      # kT*[ c ln c + s ln s + (1-c-s) ln(1-c-s) ] + (kp/2)*(1 - c - s)^2
      # 其中 c≡cg，s≡cv；为避免 log(0) 用 eps 防护
      expression = '(kT/Valpha) * ( if(cv>eps,cv,eps) * ( log(if(cv>eps,cv,eps)) - log(if(cnu0>eps,cnu0,eps)) ) + if(cg>eps,cg,eps) * ( log(if(cg>eps,cg,eps)) - log(if(cg0>eps,cg0,eps)) ) + if(1-cv-cg>eps,1-cv-cg,eps) * ( log(if(1-cv-cg>eps,1-cv-cg,eps)) - log(if(1-cnu0-cg0>eps,1-cnu0-cg0,eps)) ) )'
      derivative_order = 2
    [../]
  
    # === 气泡自由能 f^b(cv,cg)（vdW/对数型，含 Va 与 b） ===
    [./f_b]
      type = DerivativeParsedMaterial
      property_name = f_b
      coupled_variables = 'cv cg'
      material_property_names = 'kT()'
      constant_names = 'Vi b kp nQ eps'
      constant_expressions = '${Va} ${b} ${kp} 1.0 1.0e-12'
      # 常见实现：cg * [ -ln( nQ * (Va/cg - b) ) - 1 ] ；与文献一致的结构
      expression = '(kT/Vi) * if(cg>eps,cg,eps) * ( -log(nQ*(Vi/if(cg>eps,cg,eps) - b)) - 1.0 ) + 0.5*kp*(1 - cv - cg)^2'
      derivative_order = 2
    [../]
  
    # === 多相势 f^{poly}(eta,phi_i)（无 A/B/C；固定系数的多项式势） ===
    [./f_poly]
      type = DerivativeParsedMaterial
      property_name = f_poly
      coupled_variables = 'eta phi0 phi1 phi2'
      # 采用标准多项式：每个序参量的双/四项式 + 取向间二次耦合；系数固定为 1
      expression = '(phi0^4+phi1^4+phi2^4)/4 - 0.5*(phi0^2+phi1^2+phi2^2) + (eta^4/4 - 0.5*eta^2) + (phi0^2*phi1^2 + phi0^2*phi2^2 + phi1^2*phi2^2) + 0.25'
      derivative_order = 2
    [../]
  
    # === 化学耦合 f^{chem}(cv; phi_i, eta)（按需，可为0）===
    [./f_chem]
      type = DerivativeParsedMaterial
      property_name = f_chem
      coupled_variables = 'cv phi0 phi1 phi2 eta'
      constant_names = 'S_l'
      constant_expressions = '${S_l}'
      # 若不需要这项，可令 S_l=0
      expression = 'S_l*cv*(1 - (phi0^2+phi1^2+phi2^2+eta^2))'
      derivative_order = 2
    [../]
  
    # === 形函数 h(eta) & j(eta) ===
    [./h_eta]
      type = ParsedMaterial
      property_name = h_eta
      coupled_variables = 'eta'
      expression = '(1-eta)^2'
    [../]
    [./j_eta]
      type = ParsedMaterial
      property_name = j_eta
      coupled_variables = 'eta'
      expression = 'eta^2'
    [../]
  
    # === 组装 f^0，并作为相场 free_energy 的入口（与 Modules 中一致） ===
    [./f0]
      type = DerivativeParsedMaterial
      property_name = f0
      coupled_variables = 'cv cg cs eta phi0 phi1 phi2'
      material_property_names = 'h_eta(eta) j_eta(eta) f_m(cv,cg,cs) f_b(cv,cg) f_poly(eta,phi0,phi1,phi2) f_chem(cv,phi0,phi1,phi2,eta)'
      constant_names = 'omega'
      constant_expressions = '${omega}'
      expression = 'h_eta*f_m + j_eta*f_b + omega*f_poly + f_chem'
      derivative_order = 2
    [../]
  
    # === 晶界增强的成核概率（只在晶界、未成核区更容易成核） ===
    [./probability_bubble]
      type = ParsedMaterial
      property_name = P_bubble
      coupled_variables = 'bnds bubble_mask_combined'
      constant_names = 'base_prob min_mask'
      constant_expressions = '1e-3 0.05'
      expression = 'max((0.7 - bnds)*base_prob, 0) * max(1 - bubble_mask_combined, min_mask)'
      outputs = exodus
    [../]
  
    # === 离散成核 —— 对 cg（供气）与 eta（把气泡相拉到1） ===
    [./nucleation_cg]
      type = DiscreteNucleation
      property_name = Fn_cg
      op_names  = cg
      op_values = 0.95
      penalty   = 1e6
      penalty_mode = MIN
      map = map_cg
      outputs = exodus
    [../]
    [./nucleation_eta]
      type = DiscreteNucleation
      property_name = Fn_eta
      op_names  = eta
      op_values = 1.0
      penalty   = 1e6
      penalty_mode = MIN
      map = map_eta
      outputs = exodus
    [../]
  
    # === 供 PhaseField 使用的最终自由能 ===
    [./f_sum]
      type = DerivativeSumMaterial
    property_name = f_bubble_total
    coupled_variables = 'cv cg eta phi0 phi1 phi2'
      sum_materials = 'f0 Fn_cg Fn_eta'
    [../]
[]
[UserObjects]
    # 晶界 Voronoi 生成（与你原来一致）
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 2
    int_width = 0.2e3
  [../]
  
  [./grain_tracker]
    type = GrainTracker
  [../]
  
    # ===== 气泡成核：对 cg（供气） =====
  [./inserter_cg]
    type = DiscreteNucleationInserter
    hold_time  = 0.001
    probability = P_bubble     # 来自 Materials/probability_bubble
    radius     = 0.2e3         # 成核半径（nm 数值；按需要调整）
  [../]
  
  [./map_cg]
    type = DiscreteNucleationMap
    periodic = cg              # 与被成核的变量一致
    inserter = inserter_cg
  [../]
  
    # ===== 气泡成核：对 eta（把气泡序参量拉到1）=====
  [./inserter_eta]
    type = DiscreteNucleationInserter
    hold_time  = 0.001
    probability = P_bubble     # 同一概率场：晶界增强 & 禁止已成核重复成核
    radius     = 0.2e3
  [../]
  
  [./map_eta]
    type = DiscreteNucleationMap
    periodic = eta
    inserter = inserter_eta
  [../]
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

  nl_rel_tol = 1e-12 # 非线性求解的相对容差
  nl_abs_tol = 1e-7 # 非线性求解的绝对容差
  l_tol = 1e-12  # 线性求解的容差
  l_abs_tol = 1e-8 # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 100

  dt = 0.0001
  # [./Adaptivity]
  #   max_h_level = 2
  #   initial_adaptivity = 1
  #   refine_fraction = 0.9
  #   coarsen_fraction = 0.1
  # [../]
[]

[Outputs]
  exodus = true
  file_base = 'results/FB1'
[]