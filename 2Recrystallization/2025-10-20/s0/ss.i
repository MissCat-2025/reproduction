# conda activate moose && mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i test1.i
# mpirun -n 9 /home/cy/projects/xxx/xxx-opt -i UN_fisson_bubble.i
# 
# ========== 动力学方程实现说明 ==========
# 本输入文件实现以下动力学方程:
#
# 1. 空位演化 (cv):
#    ∂cν/∂t = ∇·(Mν∇(δF/δcν)) + Pν(r,t) + ξν(r,t)
#    实现位置: [Modules]/Conserved/cv (扩散项) + [Kernels]/cv_source_Pnu (源项)
#
# 2. 气体原子演化 (cg):
#    ∂cg/∂t = ∇·(Mg∇(δF/δcg)) + Pg(r,t) + ξg(r,t) - Rg(r,t)
#    实现位置: [Modules]/Conserved/cg (扩散项) + [Kernels]/cg_source_Pg/Rg (源项)
#
# 3. 固溶原子 (cs):
#    ∂cs/∂t = ∇·(Ms∇(δF/δcs))  [当前冻结: ∂cs/∂t=0]
#    实现位置: [Kernels]/cs_null (冻结)
#
# 4. 气泡序参量 (η):
#    ∂η/∂t = -Lη(δF/δη)
#    实现位置: [Modules]/Nonconserved/eta
#
# 5. 晶粒序参量 (φi):
#    ∂φi/∂t = -Lφ(δF/δφi)  [当前冻结: ∂φi/∂t=0]
#    实现位置: [Kernels]/phi*_null (冻结)
#
# 源项定义:
#   Pν(r,t) = VG (if η<0.8 and R1≤Pcasc, else 0) - 级联空位生成
#   Pg(r,t) = 2(1-η)²ΛΩfr·R3 - 裂变气体生成
#   Rg(r,t) = η²·bi·cg - 气体复合损失
#
# 初始条件:
#   cg(t=0) = 0
#   cν(t=0) = exp(-Efv/kBT)
#   cs(t=0) = 局部U原子浓度
#   η(t=0) = 0
#   φi(t=0) = Voronoi多晶结构
#本文引入了一个额外的可视化参数Ф，通过分析不同相内非守恒相场变量的值及其在扩散界面的分布特征，可视化参数可以定义为以下表达式：
# Φ=2[∑_(i=1)^r φ_i^2 ]+[∑_(j=r+1)^p φ_j^2]-η^2
# ==========================================
#辐照气泡的成核和生长
#   气泡核半径为3.2dx，过饱和空位浓度为0.2。成核概率:
#   J^*=κ_1*e^(κ_2/Δc_ν )，
#   式中，κ1、κ2为常数，值都为1e-5；Δcv为空位过饱和浓度，值为0.2。
#   成核概率引入了具有过饱和空位的圆形区域，该区域与气体原子结合形成气泡，这就是模拟中使用的成核方法。

#   一个铀原子裂变过程中气体原子（Xe和Kr）的产额为0.251（Olander，1976），



# 参数设置
# --- 基本常数 ---
JtoeV = 6.24150974e18           # J → eV 转换
length_scale = 1e9              # nm 为长度尺度
time_scale = 1.0e-6                  # ns 为时间尺度
T = 1000                        # K 温度

# --- 自由能势垒 / 势阱参数 ---
omega_SI = 1.54e8               # J/m³ (Potential height)
kp_SI = 6.4e11                  # J/m³ (Prefactor)
Va_SI = 0.04092e-27             # m³/atom (0.04092 nm³)
b_SI = 0.085e-27                # m³/atom (van der Waals constant)

# —— 化学耦合项：S_l（先用 0.05*omega_SI，数值可扫参） —— （当前未使用，已注释）
# S_l_SI = ${fparse 0.05*omega_SI}                     # J/m^3
# S_l     = ${fparse S_l_SI * JtoeV / length_scale^3}  # eV/nm^3

# --- 梯度能系数 ---
kappa_v_SI = 3.38e-8            # J/m (vacancy gradient)
kappa_eta_SI = 1.67e-9          # J/m (order parameter gradient)

# --- 迁移率与界面迁移率 ---
Mv_SI = 2.69e-26                # m⁵/(J·s) vacancy mobility
Mg_SI = 1.33e-28                # m⁵/(J·s) gas atom mobility
L_SI = 4.33e-12                  # m³/(J·s) interface mobility

# --- UN 物理参数 ---
Ev_f_eV = 7.85                  # Vacancy formation energy (eV)
# Ev_m_eV = 2.2                   # Vacancy migration energy (eV) - 未使用
# Eg_m_eV = 1.79                  # Gas atom migration energy (eV) - 未使用
# Dv_m_SI = 4.0e-7                # Vacancy diffusivity (m²/s)
# Dg_m_SI = 4.05e-7               # Gas diffusivity (m²/s)
# gamma_gb_SI = 0.70              # J/m² Grain boundary energy
# gamma_s_SI = 1.62               # J/m² Surface energy

# --- 裂变与级联参数 ---
fr_SI = 1e19                     # 裂变速率1e13 fission rate (1/cm³/s) = 1e19 (1/m³/s)
Lambda = 0.25                    # 气体原子产生份额常数 (无量纲)
b0_SI = 1e-24                    # 分辨率常数 resolution constant (cm³) = 1×10⁻¹⁸ cm³=1e-24m³
# Pcasc = 1e-1                     # 级联发生概率 (1/ns) 假设值，需根据实际调整
VG = 0.1                        # 级联空位浓度最大增加 (无量纲) 假设值

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

# --- 裂变与级联参数换算 ---
# kB_eV = 8.617333262e-5           # Boltzmann constant (eV/K)
fr = ${fparse fr_SI/length_scale^3/time_scale}                   # nm^-3·ns^-1 
b0 = ${fparse b0_SI*length_scale^3}                    # nm³
bi = ${fparse b0 * fr}               # 分辨率 bi=b0*fr (nm^3)
# cv0 = ${fparse exp(-Ev_f_eV/(kB_eV*T))}                    # 初始空位浓度 c0ν=exp(-Efv/kBT)
# cv0_hint = 1.0e-6                                           # 兼容保留：旧参考量级（未直接使用）

# --- 注释 ---
# Ev_f_eV, Ev_m_eV, Eg_m_eV 保留在 eV 单位用于计算扩散激活能。
# 如果要耦合温度依赖扩散率，可按 D = D0 * exp(-Em/kT) 转换为有效迁移率。

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 50
  ymax = 50
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
  # ===== NUCLEATION SYSTEM — OVERVIEW =====
  # 目的: 通过“概率 → 插入器 → 映射 → 惩罚势(材料) → 汇总自由能”的链路，
  #       实现对气泡变量的离散成核。
  # 数据流一览:
  #   1) 概率:  基于过饱和度 ln(cv/cv_eq) → P_bubble（对 log 加 eps 防护）
  #   2) 插入:  P_bubble → DiscreteNucleationInserter(inserter_*)
  #   3) 映射:  inserter_* → DiscreteNucleationMap(map_*)
  #   4) 惩罚:  map_* → Fn_* (Materials/DiscreteNucleation)
  #   5) 汇总:  f_bubble_total = f0 + Fn_cg + Fn_eta
  # 晶界掩膜（用于后续成核概率/界面能等）
  [./bnds]
    family = LAGRANGE
    order = FIRST
  [../]
  
  # 源项辅助变量（用于存储 P_ν, P_g, R_g 并耦合到 Kernels）
  [./P_nu_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./P_g_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./R_g_aux]
    family = MONOMIAL
    order = CONSTANT
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
  
  # 5) 源项计算：P_ν, P_g, R_g
  [./P_nu_calc]
    type = MaterialRealAux
    variable = P_nu_aux
    property = P_nu
    execute_on = 'timestep_end'
  [../]
  [./P_g_calc]
    type = MaterialRealAux
    variable = P_g_aux
    property = P_g
    execute_on = 'timestep_end'
  [../]
  [./R_g_calc]
    type = MaterialRealAux
    variable = R_g_aux
    property = R_g
    execute_on = 'timestep_end'
  [../]
[]
[Modules]
  [./PhaseField]
    # ====== 动力学方程自动实现 ======
    # Conserved: 实现 Cahn-Hilliard 方程 ∂c/∂t = ∇·(M∇(δF/δc))
    # Nonconserved: 实现 Allen-Cahn 方程 ∂η/∂t = -L(δF/δη)
    # 源项 P_ν, P_g, R_g 需在 [Kernels] 中通过 CoupledForce/BodyForce 添加
    [./Conserved]
      # 空位演化: ∂cν/∂t = ∇·(Mν∇(δF/δcν)) + Pν(基体级联生成)
      [./cv]
        free_energy = f_bubble_total
        mobility    = Mv
        kappa       = kappa_cv
        solve_type  = REVERSE_SPLIT
        args = 'cg eta phi0 phi1 phi2'
      [../]
      # 气体原子演化: ∂cg/∂t = ∇·(Mg∇(δF/δcg)) + Pg(裂变生成) - Rg(复合损失)
      [./cg]
        free_energy = f_bubble_total
        mobility    = Mg
        kappa       = kappa_cg
        solve_type  = REVERSE_SPLIT
        args = 'cv eta phi0 phi1 phi2'
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
      # 气泡序参量演化: ∂η/∂t = -Lη(δF/δη)
      [./eta]
        free_energy = f_bubble_total
        mobility    = L_eta
        kappa       = kappa_eta
        args = 'cv cg phi0 phi1 phi2'
      [../]
      # 晶粒序参量 φi 已通过 NullKernel 冻结（∂φi/∂t = 0），不演化
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
  
    # 守恒场与气泡序参量的初始值
    # 初始条件实现动力学方程的初始态：cg=0, c0ν=exp(-Efv/kBT), cs为局部U原子浓度, η=0
  [./cv_init]
    type = ConstantIC
    variable = cv
    value = 1e-6   # 从材料属性读取：cv_eq = exp(-Ev/kT)
  [../]
  [./cg_init]
    type = ConstantIC
    variable = cg
    value = 1.0e-20     # 避免 log(0)，物理上近似 0
  [../]
  [./cs_init]
    type = ConstantIC
    variable = cs
    value = 1         # cs_eq：单相 UN 的基体占位
  [../]
  [./eta_init]
    type = ConstantIC
    variable = eta
    value = 1.0e-20         # 初始无气泡
  [../]
  
    # （可选）给 cv/cg/eta 加极小随机扰动以避免完全对称：
    # [./cv_jitter] type=RandomIC variable=cv min=0.0 max=1e-8 [../]
[]
[Kernels]
    # ====== 源项内核：将 P_ν, P_g, R_g 加入动力学方程 ======
    
    # === 空位方程源项：∂cν/∂t = ... + P_ν ===
  [./cv_source_Pnu]
    type = CoupledForce
    variable = cv
    v = P_nu_aux          # 耦合 AuxVariable
    coef = 1.0
  [../]
  
    # === 气体方程源项：∂cg/∂t = ... + P_g - R_g ===
  [./cg_source_Pg]
    type = CoupledForce
    variable = cg
    v = P_g_aux           # 裂变气体生成
    coef = 1.0
  [../]
  [./cg_source_Rg]
    type = CoupledForce
    variable = cg
    v = R_g_aux           # 气体复合损失
    coef = -1.0           # 负号表示损失
  [../]
  
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
  
    # 冻结 cs：∂cs/∂t = 0 (若需要演化cs，注释掉此块并启用Conserved/cs)
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
  
    # 平衡空位浓度 cv_eq = exp( -Ev_f_eV / kT )
    [./cv_eq_mat]
      type = ParsedMaterial
      property_name = cv_eq
      material_property_names = 'kT()'
      constant_names = 'Ev_f_eV'
      constant_expressions = '${Ev_f_eV}'
      expression = 'exp(-Ev_f_eV / kT)'
    [../]
  
    # === PhaseField 动力学常数，供 Modules/PhaseField 读取 ===
    [./pfmobility]
      type = GenericConstantMaterial
    prop_names  = 'Mv Mg L_eta kappa_cv kappa_cg kappa_eta P_bubble'
    prop_values = '${Mv} ${Mg} ${L} ${kappa_v} ${kappa_v} ${kappa_eta} 1e-5'
      # 注：若暂不演化 cs，可将 Ms=0 且 kappa_cs=0；反之请给出实际值
    [../]
  
    # === 基体自由能 f^m(cv,cg,cs) ===
    # 采用理想混合熵 + 体积分数约束（kp）常见写法；与你 Word 文档思路一致
    [./f_m]
      type = DerivativeParsedMaterial
      property_name = f_m
      coupled_variables = 'cv cg cs'
      material_property_names = 'kT()'
      constant_names = 'kp Valpha cnu0 cg0'
      constant_expressions = '${kp} ${Va} 1.0e-6 1.0e-8'
      # kT*[ cv ln(cv) + cg ln(cg) + (1-cv-cg) ln(1-cv-cg) ] + (kp/2)*(1 - cv - cg)^2
      expression = '(kT/Valpha) * ( cv * ( log(cv) - log(cnu0) ) + cg * ( log(cg) - log(cg0) ) + (1-cv-cg) * ( log(1-cv-cg) - log(1-cnu0-cg0) ) )'
      derivative_order = 2
    [../]
  
    # === 气泡自由能 f^b(cv,cg)（vdW/对数型，含 Va 与 b） ===
    [./f_b]
      type = DerivativeParsedMaterial
      property_name = f_b
      coupled_variables = 'cv cg'
      material_property_names = 'kT()'
      constant_names = 'Vi b kp nQ'
      constant_expressions = '${Va} ${b} ${kp} 1.0'
      # 常见实现：cg * [ -ln( nQ * (Va/cg - b) ) - 1 ] ；与文献一致的结构
      expression = '(kT/Vi) * cg * ( -log(nQ*(Vi/cg - b)) - 1.0 ) + 0.5*kp*(1 - cv - cg)^2'
      derivative_order = 2
    [../]
  
    # === 多相势 f^{poly}(eta,phi_i)（无 A/B/C；固定系数的多项式势） ===
    [./f_poly]
      type = DerivativeParsedMaterial
      property_name = f_poly
      coupled_variables = 'eta phi0 phi1 phi2'
      constant_names = 'a_gb a_s'  # 晶界能，表面能
      constant_expressions = '1.5 1.8'
      # 采用标准多项式：每个序参量的双/四项式 + 取向间二次耦合 + 晶界/表面能耦合
      expression = '(phi0^4+phi1^4+phi2^4)/4 - 0.5*(phi0^2+phi1^2+phi2^2) + (eta^4/4 - 0.5*eta^2) + a_gb*(phi0^2*phi1^2 + phi0^2*phi2^2 + phi1^2*phi2^2) + a_s*eta^2*(phi0^2+phi1^2+phi2^2) + 0.25'
      derivative_order = 2
    [../]
      
  
    # === 化学耦合 f^{chem}(cv; phi_i, eta)（按需，可为0）===
    [./f_chem]
      type = DerivativeParsedMaterial
      property_name = f_chem
      coupled_variables = 'cs phi0 phi1 phi2 eta'
      constant_names = 'A c_alpha'
      constant_expressions = '1.0 1.0'
      # 化学耦合项: f_chem = (A/2)*(cs - c_alpha)^2
      expression = '0.5*A*(cs - c_alpha)^2'
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
  
    # ====== 源项材料：实现动力学方程中的 P_ν, P_g, R_g ======
    
    # === 气体生成源项 P_g(r,t) ===
    # 公式: P_g(r,t) = 2(1-η)²ΛΩfr R3 (基体中裂变气体生成)
    [./source_Pg]
      type = ParsedMaterial
      property_name = P_g
      coupled_variables = 'eta'
      constant_names = 'Lambda Omega fr'
      constant_expressions = '${Lambda} ${Va} ${fr}'
      # R3 为 0-1 随机数，此处简化为平均值 0.5，或用 ConditionalFunctionAux 实现随机
      expression = '2*(1-eta)^2 * Lambda * Omega * fr * 0.1'
      outputs = exodus
    [../]
    
    # === 空位级联生成源项 P_ν(r,t) ===
    # 公式: Pν = VG (if η<0.8 and R1≤Pcasc), else 0
    # 简化实现：在基体(η<0.8)中按概率生成，此处用平均值近似
    [./source_Pnu]
      type = ParsedMaterial
      property_name = P_nu
      coupled_variables = 'eta'
      constant_names = 'VG thres'
      constant_expressions = '${VG} 0.8'
      # 简化: 基体中平均生成速率 = VG*Pcasc, 气泡中为0
      expression = 'if(eta<thres, VG, 0)'
      outputs = exodus
    [../]
    
    # === 气体复合损失源项 R_g(r,t) ===
    # 公式: Rg = η² bi cg (气泡中气体的分辨率损失)
    [./source_Rg]
      type = ParsedMaterial
      property_name = R_g
      coupled_variables = 'eta cg'
      constant_names = 'bi'
      constant_expressions = '${bi}'
      expression = 'eta^2 * bi * cg'
      outputs = exodus
    [../]

    # === 成核概率场 P_bubble（基于过饱和度：ln(cv/cv_eq)；对 log(0) 用 eps 防护） ===
    # 用途: 为 DiscreteNucleationInserter 提供空间分布的成核概率。
    # 形式: P_bubble = max( k1 * ln( cv / cv_eq ), 0 )
    # 说明: 使用 if(cv>eps,cv,eps) 与 if(cv_eq>eps,cv_eq,eps) 防止 log(0)。

    # === 离散成核“惩罚势”(Materials) —— 将 map_* 注入自由能 ===
    # Fn_cg: 促使 cg 在成核点趋近 op_values；Fn_eta: 促使 eta 在成核点拉到1。
    # 注意: 两者均以大 penalty 注入，以最小化模式(MIN)生效。
    [./nucleation_cg]
      type = DiscreteNucleation
      property_name = Fn_cg
      op_names  = cg
      op_values = 0.95
      penalty   = 1e10
      penalty_mode = MIN
      map = map_cg
      outputs = exodus
    [../]
    [./nucleation_eta]
      type = DiscreteNucleation
      property_name = Fn_eta
      op_names  = eta
      op_values = 1.0
      penalty   = 1e10
      penalty_mode = MIN
      map = map_eta
      outputs = exodus
    [../]
  
    # === 供 PhaseField 使用的最终自由能 ===
    [./f_sum]
      type = DerivativeSumMaterial
      property_name = f_bubble_total
      coupled_variables = 'cv cg cs eta phi0 phi1 phi2'
      # 将基体/气泡/多相等基础能量 f0，与两类离散成核惩罚 Fn_* 一并汇总
      sum_materials = 'f0 Fn_cg Fn_eta'
    [../]
[]
[UserObjects]
    # 晶界 Voronoi 生成（与你原来一致）
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 2
    int_width = 1.6
  [../]
  
  [./grain_tracker]
    type = GrainTracker
  [../]
  
    # ===== 成核驱动（UserObjects）—— 插入器 Inserter 与 映射 Map =====
    # 作用顺序: probability(P_bubble) → inserter_* → map_* → Fn_*
    # 建议: radius 控制成核“核”的空间尺度；hold_time 控制开始插入的时间点。
    # --- 气泡成核：对 cg（供气） ---
  [./inserter_cg]
    type = DiscreteNucleationInserter
    hold_time  = 0.01
    probability = P_bubble     # 概率: 基于 ln(cv/cv_eq)
    radius     = 1.6           # 成核半径 = 3.2*dx，dx=0.5 → 1.6 nm
  [../]
  
  [./map_cg]
    type = DiscreteNucleationMap
    periodic = cg              # 与被成核的变量一致
    inserter = inserter_cg
  [../]
  
    # --- 气泡成核：对 eta（把气泡序参量拉到1） ---
  [./inserter_eta]
    type = DiscreteNucleationInserter
    hold_time  = 0.01
    probability = P_bubble     # 概率: 基于 ln(cv/cv_eq)
    radius     = 1.0
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
  # solve_type = 'NEWTON'
  # petsc_options = '-ksp_type=preonly -pc_type=lu -pc_factor_mat_solver_type=mumps'
  # line_search = bt
  # scheme = bdf2
  # solve_type = 'NEWTON'
  # petsc_options = '-ksp_type=gmres -pc_type=lu'
  line_search = bt
    #Preconditioned JFNK (default)
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = 'hypre boomeramg'
  # nl_max_its = 20

  nl_rel_tol = 1e-7 # 非线性求解的相对容差
  nl_abs_tol = 1e-7 # 非线性求解的绝对容差
  l_tol = 1e-7  # 线性求解的容差
  l_abs_tol = 1e-7 # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 5000

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