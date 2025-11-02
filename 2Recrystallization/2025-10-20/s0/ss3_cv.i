# conda activate moose && mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i ss3_cv.i
# ========== 晶粒序参量演化 ==========
# 本输入文件实现晶粒序参量演化:
#
# 晶粒序参量 (φi):
#    ∂φi/∂t = -Lφ(δF/δφi)
#    实现位置: [Modules]/Nonconserved/phi*


# 参数设置
# --- 基本常数 ---
JtoeV = 6.24150974e18           # J → eV 转换
length_scale = 1e9              # 1m=10^9nm 为长度尺度
time_scale = 1.0               # ns 为时间尺度
T = 1000                       # K 温度
k = 8.617333262e-5             # eV/K 玻尔兹曼常数

# --- 自由能势垒 / 势阱参数 ---
# omega_SI = 1.54e8              # J/m³ (Potential height)
kp_SI = 6.4e11                 # J/m³ (Prefactor)
Va_SI = 0.04092e-27            # m³/atom (0.04092 nm³)
b_SI = 0.085e-27               # m³/atom (van der Waals constant)

# --- 梯度能系数 ---
kappa_v_SI = 3.38e-8            # J/m (vacancy gradient)
# kappa_phi_SI = 1.67e-9          # J/m (order parameter gradient)
# kappa_eta_SI = 1.67e-9          # J/m (bubble order parameter gradient)
# --- 迁移率与界面迁移率 ---
Mv_SI = 2.69e-26                # m⁵/(J·s) vacancy mobility
# L_SI = 4.33e-6                 # m³/(J·s) interface mobility
# --- UN 物理参数 ---
Ev_f_eV = 7.85                  # Vacancy formation energy (eV)
cv0 = ${fparse exp(-Ev_f_eV/(k*T))}                    # 初始空位浓度 c0ν=exp(-Efv/kBT)
cg =  1e-10
cg0 = 0.251 #the yield of gas atoms (Xe and Kr) in the fission process of one uranium atom is 0.251 (Olander, 1976),
nQ = 1.0e30 #量子浓度???
# P_bubble = 0.0005 #气泡成核概率

VG = 1.0e-2                       # 级联空位浓度最大增加 (无量纲) 假设值 1e-8
eta = 0.0001
# --- 单位换算到 MOOSE 内部（eV、nm尺度） ---
# omega = ${fparse omega_SI * JtoeV / length_scale^3}
kp = ${fparse kp_SI * JtoeV / length_scale^3} # J/m³ (Prefactor)
Va = ${fparse Va_SI * length_scale^3} # m³/atom (0.04092 nm³)
b = ${fparse b_SI * length_scale^3} # m³/atom (van der Waals constant)
kappa_v = ${fparse kappa_v_SI / length_scale * JtoeV} # J/m (vacancy gradient)
# kappa_phi = ${fparse kappa_phi_SI / length_scale * JtoeV}
# kappa_eta = ${fparse kappa_eta_SI / length_scale * JtoeV}
# L = ${fparse L_SI * length_scale^3 / (JtoeV * time_scale)}
Mv = ${fparse Mv_SI * length_scale^5 / (JtoeV * time_scale)} # m⁵/(J·s) vacancy mobility


[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 100
  ymax = 100
  elem_type = QUAD4
[]

[Variables]
  # 守恒相场：空位浓度 cv、气体原子浓度 cg、固溶原子浓度 cs
  [./cv]
    family = LAGRANGE
    order = FIRST
  [../]
[]
[AuxVariables]
  # 源项辅助变量（用于存储 P_ν, P_g, R_g 并耦合到 Kernels）
  [./P_nu_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  # [./kkk]
  #     family = MONOMIAL
  #     order = CONSTANT
  # [../]
[]
[AuxKernels]
  # 源项计算：P_ν,
  [./P_nu_calc]
    type = MaterialRealAux
    variable = P_nu_aux
    property = P_nu
    execute_on = 'timestep_end'
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
        [../]
      []
    [../]
[]
[Kernels]
  # === 空位方程源项：∂cν/∂t = ... + P_ν ===
  [./cv_source_Pnu]
    type = CoupledForce
    variable = cv
    v = P_nu_aux          # 耦合 AuxVariable
    coef = 1.0
  [../]
    # [./phi1_null]
    #   type = NullKernel
    #   variable = kkk
    # [../]
[]
[ICs]
    # 守恒场
    # 初始条件实现动力学方程的初始态：cg=0, c0ν=exp(-Efv/kBT), cs为局部U原子浓度, η=0
    [./cv_init]
              type = ConstantIC
        variable = cv
        value = 1e-10         # 初始无气泡
      # variable = cv
      # min = ${fparse 1e-10}  
      # max = ${fparse 1e-2}
      # seed = 42
    [../]
      # [./eta_init]
      #   type = ConstantIC
      #   variable = kkk
      #   value = 1         # 初始无气泡
      # [../]
[]

[Materials]
      # === 空位级联生成源项 P_ν(r,t) ===
  # 公式: Pν = VG (if η<0.8 and R1≤Pcasc), else 0
  # 简化实现：在基体(η<0.8)中按概率生成，此处用平均值近似
  [./source_Pnu]
    type = ParsedMaterial
    property_name = P_nu
    coupled_variables = 'cv'
    constant_names = 'VG thres'
    constant_expressions = '${VG} 0.005'
    # 简化: 基体中平均生成速率 = VG*Pcasc, 气泡中为0
    expression = 'if(cv<thres, VG, 0)'
    outputs = exodus
  [../]
  # === PhaseField 动力学常数，供 Modules/PhaseField 读取 ===
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names  = 'kappa_cv Mv'
    prop_values = '${kappa_v} ${Mv}'
  [../]
  # === 基体自由能 f^m（简化版本，用于演示） ===
  [./f_m]
    type = DerivativeParsedMaterial
    property_name = f_m
    coupled_variables = 'cv'
    constant_names = 'k T Valpha cg cv0 cg0 eta'
    constant_expressions = '${k} ${T} ${Va} ${cg} ${cv0} ${cg0} ${eta}'
    # kT*[ cv ln(cv) + cg ln(cg) + (1-cv-cg) ln(1-cv-cg) ] + (kp/2)*(1 - cv - cg)^2
    expression = '(1-eta)^2*((k*T/Valpha) * ( cv * ( log(cv) - log(cv0) ) + cg * ( log(cg) - log(cg0) ) + (1-cv-cg) * ( log(1-cv-cg) - log(1-cv0-cg0) ) ))'
    derivative_order = 2
  [../]
  
  # === 气泡自由能 f^b（简化版本，用于演示） ===
  [./f_b]
    type =DerivativeParsedMaterial
    property_name = f_b
    coupled_variables = 'cv'
    constant_names = 'k T Vi b kp nQ cg eta'
    constant_expressions = '${k} ${T}  ${Va} ${b} ${kp} ${nQ} ${cg} ${eta}'
    # 常见实现：cg * [ -ln( nQ * (Va/cg - b) ) - 1 ] ；与文献一致的结构
    expression = '(eta)^2*((k*T/Vi) * cg * ( -log(nQ*(Vi/cg - b)) - 1.0 ) + 0.5*kp*(1 - cv - cg)^2)'
    derivative_order = 2
  [../]
  
  # === 化学耦合 f^{chem}（按需，可为0）===
  
  # === 组装总自由能 F = h(η)f^m + j(η)f^b + ωf^poly + f^chem ===
  [./f_total]
    type = DerivativeParsedMaterial
    property_name = f_total
    coupled_variables = 'cv'
    material_property_names = 'f_m(cv) f_b(cv)'
    expression = 'f_m+ f_b'
    derivative_order = 2
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
  solve_type = PJFNK
  # solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
petsc_options_value = 'asm      31                  preonly       lu           2'
  # solve_type = 'NEWTON'
  # petsc_options = '-ksp_type=preonly -pc_type=lu -pc_factor_mat_solver_type=mumps'
  # line_search = bt
  # scheme = bdf2

  nl_rel_tol = 1e-8 # 非线性求解的相对容差
  nl_abs_tol = 1e-8 # 非线性求解的绝对容差
  l_tol = 1e-8  # 线性求解的容差
  l_abs_tol = 1e-8 # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 5000

  dt = 0.0001
[]

[Outputs]
  exodus = true
  file_base = 'ss3_cv/ss3_cv'
[]
