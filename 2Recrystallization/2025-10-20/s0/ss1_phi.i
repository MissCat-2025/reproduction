# conda activate moose && mpirun -n 9 /home/yp/projects/reproduction/reproduction-opt -i ss1_simplified.i
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
time_scale = 1.0             # ns 为时间尺度.time_scale = 1.0时，time_scale=1ns
# T = 1000                        # K 温度

# --- 自由能势垒 / 势阱参数 ---
# omega_SI = 1.54e8               # J/m³ (Potential height)
# kp_SI = 6.4e11                  # J/m³ (Prefactor)

# --- 梯度能系数 ---
kappa_phi_SI = 1.67e-9          # J/m (order parameter gradient)

# --- 界面迁移率 ---
L_SI = 4.33e-6                 # m³/(J·s) interface mobility

# --- 单位换算到 MOOSE 内部（eV尺度） ---
# omega = ${fparse omega_SI * JtoeV / length_scale^3}
# kp = ${fparse kp_SI * JtoeV / length_scale^3}

kappa_phi = ${fparse kappa_phi_SI / length_scale * JtoeV}
L = ${fparse L_SI * length_scale^3 / (JtoeV * time_scale)}

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
  op_num = 5
  var_name_base = phi
  grain_num = 5
#   derivative_order = 2
[]

[Variables]
  # 多晶序参量
  [./PolycrystalVariables]
  [../]
[]

[AuxVariables]
  # 晶界掩膜
  [./bnds]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[AuxKernels]
  # 晶界掩膜计算
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    v = 'phi0 phi1 phi2 phi3 phi4'
    execute_on = 'initial timestep_end'
  [../]
[]

[Modules]
    [./PhaseField]
      [./Nonconserved]
        # 晶粒序参量演化: ∂φi/∂t = -Lφ(δF/δφi)
        [./phi0]
          free_energy = f_poly
          mobility    = L_phi
          kappa       = kappa_phi
          coupled_variables = 'phi1 phi2 phi3 phi4'
        [../]
        [./phi1]
          free_energy = f_poly
          mobility    = L_phi
          kappa       = kappa_phi
          coupled_variables = 'phi0 phi2 phi3 phi4'
        [../]
        [./phi2]
          free_energy = f_poly
          mobility    = L_phi
          kappa       = kappa_phi
          coupled_variables = 'phi0 phi1 phi3 phi4'
        [../]
        [./phi3]
          free_energy = f_poly
          mobility    = L_phi
          kappa       = kappa_phi
          coupled_variables = 'phi0 phi1 phi2 phi4'
        [../]
        [./phi4]
          free_energy = f_poly
          mobility    = L_phi
          kappa       = kappa_phi
          coupled_variables = 'phi0 phi1 phi2 phi3'
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
[]

[Materials]
  # === 多相势 f^{poly}(phi_i)（标准多项式势） ===
  [./f_poly]
    type = DerivativeParsedMaterial
    property_name = f_poly
    coupled_variables = 'phi0 phi1 phi2 phi3 phi4'
    constant_names = 'a_gb'  # 晶界能
    constant_expressions = '1.5'
    # 采用标准多项式：每个序参量的双/四项式 + 取向间二次耦合
    expression = '(phi0^4+phi1^4+phi2^4+phi3^4+phi4^4)/5 - 0.5*(phi0^2+phi1^2+phi2^2+phi3^2+phi4^2) + a_gb*(phi0^2*phi1^2 + phi0^2*phi2^2 + phi0^2*phi3^2 + phi0^2*phi4^2 + phi1^2*phi2^2 + phi1^2*phi3^2 + phi1^2*phi4^2 + phi2^2*phi3^2 + phi2^2*phi4^2 + phi3^2*phi4^2) + 0.25'
    derivative_order = 2
  [../]
  
  # === PhaseField 动力学常数，供 Modules/PhaseField 读取 ===
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names  = 'L_phi kappa_phi'
    prop_values = '${L} ${kappa_phi}'
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
[]

[Adaptivity]
  initial_steps = 2
  initial_marker = phase_marker
  marker = phase_marker
  max_h_level = 3
  [Markers]
    [phase_marker]
       type = ValueRangeMarker
       lower_bound = 0
       upper_bound = 0.99
      variable = bnds
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
  scheme = bdf2
  line_search = bt
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  nl_rel_tol = 1e-8 # 非线性求解的相对容差
  nl_abs_tol = 1e-8 # 非线性求解的绝对容差
  l_tol = 1e-8  # 线性求解的容差
  l_abs_tol = 1e-8 # 线性求解的绝对容差
  start_time = 0.0
  num_steps = 5000

  dt = 0.001
[]

[Outputs]
  exodus = true
  file_base = 'ss1_phi/ss1_phi'
[]
