# UO2蠕变模型参数 - 基于公式(25)-(27)
# 总蠕变率: dot_epsilon_cr = A1*sigma_e^n1*exp(-Q1/RT) + A2*sigma_e^n2*exp(-Q2/RT) + A3*sigma_e^n3*exp(-Q3/RT)

# 基本材料参数
grain_size = 10                    # 晶粒尺寸 Gr (μm)
fission_rate = 1.2e19             # 裂变率密度 f_dot (fissions/m^3-s)
theoretical_density = 95          # 理论密度百分比 D (%)
gas_constant = 8.314              # 气体常数 R (J/mol-K)

# 原始系数
a1_base = 0.3919                  # 基础系数
a2_fission = 1.31e-19            # 裂变相关系数
a3_density = 87.7                # 密度修正系数
a5_thermal = 2.0391e-25          # 热蠕变第二项系数
a6_density = 90.5                # 密度修正系数2
a8_irradiation = 3.7226e-35      # 辐照蠕变系数

# 激活能 (J/mol) - 需要根据氧化学计量比计算，这里使用典型值
Q1_activation = 376591.0          # 第一项激活能 (示例值)
Q2_activation = 552334.0          # 第二项激活能 (示例值)
Q3_activation = 21759.0           # 第三项激活能

# 转换后的蠕变参数
# 第一项热蠕变: (0.3919 + 1.31e-19*f_dot) / [(D-87.7)*Gr^2] * sigma_e * exp(-Q1/RT)
A1 = ${fparse (a1_base + a2_fission * fission_rate) / ((theoretical_density - a3_density) * grain_size * grain_size)}
n1 = 1.0
Q1 = ${Q1_activation}

# 第二项热蠕变: (2.0391e-25) / (D-90.5) * sigma_e^4.5 * exp(-Q2/RT)
A2 = ${fparse a5_thermal / (theoretical_density - a6_density)}
n2 = 4.5
Q2 = ${Q2_activation}

# 第三项辐照蠕变: (3.7226e-35*f_dot) * sigma_e * exp(-Q3/RT)
A3 = ${fparse a8_irradiation * fission_rate}
n3 = 1.0
Q3 = ${Q3_activation}

# # 验证计算的参数值（可以在输出中查看）
# A1_value = ${A1}                  # 第一项系数
# A2_value = ${A2}                  # 第二项系数
# A3_value = ${A3}                  # 第三项系数

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
[]

[Variables]
  [temp]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1000.0
  []
[]
[AuxVariables]
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_creep]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [vonMisesStress]
    type = RankTwoScalarAux
    variable = vonMises
    rank_two_tensor = stress
    execute_on = 'TIMESTEP_END'
    scalar_type = VonMisesStress
    # 不需要 index_i 和 index_j，因为我们使用 VonMisesStress 标量类型
  []
  [./creep_aux]
    type = MaterialRealAux
    property = effective_creep_strain
    variable = effective_creep
  [../]
[]
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_yy creep_strain_xx creep_strain_yy creep_strain_zz elastic_strain_yy'
  []
[]

[Functions]
  [top_pull]
    type = PiecewiseLinear
    x = '0 1'
    y = '1 1'
  []
[]

[Kernels]
  [heat]
    type = Diffusion
    variable = temp
  []
  [heat_ie]
    type = TimeDerivative
    variable = temp
  []
[]

[BCs]
  [u_top_pull]
    type = Pressure
    variable = disp_y
    boundary = top
    factor = -10.0e6
    function = top_pull
  []
  [u_bottom_fix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [u_yz_fix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [u_xy_fix]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []
  [temp_fix]
    type = DirichletBC
    variable = temp
    boundary = 'bottom top'
    value = 1000.0
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'power_law_creep'
    tangent_operator = elastic
    # output_properties = 'strain_increment'
    # outputs = exodus
  []
  [power_law_creep]
    type = PowerLawCreepStressUpdate
    coefficient = 1.0e-15
    n_exponent = 4
    activation_energy = 3.0e5
    temperature = temp
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'

  line_search = 'none'

  l_max_its = 20
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-5
  start_time = 0.0
  end_time = 1.0
  num_steps = 10
  dt = 0.1
[]

[Outputs]
  exodus = true
[]
