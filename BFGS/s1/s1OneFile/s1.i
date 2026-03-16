#This input uses PhaseField-Nonconserved Action to add phase field fracture bulk rate kernels
#弹性
E = 200e9
Nu = 0.33
#断裂
c0 = 3.1415
Ksi = 2
Gc = 3
Sigma = 6e7
b = 5e-5

m = 2
a1 = '${fparse 2*Ksi/c0/b*E*Nu/Sigma/Sigma}'
a2 = -0.5
a3 = 0

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 100
    xmax = 5e-3
    ymax = 2.5e-3
  []
  [./noncrack]
    type = BoundingBoxNodeSetGenerator
    new_boundary = noncrack
    bottom_left = '2.5e-3 0 0'
    top_right = '5e-3 0 0'
    input = gen
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules]
  [./PhaseField]
    [./Nonconserved]
      [./d]
        free_energy = F
        kappa = kappa_op
        mobility = b
      [../]
    [../]
  [../]
[]
[Physics]
  [./SolidMechanics]
    [./QuasiStatic]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy'
        save_in = 'resid_x resid_y'
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./resid_x]
  [../]
  [./resid_y]
  [../]
[]

[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = d
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = d
  [../]
[]

[BCs]
  [./ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 't'
  [../]
  [./yfix]
    type = DirichletBC
    variable = disp_y
    boundary = noncrack
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    variable = disp_x
    boundary = top
    value = 0
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'Gc b visco c0 gc_prop l'
    prop_values = '${Gc} ${b} 1e-1 ${c0} ${Gc} ${b} '
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'Gc visco'
    property_name = L
    expression = '1.0/(Gc * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'Gc b c0'
    property_name = kappa_op
    expression = 'Gc * b * 2 / c0'
  [../]

  [elast_tensor2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${Nu}
  []
  [./elastic]
    type = ComputeLinearElasticPFFractureStress
    c = d
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'fracture_energy'
    barrier_energy = 'barrier'
    decomposition_type = strain_spectral
  [../]
    [./degradation]
    type = DerivativeParsedMaterial
    property_name = degradation
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a3*d^2))*(1-eta)+eta
    coupled_variables = 'd'
    constant_names = 'p a1 a2 a3 eta '
    constant_expressions = '${m} ${a1} ${a2} ${a3} 1e-8'
    derivative_order = 2
  [../]
  [./CrackGeometricFunction]
    type = DerivativeParsedMaterial
    property_name = alpha
    coupled_variables = 'd'
    expression = 'Ksi*d-d*d*(1-Ksi)'
    constant_names       = 'Ksi'
    constant_expressions = '${Ksi}'
    derivative_order = 2
  [../]
  [./fracture_energy]
    type = DerivativeParsedMaterial
    property_name = fracture_energy
    coupled_variables = 'd'
    material_property_names = 'Gc b alpha c0'
    expression = 'Gc*alpha/(c0*b)'
    derivative_order = 2
    output_properties = 'fracture_energy'
    outputs = exodus
  [../]
  [./fracture_driving_energy]
    type = DerivativeParsedMaterial
    coupled_variables = d
        # expression = 'elastic_energy + fracture_energy + 1e1 * (max(0, -d)^2 + max(0, d-1)^2)'
    expression = 'elastic_energy + fracture_energy'
    material_property_names = 'elastic_energy(d) fracture_energy(d) '
    # expression = 'elastic_energy + fracture_energy + 1e10 * (max(0, -d)^2 + max(0, d-1)^2)'
    
    derivative_order = 2
    property_name = F
  [../]
  [./barrier_energy]
    type = ParsedMaterial
    property_name = barrier
    expression = 'Sigma*Sigma/E'
    constant_names       = 'Sigma E'
    constant_expressions = '${Sigma} ${E}'
  [../]
[]

[Postprocessors]
  [./resid_x]
    type = NodalSum
    variable = resid_x
    boundary = 2
  [../]
  [./resid_y]
    type = NodalSum
    variable = resid_y
    boundary = 2
  [../]
[]

# [Preconditioning]
#   [./smp]
#     type = SMP
#     full = true
#   [../]
# []

[Executioner]
  type = Transient

  # solve_type = PJFNK
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm      31                  preonly       lu           1'
  solve_type = NEWTON
  petsc_options_iname = '-pc_type   -snes_type        -snes_qn_type   -snes_qn_scale_type' 
  petsc_options_value = 'lu         qn               lbfgs           jacobian'  
  nl_rel_tol = 1e-4
  l_max_its = 10
  nl_max_its = 20

  dt = 1e-6
  dtmin = 1e-6
  num_steps = 50
[]

[Outputs]
  exodus = true
[]
