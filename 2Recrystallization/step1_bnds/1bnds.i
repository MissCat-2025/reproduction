# This simulation predicts GB migration of a 2D copper polycrystal with 100 grains represented with 8 order parameters
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains
# mpirun -n 12 /home/yp/projects/reproduction/reproduction-opt -i 1bnds.i
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 1000
  ymax = 1000
  elem_type = QUAD
[]
[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 15 # Number of order parameters used
  var_name_base = gr # Base name of grains
[]

[Modules]
  [PhaseField]
    [GrainGrowth]
    []
  []
[]

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 50 # Number of grains
    int_width = 23
    # rand_seed = 1
  []
  [grain_tracker]
    type = GrainTracker
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[AuxVariables]
  # Dependent variables
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds #计算晶界(Grain Boundary)位置，通过计算相邻晶粒之间的界面能量梯度。
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains #为每个独立的晶粒区域分配一个唯一的ID号。使用grain_tracker用户对象识别并追踪单个晶粒。
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices #显示每个区域使用的顺序参数变量（order parameter）的索引颜色。在多相场模型中，这可以显示哪个顺序参数变量主导了特定区域。
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  []
[]

[BCs]
  # Boundary Condition block
  [Periodic]
    [All]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    []
  []
[]

[Materials]
  [CuGrGr]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    wGB = 14 # Width of the diffuse GB
    GBmob0 = 2.5e-6 #m^4(Js) for copper from schonfelder1997molecular bibtex entry
    Q = 0.23 #eV for copper from schonfelder1997molecular bibtex entry
    GBenergy = 0.708 #J/m^2 from schonfelder1997molecular bibtex entry
  []
[]

[Postprocessors]
  # Scalar postprocessors
  [dt]
    # Outputs the current time step
    type = TimestepSize
  []
[]
[VectorPostprocessors]
  [./bnds_new_export]
    type = NodalValueSampler
    variable = 'bnds'
    sort_by = 'id'
    execute_on = 'final'  # 只在最后执行
  [../]
[]
[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # 指定时间积分格式为二阶后向差分格式(backward differentiation formula)。

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  l_max_its = 50 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 10 # Max number of nonlinear iterations

  end_time = 1
  dt = 1
[]

[Outputs]
  exodus = true
  [./csv_final]
    type = CSV
    file_base = 'results/final_bnds_data'
    execute_on = 'final'
  [../]
[]
