
# conda activate moose && dos2unix in593_15.i &&mpirun -n 4 /home/yp/projects/reproduction/reproduction-opt -i in593_15.i

#conda activate moose && dos2unix 2D.i&&mpirun -n 10 /home/yp/projects/reproduction/reproduction-opt -i 2D.i --mesh-only 2D.e
#《《下面数据取自[1]Mechanism study and theoretical simulation on heat split phenomenon in dual-cooled annular fuel element》》
# 双冷却环形燃料几何参数 (单位：mm)
pellet_inner_diameter = 9.900         # 芯块内直径
pellet_outer_diameter = 14.100         # 芯块外直径               # 轴向长度(m)

# 网格控制参数n_azimuthal = 512时网格尺寸为6.8e-5m
mesh_size = 8e-5 #网格尺寸即可
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
growth_factor = 1.006       # 径向增长因子
# 计算半径参数 (转换为米)
pellet_inner_radius = '${fparse pellet_inner_diameter/2*1e-3}'
pellet_outer_radius = '${fparse pellet_outer_diameter/2*1e-3}'
[Mesh]
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
    subdomain_id = 1
  []
  [rename1]
    type = RenameBoundaryGenerator
    input = pellet
    old_boundary = 'pellet_rmin pellet_rmax'
    new_boundary = 'pellet_inner pellet_outer'
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
    old_block = '1'
    new_block = 'pellet'
  []
[]



