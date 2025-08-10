# (已验证)这第一步就是测试生成网格文件KAERI_HANARO_UpperRod1.e   注意：这文件生成的是完整的3D模型，而不是1/4模型3D模型
# 语法要求：仅仅为了生成网格文件：Run with --mesh-only:
#https://mooseframework.inl.gov/source/meshgenerators/ConcentricCircleMeshGenerator.html
#conda activate moose && dos2unix Complete2DQuarter.i&&mpirun -n 10 /home/yp/projects/reproduction/reproduction-opt -i Complete2DQuarter.i --mesh-only KAERI_HANARO_UpperRod1.e
#《《下面数据取自[1]Thermomechanical Analysis and Irradiation Test of Sintered Dual-Cooled Annular pellet》》

# 双冷却环形燃料几何参数 (单位：mm)
inclad_inner_diameter = 9.0      # 内包壳内直径
inclad_outer_diameter = 10.14    # 内包壳外直径
pellet_inner_diameter = 10.271         # 芯块内直径
pellet_outer_diameter = 14.685         # 芯块外直径
outclad_inner_diameter = 14.76    # 外包壳内直径
outclad_outer_diameter = 15.9     # 外包壳外直径                   # 轴向长度(m)

# 网格控制参数n_azimuthal = 512时网格尺寸为6.8e-5m
n_radial_inner_clad = 3    # 内包壳径向单元数
mesh_size = 8e-5 #网格尺寸即可
n_azimuthal = '${fparse int(3.1415*(pellet_outer_diameter)/mesh_size*1e-3/4)*4}' #int()取整
n_radial_pellet = '${fparse int((pellet_outer_diameter-pellet_inner_diameter)/mesh_size*1e-3/2)}'
n_radial_outer_clad = 3    # 外包壳径向单元数
growth_factor = 1.006       # 径向增长因子
# 计算半径参数 (转换为米)
inner_clad_inner_radius = '${fparse inclad_inner_diameter/2*1e-3}'
inner_clad_outer_radius = '${fparse inclad_outer_diameter/2*1e-3}'
pellet_inner_radius = '${fparse pellet_inner_diameter/2*1e-3}'
pellet_outer_radius = '${fparse pellet_outer_diameter/2*1e-3}'
outer_clad_inner_radius = '${fparse outclad_inner_diameter/2*1e-3}'
outer_clad_outer_radius = '${fparse outclad_outer_diameter/2*1e-3}'

[Mesh]
  [inner_clad1]
    type = AnnularMeshGenerator
    nr = ${n_radial_inner_clad}
    nt = ${n_azimuthal}
    rmin = ${inner_clad_inner_radius}
    rmax = ${inner_clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 10
    boundary_name_prefix = 'inclad'
  []
  [inner_clad]
    type = SubdomainIDGenerator
    input = inner_clad1
    subdomain_id = 1
  []
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
    subdomain_id = 2
  []
  [outer_clad1]
    type = AnnularMeshGenerator
    nr = ${n_radial_outer_clad}
    nt = ${n_azimuthal}
    rmin = ${outer_clad_inner_radius}
    rmax = ${outer_clad_outer_radius}
    growth_r = ${growth_factor}
    boundary_id_offset = 30
    boundary_name_prefix = 'outclad'
  []
  [outer_clad]
    type = SubdomainIDGenerator
    input = outer_clad1
    subdomain_id = 3
  []
  [combine]
    type = CombinerGenerator
    inputs = 'inner_clad pellet outer_clad'
  []
  [rename1]
    type = RenameBoundaryGenerator
    input = combine
    old_boundary = 'inclad_rmin inclad_rmax pellet_rmin pellet_rmax outclad_rmin outclad_rmax'
    new_boundary = 'inclad_inner inclad_outer pellet_inner pellet_outer outclad_inner outclad_outer'
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
    old_block = '1 2 3'
    new_block = 'inclad pellet outclad'
  []
[]
