# 简化CT模型 - 实用版本
# 主要简化：
# 1. 将圆孔简化为在相应位置施加载荷
# 2. 缺口简化为矩形
# 3. 专注于裂纹扩展行为
# conda activate moose && dos2unix CT_model.i && mpirun -n 8 /home/yp/projects/reproduction/reproduction-opt -i CT_model.i --mesh-only CT_model.e
[Mesh]
  type = FileMesh
  file = CT_model.msh
[]