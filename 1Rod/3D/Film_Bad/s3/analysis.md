# s3 EqualValueBoundaryConstraint排查记录

## 修改内容
- 启用 EqualValueBoundaryConstraint 以保持顶面 disp_z 等值约束
- penalty 使用已有的 E / length_scale_paramete 尺度化惩罚系数

## 现象与机制判断
- EqualValueBoundaryConstraint 会选取顶面边界中的一个主节点，并强制该边界所有节点的 disp_z 与主节点一致
- 底面 disp_z 固定为 0 后，该约束使得厚度方向应变趋于常数，表现为广义平面应变的近似
- 该约束通过惩罚项实现，其刚度不随损伤退化，约束存在会抬升整体刚度并使 psie_active 空间分布趋于均匀，从而导致损伤沿径向弥散而非在裂纹尖端局部化

## 裂纹驱动力层面的本质问题
- 相场裂纹驱动力来自自由能对 d 的变分，当前模型中 psi = alpha*Gc/c0/l + g*(psie_active)
- 驱动力项由 g'(d)*psie_active 控制，psie_active 越高、分布越均匀，d 趋势越倾向弥散式增长
- EqualValueBoundaryConstraint 将顶面 z 位移与单一主节点强耦合，等效引入刚性“顶面等位移”约束，使 eps_zz 在厚度方向被锁定为常数并牵引整体应力场重分布
- 在热膨胀与外压共同作用下，该约束更易把拉应变能集中到外缘并形成近似均匀抬升的 psie_active 场，裂纹驱动力从边缘主导扩展而非尖端局部化

## 退化函数 g 的角色
- g 已在 Main.i 的 [degradation] 中定义，并在 [elasticity] 中通过 degradation_function = g 参与应力退化
- g 也在 Sub.i 的 psi 中参与裂纹驱动力项 g*(psie_active)
- EqualValueBoundaryConstraint 本身不支持材料退化，输入文件内无法直接用 g(d) 对惩罚项缩放
- 若需要约束随损伤减弱，需要在应用中实现自定义退化约束类并在惩罚能中引入 g(d)

## 排查顺序建议
- 对比启用与禁用约束时的 psie_active 分布，确认是否出现整体抬升或均匀化
- 监控 sigma_zz 与 eps_zz 的时空演化，判断惩罚项是否主导厚度方向响应
- 若可改为 2D 模型，改用 generalized plane strain 的标量应变方式替代顶面等位移约束，比较 psie_active 是否恢复局部化
