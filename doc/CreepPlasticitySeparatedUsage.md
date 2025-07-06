# 蠕变-塑性分离算法使用说明

## 概述

`SmallDeformationCreepPlasticitySeparatedMod` 实现了基于Fachinotti等人(2014)论文的两阶段蠕变-塑性分离算法。该算法在MOOSE框架下提供了更精确的蠕变-塑性耦合模拟。

## 算法特点

### 两阶段返回映射
1. **阶段1：非塑性试验状态**
   - 假设没有塑性增量(\(Δε^{pl} = 0\))
   - 求解蠕变方程：\(F_1(ε_{eq}^{cr,np}) = 0\)
   - 检查屈服条件

2. **阶段2：塑性+蠕变耦合**
   - 如果仍然屈服，同时考虑塑性和蠕变
   - 求解耦合方程：\(F_2(ε_{eq}^{pl}) = 0\)

### 优势
- **精确分离**：严格区分塑性和蠕变应变
- **框架兼容**：完全基于现有MOOSE返回映射框架
- **数值稳定**：使用状态机模式确保计算稳定性

## 使用方法

### 1. 输入文件配置

```bash
[Materials]
  [creep_plasticity]
    type = SmallDeformationCreepPlasticitySeparatedMod
    hardening_model = hardening_model
    coefficient = 1e-20        # 蠕变系数 A
    exponent = 5.0             # 蠕变指数 n
    epsilon_tolerance = 1e-9   # 数值稳定性参数
    stress_tolerance = 1e-6    # 应力容差
    base_name = mech
  []
  
  [hardening_model]
    type = IsotropicHardeningModel
    yield_stress = 300e6       # 初始屈服应力
    hardening_modulus = 1e9    # 硬化模量
    base_name = mech
  []
  
  [elasticity_model]
    type = ComputeSmallDeformationElasticityModelMod
    youngs_modulus = 2e11      # 杨氏模量
    poissons_ratio = 0.3       # 泊松比
    base_name = mech
  []
  
  [stress]
    type = ComputeSmallDeformationStressMod
    displacements = 'disp_x disp_y'
    plasticity_model = creep_plasticity
    base_name = mech
  []
[]
```

### 2. 参数说明

| 参数 | 说明 | 典型值 |
|------|------|--------|
| `coefficient` | 蠕变系数 A | 1e-20 ~ 1e-15 |
| `exponent` | 蠕变指数 n | 3.0 ~ 10.0 |
| `epsilon_tolerance` | 数值稳定性参数 | 1e-9 |
| `stress_tolerance` | 应力容差 | 1e-6 |

### 3. 输出变量

算法会自动生成以下状态变量：
- `mech_effective_plastic_strain` - 有效塑性应变
- `mech_effective_creep_strain` - 有效蠕变应变
- `mech_plastic_strain` - 塑性应变张量
- `mech_creep_strain` - 蠕变应变张量
- `mech_stress` - 应力张量

## 测试验证

### 运行测试
```bash
cd reproduction
./run_tests --re creep_plasticity_separated
```

### 预期结果
- 弹性阶段：只有弹性应变
- 屈服后：塑性应变开始累积
- 长期加载：蠕变效应逐渐显现
- 应变分离：塑性和蠕变应变独立更新

## 数值考虑

### 时间步长
- 建议使用较小的时间步长(dt < 100s)
- 蠕变过程需要足够的时间积累
- 避免过大的应变增量

### 收敛性
- 可能需要调整牛顿迭代参数
- 建议使用: `nl_abs_tol = 1e-8`, `nl_rel_tol = 1e-8`
- 如果收敛困难，可以减小时间步长

### 单位一致性
- 确保所有参数使用一致的单位系统
- 应力: Pa, 时间: s, 应变: 无量纲

## 扩展功能

### 复杂蠕变模型
可以通过重写 `computeCreepRate()` 函数来实现更复杂的蠕变本构：

```cpp
ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & creep_strain) override
{
  // 例如：Norton-Bailey模型
  // ε̇ = A * σ^n * ε^m
  const ADReal stress_safe = std::max(effective_stress, _stress_tolerance);
  const ADReal epsilon_safe = std::max(creep_strain, _epsilon_tolerance);
  
  return _coefficient * std::pow(stress_safe, _exponent) * std::pow(epsilon_safe, _strain_exponent);
}
```

### 温度依赖
可以添加温度相关的蠕变模型：

```cpp
// Arrhenius温度依赖
ADReal temperature_factor = std::exp(-_activation_energy / (_gas_constant * _temperature));
return _coefficient * std::pow(stress_safe, _exponent) * temperature_factor;
```

## 常见问题

### Q1: 收敛困难怎么办？
- 检查时间步长是否过大
- 调整容差参数
- 使用更好的初值猜测

### Q2: 蠕变效应不明显？
- 检查蠕变系数是否合理
- 确保应力水平足够高
- 延长模拟时间

### Q3: 塑性和蠕变应变异常？
- 检查硬化模型参数
- 验证弹性模量设置
- 确保边界条件正确

## 参考文献

1. Fachinotti, V. D., A. E. Albanesi, and A. Cardona. "Return Mapping for Creep and Plasticity Split." Computational Mechanics, 2014.

2. Simo, J. C., and T. J. R. Hughes. "Computational Inelasticity." Springer-Verlag, 1998.

## 支持与反馈

如果遇到问题或有改进建议，请通过以下方式联系：
- 提交GitHub Issue
- 查阅MOOSE文档
- 参考相关论文 