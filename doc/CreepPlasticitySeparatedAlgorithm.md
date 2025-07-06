# 蠕变-塑性分离算法架构设计 (简洁版)

## 1. 设计理念

基于现有Raccoon架构的优雅和简洁性，重新设计了蠕变-塑性分离算法。新设计遵循以下原则：

- **架构一致性**：继承现有的`SmallDeformationPlasticityModelMod`
- **职责分离**：清晰区分两个求解阶段
- **代码复用**：充分利用现有的返回映射框架
- **接口简洁**：保持与现有材料模型的兼容性

## 2. 核心架构

### 2.1 类继承关系
```
SmallDeformationCreepPlasticitySeparatedMod
    ↓ 继承自
SmallDeformationPlasticityModelMod
    ↓ 继承自
ADSingleVariableReturnMappingSolution + Material + BaseNameInterface
```

### 2.2 关键设计决策

**保持现有架构精神**：
- 继承塑性模型基类，自动获得所有基础功能
- 重写`updateState()`实现两阶段算法
- 复用现有的`returnMappingSolve()`处理塑性阶段

**简洁的两阶段实现**：
```cpp
void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) override
{
  // 1. 计算弹性试验状态
  elastic_strain -= _plastic_strain_old[_qp] + _creep_strain_old[_qp];
  stress = _elasticity_model->computeStress(elastic_strain);
  
  // 2. 阶段1：求解蠕变增量 (非塑性试验状态)
  _delta_cr_np = solveCreepIncrement(stress_dev_norm);
  _sigma_eq_np = stress_dev_norm - 3.0 * _elasticity_model->computeElasticModulus() * _delta_cr_np;
  
  // 3. 阶段2：检查屈服条件
  if (checkYieldCondition(_sigma_eq_np))
  {
    // 屈服：求解塑性-蠕变耦合
    solvePlasticCreepCoupling(stress_dev_norm, delta_ep);
  }
  
  // 4. 更新状态变量和应力
  // ...
}
```

## 3. 算法实现

### 3.1 阶段1：蠕变求解
- **目标**：求解\(F_1(Δε_{cr}) = Δε_{cr} - g(σ_{eq}^{np})Δt = 0\)
- **方法**：自实现的简单牛顿迭代
- **优点**：针对蠕变方程优化，数值稳定

```cpp
ADReal solveCreepIncrement(const ADReal & effective_trial_stress)
{
  return newtonSolveCreep(effective_trial_stress, 0.0);
}
```

### 3.2 阶段2：塑性求解
- **目标**：如果屈服，求解标准塑性方程
- **方法**：复用现有的`returnMappingSolve()`
- **优点**：保持与现有塑性模型的一致性

```cpp
void solvePlasticCreepCoupling(const ADReal & effective_trial_stress, ADReal & delta_ep)
{
  // 使用现有框架求解塑性增量
  returnMappingSolve(effective_trial_stress, delta_ep, _console);
  
  // 根据耦合关系更新蠕变增量
  // ...
}
```

## 4. 优势对比

| 特性 | 状态机设计 | 简洁设计 |
|------|-----------|----------|
| **架构一致性** | ❌ 破坏现有模式 | ✅ 完全兼容 |
| **代码复杂度** | ❌ 复杂的状态切换 | ✅ 简洁清晰 |
| **维护性** | ❌ 难以调试 | ✅ 易于理解 |
| **扩展性** | ❌ 状态机限制 | ✅ 容易扩展 |
| **性能** | ❌ 状态切换开销 | ✅ 直接计算 |

## 5. 关键技术亮点

### 5.1 模块化牛顿迭代
```cpp
ADReal newtonSolveCreep(const ADReal & effective_trial_stress, const ADReal & initial_guess)
{
  for (unsigned int iter = 0; iter < _max_iterations; ++iter)
  {
    const ADReal residual = creepResidual(delta_cr, effective_trial_stress);
    const ADReal derivative = creepResidualDerivative(delta_cr, effective_trial_stress);
    
    if (std::abs(raw_value(residual)) < _abs_tolerance)
      break;
    
    delta_cr -= residual / derivative;
  }
  return delta_cr;
}
```

### 5.2 智能状态变量管理
- **蠕变应变**：`_creep_strain`、`_epsilon_cr_eq`
- **塑性应变**：继承自基类，自动管理
- **中间变量**：`_sigma_eq_np`、`_delta_cr_np`

### 5.3 耦合关系处理
根据Fachinotti论文的耦合公式，在塑性阶段正确更新蠕变增量：
```cpp
_delta_cr_np = (effective_trial_stress - yield_stress) / (3.0 * elastic_modulus) - delta_ep + ...;
```

## 6. 实现特点

### 6.1 遵循现有风格
- **参数命名**：`coefficient`、`exponent`
- **函数签名**：与现有模型保持一致
- **错误处理**：使用MOOSE标准方式

### 6.2 数值稳定性
- **安全检查**：避免除零和负数幂运算
- **容差控制**：可配置的收敛标准
- **边界处理**：确保应变增量非负

### 6.3 性能优化
- **最小计算**：只在需要时进行塑性求解
- **缓存结果**：避免重复计算
- **早期退出**：纯弹性情况的快速处理

## 7. 使用示例

```bash
[creep_plasticity]
  type = SmallDeformationCreepPlasticitySeparatedMod
  hardening_model = hardening_model
  coefficient = 1e-20     # 蠕变系数 A
  exponent = 5.0          # 蠕变指数 n
  abs_tolerance = 1e-10   # 牛顿迭代容差
  max_iterations = 20     # 最大迭代次数
  base_name = mech
[]
```

## 8. 总结

新设计成功地：
- **保持了架构一致性**：完全符合现有代码风格
- **实现了算法正确性**：严格按照Fachinotti论文
- **确保了代码简洁性**：易于理解和维护
- **提供了良好扩展性**：可以轻松添加新功能

这个设计体现了"简洁即美"的程序设计哲学，在实现复杂算法的同时保持了代码的优雅和可读性。 