# J2Plasticity_C 残差导数推导

## 问题描述

在耦合蠕变-塑性模型中，塑性模型 `J2Plasticity_C` 需要考虑蠕变的影响。这导致残差方程变得更加复杂，需要仔细推导其导数。

## 残差方程

耦合蠕变-塑性模型的残差方程为：

$$F_3(\Delta\varepsilon_{eq}^{pl}) = \sigma_{eq}^{el} - C(\Delta\varepsilon_{eq}^{pl} N) - C(g(\sigma_Y(\varepsilon_{eq,0}^{pl} + \Delta\varepsilon_{eq}^{pl}), \varepsilon_c)\Delta t) - \sigma_Y(\varepsilon_{eq,0}^{pl} + \Delta\varepsilon_{eq}^{pl})$$

其中：
- $\sigma_{eq}^{el}$: 弹性试验应力
- $C$: 弹性刚度系数
- $N$: 流动方向张量
- $g(\sigma_Y, \varepsilon_c)$: 蠕变率函数，依赖于屈服应力和蠕变应变
- $\sigma_Y$: 屈服应力
- $\varepsilon_c$: 有效蠕变应变
- $\Delta t$: 时间步长

## 导数推导

对残差方程求导：

$$\frac{\partial F_3}{\partial(\Delta\varepsilon_{eq}^{pl})} = 0 - C - C\left(\frac{\partial g}{\partial \sigma_Y} \cdot \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})} + \frac{\partial g}{\partial \varepsilon_c} \cdot \frac{\partial \varepsilon_c}{\partial(\Delta\varepsilon_{eq}^{pl})}\right)\Delta t - \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}$$

## 关键假设

**独立性假设**: 在耦合蠕变-塑性模型中，蠕变应变 $\varepsilon_c$ 和塑性应变 $\varepsilon_p$ 是相互独立的状态变量。

在某个时刻，蠕变应变不会因为塑性应变增量的改变而立即改变，因此：

$$\frac{\partial \varepsilon_c}{\partial(\Delta\varepsilon_{eq}^{pl})} = 0$$

## 简化后的导数

基于独立性假设，导数简化为：

$$\frac{\partial F_3}{\partial(\Delta\varepsilon_{eq}^{pl})} = -C - C\frac{\partial g}{\partial \sigma_Y} \cdot \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}\Delta t - \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}$$

## 实现要点

1. **蠕变应变的正确获取**: 在调用蠕变率函数时，必须传入当前的蠕变应变值，而不是0.0
2. **导数项的完整性**: 虽然 $\frac{\partial g}{\partial \varepsilon_c} = 0$，但仍需要正确计算 $\frac{\partial g}{\partial \sigma_Y}$
3. **硬化模型的导数**: $\frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}$ 通过硬化模型计算

## 为什么不需要 computeCreepRateStrainDerivative？

在塑性模型的残差导数中，由于：
- 蠕变应变和塑性应变增量是独立的状态变量
- $\frac{\partial \varepsilon_c}{\partial(\Delta\varepsilon_{eq}^{pl})} = 0$

因此，`computeCreepRateStrainDerivative` 的贡献为零，不需要在导数计算中使用。

## 数学验证

最终的导数表达式为：

$$\frac{\partial F_3}{\partial(\Delta\varepsilon_{eq}^{pl})} = -C - C\frac{\partial g}{\partial \sigma_Y} \cdot \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}\Delta t - \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}$$

这与代码中的实现完全一致：

```cpp
return -elastic_coeff - creep_derivative - hardening_derivative;
```

其中：
- `elastic_coeff` = $C$
- `creep_derivative` = $C\frac{\partial g}{\partial \sigma_Y} \cdot \frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}\Delta t$
- `hardening_derivative` = $\frac{\partial \sigma_Y}{\partial(\Delta\varepsilon_{eq}^{pl})}$ 