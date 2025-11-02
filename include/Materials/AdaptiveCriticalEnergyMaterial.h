#pragma once

#include "Material.h"

/**
 * 自适应网格下的临界能量材料
 * 计算修正后的pellet_critical_energy = (1+grid_sizes/2/length_scale_parameter)*Gc
 * 其中grid_sizes自动计算为单元面积的平方根
 */
class AdaptiveCriticalEnergyMaterial : public Material
{
public:
  static InputParameters validParams();
  AdaptiveCriticalEnergyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// 长度尺度参数
  const Real & _length_scale_parameter;
  
  /// 断裂能
  const Real & _gc;
  
  /// 输出的修正后的临界能量
  MaterialProperty<Real> & _critical_energy;
};