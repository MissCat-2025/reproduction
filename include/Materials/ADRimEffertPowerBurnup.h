// ADRimEffertPowerBurnup.h
#pragma once

#include "ADMaterial.h"

class ADRimEffertPowerBurnup : public ADMaterial
{
public:
  static InputParameters validParams();
  ADRimEffertPowerBurnup(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;
  /// 功率历史函数
  const Function & _power_history;

  /// 芯块半径
  const Real _pellet_inner_radius;
  /// 芯块半径
  const Real _pellet_outer_radius;
  const Real _A;//控制左最大值
  const Real _B;//控制右最大值值
  const Real _C;//控制最小值
  const Real _D;//计算的燃耗终点，单位为0.01为1%燃耗（为MWd/kg除8.3）

  /// 总功率密度
  ADMaterialProperty<Real> & _total_power;

  /// 径向功率分布形状
  ADMaterialProperty<Real> & _radial_power_shape;

  //燃耗部分
  /// 初始燃料密度 (kg/m³)
  const ADMaterialProperty<Real> & _density;
  /// 声明材料属性
  ADMaterialProperty<Real> & _burnup;
  const MaterialProperty<Real> & _burnup_old;
  /// 计算功率分布因子
  ADReal powerFactor(const Real & r) const;
};