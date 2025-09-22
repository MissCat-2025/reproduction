//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "QuadraturePointMarker.h"
#include "Coupleable.h"
#include <deque>
#include <unordered_map>

class PhasePiledFractureHSMarker : public QuadraturePointMarker, public Coupleable
{
public:
  static InputParameters validParams();

  PhasePiledFractureHSMarker(const InputParameters & parameters);

  virtual void markerSetup() override;

protected:
  virtual MarkerValue computeQpMarker() override;

  // 基础变量
  const VariableValue & _von_mises;
  
  // 材料变量
  const ADMaterialProperty<Real> & _sigma0;
  
  // 基础常数参数
  Real _x1;    // 第一个d阈值
  Real _x2;    // 第二个d阈值
  Real _xmax;  // 最大d阈值
  Real _y1;    // 第一个应力阈值系数
  Real _y2;    // 第二个应力阈值系数

  // 时间步长检测参数
  const unsigned int _timeD = 3;           // d值变化检测次数
  const unsigned int _timeStress = 6;      // von Mises应力变化检测次数
  
  // 变化阈值参数
  Real _d_change_threshold;                // d变化的阈值
  Real _stress_change_threshold;           // 应力变化的阈值
  
  // 历史数据存储（基于单元ID）
  std::unordered_map<dof_id_type, std::deque<Real>> _d_history;
  std::unordered_map<dof_id_type, std::deque<Real>> _stress_history;
  
  // 辅助函数：检查d值和应力的变化趋势
  bool checkDChangeOverTime(Real current_d_value);
  bool checkStressChangeOverTime(Real current_stress_value);
};
