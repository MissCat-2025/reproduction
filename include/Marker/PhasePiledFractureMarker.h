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

class PhasePiledFractureMarker : public QuadraturePointMarker, public Coupleable
{
public:
  static InputParameters validParams();

  PhasePiledFractureMarker(const InputParameters & parameters);

  virtual void markerSetup() override;

protected:
  virtual MarkerValue computeQpMarker() override;

  // vonMises应力变量支持
  const VariableValue & _von_mises;
  // sigma0 辅助变量（如 sigma0_mesh 或 sigma0_field）
  // const VariableValue & _sigma0_var;
  
  // 判断条件参数
  Real _x1;    // 第一个d阈值
  Real _x2;    // 第二个d阈值
  Real _xmax;  // 最大d阈值
  Real _y1;    // 第一个vonMises应力阈值（保留但在计算时由sigma0*0.8替代）
  Real _y2;    // 第二个vonMises应力阈值（保留但在计算时由sigma0*0.8替代）

  // 时间步长判断参数
  unsigned int _time_steps;          // 检查的时间步长数量
  Real _d_change_threshold;          // d变化的阈值
  bool _enable_time_check;           // 是否启用时间步长检查
  bool _dChangeEro;                  // 检测是否发生d的突变错误

  // 历史数据存储（基于单元ID）
  std::unordered_map<dof_id_type, std::deque<Real>> _d_history;
  
  // 辅助函数：检查d值的变化趋势
  bool checkDChangeOverTime(Real current_d_value);
};
