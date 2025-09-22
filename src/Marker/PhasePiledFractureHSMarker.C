//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhasePiledFractureHSMarker.h"
#include "FEProblem.h"
#include "MooseEnum.h"

registerMooseObject("reproductionApp", PhasePiledFractureHSMarker);

InputParameters
PhasePiledFractureHSMarker::validParams()
{
  InputParameters params = QuadraturePointMarker::validParams();

  // vonMises应力变量参数
  params.addRequiredCoupledVar("von_mises_variable", "vonMises应力变量名");
  
  // sigma0材料变量参数
  params.addRequiredParam<MaterialPropertyName>("sigma0", "临界断裂强度材料属性名");
  
  // 判断条件参数
  params.addRequiredParam<Real>("x1", "d变量的第一个阈值");
  params.addRequiredParam<Real>("x2", "d变量的第二个阈值");
  params.addRequiredParam<Real>("xmax", "d变量的最大阈值");
  // 变化阈值参数
  params.addParam<Real>("y1", 0.5, "应力的第一个阈值系数");
  params.addParam<Real>("y2", 0.6, "应力的第二个阈值系数");
  params.addParam<unsigned int>("timeD", 5, "检查d变化的时间步长数量");
  params.addParam<unsigned int>("timeStress", 5, "检查von_mises_variable变化的时间步长数量");
  // 变化阈值参数
  params.addParam<Real>("d_change_threshold", 0.05, "判断d变化的阈值");
  params.addParam<Real>("stress_change_threshold", 3e6, "判断应力变化的阈值");

  params.addClassDescription("基于指定的上下界值和变量，结合vonMises应力和sigma0考虑，标记网格自适应单元");
  return params;
}

PhasePiledFractureHSMarker::PhasePiledFractureHSMarker(const InputParameters & parameters)
  : QuadraturePointMarker(parameters),
    Coupleable(this, false),
    _von_mises(coupledValue("von_mises_variable")),
    _sigma0(getADMaterialProperty<Real>("sigma0")),
    _x1(parameters.get<Real>("x1")),
    _x2(parameters.get<Real>("x2")),
    _xmax(parameters.get<Real>("xmax")),
    _y1(parameters.get<Real>("y1")),
    _y2(parameters.get<Real>("y2")),
    _timeD(parameters.get<unsigned int>("timeD")),
    _timeStress(parameters.get<unsigned int>("timeStress")),
    _d_change_threshold(parameters.get<Real>("d_change_threshold")),
    _stress_change_threshold(parameters.get<Real>("stress_change_threshold"))
{
  // 检查参数的合理性
  if (_x1 >= _x2 || _x2 >= _xmax)
    mooseError("阈值无效: 要求 x1 < x2 < xmax");
    
  if (_d_change_threshold <= 0)
    mooseError("d变化阈值必须为正数");
    
  if (_stress_change_threshold <= 0)
    mooseError("应力变化阈值必须为正数");
}

void
PhasePiledFractureHSMarker::markerSetup()
{
  // 清理过期的历史数据
  for (auto & pair : _d_history)
  {
    while (pair.second.size() > _timeD)
      pair.second.pop_front();
  }
  
  for (auto & pair : _stress_history)
  {
    while (pair.second.size() > _timeStress)
      pair.second.pop_front();
  }
}

bool
PhasePiledFractureHSMarker::checkDChangeOverTime(Real current_d_value)
{
  dof_id_type elem_id = _current_elem->id();
  
  // 更新当前单元的d值历史
  if (_d_history.find(elem_id) == _d_history.end())
    _d_history[elem_id] = std::deque<Real>();
    
  _d_history[elem_id].push_back(current_d_value);
  
  // 限制历史数据长度
  if (_d_history[elem_id].size() > _timeD)
    _d_history[elem_id].pop_front();
  
  // 检查是否有足够的历史数据
  if (_d_history[elem_id].size() < _timeD)
    return false;
    
  // 计算变化量
  const auto & history = _d_history[elem_id];
  Real min_d = *std::min_element(history.begin(), history.end());
  Real max_d = *std::max_element(history.begin(), history.end());
  Real max_change = max_d - min_d;

  return max_change < _d_change_threshold;
}

bool
PhasePiledFractureHSMarker::checkStressChangeOverTime(Real current_stress_value)
{
  dof_id_type elem_id = _current_elem->id();
  
  // 更新当前单元的应力值历史
  if (_stress_history.find(elem_id) == _stress_history.end())
    _stress_history[elem_id] = std::deque<Real>();
    
  _stress_history[elem_id].push_back(current_stress_value);
  
  // 限制历史数据长度
  if (_stress_history[elem_id].size() > _timeStress)
    _stress_history[elem_id].pop_front();
  
  // 检查是否有足够的历史数据
  if (_stress_history[elem_id].size() < _timeStress)
    return false;
    
  // 计算变化量
  const auto & history = _stress_history[elem_id];
  Real min_stress = *std::min_element(history.begin(), history.end());
  Real max_stress = *std::max_element(history.begin(), history.end());
  Real max_change = max_stress - min_stress;

  return max_change < _stress_change_threshold;
}

Marker::MarkerValue
PhasePiledFractureHSMarker::computeQpMarker()
{
  // 获取当前积分点的d、vonMises应力和sigma0值
  Real d_value = _u[_qp];
  Real von_mises_value = _von_mises[_qp];
  Real sigma0_value = MetaPhysicL::raw_value(_sigma0[_qp]);  // 转换为Real类型
  
  // 计算应力阈值（使用sigma0作为基准）
  Real stress_threshold1 = _y1 * sigma0_value;  // 第一个应力阈值
  Real stress_threshold2 = _y2 * sigma0_value;  // 第二个应力阈值
  
  // 检查d值和应力的时间变化
  bool d_stable = checkDChangeOverTime(d_value);
  bool stress_stable = checkStressChangeOverTime(von_mises_value);
  
  // 基于d和vonMises应力值的判断逻辑
  if (d_value <= _x1)
  {
    if (von_mises_value <= stress_threshold1)
      return COARSEN;     // 粗网格 - 相对粗化
    else
    {
      if (d_stable && stress_stable)
      {
        // 如果d在指定时间步长内变化不大，则粗化
        return COARSEN;
      }
      else
      {
        return REFINE;
      }
    }
  }
  else if (d_value <= _x2)
  {
    if (von_mises_value <= stress_threshold1)
      return COARSEN;     // 粗网格 - 相对粗化
    else
    {
      if (d_stable && stress_stable)
      {
        // 如果d在指定时间步长内变化不大，则粗化
        return COARSEN;
      }
      else
      {
        return REFINE;
      }
    }
  }
  else if (d_value <= _xmax)
  {
    if (von_mises_value < stress_threshold2)
      return COARSEN;      // 细细网格 - 相对细化
    else
    {
      if (d_stable && stress_stable)
      {
        // 如果d在指定时间步长内变化不大，则粗化
        return COARSEN;
      }
      else
      {
        return REFINE;
      }
    }
  }
  else if (d_value > _xmax)
  {
    return REFINE;       // 细网格 - 相对细化
  }
  else if (d_value > 0.98)
  {
    return COARSEN;       // 细网格 - 相对粗化
  }
  else
  {
    return COARSEN;
  }
}
