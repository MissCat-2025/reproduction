//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhasePiledFractureMarker.h"
#include "FEProblem.h"
#include "MooseEnum.h"

registerMooseObject("reproductionApp", PhasePiledFractureMarker);

InputParameters
PhasePiledFractureMarker::validParams()
{
  InputParameters params = QuadraturePointMarker::validParams();

  // vonMises应力变量参数 - 使用addCoupledVar添加耦合变量
  params.addRequiredCoupledVar("von_mises_variable", "vonMises应力变量名");
  
  // 判断条件参数
  params.addRequiredParam<Real>("x1", "d变量的第一个阈值");
  params.addRequiredParam<Real>("x2", "d变量的第二个阈值");
  params.addRequiredParam<Real>("xmax", "d变量的最大阈值");
  params.addRequiredParam<Real>("y1", "vonMises应力的第一个阈值");
  params.addRequiredParam<Real>("y2", "vonMises应力的第二个阈值");

  // 时间步长判断参数
  params.addParam<bool>("enable_time_check", false, "是否启用基于时间步长的d变化检查");
  params.addParam<unsigned int>("time_steps", 5, "检查d变化的时间步长数量");
  params.addParam<Real>("d_change_threshold", 0.05, "判断d变化的阈值");

  params.addClassDescription("基于指定的上下界值和变量，结合vonMises应力考虑，标记网格自适应单元");
  return params;
}

PhasePiledFractureMarker::PhasePiledFractureMarker(const InputParameters & parameters)
  : QuadraturePointMarker(parameters),
    Coupleable(this, false),
    _von_mises(coupledValue("von_mises_variable")),
    _x1(parameters.get<Real>("x1")),
    _x2(parameters.get<Real>("x2")),
    _xmax(parameters.get<Real>("xmax")),
    _y1(parameters.get<Real>("y1")),
    _y2(parameters.get<Real>("y2")),
    _dChangeEro(false),
    _time_steps(parameters.get<unsigned int>("time_steps")),
    _d_change_threshold(parameters.get<Real>("d_change_threshold")),
    _enable_time_check(parameters.get<bool>("enable_time_check"))
{
  // 检查参数的合理性
  if (_x1 >= _x2 || _x2 >= _xmax)
    mooseError("阈值无效: 要求 x1 < x2 < xmax");
    
  if (_enable_time_check && _time_steps < 2)
    mooseError("时间步长数量至少为2");
    
  if (_enable_time_check && _d_change_threshold <= 0)
    mooseError("d变化阈值必须为正数");
}

void
PhasePiledFractureMarker::markerSetup()
{
  // 如果启用了时间检查，清理过期的历史数据
  if (_enable_time_check)
  {
    // 遍历历史数据，清理过长的记录
    for (auto & pair : _d_history)
    {
      while (pair.second.size() > _time_steps)
        pair.second.pop_front();
    }
  }
}

bool
PhasePiledFractureMarker::checkDChangeOverTime(Real current_d_value)
{
  if (!_enable_time_check)
    return false;
  
  _dChangeEro = false;
  dof_id_type elem_id = _current_elem->id();
  
  // 更新当前单元的d值历史
  if (_d_history.find(elem_id) == _d_history.end())
    _d_history[elem_id] = std::deque<Real>();
    
  _d_history[elem_id].push_back(current_d_value);
  
  // 限制历史数据长度
  if (_d_history[elem_id].size() > _time_steps)
    _d_history[elem_id].pop_front();
  
  // 检查是否有足够的历史数据
  if (_d_history[elem_id].size() < _time_steps)
    return true;
    
  // 计算全局最大变化量（任意两步的最大差）和相邻两步最大变化
  const auto & history = _d_history[elem_id];
  Real min_d = *std::min_element(history.begin(), history.end());
  Real max_d = *std::max_element(history.begin(), history.end());
  Real max_change = max_d - min_d; // 列表内任意两步相差的最大值

  Real next_max_change = 0.0;      // 相邻两步之间的最大变化
  for (std::size_t i = 1; i < history.size(); ++i)
  {
    Real diff = history[i] - history[i - 1];
    if (diff < 0)
      diff = -diff;
    if (diff > next_max_change)
      next_max_change = diff;
  }

  // _dChangeEro 使用相邻两步最大变化进行判定
  _dChangeEro = (next_max_change > 0.5);

  // 粗化判定仍然使用任意两步最大差与阈值比较
  return max_change < _d_change_threshold;
}

Marker::MarkerValue
PhasePiledFractureMarker::computeQpMarker()
{
  // 获取当前积分点的d和vonMises应力值
  Real d_value = _u[_qp];
  Real von_mises_value = _von_mises[_qp];
  
  // 首先检查时间步长判断逻辑
  if (_enable_time_check && checkDChangeOverTime(d_value))
  {
    // 如果d在指定时间步长内变化不大，则粗化
    return COARSEN;
  }
  if(_dChangeEro)
  {
    return COARSEN;
  }
  // else
  // {
  //   return REFINE;
  // }
  // 基于d和vonMises应力值的判断逻辑（相对操作）
  // 粗网格 → COARSEN（相对当前级别粗化）
  // 细网格 → DO_NOTHING（保持当前级别）
  // 细细网格 → REFINE（相对当前级别细化）
  
  if (d_value <= _x1)
  {
    if (von_mises_value <= _y1)
      return COARSEN;     // 粗网格 - 相对粗化
    else
    {
      if (_enable_time_check && checkDChangeOverTime(d_value))
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
    if (von_mises_value <= _y1)
      return COARSEN;     // 粗网格 - 相对粗化
    else
      {
      if (_enable_time_check && checkDChangeOverTime(d_value))
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
    if (von_mises_value < _y2)
      return COARSEN;      // 细细网格 - 相对细化
    else
    {      
      if (_enable_time_check && checkDChangeOverTime(d_value))
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
    return REFINE;       // 细网格 - 相对粗化
  }
  else if (d_value > 0.98)
  {
    return COARSEN;       // 细网格 - 相对粗化
  }
  else
  {
    return COARSEN;
  }
  
    if(_dChangeEro)
  {
    return COARSEN;
  }
}
