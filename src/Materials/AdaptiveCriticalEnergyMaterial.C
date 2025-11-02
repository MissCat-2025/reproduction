#include "AdaptiveCriticalEnergyMaterial.h"
#include "MooseUtils.h"

registerMooseObject("reproductionApp", AdaptiveCriticalEnergyMaterial);

InputParameters
AdaptiveCriticalEnergyMaterial::validParams()
{
  InputParameters params = Material::validParams();
  
  params.addRequiredParam<Real>("length_scale_parameter", "相场正则化长度参数");
  params.addRequiredParam<Real>("gc", "断裂能");
  params.addClassDescription("计算自适应网格下的修正临界能量: pellet_critical_energy = (1+grid_sizes/2/length_scale_parameter)*Gc，其中grid_sizes自动计算为单元面积的平方根");
  
  return params;
}

AdaptiveCriticalEnergyMaterial::AdaptiveCriticalEnergyMaterial(const InputParameters & parameters)
  : Material(parameters),
    _length_scale_parameter(getParam<Real>("length_scale_parameter")),
    _gc(getParam<Real>("gc")),
    _critical_energy(declareProperty<Real>("pellet_critical_energy"))
{
}

void
AdaptiveCriticalEnergyMaterial::computeQpProperties()
{
  // 计算当前单元的面积
  Real element_volume = _current_elem->volume();
  
  // 对于2D问题，面积就是体积；对于3D问题，需要计算面积
  Real element_area;
  if (_current_elem->dim() == 2)
    element_area = element_volume;
  else if (_current_elem->dim() == 3)
    element_area = std::pow(element_volume, 2.0/3.0); // 假设单元是立方体
  else
    element_area = element_volume; // 1D情况
  
  // 计算网格尺寸为面积的平方根
  Real grid_size = std::sqrt(element_area);
  
  // 计算修正后的临界能量: pellet_critical_energy = (1+grid_sizes/2/length_scale_parameter)*Gc
  _critical_energy[_qp] = (1.0 + grid_size / 2.0 / _length_scale_parameter) * _gc;
}