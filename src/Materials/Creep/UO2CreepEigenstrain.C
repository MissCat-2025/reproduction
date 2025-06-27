// src/materials/UO2CreepEigenstrain.C
#include "UO2CreepEigenstrain.h"

registerADMooseObject("raccoonApp", UO2CreepEigenstrain);

InputParameters
UO2CreepEigenstrain::validParams()
{
  InputParameters params = ADComputeEigenstrainBase::validParams();
    // 添加这一行来声明vonMisesStress参数
  params.addClassDescription("计算UO2的蠕变特征应变");
  params.addParam<bool>("consider_psip_active", false, "Whether to consider psip_active");
  return params;
}

UO2CreepEigenstrain::UO2CreepEigenstrain(const InputParameters & parameters)
  : ADComputeEigenstrainBase(parameters),
    _creep_rate(getADMaterialProperty<RankTwoTensor>("creep_rate")),
    _creep_strain(declareADProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain")),
    _consider_psip_active(getParam<bool>("consider_psip_active")),
    _psip_active(declareADProperty<Real>("psip_active")),
    _psip_active_old(getMaterialPropertyOld<Real>("psip_active")),
    _stress_deviator(getADMaterialProperty<RankTwoTensor>("stress_deviator"))
{
}

void
UO2CreepEigenstrain::initQpStatefulProperties()
{
  _creep_strain[_qp].zero();
  _eigenstrain[_qp].zero();
  _psip_active[_qp] = 0.0;
}

void
UO2CreepEigenstrain::computeQpEigenstrain()
{
  // 更新累积蠕变应变
  _creep_strain[_qp] = _creep_strain_old[_qp] + _creep_rate[_qp] * _dt;
  
  // 设置特征应变
  _eigenstrain[_qp] = _creep_strain[_qp];
  if (_consider_psip_active)
  {
  // 计算当前时间步的蠕变能增量
  ADReal increment = _creep_rate[_qp].doubleContraction(_stress_deviator[_qp]) * _dt;
  
    // 累积蠕变能 = 历史值 + 当前增量
    _psip_active[_qp] = _psip_active_old[_qp] + increment;
  }
  else
  {
    _psip_active[_qp] = 0.0;
  }
}