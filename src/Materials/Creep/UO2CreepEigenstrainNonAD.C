//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "UO2CreepEigenstrainNonAD.h"

registerMooseObject("reproductionApp", UO2CreepEigenstrainNonAD);

InputParameters
UO2CreepEigenstrainNonAD::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("根据蠕变率显式积分累积蠕变本征应变。");
  params.addParam<MaterialPropertyName>("creep_rate_property", "creep_rate", "蠕变率张量材料属性名");
  params.addParam<MaterialPropertyName>("effective_creep_property",
                                        "effective_creep",
                                        "等效蠕变率材料属性名");
  return params;
}

UO2CreepEigenstrainNonAD::UO2CreepEigenstrainNonAD(const InputParameters & parameters)
  : ComputeEigenstrainBase(parameters),
    _creep_rate(getMaterialProperty<RankTwoTensor>(
        _base_name + getParam<MaterialPropertyName>("creep_rate_property"))),
    _effective_creep(
        getMaterialProperty<Real>(_base_name + getParam<MaterialPropertyName>("effective_creep_property"))),
    _creep_strain(declareProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain")),
    _effective_creep_strain(declareProperty<Real>(_base_name + "effective_creep_strain")),
    _effective_creep_strain_old(getMaterialPropertyOld<Real>(_base_name + "effective_creep_strain"))
{
}

void
UO2CreepEigenstrainNonAD::initQpStatefulProperties()
{
  _creep_strain[_qp].zero();
  _effective_creep_strain[_qp] = 0.0;
  _eigenstrain[_qp].zero();
}

void
UO2CreepEigenstrainNonAD::computeQpEigenstrain()
{
  _creep_strain[_qp] = _creep_strain_old[_qp] + _creep_rate[_qp] * _dt;
  _effective_creep_strain[_qp] = _effective_creep_strain_old[_qp] + _effective_creep[_qp] * _dt;
  _eigenstrain[_qp] = _creep_strain[_qp];
}

