//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PlasticityModel.h"
#include "ElasticityModel.h"

InputParameters
PlasticityModel::validParams()
{
  InputParameters params = Material::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params += BaseNameInterface::validParams();

  params.addRequiredParam<MaterialName>("hardening_model", "Name of the plastic hardening model");

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

PlasticityModel::PlasticityModel(const InputParameters & parameters)
  : Material(parameters),
    ADSingleVariableReturnMappingSolution(parameters),
    BaseNameInterface(parameters),
    _plastic_strain(declareADProperty<RankTwoTensor>(prependBaseName("plastic_strain"))),
    _plastic_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("plastic_strain"))),
    _ep(declareADProperty<Real>(prependBaseName("effective_plastic_strain"))),
    _ep_old(getMaterialPropertyOldByName<Real>(prependBaseName("effective_plastic_strain"))),
    _Np(declareADProperty<RankTwoTensor>(prependBaseName("plastic_flow_direction"))),
    _creep_model(nullptr)
{
}

void
PlasticityModel::initialSetup()
{
  _hardening_model = dynamic_cast<PlasticHardeningModel *>(&getMaterial("hardening_model"));
  if (!_hardening_model)
    paramError("hardening_model",
               "Plastic hardening model " + getParam<MaterialName>("hardening_model") +
                   " is not compatible with " + name());
}

void
PlasticityModel::setQp(unsigned int qp)
{
  _qp = qp;
  _hardening_model->setQp(qp);
}

void
PlasticityModel::setElasticityModel(ElasticityModel * elasticity_model)
{
  _elasticity_model = elasticity_model;
}

void
PlasticityModel::setCreepModel(CreepModel * creep_model)
{
  _creep_model = creep_model;
}

void
PlasticityModel::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
  _ep[_qp] = 0;
}

Real
PlasticityModel::getEffectivePlasticStrainOld() const
{
  return _ep_old[_qp];
}

const RankTwoTensor &
PlasticityModel::getPlasticStrainOld() const
{
  return _plastic_strain_old[_qp];
}
