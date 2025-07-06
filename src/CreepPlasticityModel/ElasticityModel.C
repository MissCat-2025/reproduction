//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ElasticityModel.h"
#include "CreepModel.h"

InputParameters
ElasticityModel::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

ElasticityModel::ElasticityModel(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _creep_model(nullptr),
    _elastic_strain(declareADProperty<RankTwoTensor>(prependBaseName("elastic_strain")))
{
}

void
ElasticityModel::setQp(unsigned int qp)
{
  _qp = qp;
  if (_creep_model)
    _creep_model->setQp(qp);
}

void
ElasticityModel::setCreepModel(
    CreepModel * creep_model)
{
  _creep_model = creep_model;
  _creep_model->setElasticityModel(this);
}

void
ElasticityModel::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
}

void
ElasticityModel::updateState(const ADRankTwoTensor & mechanical_strain,
                                             ADRankTwoTensor & stress)
{
  _elastic_strain[_qp] = mechanical_strain;
  _creep_model->updateState(stress, _elastic_strain[_qp]);
}
