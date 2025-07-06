//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PlasticityModel.h"
#include "ElasticityModel.h"
#include "CreepModel.h"

InputParameters
CreepModel::validParams()
{
  InputParameters params = Material::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params += BaseNameInterface::validParams();
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

CreepModel::CreepModel(const InputParameters & parameters)
  : Material(parameters),
    ADSingleVariableReturnMappingSolution(parameters),
    BaseNameInterface(parameters),
    _creep_strain(declareADProperty<RankTwoTensor>(prependBaseName("creep_strain"))),
    _creep_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("creep_strain"))),
    _ec(declareADProperty<Real>(prependBaseName("effective_creep_strain"))),
    _ec_old(getMaterialPropertyOldByName<Real>(prependBaseName("effective_creep_strain"))),
    _Np(declareADProperty<RankTwoTensor>(prependBaseName("flow_direction")))
{
}

// void
// CreepModel::initialSetup()
// {
// }

void
CreepModel::setQp(unsigned int qp)
{
  _qp = qp;
}

void
CreepModel::setElasticityModel(
    ElasticityModel * elasticity_model)
{
  _elasticity_model = elasticity_model;
}

void
CreepModel::setPlasticityModel(
    PlasticityModel * plasticity_model)
{
  _plasticity_model = plasticity_model;
}

void
CreepModel::initQpStatefulProperties()
{
  _creep_strain[_qp].zero();
  _ec[_qp] = 0;
}
