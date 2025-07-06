//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SmallDeformationPlasticityModelCreep.h"
#include "SmallDeformationElasticityModelCreep.h"

InputParameters
SmallDeformationPlasticityModelCreep::validParams()
{
  InputParameters params = Material::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params += BaseNameInterface::validParams();

  params.addRequiredParam<MaterialName>("hardening_model", "Name of the plastic hardening model");

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

SmallDeformationPlasticityModelCreep::SmallDeformationPlasticityModelCreep(const InputParameters & parameters)
  : Material(parameters),
    ADSingleVariableReturnMappingSolution(parameters),
    BaseNameInterface(parameters),
    _plastic_strain(declareADProperty<RankTwoTensor>(prependBaseName("plastic_strain"))),
    _plastic_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("plastic_strain"))),
    _ep(declareADProperty<Real>(prependBaseName("effective_plastic_strain"))),
    _ep_old(getMaterialPropertyOldByName<Real>(prependBaseName("effective_plastic_strain"))),
    _Np(declareADProperty<RankTwoTensor>(prependBaseName("flow_direction")))
{
}

void
SmallDeformationPlasticityModelCreep::initialSetup()
{
  _hardening_model = dynamic_cast<PlasticHardeningModel *>(&getMaterial("hardening_model"));
  if (!_hardening_model)
    paramError("hardening_model",
               "Plastic hardening model " + getParam<MaterialName>("hardening_model") +
                   " is not compatible with " + name());
}

void
SmallDeformationPlasticityModelCreep::setQp(unsigned int qp)
{
  _qp = qp;
  _hardening_model->setQp(qp);
}

void
SmallDeformationPlasticityModelCreep::setElasticityModel(
    SmallDeformationElasticityModelCreep * elasticity_model)
{
  _elasticity_model = elasticity_model;
}

void
SmallDeformationPlasticityModelCreep::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
  _ep[_qp] = 0;
}
