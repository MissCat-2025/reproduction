//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SmallDeformationPlasticityModelMod.h"
#include "SmallDeformationElasticityModelMod.h"

InputParameters
SmallDeformationPlasticityModelMod::validParams()
{
  InputParameters params = Material::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params += BaseNameInterface::validParams();

  params.addRequiredParam<MaterialName>("hardening_model", "Name of the plastic hardening model");

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

SmallDeformationPlasticityModelMod::SmallDeformationPlasticityModelMod(const InputParameters & parameters)
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
SmallDeformationPlasticityModelMod::initialSetup()
{
  _hardening_model = dynamic_cast<PlasticHardeningModel *>(&getMaterial("hardening_model"));
  if (!_hardening_model)
    paramError("hardening_model",
               "Plastic hardening model " + getParam<MaterialName>("hardening_model") +
                   " is not compatible with " + name());
}

void
SmallDeformationPlasticityModelMod::setQp(unsigned int qp)
{
  _qp = qp;
  _hardening_model->setQp(qp);
}

void
SmallDeformationPlasticityModelMod::setElasticityModel(
    SmallDeformationElasticityModelMod * elasticity_model)
{
  _elasticity_model = elasticity_model;
}

void
SmallDeformationPlasticityModelMod::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
  _ep[_qp] = 0;
}
