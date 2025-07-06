//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SmallDeformationElasticityModelMod.h"
#include "SmallDeformationPlasticityModelMod.h"

InputParameters
SmallDeformationElasticityModelMod::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

SmallDeformationElasticityModelMod::SmallDeformationElasticityModelMod(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _plasticity_model(nullptr),
    _elastic_strain(declareADProperty<RankTwoTensor>(prependBaseName("elastic_strain")))
{
}

void
SmallDeformationElasticityModelMod::setQp(unsigned int qp)
{
  _qp = qp;
  if (_plasticity_model)
    _plasticity_model->setQp(qp);
}

void
SmallDeformationElasticityModelMod::setPlasticityModel(
    SmallDeformationPlasticityModelMod * plasticity_model)
{
  _plasticity_model = plasticity_model;
  _plasticity_model->setElasticityModel(this);
}

void
SmallDeformationElasticityModelMod::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
}

void
SmallDeformationElasticityModelMod::updateState(const ADRankTwoTensor & mechanical_strain,
                                             ADRankTwoTensor & stress)
{
  _elastic_strain[_qp] = mechanical_strain;

  if (_plasticity_model)
    _plasticity_model->updateState(stress, _elastic_strain[_qp]);
  else
    stress = computeStress(_elastic_strain[_qp]);
}
