//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeCreepPlasticityDeformationStress.h"
#include "ElasticityModel.h"
#include "CreepModel.h"

registerMooseObject("reproductionApp", ComputeCreepPlasticityDeformationStress);

InputParameters
ComputeCreepPlasticityDeformationStress::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The stress calculator given an elasticity model and a plasticity "
                             "model. Small deformation is assumed.");

  params.addRequiredParam<MaterialName>("elasticity_model",
                                        "Name of the elastic stress-strain constitutive model");
  params.addParam<MaterialName>("creep_model",
                                        "Name of the creep stress-strain constitutive model");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ComputeCreepPlasticityDeformationStress::ComputeCreepPlasticityDeformationStress(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _mechanical_strain(getADMaterialProperty<RankTwoTensor>(prependBaseName("mechanical_strain"))),
    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress")))
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
ComputeCreepPlasticityDeformationStress::initialSetup()
{
  _elasticity_model =
      dynamic_cast<ElasticityModel *>(&getMaterial("elasticity_model"));
  if (!_elasticity_model)
    paramError("elasticity_model",
               "Elasticity model " + getParam<MaterialName>("elasticity_model") +
                   " is not compatible with ComputeCreepPlasticityDeformationStress");
  _creep_model =
      isParamValid("creep_model")
          ? dynamic_cast<CreepModel *>(&getMaterial("creep_model"))
          : nullptr;
  if (_creep_model)
    _elasticity_model->setCreepModel(_creep_model);
}

void
ComputeCreepPlasticityDeformationStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
}

void
ComputeCreepPlasticityDeformationStress::computeQpProperties()
{
  _elasticity_model->setQp(_qp);
  _elasticity_model->updateState(_mechanical_strain[_qp], _stress[_qp]);
} 