//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "BaseNameInterface.h"

class ElasticityModel;
class CreepModel;

/**
 * ComputeCreepPlasticityDeformationStress computes the stress under small-strain assumptions
 */
class ComputeCreepPlasticityDeformationStress : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ComputeCreepPlasticityDeformationStress(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// The elasticity model
  ElasticityModel * _elasticity_model;

  /// The creep model
  CreepModel * _creep_model;

  /// The mechanical strain excluding eigen strains from the total strain
  const ADMaterialProperty<RankTwoTensor> & _mechanical_strain;

  /// The stress
  ADMaterialProperty<RankTwoTensor> & _stress;
}; 