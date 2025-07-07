//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "BaseNameInterface.h"
// ✅ 前向声明 - 因为只需要指针
class SmallDeformationElasticityModelMod;// 前向声明
class SmallDeformationPlasticityModelMod;// 前向声明

/**
 * ComputeSmallDeformationStressMod computes the stress under small-strain assumptions
 */
class ComputeSmallDeformationStressMod : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ComputeSmallDeformationStressMod(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// The elasticity model
  SmallDeformationElasticityModelMod * _elasticity_model; // 只需要指针，前向声明足够

  /// The plasticity model
  SmallDeformationPlasticityModelMod * _plasticity_model; // 只需要指针，前向声明足够

  /// The mechanical strain excluding eigen strains from the total strain
  const ADMaterialProperty<RankTwoTensor> & _mechanical_strain;

  /// The stress
  ADMaterialProperty<RankTwoTensor> & _stress;
}; 