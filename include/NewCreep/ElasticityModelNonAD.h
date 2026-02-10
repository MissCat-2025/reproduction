#pragma once

#include "Material.h"
#include "BaseNameInterface.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class ElasticityModelNonAD : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ElasticityModelNonAD(const InputParameters & parameters);

  virtual void setQp(unsigned int qp) { _qp = qp; }

  virtual RankTwoTensor computeStress(const RankTwoTensor & elastic_strain) = 0;
  virtual RankFourTensor computeJacobianMult(const RankTwoTensor & elastic_strain) = 0;
  void resetQpProperties() override {}

protected:
  const bool _compute;
};
