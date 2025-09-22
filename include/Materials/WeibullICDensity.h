#pragma once

#include "InitialCondition.h"
#include "MooseRandom.h"

class WeibullICDensity : public InitialCondition
{
public:
  static InputParameters validParams();
  WeibullICDensity(const InputParameters & parameters);

  virtual Real value(const Point & /*p*/) override;

protected:
  const Real _scale;
  const Real _shape;
  const Real _location;
  const unsigned int _seed;
};