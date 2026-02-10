#pragma once

#include "EqualValueBoundaryConstraint.h"

class EqualValueBoundaryConstraintDegraded : public EqualValueBoundaryConstraint
{
public:
  static InputParameters validParams();

  EqualValueBoundaryConstraintDegraded(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::ConstraintType type) override;
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

private:
  Real degradationFactorPrimary() const;
  Real degradationFactorSecondary() const;

  const VariableValue * _degradation_primary;
  const VariableValue * _degradation_secondary;
  const bool _use_degradation;
  const Real _degradation_power;
  const Real _degradation_a2;
  const Real _degradation_a3;
  const Real _degradation_eta;
  const bool _use_rational_degradation;
  const VariableValue * _sigma0_primary;
  const VariableValue * _sigma0_secondary;
  const Real _youngs_modulus;
  const Real _fracture_energy;
  const Real _length_scale;
};
