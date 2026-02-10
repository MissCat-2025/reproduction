#pragma once

#include "Kernel.h"

class ComplexDiffusionKernel : public Kernel
{
public:
  static InputParameters validParams();

  ComplexDiffusionKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  Real computeD(Real u, Real T) const;
  Real computeDdu(Real u, Real T, Real D) const;
  Real computeDdT(Real u, Real T, Real D) const;

  Real computeF(Real u) const;
  Real computeFdU(Real u) const;

  Real computeQStar(Real exp_term) const;
  Real computeQStardU(Real exp_term) const;

  Real computeTempCoef(Real u, Real T, Real F, Real Q_star) const;
  Real computeTempCoefdu(Real u, Real T, Real F, Real F_du, Real Q_star, Real Q_star_du) const;
  Real computeTempCoefdT(Real temp_coef, Real T) const;

  const VariableValue & _T;
  const VariableGradient & _grad_T;
  const unsigned int _T_var;

  const Real _R;
};
