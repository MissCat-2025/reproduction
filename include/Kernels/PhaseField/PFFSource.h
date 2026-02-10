#pragma once

#include "Kernel.h"

class PFFSource : public Kernel
{
public:
  static InputParameters validParams();

  PFFSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  const MaterialProperty<Real> & _dpsi_dd;
  const MaterialProperty<Real> & _d2psi_dd2;
};
