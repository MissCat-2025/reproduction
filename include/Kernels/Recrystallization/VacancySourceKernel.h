#pragma once

#include "Kernel.h"

class VacancySourceKernel : public Kernel
{
public:
  static InputParameters validParams();
  VacancySourceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:
  const VariableValue & _eta;     // coupled order parameter η
  const Real _Pcasc;              // cascade probability per unit volume and time
  const Real _VG;                 // maximum vacancy increment per cascade
  const Real _eta_threshold;      // threshold for η (default 0.8)
  const Real _sign;               // sign to control contribution (default -1)
};