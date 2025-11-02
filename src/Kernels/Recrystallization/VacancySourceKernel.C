#include "VacancySourceKernel.h"
#include "MooseRandom.h"
registerMooseObject("reproductionApp", VacancySourceKernel);

InputParameters
VacancySourceKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("eta", "Order parameter eta used to gate cascades (matrix if eta < threshold)");
  params.addParam<Real>("eta_threshold", 0.8, "Threshold for eta above which Pv = 0 (default 0.8)");
  params.addParam<Real>("Pcasc", 0.0, "Cascade probability per unit volume and time (0..1)");
  params.addParam<Real>("VG", 0.0, "Maximum vacancy concentration increment per cascade");
  params.addParam<Real>("sign", -1.0, "Sign of source contribution (default -1 for RHS source)");
  params.addClassDescription("Stochastic irradiation vacancy source Pv per Millett et al. (2011): Pv = 0 if eta >= threshold or R1 > Pcasc; else Pv = R2 * VG.");
  return params;
}

VacancySourceKernel::VacancySourceKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _eta(coupledValue("eta")),
    _Pcasc(getParam<Real>("Pcasc")),
    _VG(getParam<Real>("VG")),
    _eta_threshold(getParam<Real>("eta_threshold")),
    _sign(getParam<Real>("sign"))
{
}

Real
VacancySourceKernel::computeQpResidual()
{
  const Real eta_qp = _eta[_qp];
  if (eta_qp >= _eta_threshold)
    return 0.0;

  const Real R1 = MooseRandom::rand();
  if (R1 > _Pcasc)
    return 0.0;

  const Real R2 = MooseRandom::rand();
  const Real Pv = R2 * _VG;
  return _sign * _test[_i][_qp] * Pv;
}