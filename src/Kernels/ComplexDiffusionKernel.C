#include "ComplexDiffusionKernel.h"

#include <cmath>

registerMooseObject("reproductionApp", ComplexDiffusionKernel);

InputParameters
ComplexDiffusionKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("temperature", "Coupled temperature variable");
  params.addParam<Real>("R", 8.314, "Gas constant");
  return params;
}

ComplexDiffusionKernel::ComplexDiffusionKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _T(coupledValue("temperature")),
    _grad_T(coupledGradient("temperature")),
    _T_var(coupled("temperature")),
    _R(getParam<Real>("R"))
{
}

Real
ComplexDiffusionKernel::computeD(Real u, Real T) const
{
  const Real inv_T = 1.0 / T;
  const Real log_term = std::log10(1.0 + 2.0 / u);
  const Real exponent = -9.386 - 4260.0 * inv_T + 0.0012 * T * u + 0.00075 * T * log_term;
  return std::pow(10.0, exponent);
}

Real
ComplexDiffusionKernel::computeDdu(Real u, Real T, Real D) const
{
  const Real ln10 = std::log(10.0);
  return D * (0.0012 * T * ln10 - 0.0015 * T / (u * (u + 2.0)));
}

Real
ComplexDiffusionKernel::computeDdT(Real u, Real T, Real D) const
{
  const Real ln10 = std::log(10.0);
  const Real inv_T = 1.0 / T;
  const Real log_term = std::log10(1.0 + 2.0 / u);
  const Real d_exponent_dT = 4260.0 * inv_T * inv_T + 0.0012 * u + 0.00075 * log_term;
  return ln10 * D * d_exponent_dT;
}

Real
ComplexDiffusionKernel::computeF(Real u) const
{
  const Real denom1 = 1.0 - 3.0 * u;
  const Real denom2 = 1.0 - 2.0 * u;
  return (2.0 + u) / (2.0 * denom1 * denom2);
}

Real
ComplexDiffusionKernel::computeFdU(Real u) const
{
  const Real denom1 = 1.0 - 3.0 * u;
  const Real denom2 = 1.0 - 2.0 * u;
  const Real d = 2.0 * denom1 * denom2;
  const Real n = 2.0 + u;
  const Real dprime = 2.0 * (-3.0 * denom2 - 2.0 * denom1);
  return (d - n * dprime) / (d * d);
}

Real
ComplexDiffusionKernel::computeQStar(Real exp_term) const
{
  return -1380.8 - 134435.5 * exp_term;
}

Real
ComplexDiffusionKernel::computeQStardU(Real exp_term) const
{
  return 134435.5 / 0.0261 * exp_term;
}

Real
ComplexDiffusionKernel::computeTempCoef(Real u, Real T, Real F, Real Q_star) const
{
  const Real inv_T = 1.0 / T;
  const Real inv_RT2 = inv_T * inv_T / _R;
  return u * Q_star * inv_RT2 / F;
}

Real
ComplexDiffusionKernel::computeTempCoefdu(
    Real u, Real T, Real F, Real F_du, Real Q_star, Real Q_star_du) const
{
  const Real inv_T = 1.0 / T;
  const Real inv_RT2 = inv_T * inv_T / _R;

  const Real term1 = (Q_star + u * Q_star_du) / F;
  const Real term2 = (u * Q_star * F_du) / (F * F);
  return inv_RT2 * (term1 - term2);
}

Real
ComplexDiffusionKernel::computeTempCoefdT(Real temp_coef, Real T) const
{
  return -2.0 * temp_coef / T;
}

Real
ComplexDiffusionKernel::computeQpResidual()
{
  const Real T = _T[_qp];
  const Real u = _u[_qp];

  const Real D = computeD(u, T);
  const Real F = computeF(u);
  const Real exp_term = std::exp(-u / 0.0261);
  const Real Q_star = computeQStar(exp_term);
  const Real temp_coef = computeTempCoef(u, T, F, Q_star);

  const RealGradient flux = D * (_grad_u[_qp] + temp_coef * _grad_T[_qp]);
  return _grad_test[_i][_qp] * flux;
}

Real
ComplexDiffusionKernel::computeQpJacobian()
{
  const Real T = _T[_qp];
  const Real u = _u[_qp];

  const Real D = computeD(u, T);
  const Real D_du = computeDdu(u, T, D);

  const Real F = computeF(u);
  const Real F_du = computeFdU(u);

  const Real exp_term = std::exp(-u / 0.0261);
  const Real Q_star = computeQStar(exp_term);
  const Real Q_star_du = computeQStardU(exp_term);

  const Real temp_coef = computeTempCoef(u, T, F, Q_star);
  const Real temp_coef_du = computeTempCoefdu(u, T, F, F_du, Q_star, Q_star_du);

  const RealGradient grad_v = _grad_u[_qp] + temp_coef * _grad_T[_qp];
  const RealGradient term_dD = (D_du * _phi[_j][_qp]) * grad_v;
  const RealGradient term_D = D * (_grad_phi[_j][_qp] + (temp_coef_du * _phi[_j][_qp]) * _grad_T[_qp]);

  return _grad_test[_i][_qp] * (term_dD + term_D);
}

Real
ComplexDiffusionKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar != _T_var)
    return 0.0;

  const Real T = _T[_qp];
  const Real u = _u[_qp];

  const Real D = computeD(u, T);
  const Real D_dT = computeDdT(u, T, D);

  const Real F = computeF(u);
  const Real exp_term = std::exp(-u / 0.0261);
  const Real Q_star = computeQStar(exp_term);

  const Real temp_coef = computeTempCoef(u, T, F, Q_star);
  const Real temp_coef_dT = computeTempCoefdT(temp_coef, T);

  const RealGradient grad_v = _grad_u[_qp] + temp_coef * _grad_T[_qp];
  const RealGradient term_dD = (D_dT * _phi[_j][_qp]) * grad_v;
  const RealGradient term_D = D * ((temp_coef_dT * _phi[_j][_qp]) * _grad_T[_qp] + temp_coef * _grad_phi[_j][_qp]);

  return _grad_test[_i][_qp] * (term_dD + term_D);
}
