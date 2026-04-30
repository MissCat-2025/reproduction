#include "GapFluxModelRossStouteContactConduction.h"

#include <cmath>

registerMooseObject("reproductionApp", GapFluxModelRossStouteContactConduction);

InputParameters
GapFluxModelRossStouteContactConduction::validParams()
{
  InputParameters params = GapFluxModelBase::validParams();
  params.addClassDescription(
      "Gap flux model for solid-solid contact conduction using a Ross-Stoute style relation.");
  params.addRequiredCoupledVar("temperature", "The name of the temperature variable");
  params.addRequiredCoupledVar("contact_pressure",
                               "The normal contact pressure lower variable");
  params.addRequiredParam<MaterialPropertyName>(
      "primary_conductivity", "Thermal conductivity of the primary surface material");
  params.addRequiredParam<MaterialPropertyName>(
      "secondary_conductivity", "Thermal conductivity of the secondary surface material");
  params.addRequiredParam<MaterialPropertyName>(
      "primary_hardness", "Meyer hardness of the primary surface material");
  params.addRequiredParam<MaterialPropertyName>(
      "secondary_hardness", "Meyer hardness of the secondary surface material");
  params.addRangeCheckedParam<Real>(
      "fuel_roughness",
      1e-12,
      "fuel_roughness>0",
      "Fuel-side surface roughness [m]");
  params.addRangeCheckedParam<Real>(
      "clad_roughness",
      1e-12,
      "clad_roughness>0",
      "Cladding-side surface roughness [m]");

  params.addParamNamesToGroup("contact_pressure", "Gap contact flux");
  params.addParamNamesToGroup(
      "primary_conductivity secondary_conductivity primary_hardness secondary_hardness",
      "Gap contact flux");
  params.addParamNamesToGroup("fuel_roughness clad_roughness", "Gap contact flux");

  return params;
}

GapFluxModelRossStouteContactConduction::GapFluxModelRossStouteContactConduction(
    const InputParameters & parameters)
  : GapFluxModelBase(parameters),
    _primary_T(adCoupledNeighborValue("temperature")),
    _secondary_T(adCoupledValue("temperature")),
    _contact_pressure(adCoupledLowerValue("contact_pressure")),
    _primary_conductivity(getNeighborADMaterialProperty<Real>("primary_conductivity")),
    _secondary_conductivity(getADMaterialProperty<Real>("secondary_conductivity")),
    _primary_hardness(getNeighborADMaterialProperty<Real>("primary_hardness")),
    _secondary_hardness(getADMaterialProperty<Real>("secondary_hardness")),
    _fuel_roughness(getParam<Real>("fuel_roughness")),
    _clad_roughness(getParam<Real>("clad_roughness")),
    _roughness_norm(std::sqrt(_fuel_roughness * _fuel_roughness +
                              _clad_roughness * _clad_roughness)),
    _roughness_factor(std::exp(5.738 - 0.528 * std::log(3.937e7 * _fuel_roughness)))
{
}

ADReal
GapFluxModelRossStouteContactConduction::computeFlux() const
{
  using std::sqrt;

  if (_contact_pressure[_qp] <= 0.0)
    return 0.0;

  const ADReal k_sum = _primary_conductivity[_qp] + _secondary_conductivity[_qp];
  const ADReal k_m = 2.0 * _primary_conductivity[_qp] * _secondary_conductivity[_qp] / k_sum;

  const ADReal P_rel = _contact_pressure[_qp] / _secondary_hardness[_qp];

  ADReal R_mult = 2.9;
  if (P_rel <= 0.0087)
    R_mult = 333.3 * P_rel;

  ADReal h_s = 0.0;
  if (P_rel < 9.0e-6)
    h_s = 0.4166 * k_m * sqrt(P_rel) / (_roughness_norm * _roughness_factor);
  else if (P_rel < 0.003)
    h_s = 0.00125 * k_m / (_roughness_norm * _roughness_factor);
  else
    h_s = 0.4166 * k_m * P_rel * R_mult / (_roughness_norm * _roughness_factor);

  return h_s * (_primary_T[_qp] - _secondary_T[_qp]);
}