#include "GapFluxModelRossStouteGasConduction.h"

#include "Function.h"
#include "MooseEnum.h"

#include <algorithm>
#include <cmath>

namespace
{
constexpr Real cgs_conductivity_per_si = 0.002390057361376673;
constexpr Real pa_to_dyne_per_cm2 = 10.0;
constexpr Real cm_to_m = 1e-2;
constexpr Real ross_stoute_constant = 5756.0;
constexpr Real helium_molar_mass_g_per_mol = 4.002602;
constexpr Real xenon_molar_mass_g_per_mol = 131.293;
constexpr Real roughness_coefficient = 1.5;
constexpr Real pressure_floor_pa = 1e-12;
} // namespace

registerMooseObject("HeatTransferApp", GapFluxModelRossStouteGasConduction);

InputParameters
GapFluxModelRossStouteGasConduction::validParams()
{
  InputParameters params = GapFluxModelBase::validParams();
  params.addClassDescription(
    "Gap flux model for gas conduction using a Ross-Stoute style conductance.");
  params.addRequiredCoupledVar("temperature", "The name of the temperature variable");
  params.addRequiredParam<FunctionName>(
      "gas_thermal_conductivity",
      "Gas thermal conductivity as a function of the average gap temperature");
  params.addRequiredParam<FunctionName>(
      "gas_pressure", "Gap gas pressure as a function of time [Pa]");
  MooseEnum gas_models("HELIUM_ONLY HELIUM_XENON", "HELIUM_ONLY");
  params.addParam<MooseEnum>(
      "gas_mixture_mode",
      gas_models,
      "Choose the gas composition model: pure helium or helium-xenon");
  params.addRangeCheckedParam<Real>(
      "xenon_mole_fraction",
      0.0,
      "xenon_mole_fraction>=0 & xenon_mole_fraction<=1",
      "Xenon mole fraction in the gap gas mixture. Required when gas_mixture_mode = "
      "HELIUM_XENON.");
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
  params.addRangeCheckedParam<Real>(
      "fuel_jump_distance",
      0.0,
      "fuel_jump_distance>=0",
      "Fuel-side jump distance [m]");
  params.addRangeCheckedParam<Real>(
      "clad_jump_distance",
      0.0,
      "clad_jump_distance>=0",
      "Cladding-side jump distance [m]");

  params.addParamNamesToGroup("gas_thermal_conductivity gas_pressure gas_mixture_mode",
                              "Gap gas flux");
  params.addParamNamesToGroup(
      "xenon_mole_fraction fuel_roughness clad_roughness fuel_jump_distance clad_jump_distance",
      "Gap gas flux");

  return params;
}

GapFluxModelRossStouteGasConduction::GapFluxModelRossStouteGasConduction(
    const InputParameters & parameters)
  : GapFluxModelBase(parameters),
    _primary_T(adCoupledNeighborValue("temperature")),
    _secondary_T(adCoupledValue("temperature")),
    _gas_thermal_conductivity_fn(getFunction("gas_thermal_conductivity")),
    _gas_pressure_fn(getFunction("gas_pressure")),
    _gas_composition_mode(getParam<MooseEnum>("gas_mixture_mode").getEnum<GasCompositionMode>()),
    _xenon_mole_fraction(getParam<Real>("xenon_mole_fraction")),
    // M_mix uses the linear mole-fraction average, while the Kennard jump term uses
    // sum_i f_i / M_i exactly as written in the correlation.
    _mixture_molar_mass(_gas_composition_mode == GasCompositionMode::HELIUM_XENON
                            ? (1.0 - _xenon_mole_fraction) * helium_molar_mass_g_per_mol +
                                  _xenon_mole_fraction * xenon_molar_mass_g_per_mol
                            : helium_molar_mass_g_per_mol),
    _inverse_molar_mass_sum(_gas_composition_mode == GasCompositionMode::HELIUM_XENON
                                ? (1.0 - _xenon_mole_fraction) / helium_molar_mass_g_per_mol +
                                      _xenon_mole_fraction / xenon_molar_mass_g_per_mol
                                : 1.0 / helium_molar_mass_g_per_mol),
    _fuel_roughness(getParam<Real>("fuel_roughness")),
    _clad_roughness(getParam<Real>("clad_roughness")),
    _fuel_jump_distance(getParam<Real>("fuel_jump_distance")),
    _clad_jump_distance(getParam<Real>("clad_jump_distance")),
    _roughness_sum(_fuel_roughness + _clad_roughness)
{
  if (_gas_composition_mode == GasCompositionMode::HELIUM_XENON &&
      !parameters.isParamSetByUser("xenon_mole_fraction"))
    paramError("xenon_mole_fraction",
               "The xenon_mole_fraction parameter must be set when gas_mixture_mode = "
               "HELIUM_XENON.");
}

ADReal
GapFluxModelRossStouteGasConduction::evaluateTemperatureFunction(const Function & fn,
                                                                 const ADReal & T) const
{
  ADReal value = fn.value(T.value());
  value.derivatives() = T.derivatives() * fn.timeDerivative(T.value());

  return value;
}

ADReal
GapFluxModelRossStouteGasConduction::computeAccommodationCoefficient(
    const ADReal & gas_temperature) const
{
  const ADReal a_he = 0.425 - 2.3e-4 * gas_temperature;

  if (_gas_composition_mode == GasCompositionMode::HELIUM_ONLY)
    return a_he;

  const ADReal a_xe = 0.749 - 2.5e-4 * gas_temperature;

  return a_he + (a_xe - a_he) * (_mixture_molar_mass - helium_molar_mass_g_per_mol) /
                       (xenon_molar_mass_g_per_mol - helium_molar_mass_g_per_mol);
}

ADReal
GapFluxModelRossStouteGasConduction::computeJumpDistance(const ADReal & gas_temperature,
                                                          const ADReal & gas_conductivity,
                                                          Real gas_pressure) const
{
  using std::sqrt;

  const ADReal accommodation = computeAccommodationCoefficient(gas_temperature);
  const Real pressure_cgs = std::max(gas_pressure, pressure_floor_pa) * pa_to_dyne_per_cm2;
  const ADReal conductivity_cgs = gas_conductivity * cgs_conductivity_per_si;

  // Ross-Stoute / Kennard correlation is published in cgs units.
  // Conversion used here:
  // - k_g: W/(m K) -> cal/(cm s K) by multiplying by 0.002390057361376673
  // - P: Pa -> dyne/cm^2 by multiplying by 10
  // - g: cm -> m by multiplying by 1e-2
  // The term (sum_i f_i / M_i)^(-1/2) is taken directly from the paper, with M_i in g/mol.
  const ADReal jump_distance_cgs = ross_stoute_constant * ((2.0 - accommodation) / accommodation) *
                                   (conductivity_cgs * sqrt(gas_temperature) / pressure_cgs) *
                                   (1.0 / sqrt(_inverse_molar_mass_sum));

  return cm_to_m * jump_distance_cgs;
}

ADReal
GapFluxModelRossStouteGasConduction::computeFlux() const
{
  const ADReal gas_temperature = 0.5 * (_primary_T[_qp] + _secondary_T[_qp]);
  const ADReal gas_conductivity =
      evaluateTemperatureFunction(_gas_thermal_conductivity_fn, gas_temperature);
  const Real gas_pressure = _gas_pressure_fn.value(_t);

  const ADReal jump_distance =
      computeJumpDistance(gas_temperature, gas_conductivity, gas_pressure);

  const ADReal denominator = _adjusted_length + roughness_coefficient * _roughness_sum +
                             _fuel_jump_distance + _clad_jump_distance + jump_distance;

  return gas_conductivity * (_primary_T[_qp] - _secondary_T[_qp]) / denominator;
}