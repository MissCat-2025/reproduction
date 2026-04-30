#pragma once

#include "GapFluxModelBase.h"

class Function;

enum class GasCompositionMode
{
  HELIUM_ONLY,
  HELIUM_XENON
};

/**
 * Gap flux model for gas conduction using a Ross-Stoute style conductance.
 */
class GapFluxModelRossStouteGasConduction : public GapFluxModelBase
{
public:
  static InputParameters validParams();

  GapFluxModelRossStouteGasConduction(const InputParameters & parameters);

  ADReal computeFlux() const override;

protected:
  ADReal evaluateTemperatureFunction(const Function & fn, const ADReal & T) const;
  ADReal computeAccommodationCoefficient(const ADReal & gas_temperature) const;
  ADReal computeJumpDistance(const ADReal & gas_temperature,
                             const ADReal & gas_conductivity,
                             Real gas_pressure) const;

  /// Primary surface temperature
  const ADVariableValue & _primary_T;
  /// Secondary surface temperature
  const ADVariableValue & _secondary_T;

  /// Gas thermal conductivity as a function of the average gap temperature
  const Function & _gas_thermal_conductivity_fn;
  /// Gas pressure as a function of time [Pa]
  const Function & _gas_pressure_fn;

  /// Pure helium or helium/xenon accommodation model
  const GasCompositionMode _gas_composition_mode;

  /// Xenon mole fraction in the gap gas mixture
  const Real _xenon_mole_fraction;

  /// Mixture molar mass used in the accommodation coefficient interpolation [g/mol]
  const Real _mixture_molar_mass;

  /// Sum of f_i / M_i used in the Kennard jump-distance term [(mol/g)]
  const Real _inverse_molar_mass_sum;

  /// Fuel-side surface roughness [m]
  const Real _fuel_roughness;
  /// Cladding-side surface roughness [m]
  const Real _clad_roughness;
  /// Fuel-side jump distance [m]
  const Real _fuel_jump_distance;
  /// Cladding-side jump distance [m]
  const Real _clad_jump_distance;

  /// Sum of the roughness terms used in the denominator [m]
  const Real _roughness_sum;
};