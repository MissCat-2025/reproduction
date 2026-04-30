#pragma once

#include "GapFluxModelBase.h"

/**
 * Gap flux model for solid-solid contact conduction using a Ross-Stoute style relation from the
 * PCI literature.
 */
class GapFluxModelRossStouteContactConduction : public GapFluxModelBase
{
public:
  static InputParameters validParams();

  GapFluxModelRossStouteContactConduction(const InputParameters & parameters);

  ADReal computeFlux() const override;

protected:
  /// Primary surface temperature
  const ADVariableValue & _primary_T;
  /// Secondary surface temperature
  const ADVariableValue & _secondary_T;
  /// Contact pressure (lower mortar variable)
  const ADVariableValue & _contact_pressure;

  /// Thermal conductivity of the primary surface material
  const ADMaterialProperty<Real> & _primary_conductivity;
  /// Thermal conductivity of the secondary surface material
  const ADMaterialProperty<Real> & _secondary_conductivity;
  /// Meyer hardness of the primary surface material
  const ADMaterialProperty<Real> & _primary_hardness;
  /// Meyer hardness of the secondary surface material
  const ADMaterialProperty<Real> & _secondary_hardness;

  /// Fuel-side surface roughness [m]
  const Real _fuel_roughness;
  /// Cladding-side surface roughness [m]
  const Real _clad_roughness;
  /// Roughness norm used in the denominator
  const Real _roughness_norm;
  /// Surface roughness factor from the Ross-Stoute style correlation
  const Real _roughness_factor;
};