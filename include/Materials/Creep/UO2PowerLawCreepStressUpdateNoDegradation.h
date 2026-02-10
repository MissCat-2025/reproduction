//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "RadialReturnCreepStressUpdateBase.h"

template <bool is_ad>
class UO2PowerLawCreepStressUpdateNoDegradationTempl
  : public RadialReturnCreepStressUpdateBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  UO2PowerLawCreepStressUpdateNoDegradationTempl(const InputParameters & parameters);

  virtual Real computeStrainEnergyRateDensity(
      const GenericMaterialProperty<RankTwoTensor, is_ad> & stress,
      const GenericMaterialProperty<RankTwoTensor, is_ad> & strain_rate) override;

  virtual bool substeppingCapabilityEnabled() override;
  virtual void resetIncrementalMaterialProperties() override;

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;

  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericReal<is_ad>>(effective_trial_stress, scalar);
  }

  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;

  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;

protected:
  virtual GenericChainedReal<is_ad>
  computeResidualAndDerivative(const GenericReal<is_ad> & effective_trial_stress,
                               const GenericChainedReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericChainedReal<is_ad>>(effective_trial_stress, scalar);
  }

  const Real _gas_constant;

  const GenericVariableValue<is_ad> & _temperature;
  const GenericVariableValue<is_ad> & _oxygen_ratio;

  const Real _fission_rate;
  const Real _theoretical_density;
  const Real _grain_size;

  const GenericReal<is_ad> _Q3;

  const Real _a1;
  const Real _a2;
  const Real _a3;
  const Real _a5;
  const Real _a6;
  const Real _a8;

  const bool _use_transition_stress;

  GenericReal<is_ad> _exp_Q1;
  GenericReal<is_ad> _exp_Q2;
  GenericReal<is_ad> _exp_Q3;

  Real _sigma_trans;
  Real _density_term1;
  Real _density_term2;
  Real _fission_term;

  MaterialProperty<RankTwoTensor> & _creep_rate;
  MaterialProperty<Real> & _effective_creep;

  usingTransientInterfaceMembers;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_qp;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_three_shear_modulus;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain_old;

private:
  template <typename ScalarType>
  ScalarType computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                     const ScalarType & scalar);
};

typedef UO2PowerLawCreepStressUpdateNoDegradationTempl<false> UO2PowerLawCreepStressUpdateNoDegradation;
typedef UO2PowerLawCreepStressUpdateNoDegradationTempl<true> ADUO2PowerLawCreepStressUpdateNoDegradation;
