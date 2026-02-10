//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "UO2PowerLawCreepStressUpdateNoDegradation.h"

registerMooseObject("reproductionApp", UO2PowerLawCreepStressUpdateNoDegradation);
registerMooseObject("reproductionApp", ADUO2PowerLawCreepStressUpdateNoDegradation);

template <bool is_ad>
InputParameters
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription("UO2蠕变模型（径向返回），不包含退化与损伤。");

  params.addParam<Real>("gas_constant", 8.314, "气体常数 (J/mol-K)");

  params.addRequiredCoupledVar("temperature", "温度 (K)");
  params.addRequiredCoupledVar("oxygen_ratio", "氧化学计量比");

  params.addRequiredParam<Real>("fission_rate", "裂变率密度 (fissions/m^3-s)");
  params.addParam<Real>("theoretical_density", 95.0, "理论密度百分比");
  params.addParam<Real>("grain_size", 10.0, "晶粒尺寸 (微米)");

  params.addParam<Real>("a1", 0.3919, "系数a1");
  params.addParam<Real>("a2", 1.31e-19, "系数a2");
  params.addParam<Real>("a3", 87.7, "系数a3");
  params.addParam<Real>("a5", 2.0391e-25, "系数a5");
  params.addParam<Real>("a6", 90.5, "系数a6");
  params.addParam<Real>("a8", 3.7226e-35, "系数a8");
  params.addParam<Real>("Q3", 21759.0, "激活能Q3 (J/mol)");

  params.addParam<bool>("use_transition_stress", false, "是否使用应力转变点");

  return params;
}

template <bool is_ad>
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::UO2PowerLawCreepStressUpdateNoDegradationTempl(
    const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _gas_constant(this->template getParam<Real>("gas_constant")),
    _temperature(this->template coupledGenericValue<is_ad>("temperature")),
    _oxygen_ratio(this->template coupledGenericValue<is_ad>("oxygen_ratio")),
    _fission_rate(this->template getParam<Real>("fission_rate")),
    _theoretical_density(this->template getParam<Real>("theoretical_density")),
    _grain_size(this->template getParam<Real>("grain_size")),
    _Q3(this->template getParam<Real>("Q3")),
    _a1(this->template getParam<Real>("a1")),
    _a2(this->template getParam<Real>("a2")),
    _a3(this->template getParam<Real>("a3")),
    _a5(this->template getParam<Real>("a5")),
    _a6(this->template getParam<Real>("a6")),
    _a8(this->template getParam<Real>("a8")),
    _use_transition_stress(this->template getParam<bool>("use_transition_stress")),
    _exp_Q1(0.0),
    _exp_Q2(0.0),
    _exp_Q3(0.0),
    _sigma_trans(0.0),
    _density_term1(0.0),
    _density_term2(0.0),
    _fission_term(0.0),
    _creep_rate(this->template declareProperty<RankTwoTensor>("creep_rate")),
    _effective_creep(this->template declareProperty<Real>("effective_creep"))
{
}

template <bool is_ad>
void
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);

  const GenericReal<is_ad> T = _temperature[_qp];
  const GenericReal<is_ad> RT = _gas_constant * T;
  const GenericReal<is_ad> inv_RT = 1.0 / RT;

  GenericReal<is_ad> x = _oxygen_ratio[_qp];

  const GenericReal<is_ad> log_x = std::log10(x);
  const GenericReal<is_ad> exp_common = std::exp(-20.0 / log_x - 8.0);
  const GenericReal<is_ad> denom = 1.0 / (exp_common + 1.0);

  const GenericReal<is_ad> Q1 = 74829.0 * denom + 301762.0;
  const GenericReal<is_ad> Q2 = 83143.0 * denom + 469191.0;

  _exp_Q1 = std::exp(-Q1 * inv_RT);
  _exp_Q2 = std::exp(-Q2 * inv_RT);
  _exp_Q3 = std::exp(-_Q3 * inv_RT);

  _sigma_trans = 1.6547e7 * std::pow(_grain_size, 0.5714);
  _density_term1 = 1.0 / ((_theoretical_density - _a3) * _grain_size * _grain_size);
  _density_term2 = 1.0 / (_theoretical_density - _a6);
  _fission_term = _a1 + _a2 * _fission_rate;
}

template <bool is_ad>
template <typename ScalarType>
ScalarType
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::computeResidualInternal(
    const GenericReal<is_ad> & effective_trial_stress, const ScalarType & scalar)
{
  const ScalarType effective_stress = effective_trial_stress - _three_shear_modulus * scalar;
  // if (MetaPhysicL::raw_value(effective_stress) <= 0.0)
  //   return -scalar;

  ScalarType creep_th1 = 0.0;
  ScalarType creep_th2 = 0.0;

  if (_use_transition_stress)
  {
    if (MetaPhysicL::raw_value(effective_stress) < _sigma_trans)
    {
      creep_th1 = _fission_term * _density_term1 * effective_stress * _exp_Q1;
      creep_th2 = 0.0;
    }
    else
    {
      creep_th1 = _fission_term * _density_term1 * _sigma_trans * _exp_Q1;
      creep_th2 = _a5 * _density_term2 * std::pow(effective_stress, 4.5) * _exp_Q2;
    }
  }
  else
  {
    creep_th1 = _fission_term * _density_term1 * effective_stress * _exp_Q1;
    creep_th2 = _a5 * _density_term2 * std::pow(effective_stress, 4.5) * _exp_Q2;
  }

  const ScalarType creep_ir = _a8 * _fission_rate * effective_stress * _exp_Q3;
  const ScalarType scalar_rate = (creep_th1 + creep_th2 + creep_ir);

  return scalar_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> effective_stress = effective_trial_stress - _three_shear_modulus * scalar;
  if (MetaPhysicL::raw_value(effective_stress) <= 0.0)
    return -1.0;

  GenericReal<is_ad> d_creep_th1_d_stress = 0.0;
  GenericReal<is_ad> d_creep_th2_d_stress = 0.0;

  if (_use_transition_stress)
  {
    if (MetaPhysicL::raw_value(effective_stress) < _sigma_trans)
    {
      d_creep_th1_d_stress = _fission_term * _density_term1 * _exp_Q1;
      d_creep_th2_d_stress = 0.0;
    }
    else
    {
      d_creep_th1_d_stress = 0.0;
      d_creep_th2_d_stress =
          4.5 * _a5 * _density_term2 * std::pow(effective_stress, 3.5) * _exp_Q2;
    }
  }
  else
  {
    d_creep_th1_d_stress = _fission_term * _density_term1 * _exp_Q1;
    d_creep_th2_d_stress =
        4.5 * _a5 * _density_term2 * std::pow(effective_stress, 3.5) * _exp_Q2;
  }

  const GenericReal<is_ad> d_creep_ir_d_stress = _a8 * _fission_rate * _exp_Q3;
  const GenericReal<is_ad> d_scalar_rate_d_stress =
      d_creep_th1_d_stress + d_creep_th2_d_stress + d_creep_ir_d_stress;

  const GenericReal<is_ad> d_scalar_rate_d_scalar = d_scalar_rate_d_stress * (-_three_shear_modulus);
  return d_scalar_rate_d_scalar * _dt - 1.0;
}

template <bool is_ad>
Real
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::computeStrainEnergyRateDensity(
    const GenericMaterialProperty<RankTwoTensor, is_ad> & stress,
    const GenericMaterialProperty<RankTwoTensor, is_ad> & strain_rate)
{
  return MetaPhysicL::raw_value(stress[_qp].doubleContraction(strain_rate[_qp]));
}

template <bool is_ad>
void
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;

  RankTwoTensor increment;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      increment(i, j) = MetaPhysicL::raw_value(plastic_strain_increment(i, j));

  _creep_rate[_qp] = increment / _dt;
  _effective_creep[_qp] =
      MetaPhysicL::raw_value(this->effectiveInelasticStrainIncrement()) / _dt;
}

template <bool is_ad>
void
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
  _creep_rate[_qp].zero();
  _effective_creep[_qp] = 0.0;
}

template <bool is_ad>
bool
UO2PowerLawCreepStressUpdateNoDegradationTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->_use_substepping != RadialReturnStressUpdateTempl<is_ad>::SubsteppingType::NONE;
}

template class UO2PowerLawCreepStressUpdateNoDegradationTempl<false>;
template class UO2PowerLawCreepStressUpdateNoDegradationTempl<true>;
template Real UO2PowerLawCreepStressUpdateNoDegradationTempl<false>::computeResidualInternal<Real>(
    const Real &, const Real &);
template ADReal
UO2PowerLawCreepStressUpdateNoDegradationTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                                      const ADReal &);
template ChainedReal
UO2PowerLawCreepStressUpdateNoDegradationTempl<false>::computeResidualInternal<ChainedReal>(
    const Real &, const ChainedReal &);
template ChainedADReal
UO2PowerLawCreepStressUpdateNoDegradationTempl<true>::computeResidualInternal<ChainedADReal>(
    const ADReal &, const ChainedADReal &);
