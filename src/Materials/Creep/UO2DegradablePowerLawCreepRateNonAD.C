//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "UO2DegradablePowerLawCreepRateNonAD.h"

#include "RankTwoScalarTools.h"

registerMooseObject("reproductionApp", UO2DegradablePowerLawCreepRateNonAD);

InputParameters
UO2DegradablePowerLawCreepRateNonAD::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("UO2蠕变率模型（显式，用旧应力），可选相场退化函数。");

  params.addRequiredCoupledVar("temperature", "温度 (K)");
  params.addParam<Real>("oxygen_ratio_value", 1.00001, "氧化学计量比（未耦合时使用）");
  params.addCoupledVar("oxygen_ratio", "氧化学计量比（可选耦合变量）");

  params.addRequiredParam<Real>("fission_rate", "裂变率密度 (fissions/m^3-s)");
  params.addParam<Real>("theoretical_density", 95.0, "理论密度百分比");
  params.addParam<Real>("grain_size", 10.0, "晶粒尺寸 (微米)");
  params.addParam<Real>("gas_constant", 8.314, "气体常数 (J/mol-K)");

  params.addParam<Real>("a1", 0.3919, "系数a1");
  params.addParam<Real>("a2", 1.31e-19, "系数a2");
  params.addParam<Real>("a3", 87.7, "系数a3");
  params.addParam<Real>("a6", 90.5, "系数a6");
  params.addParam<Real>("a7", 3.7226e-35, "系数a7");
  params.addParam<Real>("a8", 2.0391e-25, "系数a8");

  params.addParam<bool>("use_transition_stress", false, "是否使用应力转变点");

  params.addParam<MaterialPropertyName>("degradation_function", "g", "退化函数材料属性名");
  params.addParam<bool>("use_stress_degradation", true, "是否对有效应力应用退化");

  return params;
}

UO2DegradablePowerLawCreepRateNonAD::UO2DegradablePowerLawCreepRateNonAD(
    const InputParameters & parameters)
  : Material(parameters),
    _gas_constant(getParam<Real>("gas_constant")),
    _fission_rate(getParam<Real>("fission_rate")),
    _theoretical_density(getParam<Real>("theoretical_density")),
    _grain_size(getParam<Real>("grain_size")),
    _Q3(21759.0),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _a3(getParam<Real>("a3")),
    _a6(getParam<Real>("a6")),
    _a7(getParam<Real>("a7")),
    _a8(getParam<Real>("a8")),
    _use_transition_stress(getParam<bool>("use_transition_stress")),
    _temperature(coupledValue("temperature")),
    _oxygen_ratio(isCoupled("oxygen_ratio") ? &coupledValue("oxygen_ratio") : nullptr),
    _oxygen_ratio_value(getParam<Real>("oxygen_ratio_value")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress")),
    _g(getMaterialProperty<Real>(getParam<MaterialPropertyName>("degradation_function"))),
    _use_stress_degradation(getParam<bool>("use_stress_degradation")),
    _creep_rate(declareProperty<RankTwoTensor>("creep_rate")),
    _effective_creep(declareProperty<Real>("effective_creep"))
{
}

void
UO2DegradablePowerLawCreepRateNonAD::computeQpProperties()
{
  const Real T = _temperature[_qp];
  Real x = _oxygen_ratio ? (*_oxygen_ratio)[_qp] : _oxygen_ratio_value;
  if (x <= 1.0)
    x = 1.0 + 1e-12;

  const RankTwoTensor stress_old = _stress_old[_qp];
  const RankTwoTensor stress_dev = stress_old.deviatoric();
  const Real sigma_vm = RankTwoScalarTools::vonMisesStress(stress_old);

  const Real g = _use_stress_degradation ? _g[_qp] : 1.0;
  const Real stress = g * sigma_vm;

  if (stress <= 0.0 || sigma_vm <= 1e-16 || T <= 0.0)
  {
    _effective_creep[_qp] = 0.0;
    _creep_rate[_qp].zero();
    return;
  }

  const Real RT = _gas_constant * T;
  const Real inv_RT = 1.0 / RT;

  const Real log_x = std::log10(x);
  const Real exp_common = std::exp(-20.0 / log_x - 8.0);
  const Real denom = 1.0 / (exp_common + 1.0);

  const Real Q1 = 74829.0 * denom + 301762.0;
  const Real Q2 = 83143.0 * denom + 469191.0;

  const Real exp_Q1 = std::exp(-Q1 * inv_RT);
  const Real exp_Q2 = std::exp(-Q2 * inv_RT);
  const Real exp_Q3 = std::exp(-_Q3 * inv_RT);

  const Real sigma_trans = 1.6547e7 * std::pow(_grain_size, 0.5714);

  const Real density_term1 = 1.0 / ((_theoretical_density - _a3) * _grain_size * _grain_size);
  const Real density_term2 = 1.0 / (_theoretical_density - _a6);
  const Real fission_term = _a1 + _a2 * _fission_rate;

  Real creep_th1 = 0.0;
  Real creep_th2 = 0.0;
  if (_use_transition_stress)
  {
    if (stress < sigma_trans)
    {
      creep_th1 = fission_term * density_term1 * stress * exp_Q1;
      creep_th2 = 0.0;
    }
    else
    {
      creep_th1 = fission_term * density_term1 * sigma_trans * exp_Q1;
      creep_th2 = _a8 * density_term2 * std::pow(stress, 4.5) * exp_Q2;
    }
  }
  else
  {
    creep_th1 = fission_term * density_term1 * stress * exp_Q1;
    creep_th2 = _a8 * density_term2 * std::pow(stress, 4.5) * exp_Q2;
  }

  const Real creep_ir = _a7 * _fission_rate * stress * exp_Q3;
  const Real scalar_rate = creep_th1 + creep_th2 + creep_ir;

  _effective_creep[_qp] = scalar_rate;

  const RankTwoTensor flow_direction = stress_dev * (1.5 / sigma_vm);
  _creep_rate[_qp] = flow_direction * scalar_rate;
}

