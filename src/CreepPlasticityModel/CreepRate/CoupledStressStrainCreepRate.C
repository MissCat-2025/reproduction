//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "CoupledStressStrainCreepRate.h"

registerMooseObject("reproductionApp", CoupledStressStrainCreepRate);

InputParameters
CoupledStressStrainCreepRate::validParams()
{
  InputParameters params = J2CreepPlasticity::validParams();
  params.addClassDescription("Coupled stress-strain creep rate model: "
                            "dot_epsilon_c_eq = A1*sigma_eq^n1*epsilon_c_eq^m1 + "
                            "A2*sigma_eq^n2*epsilon_c_eq^m2 + A3*sigma_eq^n3*epsilon_c_eq^m3");
  
  // 材料常数 - 基于论文中的316H不锈钢在550°C的参数
  params.addParam<Real>("A1", 4.33e-32, "第一项系数");
  params.addParam<Real>("n1", 10.08, "第一项应力指数");
  params.addParam<Real>("m1", -0.633, "第一项应变指数");
  
  params.addParam<Real>("A2", 3.6e-25, "第二项系数"); 
  params.addParam<Real>("n2", 8.9, "第二项应力指数");
  params.addParam<Real>("m2", 1.7, "第二项应变指数");
  
  params.addParam<Real>("A3", 1.5e-25, "第三项系数");
  params.addParam<Real>("n3", 10.0, "第三项应力指数");
  params.addParam<Real>("m3", 3.9, "第三项应变指数");
  
  return params;
}

CoupledStressStrainCreepRate::CoupledStressStrainCreepRate(const InputParameters & parameters)
  : J2CreepPlasticity(parameters),
    _A1(getParam<Real>("A1")),
    _n1(getParam<Real>("n1")),
    _m1(getParam<Real>("m1")),
    _A2(getParam<Real>("A2")),
    _n2(getParam<Real>("n2")),
    _m2(getParam<Real>("m2")),
    _A3(getParam<Real>("A3")),
    _n3(getParam<Real>("n3")),
    _m3(getParam<Real>("m3"))
{
}

ADReal
CoupledStressStrainCreepRate::computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // Coupled stress-strain creep rate: 
  // dot_epsilon_c_eq = A1*sigma_eq^n1*epsilon_c_eq^m1 + A2*sigma_eq^n2*epsilon_c_eq^m2 + A3*sigma_eq^n3*epsilon_c_eq^m3
  
  // 避免数值问题，设置最小值
  // ADReal effective_stressMpa = effective_stress/1e6;
  ADReal sigma_eq = (effective_stress <= 1e-12) ? 1e-12 : effective_stress;
  ADReal epsilon_c_eq = (effective_creep_strain <= 1e-12) ? 1e-12 : effective_creep_strain;
  
  // 计算三个项
  ADReal term1 = _A1 * std::pow(sigma_eq, _n1) * std::pow(epsilon_c_eq, _m1);
  ADReal term2 = _A2 * std::pow(sigma_eq, _n2) * std::pow(epsilon_c_eq, _m2);
  ADReal term3 = _A3 * std::pow(sigma_eq, _n3) * std::pow(epsilon_c_eq, _m3);
  
  return term1 + term2 + term3;
}

ADReal
CoupledStressStrainCreepRate::computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // Derivative of creep rate with respect to effective stress:
  // d(dot_epsilon_c_eq)/d(sigma_eq) = A1*n1*sigma_eq^(n1-1)*epsilon_c_eq^m1 + 
  //                                   A2*n2*sigma_eq^(n2-1)*epsilon_c_eq^m2 + 
  //                                   A3*n3*sigma_eq^(n3-1)*epsilon_c_eq^m3
  // ADReal effective_stressMpa = effective_stress/1e6;
  if (effective_stress <= 1e-12)
    return 0.0;
  
  // 避免数值问题，设置最小值
  ADReal sigma_eq = (effective_stress <= 1e-12) ? 1e-12 : effective_stress;
  ADReal epsilon_c_eq = (effective_creep_strain <= 1e-12) ? 1e-12 : effective_creep_strain;
  
  // 计算三个项的导数
  ADReal dterm1_dsigma = _A1 * _n1 * std::pow(sigma_eq, _n1 - 1.0) * std::pow(epsilon_c_eq, _m1);
  ADReal dterm2_dsigma = _A2 * _n2 * std::pow(sigma_eq, _n2 - 1.0) * std::pow(epsilon_c_eq, _m2);
  ADReal dterm3_dsigma = _A3 * _n3 * std::pow(sigma_eq, _n3 - 1.0) * std::pow(epsilon_c_eq, _m3);
  
  return dterm1_dsigma + dterm2_dsigma + dterm3_dsigma;
}

ADReal
CoupledStressStrainCreepRate::computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // Derivative of creep rate with respect to effective creep strain:
  // d(dot_epsilon_c_eq)/d(epsilon_c_eq) = A1*m1*sigma_eq^n1*epsilon_c_eq^(m1-1) + 
  //                                        A2*m2*sigma_eq^n2*epsilon_c_eq^(m2-1) + 
  //                                        A3*m3*sigma_eq^n3*epsilon_c_eq^(m3-1)
  // ADReal effective_stress = effective_stress/1e6;
  if (effective_creep_strain <= 1e-12)
    return 0.0;
  
  // 避免数值问题，设置最小值
  ADReal sigma_eq = (effective_stress <= 1e-12) ? 1e-12 : effective_stress;
  ADReal epsilon_c_eq = (effective_creep_strain <= 1e-12) ? 1e-12 : effective_creep_strain;
  
  // 计算三个项的导数
  ADReal dterm1_depsilon = _A1 * _m1 * std::pow(sigma_eq, _n1) * std::pow(epsilon_c_eq, _m1 - 1.0);
  ADReal dterm2_depsilon = _A2 * _m2 * std::pow(sigma_eq, _n2) * std::pow(epsilon_c_eq, _m2 - 1.0);
  ADReal dterm3_depsilon = _A3 * _m3 * std::pow(sigma_eq, _n3) * std::pow(epsilon_c_eq, _m3 - 1.0);
  
  return dterm1_depsilon + dterm2_depsilon + dterm3_depsilon;
} 