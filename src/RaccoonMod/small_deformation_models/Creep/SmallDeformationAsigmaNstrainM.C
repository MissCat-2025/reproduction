//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SmallDeformationAsigmaNstrainM.h"

registerMooseObject("reproductionApp", SmallDeformationAsigmaNstrainM);

InputParameters
SmallDeformationAsigmaNstrainM::validParams()
{
  InputParameters params = SmallDeformationJ2PlasticityMod::validParams();
  params.addClassDescription("Small deformation creep model using the formula: "
                             "$$\\dot{\\varepsilon}_{c,eq} = A_1(\\sigma_{eq})^{n_1}(\\varepsilon_{c,eq})^{m_1} + "
                             "A_2(\\sigma_{eq})^{n_2}(\\varepsilon_{c,eq})^{m_2} + "
                             "A_3(\\sigma_{eq})^{n_3}(\\varepsilon_{c,eq})^{m_3}$$");

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

SmallDeformationAsigmaNstrainM::SmallDeformationAsigmaNstrainM(const InputParameters & parameters)
  : SmallDeformationJ2PlasticityMod(parameters),
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
SmallDeformationAsigmaNstrainM::computeResidual(const ADReal & effective_trial_stress,
                                                const ADReal & delta_ep)
{
  const ADReal dt = _dt;
  
  // 计算当前累积蠕变应变
  const ADReal epsilon_c_eq = _ep_old[_qp] + delta_ep;
  
  // 关键修正：对于负指数，需要更大的安全值
  const ADReal epsilon_c_eq_safe = std::max(epsilon_c_eq, 1.0e-3);  // 改为1e-6
  
  // 单位转换：Pa -> MPa
  const ADReal sigma_eq_MPa = effective_trial_stress / 1.0e6;
  
  // 计算三个蠕变率项（原始单位：1/h）
  ADReal creep_rate_h = 0.0;
  
  // 第一项：A1 * sigma^n1 * epsilon^m1（注意m1是负数）
  ADReal term1 = 0.0;
  if (sigma_eq_MPa > 0.0)
  {
    // 对于负指数，使用更安全的计算
    term1 = _A1 * std::pow(sigma_eq_MPa, _n1) * std::pow(epsilon_c_eq_safe, _m1);
    creep_rate_h += term1;
  }
  
  // 第二项：A2 * sigma^n2 * epsilon^m2
  ADReal term2 = 0.0;
  if (sigma_eq_MPa > 0.0 && epsilon_c_eq > 0.0)
  {
    term2 = _A2 * std::pow(sigma_eq_MPa, _n2) * std::pow(epsilon_c_eq, _m2);
    creep_rate_h += term2;
  }
  
  // 第三项：A3 * sigma^n3 * epsilon^m3
  ADReal term3 = 0.0;
  if (sigma_eq_MPa > 0.0 && epsilon_c_eq > 0.0)
  {
    term3 = _A3 * std::pow(sigma_eq_MPa, _n3) * std::pow(epsilon_c_eq, _m3);
    creep_rate_h += term3;
  }
  
  // 单位转换：1/h -> 1/s
  const ADReal creep_rate_s = creep_rate_h / 3600.0;
  
  // 计算残差
  const ADReal expected_creep_increment = creep_rate_s * dt;
  const ADReal residual = delta_ep - expected_creep_increment;
  
  return residual;
}

ADReal
SmallDeformationAsigmaNstrainM::computeDerivative(const ADReal & effective_trial_stress,
                                                  const ADReal & delta_ep)
{
  const ADReal dt = _dt;
  
  if (dt <= 0.0)
    return 1.0;

  const ADReal epsilon_c_eq = _ep_old[_qp] + delta_ep;
  
  // 与computeResidual保持一致的安全值
  const ADReal epsilon_c_eq_safe = std::max(epsilon_c_eq, 1.0e-6);
  const ADReal sigma_eq_MPa = effective_trial_stress / 1.0e6;
  
  // 计算蠕变率对蠕变应变的导数
  ADReal d_creep_rate_d_eps = 0.0;
  
  // 第一项的导数：A1 * sigma^n1 * epsilon^m1
  if (_m1 != 0.0 && sigma_eq_MPa > 0.0)
  {
    // 对于负指数，使用安全值
    d_creep_rate_d_eps += _A1 * std::pow(sigma_eq_MPa, _n1) * _m1 * std::pow(epsilon_c_eq_safe, _m1 - 1.0);
  }
  
  // 第二项的导数：A2 * sigma^n2 * epsilon^m2
  if (_m2 != 0.0 && sigma_eq_MPa > 0.0 && epsilon_c_eq > 0.0)
  {
    d_creep_rate_d_eps += _A2 * std::pow(sigma_eq_MPa, _n2) * _m2 * std::pow(epsilon_c_eq, _m2 - 1.0);
  }
  
  // 第三项的导数：A3 * sigma^n3 * epsilon^m3
  if (_m3 != 0.0 && sigma_eq_MPa > 0.0 && epsilon_c_eq > 0.0)
  {
    d_creep_rate_d_eps += _A3 * std::pow(sigma_eq_MPa, _n3) * _m3 * std::pow(epsilon_c_eq, _m3 - 1.0);
  }
  
  // 单位转换：1/h -> 1/s
  const ADReal d_creep_rate_s_d_eps = d_creep_rate_d_eps / 3600.0;
  
  // 残差的导数：∂R/∂(δεp) = 1 - ∂(蠕变率)/∂εp * dt
  return 1.0 - d_creep_rate_s_d_eps * dt;
}