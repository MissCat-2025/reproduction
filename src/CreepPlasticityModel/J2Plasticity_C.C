//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "J2Plasticity_C.h"

registerMooseObject("reproductionApp", J2Plasticity_C);

InputParameters
J2Plasticity_C::validParams()
{
  InputParameters params = PlasticityModel::validParams();
  params.addClassDescription("Small deformation $J_2$ plasticity. The plastic deformation is "
                             "updated using the additive decompsition of strain.");
  return params;
}

J2Plasticity_C::J2Plasticity_C(const InputParameters & parameters)
  : PlasticityModel(parameters)
{
}

void
J2Plasticity_C::updateState(ADRankTwoTensor & stress,
                                          ADRankTwoTensor & elastic_strain)
{
  //从Creep过来的stress、elastic_strain已经是试应力张量和试应变张量了
  // 1. 初始化塑性应变增量
  ADReal delta_ep = 0.0;
  
  // 2. 计算流动方向 (Prandtl-Reuss流动准则)
  ADRankTwoTensor stress_dev = stress.deviatoric();  // 偏应力
  ADReal stress_dev_norm = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  if (stress_dev_norm > 1e-12)
    _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;  // 归一化的流动方向
  else
    _Np[_qp].zero();
  
  returnMappingSolve(stress_dev_norm, delta_ep, _console);  // 牛顿迭代求解
    
  // 4. 更新状态变量
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + delta_ep * _Np[_qp];
  
  // 5. 更新最终应力
  elastic_strain -= delta_ep * _Np[_qp];
  stress = _elasticity_model->computeStress(elastic_strain);
}

Real
J2Plasticity_C::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                       const ADReal & delta_ep)
{
  return raw_value(
      effective_trial_stress -
      _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]));
}

ADReal
J2Plasticity_C::computeResidual(const ADReal & effective_trial_stress,
                                              const ADReal & delta_ep)
{
  // 简化的J2塑性残差计算，不再考虑蠕变耦合
  // 因为蠕变-塑性耦合已经在蠕变模型中完成
  
  // 计算当前的有效塑性应变
  ADReal ep_current = _ep_old[_qp] + delta_ep;
  
  // 计算当前的屈服应力
  ADReal yield_stress = _hardening_model->plasticEnergy(ep_current, 1);
  
  // 标准的J2塑性残差：σ_eq^el - C(Δε_eq^pl N) - σ_Y
  ADReal elastic_term = effective_trial_stress - 
                       _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]);
  
  // 返回简化的残差（不包含蠕变项）
  return elastic_term - yield_stress;
}

ADReal
J2Plasticity_C::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  // 简化的J2塑性导数计算，不再考虑蠕变耦合
  
  // 计算当前的有效塑性应变
  ADReal ep_current = _ep_old[_qp] + delta_ep;
  
  // 计算弹性系数项
  ADReal elastic_coeff = _elasticity_model->computeStress(_Np[_qp]).doubleContraction(_Np[_qp]);
  
  // 计算硬化模型的导数 ∂σ_Y/∂(Δε_eq^pl)
  ADReal hardening_derivative = _hardening_model->plasticEnergy(ep_current, 2);
  
  // 返回简化的导数（不包含蠕变项）
  return -elastic_coeff - hardening_derivative;
}
