//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "J2Plasticity_C.h"
#include "CreepModel.h"
#include "J2Creep_P.h"

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
  // F_3 (Δε_eq^pl) = σ_eq^el - C(Δε_eq^pl N) - C(g(σ_Y (ε_eq_0^pl + Δε_eq^pl))Δt) - σ_Y (ε_eq_0^pl + Δε_eq^pl)
  // 这个塑性是专门给蠕变准备的，因此需要考虑蠕变，需要计算的残差如上。
  // 其他特别注意需要调用蠕变模型的蠕变率g(x)，即computeCreepRate((σ_Y (ε_eq_0^pl + Δε_eq^pl)))
  
  // 计算当前的有效塑性应变
  ADReal ep_current = _ep_old[_qp] + delta_ep;
  
  // 计算当前的屈服应力
  ADReal yield_stress = _hardening_model->plasticEnergy(ep_current, 1);
  
  // 计算第一项：σ_eq^el - C(Δε_eq^pl N)
  ADReal elastic_term = effective_trial_stress - 
                       _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]);
  
  // 计算蠕变项：C(g(σ_Y, ε_c)Δt)
  ADReal creep_term = 0.0;
  if (_creep_model)
  {
    // 获取蠕变率 g(σ_Y, ε_c)
    J2Creep_P * j2_creep_model = dynamic_cast<J2Creep_P *>(_creep_model);
    if (j2_creep_model)
    {
      // 获取当前蠕变应变 - 从蠕变模型获取
      ADReal current_creep_strain = j2_creep_model->getEffectiveCreepStrain();
      ADReal creep_rate = j2_creep_model->getCreepRate(yield_stress, current_creep_strain);
      creep_term = _elasticity_model->computeStress(_Np[_qp]*creep_rate * _dt).doubleContraction(_Np[_qp]);
    }
  }
  
  // 计算完整的残差
  return elastic_term - creep_term - yield_stress;
}

ADReal
J2Plasticity_C::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  // 对残差 F_3 (Δε_eq^pl) = σ_eq^el - C(Δε_eq^pl N) - C(g(σ_Y (ε_eq_0^pl + Δε_eq^pl), ε_c)Δt) - σ_Y (ε_eq_0^pl + Δε_eq^pl)
  // 完整的导数推导：
  // ∂F_3/∂(Δε_eq^pl) = -C - C(∂g/∂σ_Y * ∂σ_Y/∂(Δε_eq^pl) + ∂g/∂ε_c * ∂ε_c/∂(Δε_eq^pl))Δt - ∂σ_Y/∂(Δε_eq^pl)
  // 
  // 由于蠕变应变ε_c和塑性应变增量Δε_eq^pl是独立的状态变量，所以：
  // ∂ε_c/∂(Δε_eq^pl) = 0
  // 
  // 因此简化为：
  // ∂F_3/∂(Δε_eq^pl) = -C - C(∂g/∂σ_Y * ∂σ_Y/∂(Δε_eq^pl))Δt - ∂σ_Y/∂(Δε_eq^pl)
  
  // 计算当前的有效塑性应变
  ADReal ep_current = _ep_old[_qp] + delta_ep;
  
  // 计算弹性系数项
  ADReal elastic_coeff = _elasticity_model->computeStress(_Np[_qp]).doubleContraction(_Np[_qp]);
  
  // 计算硬化模型的导数 ∂σ_Y/∂(Δε_eq^pl)
  ADReal hardening_derivative = _hardening_model->plasticEnergy(ep_current, 2);
  
  // 计算蠕变项的导数
  ADReal creep_derivative = 0.0;
  if (_creep_model)
  {
    J2Creep_P * j2_creep_model = dynamic_cast<J2Creep_P *>(_creep_model);
    if (j2_creep_model)
    {
      // 计算当前屈服应力
      ADReal yield_stress = _hardening_model->plasticEnergy(ep_current, 1);
      
      // 获取当前蠕变应变
      ADReal current_creep_strain = j2_creep_model->getEffectiveCreepStrain();
      
      // 获取蠕变率对应力的导数 ∂g/∂σ_Y
      ADReal dg_dsigma = j2_creep_model->getCreepRateStressDerivative(yield_stress, current_creep_strain);
      
      // 蠕变项的导数：C(∂g/∂σ_Y * ∂σ_Y/∂(Δε_eq^pl))Δt
      // 注意：∂g/∂ε_c * ∂ε_c/∂(Δε_eq^pl) = 0，因为蠕变应变不直接依赖于塑性应变增量
      creep_derivative = elastic_coeff * dg_dsigma * hardening_derivative * _dt;
    }
  }
  
  // 返回完整的导数
  return -elastic_coeff - creep_derivative - hardening_derivative;
}
