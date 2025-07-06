//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "J2Creep_P.h"

registerMooseObject("reproductionApp", J2Creep_P);

InputParameters
J2Creep_P::validParams()
{
  InputParameters params = CreepModel::validParams();
  params.addClassDescription("Small deformation $J_2$ creep. The creep deformation is "
                             "updated using the additive decompsition of strain.");
  return params;
}

J2Creep_P::J2Creep_P(const InputParameters & parameters)
  : CreepModel(parameters)
{
}

void
J2Creep_P::updateState(ADRankTwoTensor & stress,
                                          ADRankTwoTensor & elastic_strain)
{
  // 1. 假设没有塑性增量，计算试验状态
  ADReal delta_ec = 0;
  elastic_strain -= _creep_strain_old[_qp];  // 去掉旧的蠕变应变
  elastic_strain -= _plastic_strain_old[_qp];// 去掉旧的塑性应变
  stress = _elasticity_model->computeStress(elastic_strain);  // 计算试应力
  
  // 2. 计算流动方向 (Prandtl-Reuss流动准则)
  ADRankTwoTensor stress_dev = stress.deviatoric();  // 偏应力
  ADReal stress_dev_norm = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;  // 归一化的流动方向
  
  // 3. 检查屈服条件
  returnMappingSolve(stress_dev_norm, delta_ec, _console);  // 牛顿迭代求解
  // 5. 试应变-蠕变增量=弹性应变
  elastic_strain2 =elastic_strain- delta_ec * _Np[_qp];
  //弹性应变计算出对应的弹性应力
  stress2 = _elasticity_model->computeStress(elastic_strain2);
  // 计算流动方向
  ADRankTwoTensor stress_dev2 = stress2.deviatoric();  // 偏应力
  ADReal stress_dev_norm2 = std::sqrt(1.5 * stress_dev2.doubleContraction(stress_dev2));//有效试应力σ{tria}_eff
  _Np[_qp] = 1.5 * stress_dev2 / stress_dev_norm2;  // 流动方向

  ADReal phi = computeResidual(stress_dev_norm2, delta_ec);  // 残差 = f(σ,εp)
 
  if (phi > 0)
  {
    _plastic_model->updateState(stress, elastic_strain);
  }  // 如果屈服
  else
  {
    // 4. 更新状态变量
    _ec[_qp] = _ec_old[_qp] + delta_ec;
    _creep_strain[_qp] = _creep_strain_old[_qp] + delta_ec * _Np[_qp];
    // 5. 更新最终应力
    elastic_strain3 = elastic_strain2- delta_ec * _Np[_qp];
    stress3 = _elasticity_model->computeStress(elastic_strain3);
  }
  
}

Real
SmallDeformationJ2PlasticityMod::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                       const ADReal & delta_ep)
{
  return raw_value(
      effective_trial_stress -
      _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]));
}

ADReal
SmallDeformationJ2PlasticityMod::computeResidual(const ADReal & effective_trial_stress,
                                              const ADReal & delta_ep)
{
  return effective_trial_stress -
         _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 1);
}

ADReal
SmallDeformationJ2PlasticityMod::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  return -_elasticity_model->computeStress(_Np[_qp]).doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 2);
}
