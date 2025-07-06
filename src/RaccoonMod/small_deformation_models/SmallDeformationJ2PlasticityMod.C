//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SmallDeformationJ2PlasticityMod.h"

registerMooseObject("reproductionApp", SmallDeformationJ2PlasticityMod);

InputParameters
SmallDeformationJ2PlasticityMod::validParams()
{
  InputParameters params = SmallDeformationPlasticityModelMod::validParams();
  params.addClassDescription("Small deformation $J_2$ plasticity. The plastic deformation is "
                             "updated using the additive decompsition of strain.");
  return params;
}

SmallDeformationJ2PlasticityMod::SmallDeformationJ2PlasticityMod(const InputParameters & parameters)
  : SmallDeformationPlasticityModelMod(parameters)
{
}

void
SmallDeformationJ2PlasticityMod::updateState(ADRankTwoTensor & stress,
                                          ADRankTwoTensor & elastic_strain)
{
  // 1. 假设没有塑性增量，计算试验状态
  ADReal delta_ep = 0;
  elastic_strain -= _plastic_strain_old[_qp];  // 去掉旧的塑性应变
  stress = _elasticity_model->computeStress(elastic_strain);  // 计算试验应力
  
  // 2. 计算流动方向 (Prandtl-Reuss流动准则)
  ADRankTwoTensor stress_dev = stress.deviatoric();  // 偏应力
  ADReal stress_dev_norm = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;  // 归一化的流动方向
  
  // 3. 检查屈服条件
  ADReal phi = computeResidual(stress_dev_norm, delta_ep);  // 残差 = f(σ,εp)
  if (phi > 0)  // 如果屈服
    returnMappingSolve(stress_dev_norm, delta_ep, _console);  // 牛顿迭代求解
    
  // 4. 更新状态变量
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + delta_ep * _Np[_qp];
  
  // 5. 更新最终应力
  elastic_strain -= delta_ep * _Np[_qp];
  stress = _elasticity_model->computeStress(elastic_strain);
  _hardening_model->plasticEnergy(_ep[_qp]);
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
