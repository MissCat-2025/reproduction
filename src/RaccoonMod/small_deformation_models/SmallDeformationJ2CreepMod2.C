// //* This file is part of the RACCOON application
// //* being developed at Dolbow lab at Duke University
// //* http://dolbow.pratt.duke.edu

// #include "SmallDeformationJ2PlasticityMod.h"

// registerMooseObject("reproductionApp", SmallDeformationJ2PlasticityMod);

// InputParameters
// SmallDeformationJ2PlasticityMod::validParams()
// {
//   InputParameters params = SmallDeformationPlasticityModelMod::validParams();
//   params.addClassDescription("Small deformation $J_2$ plasticity. The plastic deformation is "
//                              "updated using the additive decompsition of strain.");
//   return params;
// }

// SmallDeformationJ2PlasticityMod::SmallDeformationJ2PlasticityMod(const InputParameters & parameters)
//   : SmallDeformationPlasticityModelMod(parameters)
// {
// }

// void
// SmallDeformationJ2PlasticityMod::updateState(ADRankTwoTensor & stress,
//                                           ADRankTwoTensor & elastic_strain)
// {
//   // 1. 假设没有塑性增量，计算试验状态
//   ADReal delta_e_cr = 0;
//   elastic_strain -= _plastic_strain_old[_qp] - _creep_strain_old[_qp];  // 去掉旧的塑性应变，旧的蠕变应变
//   stress = _elasticity_model->computeStress(elastic_strain);  // 计算试验应力
  
//   // 2. 计算流动方向 (Prandtl-Reuss流动准则)
//   ADRankTwoTensor stress_dev = stress.deviatoric();  // 偏应力
//   ADReal stress_dev_norm = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));//有效弹性试应力σ{tria}_eff
//   _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;  // 流动方向
  
//   // 1.迭代直接求解σ_eq^np=σ_eq^el-C(g(σ_eq^np )Δt)解出σ_eq^np，
//   //即使用returnMappingSolve配合即如下computeResidual1求出delta_e_cr、stress_dev_norm
//   /* ADReal
// SmallDeformationJ2PowerLawCreepMod::computeResidual1(const ADReal & effective_trial_stress,
//                                                     const ADReal & delta_ep)
// {
//   const ADReal stress_delta =
//       effective_trial_stress -
//       _elasticity_model->computeStress(delta_ep * _Np[_qp])
//           .doubleContraction(_Np[_qp]);
//   const ADReal yield_stress =
//       _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 1) +
//       _hardening_model->plasticDissipation(delta_ep, _ep_old[_qp] + delta_ep, 1);
//   const ADReal creep_rate = _coefficient * std::pow(stress_delta / yield_stress, _exponent);
//   return creep_rate * _dt - delta_ep;
// } */
  
//   returnMappingSolve(stress_dev_norm, delta_e_cr, _console);  // 牛顿迭代求解
//   // 5. 试应变-蠕变增量=弹性应变
//   elastic_strain -= delta_e_cr * _Np[_qp];
//   //弹性应变计算出对应的弹性应力
//   stress = _elasticity_model->computeStress(elastic_strain);
//   // 计算流动方向
//   ADRankTwoTensor stress_dev = stress.deviatoric();  // 偏应力
//   ADReal stress_dev_norm2 = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));//有效试应力σ{tria}_eff
//   _Np[_qp] = 1.5 * stress_dev / stress_dev_norm2;  // 流动方向

//   //此时，迭代计算出来的无塑性蠕变增量为delta_ep，有效试应力为stress_dev_norm2
//   //注意以上都是假设没有塑性应变增量


//   //如果无塑性应变时的应力小于屈服应力，则判定为蠕变+弹性
//   ADReal phi1 = computeResidual1(stress_dev_norm2, delta_e_cr);
//   //注意，此时
//   //如果无塑性应变时的应力大于屈服应力，则判定为蠕变+塑性+弹性

//   if (phi1 > 0)
//     // 如果无塑性应变时的应力大于屈服应力，则判定为蠕变+塑性+弹性
//     // returnMappingSolve2(stress_dev_norm, delta_ep, _console);  // 牛顿迭代求解
//     //需要考虑塑性应变增量
//     //需要求解的方程残差为
//     //F_3 (〖Δε〗_eq^pl )=σ_eq^el-C(〖Δε〗_eq^pl N)-C(g(σ_Y (ε_(〖eq〗_0)^pl+〖Δε〗_eq^pl ))Δt)-σ_Y (ε_(〖eq〗_0)^pl+〖Δε〗_eq^pl )
//     returnMappingSolve2(stress_dev_norm, delta_ep, _console);
//     //σ_eq=σ_eq^el-_elasticity_model(Δε_cr*_Np)-_elasticity_model(Δε_p*_Np)

//   // 4. 更新状态变量
//   _ep[_qp] = _ep_old[_qp] + delta_ep;
//   _plastic_strain[_qp] = _plastic_strain_old[_qp] + delta_ep * _Np[_qp];
  
//   // 5. 更新最终应力
//   elastic_strain -= delta_ep * _Np[_qp];
//   stress = _elasticity_model->computeStress(elastic_strain);
//   _hardening_model->plasticEnergy(_ep[_qp]);
// }

// Real
// SmallDeformationJ2PlasticityMod::computeReferenceResidual(const ADReal & effective_trial_stress,
//                                                        const ADReal & delta_ep)
// {
//   return raw_value(
//       effective_trial_stress -
//       _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]));
// }

// ADReal
// SmallDeformationJ2PlasticityMod::computeResidual(const ADReal & effective_trial_stress,
//                                               const ADReal & delta_ep)
// {
//   return effective_trial_stress -
//          _elasticity_model->computeStress(delta_ep * _Np[_qp]).doubleContraction(_Np[_qp]) -
//          _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 1);
// }

// ADReal
// SmallDeformationJ2PlasticityMod::computeDerivative(const ADReal & /*effective_trial_stress*/,
//                                                 const ADReal & delta_ep)
// {
//   return -_elasticity_model->computeStress(_Np[_qp]).doubleContraction(_Np[_qp]) -
//          _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 2);
// }
