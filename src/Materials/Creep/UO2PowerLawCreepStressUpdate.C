//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "UO2PowerLawCreepStressUpdate.h"
#include "ElasticityTensorTools.h"
#include "RaccoonUtils.h"

registerMooseObject("reproductionApp", UO2PowerLawCreepStressUpdate);

InputParameters
UO2PowerLawCreepStressUpdate::validParams()
{
  InputParameters params = ADRadialReturnCreepStressUpdateBase::validParams();
  params.addClassDescription(
      "UO2蠕变模型，考虑裂变、温度和相场断裂退化。使用复杂的蠕变率表达式。");

  // 物理常数
  params.addParam<Real>("gas_constant", 8.314, "气体常数 (J/mol-K)");
  
  // 温度相关
  params.addRequiredCoupledVar("temperature", "温度 (K)");
  params.addRequiredCoupledVar("oxygen_ratio", "氧化学计量比");
  
  // 裂变与材料特性参数
  params.addRequiredParam<Real>("fission_rate", "裂变率密度 (fissions/m^3-s)");
  params.addParam<Real>("theoretical_density", 95.0, "理论密度百分比");
  params.addParam<Real>("grain_size", 10.0, "晶粒尺寸 (微米)");
  
  // 模型系数
  params.addParam<Real>("a1", 0.3919, "系数a1");
  params.addParam<Real>("a2", 1.31e-19, "系数a2");
  params.addParam<Real>("a3", 87.7, "系数a3");
  params.addParam<Real>("a6", 90.5, "系数a6");
  params.addParam<Real>("a7", 3.7226e-35, "系数a7");
  params.addParam<Real>("a8", 2.0391e-25, "系数a8");
  
  // 相场断裂相关参数
  params.addParam<MaterialPropertyName>("degradation_function", "g", "退化函数材料属性名");

  return params;
}

UO2PowerLawCreepStressUpdate::UO2PowerLawCreepStressUpdate(
    const InputParameters & parameters)
  : ADRadialReturnCreepStressUpdateBase(parameters),
    _gas_constant(getParam<Real>("gas_constant")),
    _temperature(adCoupledValue("temperature")),
    _oxygen_ratio(adCoupledValue("oxygen_ratio")),
    _fission_rate(getParam<Real>("fission_rate")),
    _theoretical_density(getParam<Real>("theoretical_density")),
    _grain_size(getParam<Real>("grain_size")),
    _Q1(0.0),  // 将在computeStressInitialize中计算
    _Q2(0.0),  // 将在computeStressInitialize中计算
    _Q3(21759.0),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _a3(getParam<Real>("a3")),
    _a6(getParam<Real>("a6")),
    _a7(getParam<Real>("a7")),
    _a8(getParam<Real>("a8")),
    _exp_Q1(0.0),  // 将在computeStressInitialize中计算
    _exp_Q2(0.0),  // 将在computeStressInitialize中计算
    _exp_Q3(0.0),  // 将在computeStressInitialize中计算
    _creep_strain(declareADProperty<RankTwoTensor>("creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>("creep_strain")),
    _g(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("degradation_function"))),
    _density_term1(0.0),  // 将在computeStressInitialize中计算
    _density_term2(0.0),  // 将在computeStressInitialize中计算
    _fission_term(0.0),    // 将在computeStressInitialize中计算
    _deviatoric_trial_stress(declareADProperty<RankTwoTensor>("deviatoric_trial_stress"))
{
}

void
UO2PowerLawCreepStressUpdate::computeStressInitialize(
    const ADReal & effective_trial_stress, const ADRankFourTensor & elasticity_tensor)
{
  // 预计算常用值
  const ADReal T = _temperature[_qp];
  const ADReal RT = _gas_constant * T;
  const ADReal inv_RT = 1.0 / RT;
  
  // 使用氧化学计量比
  const ADReal x = _oxygen_ratio[_qp];
  const ADReal log_x = std::log10(x);
  const ADReal exp_common = std::exp(-20.0 / log_x - 8.0);
  const ADReal denom = 1.0 / (exp_common + 1.0);
  
  // 计算激活能
  _Q1 = 74829.0 * denom + 301762.0;
  _Q2 = 83143.0 * denom + 469191.0;
  
  // 预计算指数项
  _exp_Q1 = std::exp(-_Q1 * inv_RT);
  _exp_Q2 = std::exp(-_Q2 * inv_RT);
  _exp_Q3 = std::exp(-_Q3 * inv_RT);
  
  // 预计算密度和晶粒尺寸相关项
  _density_term1 = 1.0 / ((_theoretical_density - _a3) * _grain_size * _grain_size);
  _density_term2 = 1.0 / (_theoretical_density - _a6);
  
  // 预计算裂变率项
  _fission_term = _a1 + _a2 * _fission_rate;

  // 1) 先调用基类，原本会做： _three_shear_modulus = 3*μ
  ADRadialReturnCreepStressUpdateBase::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);
    // 保存旧的 three_shear_modulus（即 3*μ）
  ADReal old_three_shear = _three_shear_modulus;
  // 计算并赋值新的 GL
  ADReal GL = compute_three_shear_modulus_New(elasticity_tensor);
  _three_shear_modulus = GL;
}


// —— 重载 updateState：先保存偏应力，再调用基类 —— 
// 2）修正 updateState 签名，去掉 template<…>，确保和头文件一致：
void
UO2PowerLawCreepStressUpdate::updateState(
    GenericRankTwoTensor<true> & strain_increment,
    GenericRankTwoTensor<true> & inelastic_strain_increment,
    const GenericRankTwoTensor<true> & rotation_increment,
    GenericRankTwoTensor<true> & stress_new,
    const RankTwoTensor & stress_old,
    const GenericRankFourTensor<true> & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old,
    bool compute_full_tangent_operator,
    RankFourTensor & tangent_operator)
{
  // 1) 保存当前的弹性偏应力
  _deviatoric_trial_stress[_qp] = stress_new.deviatoric();

  // 2) 调用父类实现，完成普通的 RadialReturn 逻辑
  ADRadialReturnCreepStressUpdateBase::updateState(strain_increment,
                                                   inelastic_strain_increment,
                                                   rotation_increment,
                                                   stress_new,
                                                   stress_old,
                                                   elasticity_tensor,
                                                   elastic_strain_old,
                                                   compute_full_tangent_operator,
                                                   tangent_operator);
}

// 1）在文件末尾添加这一段：GL 的具体计算
ADReal
UO2PowerLawCreepStressUpdate::compute_three_shear_modulus_New(const GenericRankFourTensor<true> & elasticity_tensor)
{
  // 取上次 updateState 保存的偏应力
  ADRankTwoTensor dev = _deviatoric_trial_stress[_qp];
  // 计算 ||dev||^2 并防止为零
  ADReal dev2 = dev.doubleContraction(dev);
  Real  val2 = MetaPhysicL::raw_value(dev2);
  if (MooseUtils::absoluteFuzzyEqual(val2, 0.0))
    dev2 = libMesh::TOLERANCE * libMesh::TOLERANCE;
  // 规范化因子 sqrt(1.5*dev2)
  ADReal norm = std::sqrt(1.5 * dev2);
  // 流向张量 Np
  ADRankTwoTensor Np = (1.5 * dev) / norm;

  // 谱分解并取正态部分
  ADRankTwoTensor Np_p = RaccoonUtils::spectralDecomposition(Np);
  ADReal alpha = Np_p.doubleContraction(Np_p);

  // 从弹性张量拿 μ
  ADReal mu = ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);
  // 退化函数 g
  ADReal gval = _g[_qp];

  // 返回 GL
  return 2.0 * mu * (1.5 + (gval - 1.0) * alpha);
}

ADReal
UO2PowerLawCreepStressUpdate::computeResidual(
    const ADReal & effective_trial_stress, const ADReal & scalar)
{
  // 计算径向返回迭代中的有效应力 (sigma_eff_rr = effective_trial_stress - G_L * scalar)
  ADReal stress_rr = effective_trial_stress - _three_shear_modulus * scalar;
  
  // 计算实际驱动蠕变本构的应力 (stress_creep_driving = g * sigma_eff_rr)
  ADReal stress_creep_driving = _g[_qp] * stress_rr;

  // 计算各分量蠕变率，注意这里使用的是 stress_creep_driving
  ADReal creep_th1 = _fission_term * _density_term1 * stress_creep_driving * _exp_Q1;
  ADReal creep_th2 = _a8 * _density_term2 * std::pow(stress_creep_driving, 4.5) * _exp_Q2;
  
  // 辐照蠕变，注意这里使用的是 stress_creep_driving
  ADReal creep_ir = _a7 * _fission_rate * stress_creep_driving * _exp_Q3;
  
  // 计算总蠕变率
  ADReal scalar_rate = (creep_th1 + creep_th2 + creep_ir);
  
  // 计算残差
  return scalar_rate * _dt - scalar;
}

ADReal
UO2PowerLawCreepStressUpdate::computeDerivative(
    const ADReal & effective_trial_stress, const ADReal & scalar)
{
  // 计算径向返回迭代中的有效应力 (sigma_eff_rr = effective_trial_stress - G_L * scalar)
  ADReal stress_rr = effective_trial_stress - _three_shear_modulus * scalar;

  // 计算实际驱动蠕变本构的应力 (stress_creep_driving = g * sigma_eff_rr)
  ADReal stress_creep_driving = _g[_qp] * stress_rr;
  
  // 计算蠕变本构中各项对 stress_creep_driving 的导数
  ADReal d_creep_th1_d_scd = _fission_term * _density_term1 * _exp_Q1;
  ADReal d_creep_th2_d_scd = 0.0; // 初始化
  if (MetaPhysicL::raw_value(stress_creep_driving) > 0 || MetaPhysicL::raw_value(stress_creep_driving) < 0 ) // 避免对0取pow(..., 3.5)
      d_creep_th2_d_scd = 4.5 * _a8 * _density_term2 * std::pow(stress_creep_driving, 3.5) * _exp_Q2;
  
  ADReal d_creep_ir_d_scd = _a7 * _fission_rate * _exp_Q3;
  
  // 计算 stress_creep_driving 对 scalar 的导数
  // stress_creep_driving = _g[_qp] * (effective_trial_stress - _three_shear_modulus * scalar)
  // d(stress_creep_driving)/d(scalar) = _g[_qp] * (-_three_shear_modulus)
  ADReal d_stress_creep_driving_d_scalar = -_g[_qp] * _three_shear_modulus;
  
  // 应用链式法则计算总蠕变率对 scalar 的导数
  ADReal d_scalar_rate_d_scalar = (d_creep_th1_d_scd + d_creep_th2_d_scd + d_creep_ir_d_scd) 
                                * d_stress_creep_driving_d_scalar;
  
  // 残差的导数
  return d_scalar_rate_d_scalar * _dt - 1.0;
}

Real
UO2PowerLawCreepStressUpdate::computeStrainEnergyRateDensity(
    const ADMaterialProperty<RankTwoTensor> & stress,
    const ADMaterialProperty<RankTwoTensor> & strain_rate)
{
  // 简化的应变能率密度计算
  return MetaPhysicL::raw_value(stress[_qp].doubleContraction(strain_rate[_qp]));
}