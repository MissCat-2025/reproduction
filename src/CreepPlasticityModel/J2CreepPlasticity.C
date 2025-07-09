//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "J2CreepPlasticity.h"
#include "PlasticityModel.h"
#include "PlasticHardeningModel.h"
#include "ElasticityModel.h"
registerMooseObject("reproductionApp", J2CreepPlasticity);

InputParameters
J2CreepPlasticity::validParams()
{
  InputParameters params = CreepModel::validParams();
  params.addClassDescription("J2 creep model with plasticity check. Implements the algorithm "
                             "for separated creep and plasticity with radial return mapping.");
  //塑性 (可选)
  params.addParam<MaterialName>("plasticity_model", "Name of the plasticity model (optional)");
  params.addParam<MaterialName>("hardening_model", "Name of the plastic hardening model (required when plasticity_model is provided)");
  
  // 牛顿迭代参数
  params.addParam<int>("max_iter", 1000, "Maximum number of iterations for Newton-Raphson solve");
  params.addParam<Real>("tolerance", 1e-10, "Tolerance for Newton-Raphson convergence");
    // 蠕变输出属性参数
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "creep_energy_density",
      "psic",
      "Name of the creep energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "gc", "The degradation function for creep");
  
  return params;
}

J2CreepPlasticity::J2CreepPlasticity(const InputParameters & parameters) : 
  CreepModel(parameters),
  DerivativeMaterialPropertyNameInterface(),
  _plasticity_model(nullptr),
  _hardening_model(nullptr),
  
  // 初始化蠕变输出属性
  _effective_creep_strain(declareADProperty<Real>("effective_creep_strain")),
  _d_name(getVar("phase_field", 0)->name()),
    // 修复参数获取 - 使用正确的参数名
  _psic_name(prependBaseName("creep_energy_density", true)),
  _psic(declareADProperty<Real>(_psic_name)),
  _psic_active(declareADProperty<Real>(_psic_name + "_active")),
  _psic_active_old(getMaterialPropertyOld<Real>(_psic_name + "_active")),
  _dpsic_dd(declareADProperty<Real>(derivativePropertyName(_psic_name, {_d_name}))),
  
  // 获取降解函数及其导数
  _gc(getADMaterialProperty<Real>(prependBaseName("degradation_function", true))),
  _dgc_dd(getADMaterialProperty<Real>(derivativePropertyName(prependBaseName("degradation_function", true), {_d_name})))
{
}

void
J2CreepPlasticity::initialSetup()
{
  // 塑性模型是可选的
  if (isParamValid("plasticity_model"))
  {
    _plasticity_model = dynamic_cast<PlasticityModel *>(&getMaterial("plasticity_model"));
    if (!_plasticity_model)
      paramError("plasticity_model",
                 "Plasticity model " + getParam<MaterialName>("plasticity_model") +
                     " is not compatible with " + name());
    
    // 如果提供了塑性模型，则需要硬化模型
    if (!isParamValid("hardening_model"))
      paramError("hardening_model",
                 "Hardening model must be provided when plasticity_model is specified");
    
    _hardening_model = dynamic_cast<PlasticHardeningModel *>(&getMaterial("hardening_model"));
    if (!_hardening_model)
      paramError("hardening_model",
                 "Hardening model " + getParam<MaterialName>("hardening_model") +
                     " is not compatible with " + name());
  }
}

void
J2CreepPlasticity::setQp(unsigned int qp)
{
  // 调用基类的setQp方法
  CreepModel::setQp(qp);
  
  // 设置塑性模型的qp
  if (_plasticity_model)
    _plasticity_model->setQp(qp);
  
  // 设置硬化模型的qp
  if (_hardening_model)
    _hardening_model->setQp(qp);
}

void
J2CreepPlasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain)
{
  // 设置塑性模型的蠕变模型指针（用于蠕变-塑性耦合计算）
  if (_plasticity_model)
  {
    _plasticity_model->setCreepModel(this);
  }
    
  // Step 1: 计算弹性试验状态（移除旧的蠕变应变和塑性应变）
  ADRankTwoTensor elastic_trial_strain = elastic_strain;
  elastic_trial_strain -= _creep_strain_old[_qp];  // 减去旧的蠕变应变张量
  
  // 如果有塑性模型，也要减去旧的塑性应变张量
  if (_plasticity_model)
    elastic_trial_strain -= _plasticity_model->getPlasticStrainOld();
  
  // Step 2: 计算非塑性试验状态（假设没有塑性增量，只考虑蠕变）
  ADReal delta_ec = 0.0;
  
  // 计算试验应力
  ADRankTwoTensor stress_trial = _elasticity_model->computeStress(elastic_trial_strain);
  
  // 计算偏应力和有效应力
  ADRankTwoTensor stress_dev = stress_trial.deviatoric();
  ADReal effective_stress_trial = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 计算流动方向
  if (effective_stress_trial > 1e-12)
    _Nc[_qp] = 1.5 * stress_dev / effective_stress_trial;
  else
    _Nc[_qp].zero();
  
  // 求解非塑性蠕变增量（使用径向返回映射）
  if (effective_stress_trial > 1e-12)
  {
    returnMappingSolve(effective_stress_trial, delta_ec, _console);
  }
  
  // 更新蠕变应变和有效蠕变应变
  ADRankTwoTensor creep_strain_np = _creep_strain_old[_qp] + delta_ec * _Nc[_qp];
  ADReal ec_np = _ec_old[_qp] + delta_ec;
  
  // 计算非塑性状态下的弹性应变和应力
  ADRankTwoTensor elastic_strain_np = elastic_trial_strain - delta_ec * _Nc[_qp];
  ADRankTwoTensor stress_np = _elasticity_model->computeStress(elastic_strain_np);
  
  // 重新计算偏应力和有效应力
  stress_dev = stress_np.deviatoric();
  ADReal effective_stress_np = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 更新流动方向
  if (effective_stress_np > 1e-12)
    _Nc[_qp] = 1.5 * stress_dev / effective_stress_np;
  else
    _Nc[_qp].zero();
  
  // Step 3: 检查屈服条件
  if (!_plasticity_model || !_hardening_model || f_no_plastic_strain(effective_stress_np))
  {
    // 情况1：没有塑性模型或硬化模型，或不屈服，只有蠕变
    _creep_strain[_qp] = creep_strain_np;
    _ec[_qp] = ec_np;
    elastic_strain = elastic_strain_np;
    stress = stress_np;
    
    // 如果有硬化模型但不屈服，仍需要更新塑性应变能（使用旧的塑性应变）
    if (_plasticity_model && _hardening_model)
    {
      Real ep_current = _plasticity_model->getEffectivePlasticStrainOld();
      _hardening_model->plasticEnergy(ep_current, 0);
    }
  }  
  else
  {
    // 情况2：有塑性模型且屈服，需要塑性-蠕变耦合求解
    // 根据论文算法，求解F3方程：
    // F3(Δε_eq^pl) = σ_eq^el - C(Δε_eq^pl N) - C(g(σ_Y(ε_eq0^pl + Δε_eq^pl))Δt) - σ_Y(ε_eq0^pl + Δε_eq^pl)
    
    solveCreepPlasticityCoupled(stress, elastic_strain, elastic_trial_strain, 
                                effective_stress_trial, delta_ec, creep_strain_np, ec_np);
    // 注意：塑性应变能在solveCreepPlasticityCoupled中已经更新
  }
    
  // 更新输出属性
  // 输出有效蠕变应变（标量）
  _effective_creep_strain[_qp] = _ec[_qp];
  
  // 计算蠕变能量密度：基于有效应力与有效蠕变应变的积分
  // 对于J2蠕变，蠕变能密度 W_c = ∫ σ_eq * dε_eq^c
  // 在增量形式下，使用梯形积分规则：W_c = W_c_old + 0.5 * (σ_eq_old + σ_eq) * Δε_eq^c
  
  // 计算当前的有效应力
  stress_dev = stress.deviatoric();
  ADReal effective_stress = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 计算蠕变应变增量
  ADReal delta_ec_total = _ec[_qp] - _ec_old[_qp];
  
  // 获取旧的蠕变能量密度
  ADReal creep_energy_old = _psic_active_old[_qp];
  
  // 蠕变能量密度的增量计算
  // 对于小增量，使用简化的矩形积分：W_c += σ_eq * Δε_eq^c
  // 对于更精确的计算，可以使用梯形积分或其他数值积分方法
  ADReal creep_energy_increment = effective_stress * delta_ec_total;
  
  // 确保蠕变能量密度非负（物理要求）
  ADReal creep_energy = creep_energy_old + std::max(0.0, raw_value(creep_energy_increment));
  
  // 设置蠕变能量密度
  _psic_active[_qp] = creep_energy;
  _psic[_qp] = _gc[_qp] * _psic_active[_qp];
  _dpsic_dd[_qp] = _dgc_dd[_qp] * _psic_active[_qp];
}

bool
J2CreepPlasticity::f_no_plastic_strain(const ADReal & effective_stress)
{
  // 如果没有塑性模型，只有蠕变
  if (!_plasticity_model || !_hardening_model)
    return true;
  
  // 通过自己的硬化模型获取屈服应力
  // 使用上一步的塑性应变进行硬化计算
  Real ep_old_value = _plasticity_model->getEffectivePlasticStrainOld();
  
  // 通过硬化模型获取当前屈服应力
  ADReal yield_stress = _hardening_model->plasticEnergy(ep_old_value, 1);
  
  // 屈服函数: f = σ_eq - σ_Y
  ADReal f = effective_stress - yield_stress;
  
  // 返回 f <= 0 表示不屈服（只有蠕变）
  // 返回 f > 0 表示屈服（需要考虑塑性）
  return f <= 0.0;
}

ADReal
J2CreepPlasticity::computeResidual(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 根据文档公式:
  // F1(ε_eq^(cr,np)) = ε_eq^(cr,np) - ε_eq0^cr - g(σ_eq^np, ε_eq^(cr,np)) * Δt = 0
  // 其中 σ_eq^np = σ_eq^el - C * (ε_eq^(cr,np) - ε_eq0^cr)
  
  ADReal ec_current = _ec_old[_qp] + delta_ec;
  
  // 计算非塑性有效应力
  ADReal C2 = _elasticity_model->computeStress(_Nc[_qp]).doubleContraction(_Nc[_qp]);
  ADReal effective_stress_np = effective_trial_stress - C2 * delta_ec;
  
  // // DEBUG_OUTPUT: 输出计算过程
  // if (_qp == 0)
  // {
  //   _console << "DEBUG_RESIDUAL: effective_trial_stress = " << raw_value(effective_trial_stress) << std::endl;
  //   _console << "DEBUG_RESIDUAL: delta_ec = " << raw_value(delta_ec) << std::endl;
  //   _console << "DEBUG_RESIDUAL: ec_current = " << raw_value(ec_current) << std::endl;
  //   _console << "DEBUG_RESIDUAL: C2 = " << raw_value(C2) << std::endl;
  //   _console << "DEBUG_RESIDUAL: effective_stress_np = " << raw_value(effective_stress_np) << std::endl;
  //   _console << "DEBUG_RESIDUAL: _dt = " << _dt << std::endl;
  // }
  
  // 计算蠕变率
  ADReal creep_rate = computeCreepRate(effective_stress_np, ec_current);
  
  // // DEBUG_OUTPUT: 输出蠕变率
  // if (_qp == 0)
  // {
  //   _console << "DEBUG_RESIDUAL: creep_rate = " << raw_value(creep_rate) << std::endl;
  //   _console << "DEBUG_RESIDUAL: creep_rate * _dt = " << raw_value(creep_rate * _dt) << std::endl;
  // }
  
  // 返回残差
  ADReal residual = delta_ec - creep_rate * _dt;
  
  // // DEBUG_OUTPUT: 输出残差
  // if (_qp == 0)
  // {
  //   _console << "DEBUG_RESIDUAL: residual = " << raw_value(residual) << std::endl;
  // }
  
  return residual;
}

ADReal
J2CreepPlasticity::computeDerivative(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 残差对 delta_ec 的导数
  ADReal ec_current = _ec_old[_qp] + delta_ec;
  
  // 计算非塑性有效应力
  ADReal C = _elasticity_model->computeStress(_Nc[_qp]).doubleContraction(_Nc[_qp]);
  ADReal effective_stress_np = effective_trial_stress - C * delta_ec;
  
  // 计算蠕变率的导数
  ADReal dg_dsigma = computeCreepRateStressDerivative(effective_stress_np, ec_current);
  ADReal dg_dec = computeCreepRateStrainDerivative(effective_stress_np, ec_current);
  
  // 返回导数: ∂F/∂(δεc) = 1 - Δt * (∂g/∂σ * (-C) + ∂g/∂εc)
  return 1.0 - _dt * (dg_dsigma * (-C) + dg_dec);
}

Real
J2CreepPlasticity::computeReferenceResidual(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 计算参考残差（用于收敛判断）
  ADReal residual = computeResidual(effective_trial_stress, delta_ec);
  return raw_value(residual);
}

// 蠕变率计算的默认实现（可以被子类覆盖）
ADReal
J2CreepPlasticity::computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 默认实现：简单的幂律蠕变 g(σ) = A * σ^n
  // 子类应该覆盖这个函数以实现具体的蠕变模型
  if (effective_stress <= 1e-12)
    return 0.0;
  
  // 临时使用简单的线性关系
  ADReal A = 1.0e-6; // 蠕变系数，应该从输入参数获取
  return A * effective_stress;
}

ADReal
J2CreepPlasticity::computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 蠕变率对应力的导数
  // 默认实现：对应简单线性关系的导数
  if (effective_stress <= 1e-12)
    return 0.0;
  
  ADReal A = 1.0e-6; // 蠕变系数
  return A;
}

ADReal
J2CreepPlasticity::computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 蠕变率对蠕变应变的导数
  // 默认实现：对于不依赖蠕变应变的模型，导数为0
  return 0.0;
}

void
J2CreepPlasticity::solveCreepPlasticityCoupled(ADRankTwoTensor & stress, 
                                       ADRankTwoTensor & elastic_strain,
                                       const ADRankTwoTensor & elastic_trial_strain,
                                       const ADReal & effective_stress_trial,
                                       const ADReal & delta_ec_np,
                                       const ADRankTwoTensor & creep_strain_np,
                                       const ADReal & ec_np)
{
  // 根据论文算法求解F3方程：
  // F3(Δε_eq^pl) = σ_eq^el - C(Δε_eq^pl N) - C(g(σ_Y(ε_eq0^pl + Δε_eq^pl))Δt) - σ_Y(ε_eq0^pl + Δε_eq^pl) = 0
  
  if (!_plasticity_model || !_hardening_model)
  {
    // 如果没有塑性模型或硬化模型，回退到只有蠕变
    _creep_strain[_qp] = creep_strain_np;
    _ec[_qp] = ec_np;
    elastic_strain = elastic_trial_strain - delta_ec_np * _Nc[_qp];
    stress = _elasticity_model->computeStress(elastic_strain);
    return;
  }
  
  // 获取弹性系数 C = 2μ * 3/2 = 3μ (对于J2塑性)
  ADReal C = _elasticity_model->computeStress(_Nc[_qp]).doubleContraction(_Nc[_qp]);
  
  
  // 获取上一步的塑性应变
  Real ep_old = _plasticity_model->getEffectivePlasticStrainOld();
  
  // 初始化塑性应变增量
  ADReal delta_ep = 0.0;
  
  // 牛顿迭代求解F3方程
  int max_iter = getParam<int>("max_iter");
  Real tolerance = getParam<Real>("tolerance");
  
  for (int iter = 0; iter < max_iter; iter++)
  {
    // 计算当前有效塑性应变
    ADReal ep_current = ep_old + delta_ep;
    
    // 计算当前屈服应力
    ADReal yield_stress = _hardening_model->plasticEnergy(ep_current, 1);
    
    // 计算蠕变率 g(σ_Y, ε_c)
    // 注意：这里使用屈服应力作为蠕变的驱动应力
    ADReal creep_rate = computeCreepRate(yield_stress, ec_np);
    
    // 计算F3残差
    ADReal F3 = effective_stress_trial - _elasticity_model->computeStress(delta_ep*_Nc[_qp]).doubleContraction(_Nc[_qp]) 
    - _elasticity_model->computeStress(creep_rate * _dt*_Nc[_qp]).doubleContraction(_Nc[_qp]) - yield_stress;
    
    // 检查收敛
    if (std::abs(raw_value(F3)) < tolerance)
      break;
    
    // 计算导数 dF3/d(Δε_eq^pl)
    ADReal hardening_derivative = _hardening_model->plasticEnergy(ep_current, 2);  // dσ_Y/dε_eq^pl
    ADReal creep_stress_derivative = computeCreepRateStressDerivative(yield_stress, ec_np);  // dg/dσ_Y
    
    // dF3/d(Δε_eq^pl) = -C - C * (dg/dσ_Y * dσ_Y/dε_eq^pl) * Δt - dσ_Y/dε_eq^pl
    ADReal dF3_ddelta_ep = -C - C * creep_stress_derivative * hardening_derivative * _dt - hardening_derivative;
    
    // 牛顿迭代更新
    delta_ep -= F3 / dF3_ddelta_ep;
    
    // 防止负值
    if (delta_ep < 0.0)
      delta_ep = 0.0;
  }
  
  // 使用求解得到的塑性应变增量更新状态
  
  // 1. 更新最终的塑性应变
  ADReal ep_final = ep_old + delta_ep;
  
  // 2. 计算最终的屈服应力和蠕变率
  ADReal yield_stress_final = _hardening_model->plasticEnergy(ep_final, 1);
  ADReal creep_rate_final = computeCreepRate(yield_stress_final, ec_np);
  ADReal delta_ec_final = creep_rate_final * _dt;
  
  // 3. 更新蠕变应变
  _creep_strain[_qp] = _creep_strain_old[_qp] + delta_ec_final * _Nc[_qp];
  _ec[_qp] = _ec_old[_qp] + delta_ec_final;
  
  // 4. 更新弹性应变（减去塑性应变增量和蠕变应变增量）
  elastic_strain = elastic_trial_strain - delta_ep * _Nc[_qp] - delta_ec_final * _Nc[_qp];
  
  // 5. 计算最终应力
  stress = _elasticity_model->computeStress(elastic_strain);
  
  // 6. 手动更新塑性模型的状态变量（避免调用updateState）
  updatePlasticityModelState(delta_ep);
  
  // 7. 更新塑性应变能
  // 使用最终的塑性应变值来计算塑性应变能
  _hardening_model->plasticEnergy(ep_final, 0);
}

void
J2CreepPlasticity::updatePlasticityModelState(const ADReal & delta_ep)
{
  // 手动更新塑性模型的状态变量，避免循环调用
  if (_plasticity_model)
  {
    // 设置塑性模型的当前quadrature point
    _plasticity_model->setQp(_qp);
    
    // 使用新的接口直接设置塑性应变状态
    _plasticity_model->setPlasticStrainState(delta_ep, _Nc[_qp]);
  }
}
