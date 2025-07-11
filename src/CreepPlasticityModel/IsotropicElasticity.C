// IsotropicElasticity.C 修正版
#include "IsotropicElasticity.h"
#include "RaccoonUtils.h"

registerMooseObject("reproductionApp", IsotropicElasticity);

InputParameters
IsotropicElasticity::validParams()
{
  InputParameters params = ElasticityModel::validParams();
  params.addClassDescription("Isotropic elasticity with threshold energy option using E and nu as inputs.");
  
  params.addRequiredParam<MaterialPropertyName>("youngs_modulus", "Young's modulus E");
  params.addRequiredParam<MaterialPropertyName>("poissons_ratio", "Poisson's ratio nu");
  
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  
  MooseEnum decomp_options("NONE SPECTRAL VOLDEV MAXPRINCIPAL", "NONE");
  params.addParam<MooseEnum>("decomposition", decomp_options, "The decomposition method");
  
  params.addParam<bool>("use_threshold", false, "Whether to use threshold energy");
  params.addParam<MaterialPropertyName>("tensile_strength", "sigma0", "Critical stress for threshold");
  
  params.addParam<bool>("use_history_max", false, "Whether to use history maximum variable H for damage irreversibility");
  
  return params;
}

IsotropicElasticity::IsotropicElasticity(
    const InputParameters & parameters)
  : ElasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _youngs_modulus(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("youngs_modulus"))),
    _poissons_ratio(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("poissons_ratio"))),
    
    _d_name(getVar("phase_field", 0)->name()),
    
    // 修复参数获取 - 使用正确的参数名
    _psie_name(getParam<MaterialPropertyName>("strain_energy_density")),
    _psie(declareADProperty<Real>(_psie_name)),
    _psie_active(declareADProperty<Real>(_psie_name + "_active")),
    _dpsie_dd(declareADProperty<Real>(derivativePropertyName(_psie_name, {_d_name}))),
    
    // 修复参数获取 - 使用正确的参数名
    _g_name(getParam<MaterialPropertyName>("degradation_function")),
    _g(getADMaterialProperty<Real>(_g_name)),
    _dg_dd(getADMaterialProperty<Real>(derivativePropertyName(_g_name, {_d_name}))),
    
    _use_threshold(getParam<bool>("use_threshold")),
    _sigma_c(_use_threshold ? getADMaterialProperty<Real>(getParam<MaterialPropertyName>("tensile_strength")) 
                            : declareADProperty<Real>("dummy_sigma_c")),
    
    _use_history_max(getParam<bool>("use_history_max")),
    _history_max(_use_history_max ? declareProperty<Real>("history_max") : declareProperty<Real>("dummy_history_max")),
    _history_max_old(_use_history_max ? getMaterialPropertyOld<Real>("history_max") : getMaterialPropertyOld<Real>("dummy_history_max")),
    
    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>())
{
}

void
IsotropicElasticity::updateHistoryMax(ADReal psie_active_current)
{
  if (_use_history_max)
  {
    // 历史最大变量：H = max(H_old, psie_active_current)
    Real psie_active_value = MetaPhysicL::raw_value(psie_active_current);
    Real history_max_old_value = _history_max_old[_qp];
    _history_max[_qp] = std::max(history_max_old_value, psie_active_value);
  }
}

ADReal
IsotropicElasticity::applyThreshold(ADReal psie_active_raw)
{
  ADReal psie_active_final = psie_active_raw;
  
  // 首先应用历史最大变量
  if (_use_history_max)
  {
    // 使用历史最大值
    psie_active_final = _history_max[_qp];
  }
  
  // 然后应用阈值：max(psie_active, 0.5*sigma_c^2/E)
  if (_use_threshold)
  {
    ADReal threshold = 0.5 * _sigma_c[_qp] * _sigma_c[_qp] / _youngs_modulus[_qp];
    psie_active_final = std::max(psie_active_final, threshold);
  }
  
  return psie_active_final;
}

ADRankTwoTensor
IsotropicElasticity::computeStress(const ADRankTwoTensor & strain)
{
  switch (_decomposition)
  {
    case Decomposition::none:
      return computeStressNoDecomposition(strain);
    case Decomposition::spectral:
      return computeStressSpectralDecomposition(strain);
    case Decomposition::voldev:
      return computeStressVolDevDecomposition(strain);
    case Decomposition::maxprincipal:
      return computeStressMaxPrincipalDecomposition(strain);
    default:
      mooseError("Unknown decomposition type in IsotropicElasticity");
  }
  return ADRankTwoTensor(); // 不会执行到这里
}

ADReal
IsotropicElasticity::computeThreeShearModulus()
{
  // 计算3倍剪切模量：3G = 3 * E / (2 * (1 + nu))
  const ADReal nu = _poissons_ratio[_qp];
  const ADReal E = _youngs_modulus[_qp];
  const ADReal three_G = 3.0 * E / (2.0 * (1.0 + nu));
  
  return three_G;
}

ADRankTwoTensor
IsotropicElasticity::computeStressIntact(const ADRankTwoTensor & strain)
{
  // 计算拉梅常数
  const ADReal nu = _poissons_ratio[_qp];
  const ADReal E = _youngs_modulus[_qp];
  const ADReal K = E  / (3.0 *  (1.0 - 2.0 * nu));
  const ADReal G = E / (2.0 * (1.0 + nu));
  
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankTwoTensor stress_intact = K * strain.trace() * I2 + 2.0 * G * strain.deviatoric();
  return stress_intact;
}



ADRankTwoTensor
IsotropicElasticity::computeStressNoDecomposition(const ADRankTwoTensor & strain)
{
  // 计算拉梅常数
  const ADReal nu = _poissons_ratio[_qp];
  const ADReal E = _youngs_modulus[_qp];
  const ADReal K = E  / (3.0 *  (1.0 - 2.0 * nu));
  const ADReal G = E / (2.0 * (1.0 + nu));
  
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankTwoTensor stress_intact = K * strain.trace() * I2 + 2.0 * G * strain.deviatoric();
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  ADReal psie_active_raw = 0.5 * stress_intact.doubleContraction(strain);
  
  // 更新历史最大变量
  updateHistoryMax(psie_active_raw);
  
  // 应用阈值（考虑历史最大变量）
  _psie_active[_qp] = applyThreshold(psie_active_raw);
  
  _psie[_qp] = _g[_qp] * _psie_active[_qp];
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

ADRankTwoTensor
IsotropicElasticity::computeStressSpectralDecomposition(const ADRankTwoTensor & strain)
{
  // 计算拉梅常数
  const ADReal nu = _poissons_ratio[_qp];
  const ADReal E = _youngs_modulus[_qp];
  const ADReal K = E  / (3.0 *  (1.0 - 2.0 * nu));
  const ADReal G = E / (2.0 * (1.0 + nu));
  const ADReal lambda = K - 2 * G / LIBMESH_DIM;
  
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADReal strain_tr = strain.trace();
  ADReal strain_tr_pos = RaccoonUtils::Macaulay(strain_tr);

  // 谱分解
  ADRankTwoTensor strain_pos = RaccoonUtils::spectralDecomposition(strain);

  // 应力计算
  ADRankTwoTensor stress_intact = K * strain.trace() * I2 + 2.0 * G * strain.deviatoric();
  ADRankTwoTensor stress_pos = lambda * strain_tr_pos * I2 + 2.0 * G * strain_pos;
  ADRankTwoTensor stress_neg = stress_intact - stress_pos;
  ADRankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // 能量密度计算
  ADReal psie_intact = 0.5 * lambda * strain_tr * strain_tr + G * strain.doubleContraction(strain);
  ADReal psie_active_raw = 0.5 * lambda * strain_tr_pos * strain_tr_pos +
                          G * strain_pos.doubleContraction(strain_pos);
  
  // 更新历史最大变量
  updateHistoryMax(psie_active_raw);
  
  // 应用阈值（考虑历史最大变量）
  _psie_active[_qp] = applyThreshold(psie_active_raw);
  
  ADReal psie_inactive = psie_intact - psie_active_raw;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

ADRankTwoTensor
IsotropicElasticity::computeStressVolDevDecomposition(const ADRankTwoTensor & strain)
{
  // 计算拉梅常数
  const ADReal nu = _poissons_ratio[_qp];
  const ADReal E = _youngs_modulus[_qp];
  const ADReal G = E / (2.0 * (1.0 + nu));
  const ADReal K = E / (3.0 * (1.0 - 2.0 * nu));
  
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // 体积-偏差分解
  ADReal strain_tr = strain.trace();
  ADReal strain_tr_pos = RaccoonUtils::Macaulay(strain_tr);
  ADReal strain_tr_neg = strain_tr - strain_tr_pos;
  ADRankTwoTensor strain_dev = strain.deviatoric();

  // 应力计算
  ADRankTwoTensor stress_intact = K * strain.trace() * I2 + 2.0 * G * strain.deviatoric();
  ADRankTwoTensor stress_neg = K * strain_tr_neg * I2;
  ADRankTwoTensor stress_pos = stress_intact - stress_neg;
  ADRankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // 能量密度计算
  ADReal psie_intact = 0.5 * K * strain_tr * strain_tr + G * strain_dev.doubleContraction(strain_dev);
  ADReal psie_inactive = 0.5 * K * strain_tr_neg * strain_tr_neg;
  ADReal psie_active_raw = psie_intact - psie_inactive;
  
  // 更新历史最大变量
  updateHistoryMax(psie_active_raw);
  
  // 应用阈值（考虑历史最大变量）
  _psie_active[_qp] = applyThreshold(psie_active_raw);
  
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

ADRankTwoTensor
IsotropicElasticity::computeStressMaxPrincipalDecomposition(const ADRankTwoTensor & strain)
{
  // 基于 SmallDeformationH.C 的最大主应力分解算法
  const ADReal nu = _poissons_ratio[_qp];
  const ADReal E = _youngs_modulus[_qp];
  const ADReal K = E / (3.0 * (1.0 - 2.0 * nu));
  const ADReal G = E / (2.0 * (1.0 + nu));
  
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankTwoTensor stress_intact = K * strain.trace() * I2 + 2.0 * G * strain.deviatoric();
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  // 使用 symmetricEigenvalues 获取特征值，避免与 Point 不兼容问题
  std::vector<ADReal> eigenvals(LIBMESH_DIM);
  stress_intact.symmetricEigenvalues(eigenvals);

  // 找出最大特征值
  ADReal sigma_bar_eq = eigenvals[0];
  for (unsigned int i = 1; i < LIBMESH_DIM; i++)
    if (eigenvals[i] > sigma_bar_eq)
      sigma_bar_eq = eigenvals[i];
      
  // 应用 Macaulay 括号（只取正部分）
  sigma_bar_eq = RaccoonUtils::Macaulay(sigma_bar_eq);

  // 计算损伤驱动力：Y_bar = 0.5 * sigma_bar_eq^2 / E
  ADReal Y_bar = 0.5 * sigma_bar_eq * sigma_bar_eq / E;
  
  // 更新历史最大变量
  updateHistoryMax(Y_bar);
  
  // 应用阈值（考虑历史最大变量）
  _psie_active[_qp] = applyThreshold(Y_bar);
  
  _psie[_qp] = _g[_qp] * _psie_active[_qp];
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}