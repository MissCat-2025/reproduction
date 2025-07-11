// IsotropicElasticity.h
#pragma once

#include "ElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * 各向同性弹性材料，可选阈值能量，直接使用杨氏模量和泊松比作为输入参数
 * 支持历史最大变量H以确保损伤不可逆性
 */
class IsotropicElasticity : public ElasticityModel,
                                                    public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();
  IsotropicElasticity(const InputParameters & parameters);
  
  virtual ADRankTwoTensor computeStress(const ADRankTwoTensor & strain) override;
  
  /// 计算3倍剪切模量：3G = 3 * E / (2 * (1 + nu))
  virtual ADReal computeThreeShearModulus();

  ADRankTwoTensor computeStressIntact(const ADRankTwoTensor & strain);

protected:
  // 弹性参数
  const ADMaterialProperty<Real> & _youngs_modulus;
  const ADMaterialProperty<Real> & _poissons_ratio;
  
  // 相场变量名
  const VariableName _d_name;
  
  // 能量密度及其导数
  const MaterialPropertyName _psie_name;
  ADMaterialProperty<Real> & _psie;
  ADMaterialProperty<Real> & _psie_active;
  ADMaterialProperty<Real> & _dpsie_dd;
  
  // 退化函数及其导数
  const MaterialPropertyName _g_name;
  const ADMaterialProperty<Real> & _g;
  const ADMaterialProperty<Real> & _dg_dd;
  
  // 阈值参数
  const bool _use_threshold;
  const ADMaterialProperty<Real> & _sigma_c;
  
  // 历史最大变量H参数
  const bool _use_history_max;
  MaterialProperty<Real> & _history_max;
  const MaterialProperty<Real> & _history_max_old;
  
  // 分解方法
  enum class Decomposition { none, spectral, voldev, maxprincipal };
  const Decomposition _decomposition;
  
  // 计算方法
  ADRankTwoTensor computeStressNoDecomposition(const ADRankTwoTensor & strain);
  ADRankTwoTensor computeStressSpectralDecomposition(const ADRankTwoTensor & strain);
  ADRankTwoTensor computeStressVolDevDecomposition(const ADRankTwoTensor & strain);
  ADRankTwoTensor computeStressMaxPrincipalDecomposition(const ADRankTwoTensor & strain);
  
  // 应用阈值
  ADReal applyThreshold(ADReal psie_active_raw);
  
  // 更新历史最大变量
  void updateHistoryMax(ADReal psie_active_current);
};