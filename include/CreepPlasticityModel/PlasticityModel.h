//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "ADSingleVariableReturnMappingSolution.h"
#include "BaseNameInterface.h"
#include "PlasticHardeningModel.h"

class ElasticityModel;
class CreepModel;

class PlasticityModel : public Material,
                                        public ADSingleVariableReturnMappingSolution,
                                        public BaseNameInterface
{
public:
  static InputParameters validParams();

  PlasticityModel(const InputParameters & parameters);

  virtual void initialSetup() override;

  /// Set the current quadrature point
  virtual void setQp(unsigned int qp);

  /// Set the associated elasticity model
  virtual void setElasticityModel(ElasticityModel * elasticity_model);

  /// Set the associated creep model
  virtual void setCreepModel(CreepModel * creep_model);

  /**
   * Update the stress and elastic strain if need to following the specified plastic flow
   * @param stress         The stress
   * @param elastic_strain The elastic strain
   */
  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) = 0;

  /**
   * Get the effective plastic strain from the previous time step
   * @return The effective plastic strain at the current quadrature point
   */
  virtual Real getEffectivePlasticStrainOld() const;

  /**
   * Get the plastic strain tensor from the previous time step
   * @return The plastic strain tensor at the current quadrature point
   */
  virtual const RankTwoTensor & getPlasticStrainOld() const;

  /**
   * Manually set the plastic strain state (used by creep-plasticity coupling to avoid circular calls)
   * @param delta_ep The increment of effective plastic strain
   * @param flow_direction The plastic flow direction
   */
  virtual void setPlasticStrainState(const ADReal & delta_ep, const ADRankTwoTensor & flow_direction);

  // @{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are
  // unused in all inheriting classes and should not be overwritten.
  void resetQpProperties() final {}
  void resetProperties() final {}
  // @}

protected:
  virtual void initQpStatefulProperties() override;

  /// The elasticity model
  ElasticityModel * _elasticity_model;

  /// The creep model
  CreepModel * _creep_model;

  /// The plastic strain
  ADMaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// The (scalar) effective plastic strain
  ADMaterialProperty<Real> & _ep;
  const MaterialProperty<Real> & _ep_old;

  /// The flow direction
  ADMaterialProperty<RankTwoTensor> & _Np;

  /// The plastic hardening model
  PlasticHardeningModel * _hardening_model;
};
