//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "ADSingleVariableReturnMappingSolution.h"
#include "BaseNameInterface.h"

class ElasticityModel;
class PlasticityModel;
class CreepModel : public Material,
                                        public ADSingleVariableReturnMappingSolution,
                                        public BaseNameInterface
{
public:
  static InputParameters validParams();

  CreepModel(const InputParameters & parameters);

  // virtual void initialSetup() override;

  /// Set the current quadrature point
  virtual void setQp(unsigned int qp);

  /// Set the associated elasticity model
  virtual void setElasticityModel(ElasticityModel * elasticity_model);

  /// Set the associated plasticity model
  virtual void setPlasticityModel(PlasticityModel * plasticity_model);

  /**
   * Update the stress and elastic strain if need to following the specified plastic flow
   * @param stress         The stress
   * @param elastic_strain The elastic strain
   */
  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) = 0;

  // @{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are
  // unused in all inheriting classes and should not be overwritten.
  void resetQpProperties() final {}
  void resetProperties() final {}
  // @}

protected:
  virtual void initQpStatefulProperties() override;

  /// The elasticity model
  ElasticityModel * _elasticity_model;
  /// The plasticity model
  PlasticityModel * _plasticity_model;
  /// The creep strain
  ADMaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;

  /// The (scalar) effective creep strain
  ADMaterialProperty<Real> & _ec;
  const MaterialProperty<Real> & _ec_old;

  /// The flow direction
  ADMaterialProperty<RankTwoTensor> & _Nc;
};
