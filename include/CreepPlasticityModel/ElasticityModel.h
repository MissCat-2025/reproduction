//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "BaseNameInterface.h"

class CreepModel;
class ElasticityModel : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ElasticityModel(const InputParameters & parameters);

  /// Set the current quadrature point
  virtual void setQp(unsigned int qp);


  virtual void setCreepModel(CreepModel * creep_model);

  /**
   * Compute the stress given the mechanical strain. Also performs the plasticity update, if any.
   * @param mechanical_strain The mechanical strain
   * @param stress            The stress
   */
  virtual void updateState(const ADRankTwoTensor & mechanical_strain, ADRankTwoTensor & stress);

  /**
   * The stress-strain relation
   * @param strain The given strain
   * @return The computed stress given the strain and the constitutive relation
   */
  virtual ADRankTwoTensor computeStress(const ADRankTwoTensor & strain) = 0;

  // @{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are
  // unused in all inheriting classes and should not be overwritten.
  void resetQpProperties() final {}
  void resetProperties() final {}
  // @}

protected:
  virtual void initQpStatefulProperties() override;

  /// The optional creep model
  CreepModel * _creep_model;

  /// The elastic strain
  ADMaterialProperty<RankTwoTensor> & _elastic_strain;
};
