#include "ComputeCreepPlasticityDeformationStressNew.h"

#include "ElasticityModelNonAD.h"

#include <algorithm>

registerMooseObject("reproductionApp", ComputeCreepPlasticityDeformationStressNew);

InputParameters
ComputeCreepPlasticityDeformationStressNew::validParams()
{
  InputParameters params = ComputeMultipleInelasticStressBase::validParams();
  params.addClassDescription(
      "Compute stress using ComputeMultipleInelasticStress and drive an optional "
      "ElasticityModel to compute energy-based outputs from the elastic strain.");

  params.addParam<std::vector<MaterialName>>(
      "inelastic_models",
      {},
      "The material objects to use to calculate stress and inelastic strains. "
      "Note: specify creep models first and plasticity models second. "
      "Leave empty to compute purely elastic stress.");

  params.addParam<MaterialName>(
      "elasticity_model",
      "Optional ElasticityModel to evaluate using the computed elastic strain.");

  params.addParam<bool>(
      "use_elasticity_model_for_stress",
      false,
      "If true, replaces the main mechanics stress with the stress returned by the provided "
      "elasticity_model (computed from the elastic strain). This is useful for tension/compression "
      "split models that should also drive the equilibrium stress.");

  params.addParam<bool>("apply_damage_degradation",
                        false,
                        "If true, multiplies the computed stress (and Jacobian multiplier, when "
                        "available) by a degradation function g(d) evaluated as a material property. "
                        "This enables damage to affect the main mechanics equilibrium.");
  params.addParam<MaterialPropertyName>(
      "degradation_function", "g", "Material property name for the degradation function g(d)");
  params.addParam<Real>("min_degradation",
                        1e-6,
                        "Lower bound applied to g(d) when degrading the main mechanics stress. "
                        "Helps prevent singular/ill-conditioned systems as d approaches 1.");

  params.addParam<bool>("debug", false, "Enable ComputeCreepPlasticityDeformationStressNew debug output");
  params.addParam<int>("debug_qp", -1, "Only print debug for this qp (-1 prints all qps)");
  params.addParam<unsigned int>("debug_step_interval",
                                1,
                                "Print debug every N time steps (only when debug=true)");

  return params;
}

ComputeCreepPlasticityDeformationStressNew::ComputeCreepPlasticityDeformationStressNew(
    const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
    _has_elasticity_model(parameters.isParamSetByUser("elasticity_model")),
    _elasticity_model(nullptr),
    _use_elasticity_model_for_stress(getParam<bool>("use_elasticity_model_for_stress")),
    _apply_damage_degradation(getParam<bool>("apply_damage_degradation")),
    _g_name(getParam<MaterialPropertyName>("degradation_function")),
    _g(nullptr),
    _min_degradation(getParam<Real>("min_degradation")),
    _debug(getParam<bool>("debug")),
    _debug_qp(getParam<int>("debug_qp")),
    _debug_step_interval(getParam<unsigned int>("debug_step_interval"))
{
  if (_apply_damage_degradation)
    _g = &getMaterialProperty<Real>(_g_name);
}

void
ComputeCreepPlasticityDeformationStressNew::initialSetup()
{
  ComputeMultipleInelasticStress::initialSetup();

  if (_has_elasticity_model)
  {
    _elasticity_model = dynamic_cast<ElasticityModelNonAD *>(&getMaterial("elasticity_model"));
    if (!_elasticity_model)
      paramError("elasticity_model",
                 "Elasticity model " + getParam<MaterialName>("elasticity_model") +
                     " is not compatible with ComputeCreepPlasticityDeformationStressNew");
  }
}

void
ComputeCreepPlasticityDeformationStressNew::computeQpStress()
{
  ComputeMultipleInelasticStress::computeQpStress();

  if (_has_elasticity_model)
  {
    _elasticity_model->setQp(_qp);
    if (_use_elasticity_model_for_stress)
    {
      _stress[_qp] = _elasticity_model->computeStress(_elastic_strain[_qp]);
      if (_fe_problem.currentlyComputingJacobian())
        _Jacobian_mult[_qp] = _elasticity_model->computeJacobianMult(_elastic_strain[_qp]);
    }
    else
      _elasticity_model->computeStress(_elastic_strain[_qp]);
  }

  if (_apply_damage_degradation)
  {
    const Real g_eff = std::max((*_g)[_qp], _min_degradation);
    _stress[_qp] *= g_eff;
    if (_fe_problem.currentlyComputingJacobian())
      _Jacobian_mult[_qp] *= g_eff;
  }

  if (_debug && (_debug_qp < 0 || static_cast<int>(_qp) == _debug_qp) &&
      (_debug_step_interval == 0 || (_t_step % _debug_step_interval == 0)))
  {
    const auto & e = _elastic_strain[_qp];
    const auto & inel = _inelastic_strain[_qp];
    const auto & de = _strain_increment[_qp];
    const auto & eo = _elastic_strain_old[_qp];
    const auto & inelo = _inelastic_strain_old[_qp];
    const auto & s = _stress[_qp];

    _console << "ComputeCreepPlasticityDeformationStressNew qp=" << _qp << " t_step=" << _t_step
             << " t=" << _t << " dt=" << _fe_problem.dt()
             << " de_00=" << de(0, 0) << " de_11=" << de(1, 1) << " de_22=" << de(2, 2)
             << " de_01=" << de(0, 1) << " de_12=" << de(1, 2) << " de_02=" << de(0, 2)
             << " e_old_00=" << eo(0, 0) << " e_old_11=" << eo(1, 1) << " e_old_22=" << eo(2, 2)
             << " e_old_01=" << eo(0, 1) << " e_old_12=" << eo(1, 2) << " e_old_02=" << eo(0, 2)
             << " inel_old_00=" << inelo(0, 0) << " inel_old_11=" << inelo(1, 1)
             << " inel_old_22=" << inelo(2, 2) << " inel_old_01=" << inelo(0, 1)
             << " inel_old_12=" << inelo(1, 2) << " inel_old_02=" << inelo(0, 2)
             << " e_00=" << e(0, 0) << " e_11=" << e(1, 1) << " e_22=" << e(2, 2)
             << " e_01=" << e(0, 1) << " e_12=" << e(1, 2) << " e_02=" << e(0, 2)
             << " inel_00=" << inel(0, 0) << " inel_11=" << inel(1, 1) << " inel_22=" << inel(2, 2)
             << " inel_01=" << inel(0, 1) << " inel_12=" << inel(1, 2) << " inel_02=" << inel(0, 2)
             << " stress_00=" << s(0, 0) << " stress_11=" << s(1, 1) << " stress_22=" << s(2, 2)
             << " stress_01=" << s(0, 1) << " stress_12=" << s(1, 2) << " stress_02=" << s(0, 2)
             << std::endl;
  }
}
