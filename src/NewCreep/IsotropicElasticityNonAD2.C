// #include "IsotropicElasticityNonAD.h"

// #include <algorithm>
// #include <cmath>

// registerMooseObject("reproductionApp", IsotropicElasticityNonAD);

// namespace
// {
// void
// computeIsotropicElasticConstants(const Real E,
//                                  const Real nu,
//                                  const unsigned int dim,
//                                  const MooseEnum & kinematic_assumption,
//                                  Real & G,
//                                  Real & lambda,
//                                  Real & K)
// {
//   G = E / (2.0 * (1.0 + nu));

//   if (dim == 2)
//   {
//     if (kinematic_assumption == "PLANE_STRESS")
//       lambda = E * nu / (1.0 - nu * nu);
//     else
//       lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

//     K = lambda + (2.0 * G) / dim;
//   }
//   else
//   {
//     lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
//     K = lambda + (2.0 * G) / dim;
//   }
// }

// inline Real
// macaulay(const Real x)
// {
//   return x > 0.0 ? x : 0.0;
// }

// RankTwoTensor
// spectralPositivePart(const RankTwoTensor & r2t)
// {
//   RankTwoTensor eigvecs;
//   std::vector<Real> eigvals(RankTwoTensor::N);
//   r2t.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

//   for (auto & x : eigvals)
//     x = macaulay(x);

//   RankTwoTensor eigvals_pos;
//   eigvals_pos.fillFromInputVector(eigvals);
//   return eigvecs * eigvals_pos * eigvecs.transpose();
// }
// }

// InputParameters
// IsotropicElasticityNonAD::validParams()
// {
//   InputParameters params = ElasticityModelNonAD::validParams();
//   params.addClassDescription(
//       "Non-AD isotropic elasticity utility for energy outputs (psie, psie_active) with optional "
//       "damage degradation function.");

//   params.addRequiredParam<Real>("youngs_modulus", "Young's modulus E");
//   params.addRequiredParam<Real>("poissons_ratio", "Poisson's ratio nu");

//   params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
//   params.addParam<MaterialPropertyName>(
//       "strain_energy_density",
//       "psie",
//       "Name of the strain energy density computed by this material model");
//   params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");

//   MooseEnum decomp_options("NONE SPECTRAL VOLDEV MAXPRINCIPAL", "NONE");
//   params.addParam<MooseEnum>("decomposition", decomp_options, "The decomposition method");

//   MooseEnum kinematic_options("PLANE_STRAIN PLANE_STRESS", "PLANE_STRAIN");
//   params.addParam<MooseEnum>(
//       "kinematic_assumption", kinematic_options, "Kinematic assumption for 2D problems");

//   params.addParam<bool>("use_threshold", false, "Whether to use threshold energy");
//   params.addParam<MaterialPropertyName>(
//       "tensile_strength", "sigma0", "Critical stress for threshold");

//   params.addParam<bool>("use_history_max",
//                         false,
//                         "Whether to use history maximum variable H for damage irreversibility");
  
//   params.addParam<bool>("degrade_out_of_plane_strain",
//                         true,
//                         "If false, the out-of-plane (zz) stress component will NOT be degraded by damage. "
//                         "Useful for Generalized Plane Strain to avoid singularity.");

//   params.addParam<bool>("debug", false, "Enable IsotropicElasticityNonAD debug output");
//   params.addParam<int>("debug_qp", -1, "Only print debug for this qp (-1 prints all qps)");
//   params.addParam<unsigned int>(
//       "debug_step_interval", 1, "Print debug every N time steps (only when debug=true)");

//   return params;
// }

// IsotropicElasticityNonAD::IsotropicElasticityNonAD(const InputParameters & parameters)
//   : ElasticityModelNonAD(parameters),
//     DerivativeMaterialPropertyNameInterface(),
//     _youngs_modulus(getParam<Real>("youngs_modulus")),
//     _poissons_ratio(getParam<Real>("poissons_ratio")),
//     _phase_field(coupledValue("phase_field")),
//     _d_name(getVar("phase_field", 0)->name()),
//     _g_name(getParam<MaterialPropertyName>("degradation_function")),
//     _g(getMaterialProperty<Real>(_g_name)),
//     _dg_dd(getMaterialProperty<Real>(derivativePropertyNameFirst(_g_name, _d_name))),
//     _psie_active(declareProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density") + "_active")),
//     _psie(declareProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density"))),
//     _dpsie_dd(declareProperty<Real>("dpsie_dd")),
//     _kinematic_assumption(getParam<MooseEnum>("kinematic_assumption")),
//     _decomposition(getParam<MooseEnum>("decomposition")),
//     _use_threshold(getParam<bool>("use_threshold")),
//     _use_history_max(getParam<bool>("use_history_max")),
//     _tensile_strength(_use_threshold ? getMaterialProperty<Real>(getParam<MaterialPropertyName>("tensile_strength"))
//                                      : declareProperty<Real>("dummy_sigma_c")),
//     _history_max(_use_history_max ? &declareProperty<Real>("history_max") : nullptr),
//     _history_max_old(_use_history_max ? &getMaterialPropertyOld<Real>("history_max") : nullptr),
//     _debug(getParam<bool>("debug")),
//     _debug_qp(getParam<int>("debug_qp")),
//     _debug_step_interval(getParam<unsigned int>("debug_step_interval")),
//     _degrade_out_of_plane_strain(getParam<bool>("degrade_out_of_plane_strain"))
// {
// }

// void
// IsotropicElasticityNonAD::initQpStatefulProperties()
// {
//   if (_use_history_max)
//     (*_history_max)[_qp] = 0.0;
// }

// void
// IsotropicElasticityNonAD::updateHistoryMax(const Real Y_bar)
// {
//   if (_use_history_max)
//     (*_history_max)[_qp] = std::max((*_history_max_old)[_qp], Y_bar);
// }

// Real
// IsotropicElasticityNonAD::applyThreshold(const Real psie_active_raw) const
// {
//   Real psie_active_final = psie_active_raw;

//   if (_use_history_max)
//     psie_active_final = (*_history_max)[_qp];

//   if (_use_threshold)
//   {
//     const Real threshold = 0.5 * _tensile_strength[_qp] * _tensile_strength[_qp] / _youngs_modulus;
//     psie_active_final = std::max(psie_active_final, threshold);
//   }

//   return psie_active_final;
// }

// RankTwoTensor
// IsotropicElasticityNonAD::computeStress(const RankTwoTensor & elastic_strain)
// {
//   if (_decomposition == "NONE")
//     return computeStressNone(elastic_strain);
//   if (_decomposition == "SPECTRAL")
//     return computeStressSpectralDecomposition(elastic_strain);
//   if (_decomposition == "VOLDEV")
//     return computeStressVolDevDecomposition(elastic_strain);
//   if (_decomposition == "MAXPRINCIPAL")
//     return computeStressMaxPrincipalDecomposition(elastic_strain);

//   mooseError("Unknown decomposition type in IsotropicElasticityNonAD");
// }

// RankFourTensor
// IsotropicElasticityNonAD::computeJacobianMult(const RankTwoTensor & elastic_strain)
// {
//   Real G, lambda, K;
//   computeIsotropicElasticConstants(
//       _youngs_modulus, _poissons_ratio, _mesh.dimension(), _kinematic_assumption, G, lambda, K);

//   const RankTwoTensor I2(RankTwoTensor::initIdentity);
//   const RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
//   const RankFourTensor I2I2 = I2.outerProduct(I2);
//   const RankFourTensor C_intact = lambda * I2I2 + 2.0 * G * I4sym;

//   if (_decomposition == "NONE" || _decomposition == "MAXPRINCIPAL")
//   {
//     RankFourTensor C = _g[_qp] * C_intact;
//     if (!_degrade_out_of_plane_strain && _kinematic_assumption == "PLANE_STRAIN")
//       C(2, 2, 2, 2) = C_intact(2, 2, 2, 2); // Restore C3333
//     return C;
//   }

//   if (_decomposition == "VOLDEV")
//   {
//     const Real strain_tr = elastic_strain.trace();
//     const Real strain_tr_neg = std::min(strain_tr, 0.0);

//     RankFourTensor Jacobian_neg;
//     if (strain_tr_neg < 0.0)
//     {
//       const Real Kvol = lambda + 2.0 * G / 3.0;
//       Jacobian_neg = Kvol * I2I2;
//     }

//     const RankFourTensor Jacobian_pos = C_intact - Jacobian_neg;
//     RankFourTensor Jacobian = _g[_qp] * Jacobian_pos + Jacobian_neg;
    
//     if (!_degrade_out_of_plane_strain && _kinematic_assumption == "PLANE_STRAIN")
//       Jacobian(2, 2, 2, 2) = C_intact(2, 2, 2, 2); // Restore C3333

//     return Jacobian;
//   }

//   if (_decomposition == "SPECTRAL")
//   {
//     std::vector<Real> eigval;
//     RankTwoTensor eigvec;
//     const RankFourTensor Ppos =
//         elastic_strain.positiveProjectionEigenDecomposition(eigval, eigvec);

//     RankFourTensor Jacobian = (I4sym - (1.0 - _g[_qp]) * Ppos) * C_intact;

//     if (!_degrade_out_of_plane_strain && _kinematic_assumption == "PLANE_STRAIN")
//       Jacobian(2, 2, 2, 2) = C_intact(2, 2, 2, 2); // Restore C3333

//     return Jacobian;
//   }

//   mooseError("Unknown decomposition type in IsotropicElasticityNonAD");
// }

// RankTwoTensor
// IsotropicElasticityNonAD::computeStressNone(const RankTwoTensor & strain)
// {
//   Real G, lambda, K;
//   computeIsotropicElasticConstants(
//       _youngs_modulus, _poissons_ratio, _mesh.dimension(), _kinematic_assumption, G, lambda, K);

//   const RankTwoTensor I2(RankTwoTensor::initIdentity);
//   const RankTwoTensor stress_intact = lambda * strain.trace() * I2 + 2.0 * G * strain;
  
//   RankTwoTensor stress = _g[_qp] * stress_intact;

//   if (!_degrade_out_of_plane_strain && _kinematic_assumption == "PLANE_STRAIN")
//   {
//       // Restore the zz component to its intact value
//       stress(2, 2) = stress_intact(2, 2);
//   }

//   const Real strain_tr = strain.trace();
//   const Real psie_active_raw =
//       0.5 * lambda * strain_tr * strain_tr + G * strain.doubleContraction(strain);

//   updateHistoryMax(psie_active_raw);
//   _psie_active[_qp] = applyThreshold(psie_active_raw);
//   _psie[_qp] = _g[_qp] * _psie_active[_qp];
//   _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

//   if (_debug &&
//       (_debug_qp < 0 || static_cast<int>(_qp) == _debug_qp) &&
//       (_debug_step_interval == 0 || (_t_step % _debug_step_interval == 0)))
//   {
//     _console << "IsotropicElasticityNonAD::NONE qp=" << _qp << " t_step=" << _t_step
//              << " g=" << _g[_qp] << " dg_dd=" << _dg_dd[_qp] << " strain_tr=" << strain_tr
//              << " psie_active_raw=" << psie_active_raw << " psie_active=" << _psie_active[_qp]
//              << " psie=" << _psie[_qp] << " dpsie_dd=" << _dpsie_dd[_qp]
//              << " strain_00=" << strain(0, 0) << " strain_11=" << strain(1, 1)
//              << " strain_22=" << strain(2, 2) << " strain_01=" << strain(0, 1)
//              << " strain_12=" << strain(1, 2) << " strain_02=" << strain(0, 2)
//              << " stress_00=" << stress(0, 0) << " stress_11=" << stress(1, 1)
//              << " stress_22=" << stress(2, 2) << " stress_01=" << stress(0, 1)
//              << " stress_12=" << stress(1, 2) << " stress_02=" << stress(0, 2) << std::endl;
//   }

//   return stress;
// }

// RankTwoTensor
// IsotropicElasticityNonAD::computeStressSpectralDecomposition(const RankTwoTensor & strain)
// {
//   Real G, lambda, K;
//   computeIsotropicElasticConstants(
//       _youngs_modulus, _poissons_ratio, _mesh.dimension(), _kinematic_assumption, G, lambda, K);

//   const RankTwoTensor I2(RankTwoTensor::initIdentity);
//   const Real strain_tr = strain.trace();
//   const Real strain_tr_pos = macaulay(strain_tr);

//   const RankTwoTensor strain_pos = spectralPositivePart(strain);

//   const RankTwoTensor stress_intact = lambda * strain_tr * I2 + 2.0 * G * strain;
//   const RankTwoTensor stress_pos = lambda * strain_tr_pos * I2 + 2.0 * G * strain_pos;
//   const RankTwoTensor stress_neg = stress_intact - stress_pos;
//   RankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

//   if (!_degrade_out_of_plane_strain && _kinematic_assumption == "PLANE_STRAIN")
//   {
//       // Restore the zz component to its intact value
//       stress(2, 2) = stress_intact(2, 2);
//   }

//   const Real psie_intact =
//       0.5 * lambda * strain_tr * strain_tr + G * strain.doubleContraction(strain);
//   const Real psie_active_raw =
//       0.5 * lambda * strain_tr_pos * strain_tr_pos + G * strain_pos.doubleContraction(strain_pos);

//   updateHistoryMax(psie_active_raw);
//   _psie_active[_qp] = applyThreshold(psie_active_raw);

//   const Real psie_inactive = psie_intact - psie_active_raw;
//   _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
//   _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

//   if (_debug &&
//       (_debug_qp < 0 || static_cast<int>(_qp) == _debug_qp) &&
//       (_debug_step_interval == 0 || (_t_step % _debug_step_interval == 0)))
//   {
//     _console << "IsotropicElasticityNonAD::SPECTRAL qp=" << _qp << " t_step=" << _t_step
//              << " g=" << _g[_qp] << " dg_dd=" << _dg_dd[_qp] << " strain_tr=" << strain_tr
//              << " strain_tr_pos=" << strain_tr_pos << " psie_intact=" << psie_intact
//              << " psie_active_raw=" << psie_active_raw << " psie_active=" << _psie_active[_qp]
//              << " psie=" << _psie[_qp] << " dpsie_dd=" << _dpsie_dd[_qp]
//              << " strain_00=" << strain(0, 0) << " strain_11=" << strain(1, 1)
//              << " strain_22=" << strain(2, 2) << " strain_01=" << strain(0, 1)
//              << " strain_12=" << strain(1, 2) << " strain_02=" << strain(0, 2)
//              << " stress_00=" << stress(0, 0) << " stress_11=" << stress(1, 1)
//              << " stress_22=" << stress(2, 2) << " stress_01=" << stress(0, 1)
//              << " stress_12=" << stress(1, 2) << " stress_02=" << stress(0, 2) << std::endl;
//   }

//   return stress;
// }

// RankTwoTensor
// IsotropicElasticityNonAD::computeStressVolDevDecomposition(const RankTwoTensor & strain)
// {
//   Real G, lambda, K;
//   computeIsotropicElasticConstants(
//       _youngs_modulus, _poissons_ratio, _mesh.dimension(), _kinematic_assumption, G, lambda, K);

//   const RankTwoTensor I2(RankTwoTensor::initIdentity);
//   const Real Kvol = lambda + 2.0 * G / 3.0;

//   const Real strain_tr = strain.trace();
//   const Real strain_tr_pos = macaulay(strain_tr);
//   const Real strain_tr_neg = strain_tr - strain_tr_pos;

//   const RankTwoTensor stress_intact = lambda * strain_tr * I2 + 2.0 * G * strain;
//   const RankTwoTensor stress_neg = Kvol * strain_tr_neg * I2;
//   const RankTwoTensor stress_pos = stress_intact - stress_neg;
//   RankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

//   if (!_degrade_out_of_plane_strain && _kinematic_assumption == "PLANE_STRAIN")
//   {
//       // Restore the zz component to its intact value
//       stress(2, 2) = stress_intact(2, 2);
//   }

//   const Real psie_intact =
//       0.5 * lambda * strain_tr * strain_tr + G * strain.doubleContraction(strain);
//   const Real psie_inactive = 0.5 * Kvol * strain_tr_neg * strain_tr_neg;
//   const Real psie_active_raw = psie_intact - psie_inactive;

//   updateHistoryMax(psie_active_raw);
//   _psie_active[_qp] = applyThreshold(psie_active_raw);
//   _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
//   _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

//   if (_debug &&
//       (_debug_qp < 0 || static_cast<int>(_qp) == _debug_qp) &&
//       (_debug_step_interval == 0 || (_t_step % _debug_step_interval == 0)))
//   {
//     _console << "IsotropicElasticityNonAD::VOLDEV qp=" << _qp << " t_step=" << _t_step
//              << " g=" << _g[_qp] << " dg_dd=" << _dg_dd[_qp] << " strain_tr=" << strain_tr
//              << " strain_tr_pos=" << strain_tr_pos << " strain_tr_neg=" << strain_tr_neg
//              << " psie_intact=" << psie_intact << " psie_inactive=" << psie_inactive
//              << " psie_active_raw=" << psie_active_raw << " psie_active=" << _psie_active[_qp]
//              << " psie=" << _psie[_qp] << " dpsie_dd=" << _dpsie_dd[_qp]
//              << " strain_00=" << strain(0, 0) << " strain_11=" << strain(1, 1)
//              << " strain_22=" << strain(2, 2) << " strain_01=" << strain(0, 1)
//              << " strain_12=" << strain(1, 2) << " strain_02=" << strain(0, 2)
//              << " stress_00=" << stress(0, 0) << " stress_11=" << stress(1, 1)
//              << " stress_22=" << stress(2, 2) << " stress_01=" << stress(0, 1)
//              << " stress_12=" << stress(1, 2) << " stress_02=" << stress(0, 2) << std::endl;
//   }

//   return stress;
// }

// RankTwoTensor
// IsotropicElasticityNonAD::computeStressMaxPrincipalDecomposition(const RankTwoTensor & strain)
// {
//   Real G, lambda, K;
//   computeIsotropicElasticConstants(
//       _youngs_modulus, _poissons_ratio, _mesh.dimension(), _kinematic_assumption, G, lambda, K);

//   const RankTwoTensor I2(RankTwoTensor::initIdentity);
//   const RankTwoTensor stress_intact = lambda * strain.trace() * I2 + 2.0 * G * strain;
//   const RankTwoTensor stress = _g[_qp] * stress_intact;

//   std::vector<Real> eigenvals(RankTwoTensor::N);
//   stress_intact.symmetricEigenvalues(eigenvals);

//   Real sigma_bar_eq = eigenvals[0];
//   for (unsigned int i = 1; i < RankTwoTensor::N; ++i)
//     sigma_bar_eq = std::max(sigma_bar_eq, eigenvals[i]);

//   sigma_bar_eq = macaulay(sigma_bar_eq);
//   const Real Y_bar = 0.5 * sigma_bar_eq * sigma_bar_eq / _youngs_modulus;

//   updateHistoryMax(Y_bar);
//   _psie_active[_qp] = applyThreshold(Y_bar);
//   _psie[_qp] = _g[_qp] * _psie_active[_qp];
//   _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

//   if (_debug &&
//       (_debug_qp < 0 || static_cast<int>(_qp) == _debug_qp) &&
//       (_debug_step_interval == 0 || (_t_step % _debug_step_interval == 0)))
//   {
//     _console << "IsotropicElasticityNonAD::MAXPRINCIPAL qp=" << _qp << " t_step=" << _t_step
//              << " g=" << _g[_qp] << " dg_dd=" << _dg_dd[_qp] << " sigma_bar_eq=" << sigma_bar_eq
//              << " Y_bar=" << Y_bar << " psie_active=" << _psie_active[_qp] << " psie=" << _psie[_qp]
//              << " dpsie_dd=" << _dpsie_dd[_qp] << " strain_00=" << strain(0, 0)
//              << " strain_11=" << strain(1, 1) << " strain_22=" << strain(2, 2)
//              << " strain_01=" << strain(0, 1) << " strain_12=" << strain(1, 2)
//              << " strain_02=" << strain(0, 2) << " stress_00=" << stress(0, 0)
//              << " stress_11=" << stress(1, 1) << " stress_22=" << stress(2, 2)
//              << " stress_01=" << stress(0, 1) << " stress_12=" << stress(1, 2)
//              << " stress_02=" << stress(0, 2) << std::endl;
//   }

//   return stress;
// }
