// //* This file is part of the MOOSE framework
// //* https://mooseframework.inl.gov
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html

// #include "GapConductanceConstant.h"

// registerMooseObject("reproductionApp", GapConductanceConstant);

// InputParameters
// GapConductanceConstant::validParams()
// {
//   InputParameters params = Material::validParams();
//   params += GapConductanceConstant::actionParameters();
//   // We can't just make it required in the first place because then it would always
//   // be required in the Action, even if this model isn't used.
//   params.makeParamRequired<Real>("gap_conductance");
//   params.addClassDescription("Material to compute a constant, prescribed gap conductance");
//   // 新增：函数支持

//   return params;
// }

// InputParameters
// GapConductanceConstant::actionParameters()
// {
//   InputParameters params = emptyInputParameters();
//   params.addParam<std::string>(
//       "appended_property_name", "", "Name appended to material properties to make them unique");
//   params.addParam<Real>("gap_conductance", 0.0, "Gap conductance");
//   params.addParamNamesToGroup("gap_conductance", "Gap conductivity");
//   params.addParam<FunctionName>(
//     "gap_conductance_function",
//     "Gap conductance as a function. If provided, this will override the prescribed_gap_conductance.");
//   params.addCoupledVar("gap_conductance_function_variable",
//                      "Variable to be used in the gap_conductance_function in place of time");

//   return params;
// }

// GapConductanceConstant::GapConductanceConstant(const InputParameters & parameters)
//   : Material(parameters),
//     _prescribed_gap_conductance(getParam<Real>("gap_conductance")),
//     _appended_property_name(getParam<std::string>("appended_property_name")),
//     _gap_conductance(declareProperty<Real>("gap_conductance" + _appended_property_name)),
//     _gap_conductance_dT(declareProperty<Real>("gap_conductance" + _appended_property_name + "_dT")),
//     // 新增：函数支持
//     _gap_conductance_function(isParamValid("gap_conductance_function")
//                                   ? &getFunction("gap_conductance_function")
//                                   : nullptr),
//     _gap_conductance_function_variable(isCoupled("gap_conductance_function_variable")
//                                            ? &coupledValue("gap_conductance_function_variable")
//                                            : nullptr)
// {
//   // if (!params.isParamSetByUser("gap_conductance"))
//   //   mooseError("gap_conductance must be specified");
// }


// void
// GapConductanceConstant::computeQpProperties()
// {
//   Real conductance_value = _prescribed_gap_conductance;
  
//   // 如果提供了函数，使用函数值
//   if (_gap_conductance_function)
//   {
//     if (_gap_conductance_function_variable)
//       conductance_value = _gap_conductance_function->value(
//           (*_gap_conductance_function_variable)[_qp], _q_point[_qp]);
//     else
//       conductance_value = _gap_conductance_function->value(_t, _q_point[_qp]);
//   }
  
//   _gap_conductance[_qp] = conductance_value;
//   _gap_conductance_dT[_qp] = 0.0; // 假设函数不依赖于温度
// }
