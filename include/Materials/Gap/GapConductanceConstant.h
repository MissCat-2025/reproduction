// //* This file is part of the MOOSE framework
// //* https://mooseframework.inl.gov
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html

// #pragma once

// #include "Material.h"
// #include "Function.h" // 添加函数支持

// class GapConductanceConstant : public Material
// {
// public:
//   static InputParameters validParams();

//   GapConductanceConstant(const InputParameters & parameters);

//   static InputParameters actionParameters();

//   virtual void computeQpProperties() override;

// protected:
//   const Real & _prescribed_gap_conductance;
//   const std::string _appended_property_name;
//   MaterialProperty<Real> & _gap_conductance;
//   MaterialProperty<Real> & _gap_conductance_dT;
  
//   // 新增：函数支持
//   const Function * const _gap_conductance_function;
//   const VariableValue * const _gap_conductance_function_variable;
// };