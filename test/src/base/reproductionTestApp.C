//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "reproductionTestApp.h"
#include "reproductionApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
reproductionTestApp::validParams()
{
  InputParameters params = reproductionApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

reproductionTestApp::reproductionTestApp(InputParameters parameters) : MooseApp(parameters)
{
  reproductionTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

reproductionTestApp::~reproductionTestApp() {}

void
reproductionTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  reproductionApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"reproductionTestApp"});
    Registry::registerActionsTo(af, {"reproductionTestApp"});
  }
}

void
reproductionTestApp::registerApps()
{
  registerApp(reproductionApp);
  registerApp(reproductionTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
reproductionTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  reproductionTestApp::registerAll(f, af, s);
}
extern "C" void
reproductionTestApp__registerApps()
{
  reproductionTestApp::registerApps();
}
