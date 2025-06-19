#include "reproductionApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
reproductionApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

reproductionApp::reproductionApp(InputParameters parameters) : MooseApp(parameters)
{
  reproductionApp::registerAll(_factory, _action_factory, _syntax);
}

reproductionApp::~reproductionApp() {}

void
reproductionApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<reproductionApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"reproductionApp"});
  Registry::registerActionsTo(af, {"reproductionApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
reproductionApp::registerApps()
{
  registerApp(reproductionApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
reproductionApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  reproductionApp::registerAll(f, af, s);
}
extern "C" void
reproductionApp__registerApps()
{
  reproductionApp::registerApps();
}
