#include "ElasticityModelNonAD.h"

InputParameters
ElasticityModelNonAD::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addParam<bool>("compute", false, "The base class does not compute by itself");
  params.suppressParameter<bool>("compute");
  return params;
}

ElasticityModelNonAD::ElasticityModelNonAD(const InputParameters & parameters)
  : Material(parameters), BaseNameInterface(parameters), _compute(getParam<bool>("compute"))
{
}
