//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeSmallDeformationStressMod.h"
#include "SmallDeformationElasticityModelMod.h"
#include "SmallDeformationPlasticityModelMod.h"
/*这正是前向声明 + 完整包含的经典用法：
头文件（.h）：只用前向声明（只需要指针）
源文件（.C）：包含完整头文件（需要调用方法）
可以调用 _elasticity_model->setQp(_qp)
可以调用 _elasticity_model->updateState(...)
可以调用 _elasticity_model->setPlasticityModel(...)
*/
registerMooseObject("reproductionApp", ComputeSmallDeformationStressMod);

InputParameters
ComputeSmallDeformationStressMod::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The stress calculator given an elasticity model and a plasticity "
                             "model. Small deformation is assumed.");

  params.addRequiredParam<MaterialName>("elasticity_model",
                                        "Name of the elastic stress-strain constitutive model");
  params.addParam<MaterialName>("plasticity_model", "Name of the plasticity model");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ComputeSmallDeformationStressMod::ComputeSmallDeformationStressMod(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _mechanical_strain(getADMaterialProperty<RankTwoTensor>(prependBaseName("mechanical_strain"))),
    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress")))
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
ComputeSmallDeformationStressMod::initialSetup()
{

  _elasticity_model =
      dynamic_cast<SmallDeformationElasticityModelMod *>(&getMaterial("elasticity_model"));
  if (!_elasticity_model)
    paramError("elasticity_model",
               "Elasticity model " + getParam<MaterialName>("elasticity_model") +
                   " is not compatible with ComputeSmallDeformationStressMod");
  //条件 ? 真值 : 假值
  /*
  检查输入参数中是否提供了"plasticity_model"
  isParamValid("plasticity_model")
  
  getMaterial("plasticity_model")：从MOOSE框架获取名为"plasticity_model"的材料对象
  &：取地址，得到指针
  dynamic_cast<SmallDeformationPlasticityModel *>：安全地将指针转换为特定类型

  设置为空指针nullptr
  */
  _plasticity_model =
      isParamValid("plasticity_model")
          ? dynamic_cast<SmallDeformationPlasticityModelMod *>(&getMaterial("plasticity_model"))
          : nullptr;
  if (_plasticity_model)
    _elasticity_model->setPlasticityModel(_plasticity_model);
}

void
ComputeSmallDeformationStressMod::initQpStatefulProperties()
{
  _stress[_qp].zero();
}

void
ComputeSmallDeformationStressMod::computeQpProperties()
{
  _elasticity_model->setQp(_qp);
  _elasticity_model->updateState(_mechanical_strain[_qp], _stress[_qp]);
} 