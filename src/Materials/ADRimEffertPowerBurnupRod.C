// ADRimEffertPowerBurnupRod.C
#include "ADRimEffertPowerBurnupRod.h"
#include "Function.h"  // 或者 #include "FunctionInterface.h"
registerMooseObject("reproductionApp", ADRimEffertPowerBurnupRod);

InputParameters
ADRimEffertPowerBurnupRod::validParams()
{
  /*
  [1] Deng Y, Wu Y, Zhang D, et al. Development of a thermal–mechanical behavior coupling analysis code for a dual-cooled annular fuel element in PWRs[J]. Nuclear Engineering and Design, 2016, 301: 353-365.
    公式:
    \mathrm{term1}\ =\ 3.5\ Bu^{0.55}\exp{\left(-180\ \ \left(r\ -\ r_i\right)^{0.56}\right)}
\mathrm{term2}\ =\ 4.5\ Bu^{0.55}\exp{\left(-115\ \ \left(r_o\ -\ r\right)^{0.5}\right)}
\mathrm{term3}=\left(-15.7Bu^2+3.5Bu\right)-1
改编
    */
  InputParameters params = ADMaterial::validParams();
  params.addRequiredParam<FunctionName>("power_history", "功率历史函数");
  params.addRequiredParam<Real>("pellet_outer_radius", "棒状芯块半径 (m)");
  params.addParam<Real>("A", 4.0, "径向功率外缘峰值系数");
  params.addParam<Real>("B", -69.84627, "径向指数衰减系数 p2");
  params.addParam<Real>("C", -4.7, "燃耗二次项系数");
  params.addParam<Real>("D", 3.0, "燃耗一次项系数");
  params.addClassDescription("棒状燃料的径向功率+燃耗边缘效应 AD 材料");
  return params;
}

ADRimEffertPowerBurnupRod::ADRimEffertPowerBurnupRod(const InputParameters & parameters)
  : ADMaterial(parameters),
    _power_history(getFunctionByName("power_history")),
    _pellet_outer_radius(getParam<Real>("pellet_outer_radius")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _C(getParam<Real>("C")),
    _D(getParam<Real>("D")),
    _total_power(declareADProperty<Real>("total_power")),
    _radial_power_shape(declareADProperty<Real>("radial_power_shape")),
    _density(getADMaterialProperty<Real>("density")),
    _burnup(declareADProperty<Real>("burnup")),
    _burnup_old(getMaterialPropertyOld<Real>("burnup"))
{
}
void
ADRimEffertPowerBurnupRod::initQpStatefulProperties()
{
  // 初始化燃耗为0
  _burnup[_qp] = 1e-12;
}
ADReal
ADRimEffertPowerBurnupRod::powerFactor(const Real & r) const
{
  const ADReal Bu = _burnup_old[_qp];
  ADReal term1 =
      _A * std::pow(Bu, 0.5) * std::exp(_B * std::pow((_pellet_outer_radius - r), 0.51));
  ADReal term2 = (_C * Bu * Bu + _D * Bu) - 0.925;
  return term1 - term2;
}

void
ADRimEffertPowerBurnupRod::computeQpProperties()
{
  // 计算相对半径
  const Point & p = _q_point[_qp];
  Real r = std::sqrt(p(0) * p(0) + p(1) * p(1));
  // 计算功率分布形状
  _radial_power_shape[_qp] = powerFactor(r);
  
  // 计算当前时刻的基础功率密度 (W/m³)
  Point curr_point(_q_point[_qp]);
  const Real power_base = _power_history.value(_t, curr_point);
  // 计算总功率密度
  _total_power[_qp] = power_base * _radial_power_shape[_qp];

  //[1] Hales J D, Williamson R L, Novascone S R, et al. BISON Theory Manual The Equations behind Nuclear Fuel Analysis: INL/EXT--13-29930, 1374503[R]. 2016: INL/EXT--13-29930, 1374503.

  //燃耗部分定义
  const Real N_av = 6.022e23;    // 阿伏伽德罗常数
  const Real M_w = 270.0;        // UO2分子量 (g/mol)
  const Real alpha = 3.2845e-11; // 每次裂变释放的能量 (J/fission)
  // 计算初始重金属原子密度 (atoms/m³)
  ADReal N_f0 = _density[_qp] * N_av / M_w * 1000.0;
  const Real _dt = _fe_problem.dt();
  // 更新总燃耗
  _burnup[_qp] = _burnup_old[_qp] + _total_power[_qp] / alpha * _dt / N_f0;
}
