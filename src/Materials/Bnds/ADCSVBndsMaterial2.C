#include "ADCSVBndsMaterial2.h"
#include "MooseUtils.h"
#include <fstream>
#include <sstream>
#include <algorithm>

registerMooseObject("reproductionApp", ADCSVBndsMaterial2);

InputParameters
ADCSVBndsMaterial2::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredParam<FileName>("csv_file", "CSV文件路径");
  params.addParam<std::string>("data_column", "bnds_new", "要读取的数据列名");
  params.addParam<std::string>("x_column", "x", "x坐标列名");
  params.addParam<std::string>("y_column", "y", "y坐标列名");
  params.addParam<std::string>("property_name", "bnds", "输出的材料属性名");
  params.addParam<bool>("verbose", false, "是否输出详细信息");
  params.addClassDescription("从CSV文件读取晶界数据的AD材料");
  return params;
}

ADCSVBndsMaterial2::ADCSVBndsMaterial2(const InputParameters & parameters)
  : ADMaterial(parameters),
    _csv_file(getParam<FileName>("csv_file")),
    _data_column(getParam<std::string>("data_column")),
    _x_column(getParam<std::string>("x_column")),
    _y_column(getParam<std::string>("y_column")),
    _property_name(getParam<std::string>("property_name")),
    _verbose(getParam<bool>("verbose")),
    _csv_data(),
    _data_col_idx(-1),
    _x_col_idx(-1),
    _y_col_idx(-1),
    _bnds(declareADProperty<Real>(_property_name))
{
  // 读取CSV文件
  readCSVFile();
  
  // 验证必要的列索引
  if (_data_col_idx < 0)
    mooseError("错误：在CSV文件中找不到数据列 '", _data_column, "'");
  if (_x_col_idx < 0)
    mooseError("错误：在CSV文件中找不到X坐标列 '", _x_column, "'");
  if (_y_col_idx < 0)
    mooseError("错误：在CSV文件中找不到Y坐标列 '", _y_column, "'");
    
  if (processor_id() == 0)
  {
    _console << "成功读取 " << _csv_data.size() << " 行CSV数据" << std::endl;
    if (_verbose)
    {
      _console << "前5个数据点:" << std::endl;
      for (size_t i = 0; i < std::min(static_cast<size_t>(5), _csv_data.size()); ++i)
      {
        const auto& row = _csv_data[i];
        _console << "  [" << i << "] x=" << row[_x_col_idx] 
                 << ", y=" << row[_y_col_idx] 
                 << ", " << _data_column << "=" << row[_data_col_idx] << std::endl;
      }
    }
  }
}

void
ADCSVBndsMaterial2::initQpStatefulProperties()
{
  _bnds[_qp] = 0.0;
}

void
ADCSVBndsMaterial2::computeQpProperties()
{
  // 获取当前积分点的坐标
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);
  
  // 查找对应的bnds值
  _bnds[_qp] = findBndsValue(x, y);
}

void
ADCSVBndsMaterial2::readCSVFile()
{
  std::ifstream file(_csv_file);
  if (!file.is_open())
    mooseError("无法打开CSV文件: ", _csv_file);

  std::string line;
  std::vector<std::string> headers;
  bool first_line = true;
  
  while (std::getline(file, line))
  {
    // 移除行末字符
    while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
      line.pop_back();
      
    if (line.empty())
      continue;
    
    // 处理标题行
    if (first_line)
    {
      std::stringstream ss(line);
      std::string header;
      while (std::getline(ss, header, ','))
      {
        // 清理header
        header.erase(0, header.find_first_not_of(" \t\""));
        header.erase(header.find_last_not_of(" \t\"") + 1);
        headers.push_back(header);
      }
      
      // 查找列索引
      _data_col_idx = findColumnIndex(headers, _data_column);
      _x_col_idx = findColumnIndex(headers, _x_column);
      _y_col_idx = findColumnIndex(headers, _y_column);
      
      first_line = false;
      continue;
    }
    
    // 解析数据行
    std::stringstream ss(line);
    std::string token;
    std::vector<Real> row;
    
    try
    {
      while (std::getline(ss, token, ','))
      {
        // 清理token
        token.erase(0, token.find_first_not_of(" \t\""));
        token.erase(token.find_last_not_of(" \t\"") + 1);
        
        if (!token.empty())
        {
          Real value = MooseUtils::convert<Real>(token);
          row.push_back(value);
        }
        else
        {
          row.push_back(0.0);
        }
      }
      
      // 检查行完整性 - 只要有足够的列就添加
      if (row.size() > static_cast<size_t>(std::max({_data_col_idx, _x_col_idx, _y_col_idx})))
      {
        _csv_data.push_back(row);
      }
    }
    catch (const std::exception& e)
    {
      if (_verbose && processor_id() == 0)
        _console << "跳过错误行: " << line << " 错误: " << e.what() << std::endl;
    }
  }
  
  file.close();
  
  if (_csv_data.empty())
    mooseError("CSV文件中没有有效数据: ", _csv_file);
}

int
ADCSVBndsMaterial2::findColumnIndex(const std::vector<std::string>& headers, 
                                   const std::string& column_name) const
{
  for (size_t i = 0; i < headers.size(); ++i)
  {
    if (headers[i] == column_name)
      return static_cast<int>(i);
  }
  
  if (processor_id() == 0)
  {
    _console << "错误：找不到列 '" << column_name << "'" << std::endl;
    _console << "可用列名:";
    for (const auto& header : headers)
      _console << " '" << header << "'";
    _console << std::endl;
  }
  
  return -1;
}

Real
ADCSVBndsMaterial2::findBndsValue(Real x, Real y) const
{
  if (_csv_data.empty())
    return 0.0;
  
  Real min_distance = std::numeric_limits<Real>::max();
  Real best_value = 0.0;
  
  // 查找最近的数据点 - 参考ADCSVBndsMaterial.C的做法
  for (const auto & row : _csv_data)
  {
    if (row.size() <= static_cast<size_t>(std::max({_data_col_idx, _x_col_idx, _y_col_idx})))
      continue;
      
    Real csv_x = row[_x_col_idx];
    Real csv_y = row[_y_col_idx];
    Real data_val = row[_data_col_idx];
    
    Real distance = std::sqrt((x - csv_x) * (x - csv_x) + (y - csv_y) * (y - csv_y));
    
    if (distance < min_distance)
    {
      min_distance = distance;
      best_value = data_val;
    }
  }
  
  // 验证：确保返回值只能是0.1或0.99
  // if ((best_value - 0.5) > 0)
  //   return 1;
  // else
  //   return 0.1;
  return best_value;
}