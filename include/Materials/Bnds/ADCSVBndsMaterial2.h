#pragma once

#include "ADMaterial.h"

class ADCSVBndsMaterial2 : public ADMaterial
{
public:
  static InputParameters validParams();
  ADCSVBndsMaterial2(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

private:
  void readCSVFile();
  int findColumnIndex(const std::vector<std::string>& headers, 
                     const std::string& column_name) const;
  Real findBndsValue(Real x, Real y) const;

  const FileName _csv_file;
  const std::string _data_column;
  const std::string _x_column;
  const std::string _y_column;
  const std::string _property_name;
  const bool _verbose;

  std::vector<std::vector<Real>> _csv_data;
  int _data_col_idx;
  int _x_col_idx;
  int _y_col_idx;
  
  ADMaterialProperty<Real> & _bnds;
};