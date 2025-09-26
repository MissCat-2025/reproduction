#!/bin/bash

echo "=== 数据导出工具 (v2) ==="
echo "默认数据目录: $(pwd)/gap_conductance1/"
echo "默认.e文件:   $(pwd)/gap_conductance1/2D.e"

# 优先使用 pvpython（ParaView 自带解释器）
if command -v pvpython &> /dev/null; then
  echo "使用 pvpython 运行 export_point_data_v2.py ..."
  pvpython export_point_data_v2.py | cat
  echo "\n✅ v2 导出完成！"
  exit 0
fi

# 否则尝试激活 conda 中的 paraview_post 环境
if command -v conda &> /dev/null; then
  if conda env list | grep -q "paraview_post"; then
    echo "激活 ParaView 环境 paraview_post ..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate paraview_post
  fi
fi

# 回退使用系统 python（要求已能 import paraview.simple）
echo "使用 python 运行 export_point_data_v2.py ..."
python export_point_data_v2.py | cat

echo "\n✅ v2 导出完成！"


