#!/bin/bash

echo "=== 数据导出工具 ==="
echo "找到默认.e文件: $(pwd)/gap_conductance1/"
echo "找到默认.e文件: $(pwd)/gap_conductance1/2D.e"

# 激活ParaView环境（如果存在）
if command -v conda &> /dev/null; then
    # 尝试激活paraview_post环境
    if conda env list | grep -q "paraview_post"; then
        echo "激活ParaView环境..."
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate paraview_post
    fi
fi

echo ""

# 直接运行Python脚本
python export_point_data.py

echo ""
echo "✅ 脚本执行完成！"
