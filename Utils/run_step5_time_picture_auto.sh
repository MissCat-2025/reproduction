#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ENV_NAME="${TIME_PICTURE_CONDA_ENV:-paraview_post}"

echo "=== Step5 自动环境检查 ==="

if ! command -v conda &> /dev/null; then
    echo "错误：未找到 conda，无法自动安装 Step5 所需环境"
    exit 1
fi

ENV_EXISTS=$(conda env list | awk -v env="$ENV_NAME" '$1 == env { found=1 } END { print found ? 1 : 0 }')

if [ "$ENV_EXISTS" -eq 0 ]; then
    echo "创建环境: $ENV_NAME"
    conda create -y -n "$ENV_NAME" -c conda-forge python Pillow
else
    echo "更新环境: $ENV_NAME"
    conda install -y -n "$ENV_NAME" -c conda-forge Pillow
fi

echo "使用环境运行 Step5..."
cd "$SCRIPT_DIR"
conda run -n "$ENV_NAME" python "$SCRIPT_DIR/step5_time_picture.py" "$@"
