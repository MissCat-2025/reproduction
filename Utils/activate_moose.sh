#!/bin/bash

# 激活 MOOSE 运行环境，并把参数转给 Python 脚本
source $(conda info --base)/etc/profile.d/conda.sh
conda activate moose

# 双重确认当前环境名称
if [ "$CONDA_DEFAULT_ENV" != "moose" ]; then
    echo "❌ MOOSE环境激活失败！"
    exit 1
fi

# 将当前脚本后面的参数交给 Python 执行
exec python "$@"
