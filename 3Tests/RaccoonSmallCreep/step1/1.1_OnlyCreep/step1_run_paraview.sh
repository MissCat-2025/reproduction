#!/bin/bash

# 获取脚本所在目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "=== ParaView专用环境配置 ==="

# 1. 检查conda安装
if ! command -v conda &> /dev/null; then
    echo "安装Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"
    rm miniconda.sh
    export PATH="$HOME/miniconda/bin:$PATH"
    source "$HOME/miniconda/etc/profile.d/conda.sh"
fi

# 2. 初始化conda
source "$HOME/miniconda/etc/profile.d/conda.sh" 2> /dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"

# 3. 环境配置
ENV_NAME="paraview_post"
ENV_EXISTS=$(conda env list | grep -w "$ENV_NAME" | wc -l)

# 添加强制重建参数处理
FORCE_RECREATE=false
if [[ $1 == "--force" ]]; then
    FORCE_RECREATE=true
    echo "强制重建环境..."
    conda deactivate 2> /dev/null
    conda env remove -n "$ENV_NAME" -y
    ENV_EXISTS=0
fi

# 检查包是否正确安装的函数
check_packages() {
    echo "检查必要包是否已安装..."
    
    # 检查numpy是否存在
    if ! conda run -n "$ENV_NAME" python -c "import numpy; print('numpy版本:', numpy.__version__)" 2>/dev/null; then
        echo "❌ numpy未正确安装"
        return 1
    fi
    
    # 检查pandas是否存在
    if ! conda run -n "$ENV_NAME" python -c "import pandas; print('pandas版本:', pandas.__version__)" 2>/dev/null; then
        echo "❌ pandas未正确安装"
        return 1
    fi
    
    # 检查paraview是否存在
    if ! conda run -n "$ENV_NAME" python -c "import paraview; print('paraview可用')" 2>/dev/null; then
        echo "❌ paraview未正确安装"
        return 1
    fi
    
    # 检查openpyxl是否存在
    if ! conda run -n "$ENV_NAME" python -c "import openpyxl; print('openpyxl版本:', openpyxl.__version__)" 2>/dev/null; then
        echo "❌ openpyxl未正确安装"
        return 1
    fi
    
    echo "✅ 所有必要包都已正确安装"
    return 0
}

# 安装包的函数
install_packages() {
    echo "开始安装必要包..."
    
    # 首先尝试使用conda安装基础包
    echo "安装基础包..."
    conda install -n "$ENV_NAME" -y -c conda-forge \
        numpy \
        pandas \
        matplotlib \
        scipy \
        openpyxl \
        xlsxwriter \
        h5py
    
    if [ $? -ne 0 ]; then
        echo "❌ 基础包安装失败"
        return 1
    fi
    
    # 安装ParaView和VTK
    echo "安装ParaView..."
    conda install -n "$ENV_NAME" -y -c conda-forge paraview
    
    if [ $? -ne 0 ]; then
        echo "❌ ParaView安装失败"
        return 1
    fi
    
    echo "✅ 所有包安装完成"
    return 0
}

# 创建或检查环境
if [ "$ENV_EXISTS" -eq 0 ]; then
    echo "创建新环境: $ENV_NAME..."
    # 使用Python 3.11以获得更好的ParaView兼容性
    conda create -y -n "$ENV_NAME" python=3.11
    
    if [ $? -ne 0 ]; then
        echo "❌ 环境创建失败"
        exit 1
    fi
    
    # 安装包
    install_packages
    if [ $? -ne 0 ]; then
        echo "❌ 包安装失败"
        exit 1
    fi
    
    echo "✅ 环境创建完成"
else
    echo "使用现有环境: $ENV_NAME"
    
    # 检查Python版本
    PYTHON_VERSION=$(conda run -n "$ENV_NAME" python --version 2>&1 | grep -o "3\.[0-9]*")
    if [[ "$PYTHON_VERSION" != "3.11" ]]; then
        echo "检测到Python版本为$PYTHON_VERSION，需要Python 3.11才能兼容ParaView"
        echo "强制重建环境..."
        conda deactivate 2> /dev/null
        conda env remove -n "$ENV_NAME" -y
        conda create -y -n "$ENV_NAME" python=3.11
        
        if [ $? -ne 0 ]; then
            echo "❌ 环境重建失败"
            exit 1
        fi
        
        # 安装包
        install_packages
        if [ $? -ne 0 ]; then
            echo "❌ 包安装失败"
            exit 1
        fi
    else
        # 检查包是否正确安装
        if ! check_packages; then
            echo "包不完整，重新安装..."
            install_packages
            if [ $? -ne 0 ]; then
                echo "❌ 包重新安装失败"
                exit 1
            fi
        fi
    fi
fi

# 4. 设置环境变量
echo "配置环境变量..."
ENV_PREFIX=$(conda info --envs | grep "$ENV_NAME" | awk '{print $2}')
export PV_PYTHON_PATH="$ENV_PREFIX/lib/python3.11/site-packages"
export PV_LIB_PATH="$ENV_PREFIX/lib"
export PYTHONPATH="$PV_PYTHON_PATH:$PYTHONPATH"
export LD_LIBRARY_PATH="$PV_LIB_PATH:$LD_LIBRARY_PATH"

# 5. 最终检查
echo "最终检查环境状态..."
echo "环境路径: $ENV_PREFIX"
echo "Python路径: $(conda run -n "$ENV_NAME" which python)"

# 检查关键包
echo "检查关键包..."
conda run -n "$ENV_NAME" python -c "
import sys
print('Python版本:', sys.version)
try:
    import numpy
    print('✅ numpy:', numpy.__version__)
except ImportError as e:
    print('❌ numpy未安装:', e)
    
try:
    import pandas
    print('✅ pandas:', pandas.__version__)
except ImportError as e:
    print('❌ pandas未安装:', e)
    
try:
    import paraview
    print('✅ paraview可用')
except ImportError as e:
    print('❌ paraview未安装:', e)
    
try:
    import openpyxl
    print('✅ openpyxl:', openpyxl.__version__)
except ImportError as e:
    print('❌ openpyxl未安装:', e)
"

# 6. 运行处理脚本
echo "执行处理脚本..."
cd "$SCRIPT_DIR"

# 确保在正确的conda环境中运行Python脚本
conda run -n "$ENV_NAME" python "$SCRIPT_DIR/paraview_processor.py"

RESULT=$?
if [ $RESULT -eq 0 ]; then
    echo "✅ 处理完成"
else
    echo "❌ 处理失败，错误代码: $RESULT"
fi

exit $RESULT