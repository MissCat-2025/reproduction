#!/bin/bash

# RACCOON 安装脚本 - 简化版
# 功能：检测 -> 下载 -> 编译

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# 配置
RACCOON_DIR="/home/yp/projects/raccoon"
GITHUB_URL="https://github.com/hugary1995/raccoon.git"
MIRROR_URL="https://gitee.com/mirrors/raccoon.git"
TIMEOUT=180  # 3分钟超时

# 打印消息
print_msg() {
    echo -e "${1}${2}${NC}"
}

# 检查 RACCOON 是否存在
check_raccoon() {
    print_msg $BLUE "检查 RACCOON 安装状态..."
    
    if [ -d "$RACCOON_DIR" ] && [ -f "$RACCOON_DIR/Makefile" ]; then
        print_msg $GREEN "RACCOON 已存在"
        return 0
    else
        print_msg $YELLOW "RACCOON 不存在，需要下载"
        return 1
    fi
}

# 带超时的 Git 克隆
clone_with_timeout() {
    local url=$1
    local timeout=$2
    
    print_msg $BLUE "从 $url 克隆 RACCOON..."
    
    # 删除已存在的目录
    [ -d "$RACCOON_DIR" ] && rm -rf "$RACCOON_DIR"
    
    # 使用 timeout 命令限制克隆时间
    if timeout $timeout git clone --recursive "$url" "$RACCOON_DIR"; then
        print_msg $GREEN "克隆成功"
        return 0
    else
        print_msg $RED "克隆失败或超时"
        return 1
    fi
}

# 下载 RACCOON
download_raccoon() {
    print_msg $BLUE "开始下载 RACCOON..."
    
    # 先尝试官方源
    if clone_with_timeout "$GITHUB_URL" $TIMEOUT; then
        return 0
    fi
    
    # 官方源失败，尝试镜像源
    print_msg $YELLOW "官方源超时，尝试镜像源..."
    if clone_with_timeout "$MIRROR_URL" $TIMEOUT; then
        return 0
    fi
    
    print_msg $RED "所有源都失败"
    return 1
}

# 激活 conda 环境
activate_moose_env() {
    print_msg $BLUE "激活 moose 环境..."
    
    # 初始化 conda
    eval "$(conda shell.bash hook)"
    
    # 激活 moose 环境
    if conda activate moose; then
        print_msg $GREEN "moose 环境激活成功"
        return 0
    else
        print_msg $RED "moose 环境激活失败"
        return 1
    fi
}

# 编译 RACCOON
compile_raccoon() {
    print_msg $BLUE "开始编译 RACCOON..."
    
    cd "$RACCOON_DIR"
    
    if make -j 8; then
        print_msg $GREEN "编译成功！"
        
        # 显示编译结果
        if [ -f "raccoon-opt" ]; then
            print_msg $GREEN "生成可执行文件: raccoon-opt"
        fi
        if [ -f "raccoon-dbg" ]; then
            print_msg $GREEN "生成可执行文件: raccoon-dbg"
        fi
        
        return 0
    else
        print_msg $RED "编译失败"
        return 1
    fi
}

# 主函数
main() {
    print_msg $BLUE "========================================"
    print_msg $BLUE "RACCOON 安装脚本 - 简化版"
    print_msg $BLUE "========================================"
    
    # 1. 检查 RACCOON
    if ! check_raccoon; then
        # 2. 下载 RACCOON
        if ! download_raccoon; then
            print_msg $RED "下载失败，退出"
            exit 1
        fi
    fi
    
    # 3. 激活环境
    if ! activate_moose_env; then
        print_msg $RED "环境激活失败，退出"
        exit 1
    fi
    
    # 4. 编译
    if ! compile_raccoon; then
        print_msg $RED "编译失败，退出"
        exit 1
    fi
    
    print_msg $GREEN "========================================"
    print_msg $GREEN "RACCOON 安装完成！"
    print_msg $GREEN "位置: $RACCOON_DIR"
    print_msg $GREEN "========================================"
}

# 运行主函数
main "$@"