#!/bin/bash

# MakeReproduction.sh - 一键解决reproduction编译问题
# 专门处理权限问题、环境激活等

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# 检查当前目录
check_directory() {
    if [ ! -f "Makefile" ]; then
        print_error "请在reproduction项目根目录下运行此脚本！"
        exit 1
    fi
    print_success "当前目录: $(pwd)"
}

# 修复换行符问题
fix_line_endings() {
    print_info "修复换行符(.i/.py/.sh)为UNIX格式..."

    # 目标根目录数组
    local roots=(".. /moose" ".. /raccoon" ".")

    # 逐个根目录处理 .i/.py/.sh 文件的 CRLF
    if [ -d "../moose" ]; then
        find ../moose -type f \( -name "*.i" -o -name "*.py" -o -name "*.sh" \) -exec sed -i 's/\r$//' {} \; 2>/dev/null || true
    fi
    if [ -d "../raccoon" ]; then
        find ../raccoon -type f \( -name "*.i" -o -name "*.py" -o -name "*.sh" \) -exec sed -i 's/\r$//' {} \; 2>/dev/null || true
    fi
    find . -type f \( -name "*.i" -o -name "*.py" -o -name "*.sh" \) -exec sed -i 's/\r$//' {} \; 2>/dev/null || true

    print_success "换行符修复完成"
}

# 修复MOOSE版本头文件
fix_moose_revision() {
    print_info "检查并修复MOOSE版本头文件..."
    
    local revision_file="../moose/framework/include/base/MooseRevision.h"
    
    if [ ! -f "$revision_file" ]; then
        print_warning "MooseRevision.h文件不存在，正在创建..."
        
        # 创建目录（如果不存在）
        mkdir -p "$(dirname "$revision_file")"
        
        # 创建简单的MooseRevision.h文件
        cat > "$revision_file" << 'EOF'
#ifndef MOOSEREVISION_H
#define MOOSEREVISION_H

#define MOOSE_REVISION "unknown"
#define MOOSE_REVISION_SHORT "unknown"
#define MOOSE_REVISION_LONG "unknown"
#define MOOSE_VERSION "unknown"

#endif // MOOSEREVISION_H
EOF
        
        print_success "MooseRevision.h文件已创建"
    else
        print_success "MooseRevision.h文件已存在"
    fi
}

# 检查并编译RACCOON
check_and_build_raccoon() {
    print_info "检查RACCOON库..."
    
    local raccoon_lib="../raccoon/lib/libraccoon-opt.la"
    
    if [ ! -f "$raccoon_lib" ]; then
        print_warning "RACCOON库不存在，正在编译RACCOON..."
        
        # 进入RACCOON目录并编译
        cd ../raccoon
        if make -j8; then
            print_success "RACCOON编译成功"
        else
            print_error "RACCOON编译失败"
            return 1
        fi
        cd ../reproduction
    else
        print_success "RACCOON库已存在"
    fi
}

# 重新编译MOOSE框架和模块
rebuild_moose_modules() {
    print_info "重新编译MOOSE框架和模块..."
    
    # 清理MOOSE构建目录
    if [ -d "../moose/framework/build" ]; then
        print_info "清理MOOSE框架构建目录..."
        rm -rf ../moose/framework/build/*
    fi
    
    # 清理MOOSE模块构建目录
    for module_dir in ../moose/modules/*/build; do
        if [ -d "$module_dir" ]; then
            print_info "清理模块构建目录: $(dirname "$module_dir")"
            rm -rf "$module_dir"/*
        fi
    done
    
    # 重新编译MOOSE框架
    print_info "重新编译MOOSE框架..."
    cd ../moose/framework
    if make -j8; then
        print_success "MOOSE框架编译成功"
    else
        print_error "MOOSE框架编译失败"
        cd ../../reproduction
        return 1
    fi
    
    # 重新编译MOOSE模块
    print_info "重新编译MOOSE模块..."
    cd ../modules
    if make -j8; then
        print_success "MOOSE模块编译成功"
    else
        print_error "MOOSE模块编译失败"
        cd ../../reproduction
        return 1
    fi
    
    cd ../../reproduction
    print_success "MOOSE重新编译完成"
}

# 修复权限问题（递归激活脚本）
fix_permissions() {
    print_info "递归激活脚本执行权限(.py/.sh/.pl)..."

    # moose
    if [ -d "../moose" ]; then
        find ../moose -type f \( -name "*.py" -o -name "*.sh" -o -name "*.pl" \) -exec chmod +x {} \; 2>/dev/null || true
    fi
    # raccoon
    if [ -d "../raccoon" ]; then
        find ../raccoon -type f \( -name "*.py" -o -name "*.sh" -o -name "*.pl" \) -exec chmod +x {} \; 2>/dev/null || true
    fi
    # reproduction
    find . -type f \( -name "*.py" -o -name "*.sh" -o -name "*.pl" \) -exec chmod +x {} \; 2>/dev/null || true

    print_success "脚本权限激活完成"
}

# 清理构建目录
clean_build() {
    print_info "清理构建目录..."
    if [ -d "build" ]; then
        rm -rf build/*
        print_success "构建目录已清理"
    else
        mkdir -p build
        print_success "创建构建目录"
    fi
}

# 运行编译
run_compilation() {
    local num_jobs=${1:-8}
    print_info "使用 $num_jobs 个并行任务进行编译..."
    
    if make -j$num_jobs; then
        print_success "编译成功！"
        return 0
    else
        print_error "编译失败！"
        return 1
    fi
}

# 显示帮助
show_help() {
    echo "MakeReproduction.sh - 一键解决reproduction编译问题"
    echo ""
    echo "使用方法:"
    echo "  ./MakeReproduction.sh [选项]"
    echo ""
    echo "选项:"
    echo "  -j N        使用N个并行任务编译 (默认: 8)"
    echo "  -c          仅清理构建目录"
    echo "  -p          仅修复权限问题"
    echo "  -x          仅修复换行符与权限(不编译)"
    echo "  -r          重新编译MOOSE框架和模块"
    echo "  -h          显示帮助信息"
    echo ""
    echo "注意: 请确保已激活moose conda环境"
    echo "运行: conda activate moose"
}

# 主函数
main() {
    local num_jobs=8
    local clean_only=false
    local permission_only=false
    local rebuild_moose=false
    local fix_only=false
    
    # 解析参数
    while [[ $# -gt 0 ]]; do
        case $1 in
            -j) num_jobs="$2"; shift 2 ;;
            -c) clean_only=true; shift ;;
            -p) permission_only=true; shift ;;
            -x) fix_only=true; shift ;;
            -r) rebuild_moose=true; shift ;;
            -h|--help) show_help; exit 0 ;;
            *) print_error "未知选项: $1"; show_help; exit 1 ;;
        esac
    done
    
    print_info "开始执行MakeReproduction.sh..."
    print_info "=========================================="
    
    check_directory
    
    if [ "$clean_only" = true ]; then
        clean_build
        print_success "清理完成！"
        exit 0
    fi
    
    if [ "$permission_only" = true ]; then
        fix_permissions
        print_success "权限修复完成！"
        exit 0
    fi

    if [ "$fix_only" = true ]; then
        fix_line_endings
        fix_permissions
        print_success "换行与权限修复完成！"
        exit 0
    fi
    
    # 完整流程
    fix_line_endings
    fix_moose_revision
    
    # 如果需要重新编译MOOSE
    if [ "$rebuild_moose" = true ]; then
        rebuild_moose_modules
    fi
    
    check_and_build_raccoon
    fix_permissions
    clean_build
    
    print_info "=========================================="
    print_info "开始编译..."
    
    if run_compilation $num_jobs; then
        print_success "=========================================="
        print_success "编译成功！reproduction项目已准备就绪"
        print_success "可执行文件: reproduction-opt"
        print_success "=========================================="
    else
        print_error "=========================================="
        print_error "编译失败！"
        print_error "解决方案："
        print_error "1. 确保已激活moose环境: conda activate moose"
        print_error "2. 运行权限修复: ./MakeReproduction.sh -p"
        print_error "3. 清理构建目录: ./MakeReproduction.sh -c"
        print_error "=========================================="
        exit 1
    fi
}

main "$@"