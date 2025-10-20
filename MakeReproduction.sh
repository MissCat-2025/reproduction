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
    print_info "修复Python脚本换行符问题..."
    
    # 修复MOOSE框架脚本换行符
    if [ -d "../moose/framework/scripts" ]; then
        find ../moose/framework/scripts -name "*.py" -exec sed -i 's/\r$//' {} \; 2>/dev/null || true
        print_success "MOOSE框架脚本换行符已修复"
    fi
    
    # 修复MOOSE模块脚本换行符
    for module_dir in ../moose/modules/*/scripts; do
        if [ -d "$module_dir" ]; then
            find "$module_dir" -name "*.py" -exec sed -i 's/\r$//' {} \; 2>/dev/null || true
        fi
    done
    print_success "MOOSE模块脚本换行符已修复"
    
    # 修复RACCOON脚本换行符
    if [ -d "../raccoon/scripts" ]; then
        find ../raccoon/scripts -name "*.py" -exec sed -i 's/\r$//' {} \; 2>/dev/null || true
        print_success "RACCOON脚本换行符已修复"
    fi
    
    # 修复当前项目脚本换行符
    if [ -d "scripts" ]; then
        find scripts -name "*.py" -exec sed -i 's/\r$//' {} \; 2>/dev/null || true
        print_success "当前项目脚本换行符已修复"
    fi
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

# 修复权限问题
fix_permissions() {
    print_info "修复MOOSE框架脚本权限..."
    
    # 修复MOOSE框架脚本权限
    if [ -d "../moose/framework/scripts" ]; then
        chmod +x ../moose/framework/scripts/*.py 2>/dev/null || true
        print_success "MOOSE框架脚本权限已修复"
    fi
    
    # 修复MOOSE模块脚本权限
    for module_dir in ../moose/modules/*/scripts; do
        if [ -d "$module_dir" ]; then
            chmod +x "$module_dir"/*.py 2>/dev/null || true
        fi
    done
    print_success "MOOSE模块脚本权限已修复"
    
    # 修复RACCOON脚本权限
    if [ -d "../raccoon/scripts" ]; then
        chmod +x ../raccoon/scripts/*.py 2>/dev/null || true
        print_success "RACCOON脚本权限已修复"
    fi
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
    
    # 解析参数
    while [[ $# -gt 0 ]]; do
        case $1 in
            -j) num_jobs="$2"; shift 2 ;;
            -c) clean_only=true; shift ;;
            -p) permission_only=true; shift ;;
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
    
    # 完整流程
    fix_line_endings
    fix_moose_revision
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