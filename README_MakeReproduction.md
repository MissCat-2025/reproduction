# MakeReproduction.sh 使用说明

## 概述
`MakeReproduction.sh` 是一个专门用于解决 reproduction 项目与 MOOSE、RACCOON 之间连接关系的一键脚本。它主要解决编译过程中的权限问题、环境配置问题等。

## 功能特性
- ✅ 自动修复 MOOSE 框架脚本权限问题
- ✅ 自动修复 MOOSE 模块脚本权限问题  
- ✅ 自动修复 RACCOON 脚本权限问题
- ✅ 自动清理构建目录
- ✅ 支持并行编译
- ✅ 彩色输出，便于查看执行状态
- ✅ 详细的错误提示和解决方案

## 使用方法

### 基本用法
```bash
# 激活moose环境
conda activate moose

# 进入reproduction项目目录
cd /home/yp/projects/reproduction

# 运行脚本（使用默认8个并行任务）
./MakeReproduction.sh
```

### 命令行选项
```bash
# 使用12个并行任务编译
./MakeReproduction.sh -j 12

# 仅清理构建目录
./MakeReproduction.sh -c

# 仅修复权限问题
./MakeReproduction.sh -p

# 显示帮助信息
./MakeReproduction.sh -h
```

## 解决的问题

### 1. 权限问题
**问题描述**：
```
make: /home/yp/projects/moose/framework/scripts/write_appresource_file.py: Permission denied
/bin/sh: 1: /home/yp/projects/moose/framework/scripts/get_repo_revision.py: Permission denied
```

**解决方案**：
脚本会自动执行：
```bash
chmod +x ../moose/framework/scripts/*.py
chmod +x ../moose/modules/*/scripts/*.py
chmod +x ../raccoon/scripts/*.py
```

### 2. 构建目录问题
**问题描述**：
构建目录可能包含过期的文件，导致编译失败。

**解决方案**：
脚本会自动清理构建目录：
```bash
rm -rf build/*
```

### 3. 环境配置问题
**问题描述**：
需要确保在正确的 conda 环境中运行编译。

**解决方案**：
脚本会提示用户激活 moose 环境：
```bash
conda activate moose
```

## 项目结构要求
脚本假设以下项目结构：
```
/home/yp/projects/
├── moose/          # MOOSE框架
├── raccoon/        # RACCOON项目
└── reproduction/   # 当前项目
    ├── Makefile
    ├── MakeReproduction.sh
    └── build/
```

## 常见问题

### Q: 脚本提示"请在reproduction项目根目录下运行"
**A**: 确保在包含 Makefile 的目录下运行脚本。

### Q: 编译失败，提示找不到编译器
**A**: 确保已激活 moose conda 环境：
```bash
conda activate moose
```

### Q: 权限修复后仍然有权限问题
**A**: 手动检查文件权限：
```bash
ls -la ../moose/framework/scripts/*.py
```

### Q: 想要重新开始编译
**A**: 清理构建目录：
```bash
./MakeReproduction.sh -c
```

## 输出示例

### 成功编译
```
[INFO] 开始执行MakeReproduction.sh...
[INFO] ==========================================
[SUCCESS] 当前目录: /home/yp/projects/reproduction
[INFO] 修复MOOSE框架脚本权限...
[SUCCESS] MOOSE框架脚本权限已修复
[SUCCESS] MOOSE模块脚本权限已修复
[SUCCESS] RACCOON脚本权限已修复
[INFO] 清理构建目录...
[SUCCESS] 构建目录已清理
[INFO] ==========================================
[INFO] 开始编译...
[INFO] 使用 8 个并行任务进行编译...
[SUCCESS] 编译成功！
[SUCCESS] ==========================================
[SUCCESS] 编译成功！reproduction项目已准备就绪
[SUCCESS] 可执行文件: build/reproduction-opt
[SUCCESS] ==========================================
```

### 编译失败
```
[ERROR] 编译失败！
[ERROR] ==========================================
[ERROR] 编译失败！
[ERROR] 解决方案：
[ERROR] 1. 确保已激活moose环境: conda activate moose
[ERROR] 2. 运行权限修复: ./MakeReproduction.sh -p
[ERROR] 3. 清理构建目录: ./MakeReproduction.sh -c
[ERROR] ==========================================
```

## 技术细节

### 权限修复范围
- `../moose/framework/scripts/*.py`
- `../moose/modules/*/scripts/*.py`
- `../raccoon/scripts/*.py`

### 构建目录清理
- 删除 `build/` 目录下的所有文件
- 如果 `build/` 目录不存在，则创建它

### 并行编译
- 默认使用 8 个并行任务
- 可通过 `-j N` 参数调整并行任务数

## 维护说明
如果遇到新的编译问题，可以：
1. 在脚本中添加新的权限修复逻辑
2. 添加新的环境检查
3. 更新错误提示信息

脚本设计为模块化，易于扩展和维护。
