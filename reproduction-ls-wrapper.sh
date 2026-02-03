#!/usr/bin/env bash

# 实际可执行文件
REAL_EXEC="/home/yp/projects/reproduction/reproduction-opt"
LOG_FILE="/home/yp/projects/reproduction/moose_ls_debug.log"

# 记录时间戳
echo "--- Call at $(date) ---" >> "$LOG_FILE"
# 记录接收到的所有参数
echo "Args received: $@" >> "$LOG_FILE"

# 直接透传参数执行
exec "$REAL_EXEC" "$@"
