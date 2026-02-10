import os
import glob
import subprocess
import time
import sys
import re
import json
from datetime import datetime

#######################
# 用户配置参数
#######################

# 路径配置
# 基础配置
# base_dir = '/home/yp/projects/raccoon/FuelFracture/RodFuel/Liwei2021/MaterialParametersVerification/step4.3_ThermalCreepFractureReturnMapQuarter'

# 修改为使用脚本所在路径：
BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # 获取脚本所在目录作为基础目录
OUTPUT_DIR = os.path.join(BASE_DIR, 'parameter_studies')           # 参数研究输出目录
MOOSE_APP = "/home/yp/projects/reproduction/reproduction-opt"               # MOOSE可执行文件路径

# 运行配置
MPI_PROCESSES = 10       # MPI进程数
TIMEOUT = 36000           # 单个案例超时时间（秒）
CONDA_ENV = 'moose'      # Conda环境名称

# 输出配置
LOG_FILE = 'run.log'     # 运行日志文件名+
PROGRESS_FILE = '.run_progress.json'  # 进度文件名

# 文件匹配模式
MAIN_PATTERN = "case_*/[!main_]*.i"    # 主程序文件匹配模式
SINGLE_PATTERN = "case_*/[!main_]*.i"  # 单程序文件匹配模式
SUB_PATTERN = "sub_*.i"             # 子程序文件匹配模式
# MAIN_PATTERN = "case_*/main_*.i"    # 主程序文件匹配模式
# SINGLE_PATTERN = "case_*/[!main_]*.i"  # 单程序文件匹配模式
# SUB_PATTERN = "sub_*.i"             # 子程序文件匹配模式
#######################
# 程序代码
#######################

def activate_and_run():
    """激活MOOSE环境并重新运行此脚本"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    activate_script = os.path.join(script_dir, 'activate_moose.sh')
    
    if not os.path.exists(activate_script):
        with open(activate_script, 'w') as f:
            f.write(f'''#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate {CONDA_ENV}
if [ "$CONDA_DEFAULT_ENV" != "{CONDA_ENV}" ]; then
    echo "❌ MOOSE环境激活失败！"
    exit 1
fi
exec python "$@"
''')
        os.chmod(activate_script, 0o755)

    script_path = os.path.abspath(__file__)
    try:
        print("正在激活MOOSE环境...")
        os.execv('/bin/bash', ['/bin/bash', activate_script, script_path])
    except Exception as e:
        print(f"环境激活失败: {str(e)}")
        sys.exit(1)

def check_environment():
    """检查当前环境"""
    issues = []
    
    # 检查是否在MOOSE环境中
    current_env = os.environ.get('CONDA_DEFAULT_ENV', '')
    if current_env != CONDA_ENV:
        activate_and_run()
        return []
    
    # 检查MOOSE可执行文件
    if not os.path.exists(MOOSE_APP):
        issues.append(f"⚠ MOOSE可执行文件不存在: {MOOSE_APP}")
    
    # 检查mpirun命令
    try:
        subprocess.run(['which', 'mpirun'], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        issues.append("⚠ 未找到mpirun命令，请确保已安装MPI")
        issues.append("  Ubuntu/Debian: sudo apt-get install mpich")
        issues.append(f"  或在{CONDA_ENV}环境中: conda install mpich")
    
    return issues

def find_input_files():
    """查找所有输入文件，智能检测单程序和多程序模式"""
    cases = []
    
    # 查找所有case目录
    case_dirs = glob.glob(os.path.join(OUTPUT_DIR, "case_*"))
    
    for case_dir in case_dirs:
        # 查找该case目录下的所有.i文件
        i_files = glob.glob(os.path.join(case_dir, "*.i"))
        
        # 过滤掉子程序文件
        main_files = [f for f in i_files if not os.path.basename(f).startswith('sub_')]
        
        # 如果找到文件，添加到列表中
        if main_files:
            cases.extend(main_files)
    
    # 按case编号排序
    def get_case_number(file_path):
        match = re.search(r'case_(\d+)', file_path)
        return int(match.group(1)) if match else float('inf')
    
    return sorted(cases, key=get_case_number)

def save_progress(completed_cases):
    """保存运行进度"""
    progress_file = os.path.join(OUTPUT_DIR, PROGRESS_FILE)
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        with open(progress_file, 'w') as f:
            json.dump(completed_cases, f)
    except Exception as e:
        print(f"警告：无法保存进度信息: {str(e)}")

def load_progress():
    """加载运行进度"""
    progress_file = os.path.join(OUTPUT_DIR, PROGRESS_FILE)
    if os.path.exists(progress_file):
        try:
            with open(progress_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"警告：无法加载进度信息: {str(e)}")
    return []

def check_convergence(log_path):
    """检查运行日志中是否有收敛问题"""
    try:
        with open(log_path, 'r') as f:
            content = f.read()
            # 检查是否有收敛失败的标志
            if "Solve Did NOT Converge!" in content or "Solve Failed!" in content:
                return False, "收敛失败"
            # 检查是否有其他严重错误
            if "*** ERROR ***" in content:
                return False, "运行错误"
            # 检查是否正常完成
            if "Finished Executing" in content:
                return True, "运行完成"
    except Exception as e:
        return None, f"无法读取日志: {str(e)}"
    return None, "状态未知"

def run_case(input_path, is_first_case=False):
    """执行单个案例"""
    case_dir = os.path.dirname(input_path)
    input_name = os.path.basename(input_path)
    log_path = os.path.join(case_dir, LOG_FILE)
    
    # 预检查
    print(f"\n🔍 预检查案例目录: {case_dir}")
    print(f"   输入文件存在: {os.path.exists(os.path.join(case_dir, input_name))}")
    
    # 检查是否为多程序模式
    is_multiapp = input_name.startswith('main_')
    if is_multiapp:
        sub_pattern = os.path.join(case_dir, SUB_PATTERN)
        has_sub = bool(glob.glob(sub_pattern))
        print(f"   模式: MultiApp (子程序{'存在' if has_sub else '不存在'})")
    else:
        print("   模式: SingleApp")
    print(f"   MOOSE可执行文件权限: {oct(os.stat(MOOSE_APP).st_mode)[-3:]}")

    # 检查上次运行状态
    if os.path.exists(log_path):
        converged, message = check_convergence(log_path)
        if converged is False:
            print(f"\n⚠ 上次运行{message}，跳过此案例")
            return {
                'status': 'skipped',
                'reason': message,
                'log': log_path
            }

    # 构建命令
    cmd = ["mpirun", "-n", str(MPI_PROCESSES), MOOSE_APP, "-i", input_name]
    
    # 如果是第一个案例，检查是否存在checkpoint文件夹
    if is_first_case:
        # 检查checkpoint文件夹
        checkpoint_pattern = os.path.join(case_dir, "*_my_checkpoint_cp")
        checkpoint_folders = glob.glob(checkpoint_pattern)
        if checkpoint_folders:
            cmd.append("--recover")
            print(f"\n💡 发现checkpoint文件夹: {os.path.basename(checkpoint_folders[0])}")
            print(f"   将从上次中断处恢复运行")
    
    print(f"\n▶ 开始执行案例: {os.path.relpath(input_path, BASE_DIR)}")
    
    try:
        with open(log_path, 'a' if is_first_case else 'w') as log_file:
            # 写入日志头
            log_file.write(f"\n=== {'恢复' if is_first_case else '开始'}执行 {datetime.now().isoformat()} ===\n")
            log_file.write(f"模式: {'MultiApp' if is_multiapp else 'SingleApp'}\n")
            log_file.write(f"命令: {' '.join(cmd)}\n")
            log_file.write(f"工作目录: {case_dir}\n\n")
            log_file.flush()

            # 执行命令
            start_time = time.time()
            process = subprocess.Popen(
                cmd,
                cwd=case_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
            )

            # 实时输出
            while True:
                output = process.stdout.readline()
                if output:
                    print(f"[{datetime.now().strftime('%H:%M:%S')}] {output.strip()}")
                    log_file.write(f"[{datetime.now().isoformat()}] {output}")
                    log_file.flush()
                if process.poll() is not None and output == '':
                    break

            # 记录结果
            elapsed = time.time() - start_time
            log_file.write(f"\n=== 运行结束 ===\n")
            log_file.write(f"返回码: {process.returncode}\n")
            log_file.write(f"耗时: {elapsed:.1f}s\n")
            
            # 检查运行结果
            converged, message = check_convergence(log_path)
            if converged is False:
                return {
                    'status': 'failed',
                    'reason': message,
                    'time': round(elapsed, 1),
                    'log': log_path,
                    'recovered': is_first_case
                }
            
            return {
                'status': 'success' if process.returncode == 0 else 'failed',
                'time': round(elapsed, 1),
                'log': log_path,
                'recovered': is_first_case
            }
            
    except Exception as e:
        error_msg = f"严重错误: {str(e)}"
        print(error_msg)
        return {
            'status': 'error',
            'error': error_msg,
            'log': log_path,
            'recovered': is_first_case
        }

def main():
    # 检查环境
    issues = check_environment()
    if issues:
        print("\n环境检查发现以下问题：")
        for issue in issues:
            print(issue)
        sys.exit(1)

    # 加载进度
    completed_cases = load_progress()
    if completed_cases:
        print(f"\n发现 {len(completed_cases)} 个已完成的案例")

    # 查找待运行案例
    cases = find_input_files()
    if not cases:
        print("未找到可执行案例！")
        return
    
    # 过滤已完成案例
    cases_to_run = [case for case in cases 
                    if os.path.relpath(case, BASE_DIR) not in completed_cases]
    
    if len(cases) != len(cases_to_run):
        for case in cases:
            if os.path.relpath(case, BASE_DIR) in completed_cases:
                print(f"跳过已完成的案例: {os.path.basename(case)}")
    
    print(f"\n找到 {len(cases_to_run)} 个待执行案例")
    
    # 执行案例
    results = []
    try:
        for idx, case in enumerate(cases_to_run):
            print(f"\n=== 进度 [{idx+1}/{len(cases_to_run)}] ===")
            result = run_case(case, is_first_case=(idx == 0))
            results.append(result)
            
            if result['status'] == 'success':
                print(f"✔ 成功完成！耗时 {result['time']} 秒")
                completed_cases.append(os.path.relpath(case, BASE_DIR))
                save_progress(completed_cases)
            elif result['status'] == 'skipped':
                print(f"⏭ 跳过案例！原因: {result['reason']}")
            else:
                print(f"✖ 执行失败！日志路径: {result['log']}")
                if 'reason' in result:
                    print(f"   原因: {result['reason']}")
    except KeyboardInterrupt:
        print("\n\n检测到用户中断，保存进度...")
        save_progress(completed_cases)
        print("进度已保存，下次运行时将从中断处继续")
        sys.exit(1)
    
    # 生成报告
    success_count = sum(1 for r in results if r['status'] == 'success')
    recovered_count = sum(1 for r in results if r.get('recovered', False))
    print(f"\n执行完成：成功 {success_count}/{len(cases_to_run)} 个案例")
    if recovered_count > 0:
        print(f"其中 {recovered_count} 个案例是从中断处恢复运行的")
    print(f"详细日志请查看各案例目录下的 {LOG_FILE} 文件")

    # 清理进度文件
    progress_file = os.path.join(OUTPUT_DIR, PROGRESS_FILE)
    if os.path.exists(progress_file):
        os.remove(progress_file)

if __name__ == '__main__':
    main()