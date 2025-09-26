import os
import re
import shutil
import itertools
from datetime import datetime

# 基础配置（与原脚本一致）
# base_dir = '/home/yp/projects/raccoon/FuelFracture/RodFuel/Liwei2021/MaterialParametersVerification/step4.3_ThermalCreepFractureReturnMapQuarter'

# 修改为使用脚本所在路径：
base_dir = os.path.dirname(os.path.abspath(__file__))  # 获取脚本所在目录作为基础目录
# 主程序模板文件与子程序模板文件相同，则是单程序模式
template_main = os.path.join(base_dir, '2D.i')
template_sub = os.path.join(base_dir, '2D.i')
output_dir = os.path.join(base_dir, 'parameter_studies_series')

# Checkpoint配置，加入存档功能
checkpoint_config = '''
  [my_checkpoint]
    type = Checkpoint
    time_step_interval = 5    # 每5个时间步保存
    num_files = 2            # 保留最近4个检查点
    wall_time_interval = 600 # 每10分钟保存一次（秒）
  []'''

# 参数序列定义（“对齐配对”生成）：
# 示例：
# parameter_matrix = {
#   'initial_T': [583.15, 600, 620, 640],
#   'linePower': [1, 2, 3, 4]
# }
# 仅生成 4 个案例，而非笛卡尔积 16 个
parameter_matrix = {
    'initial_T': [552.6, 558.2, 570.7, 582.7, 589.5],
    'initial_T_out': [553.7, 560.8, 582.8, 602.5, 611.8],
    'linePower': [10, 60, 90, 70, 10]
}

def add_checkpoint_to_outputs(content):
    """在[Outputs]块中添加checkpoint配置"""
    outputs_match = re.search(r'\[Outputs\](.*?)\[\]', content, re.DOTALL)
    if outputs_match:
        outputs_block = outputs_match.group(1)
        if 'type = Checkpoint' not in outputs_block:
            new_outputs = f'[Outputs]{checkpoint_config}{outputs_block}[]'
            content = content.replace(outputs_match.group(0), new_outputs)
    else:
        content += f'\n[Outputs]{checkpoint_config}\n[]'
    return content

def extract_end_time(template_file):
    """从模板文件中提取end_time值（兼容大小写 EndTime / end_time）"""
    try:
        with open(template_file, 'r') as f:
            content = f.read()
            match = re.search(r'(?:end_time|EndTime)\s*=\s*([\d\.eE+-]+)', content)
            if match:
                return float(match.group(1))
    except Exception as e:
        print(f"警告：无法从模板文件中提取end_time: {str(e)}")
    return None

def format_scientific(value):
    """将数值格式化为科学计数法字符串"""
    if isinstance(value, float) and (abs(value) >= 1e4 or abs(value) < 1e-3):
        return f"{value:.2e}".replace('e-0', 'e-').replace('e+0', 'e+')
    return str(value)

def generate_case_name(params):
    """生成包含所有参数的短名称（按键名排序以保证稳定）"""
    items = sorted(params.items(), key=lambda kv: kv[0])
    return '_'.join([f"{k[:2]}{format_scientific(v).replace('.','_')}" for k, v in items])

def replace_parameters(content, params):
    """动态替换所有参数（数值）"""
    for param, value in params.items():
        pattern = rf'(\s*){re.escape(param)}\s*=\s*[\d\.eE+-]+(.*?)(\n)'
        replacement = f'\\1{param} = {format_scientific(value)}\\2\\3'
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

    # 特殊处理MultiApp的input_files参数
    subapp_filename = f"sub_{generate_case_name(params)}.i"
    content = re.sub(
        r'(input_files\s*=\s*)\'\'\S+\.i\'\'',
        f"\\1'{subapp_filename}'",
        content
    )
    return content

def generate_header(params):
    """生成包含参数信息的注释头"""
    header = "# === 参数研究案例（对齐配对） ===\n"
    end_time = extract_end_time(template_main)
    if end_time is not None:
        header += f"# end_time = {format_scientific(end_time)}\n"
    for k, v in sorted(params.items(), key=lambda kv: kv[0]):
        header += f"# {k}: {format_scientific(v)}\n"
    header += f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    return header

def rename_and_cleanup_files(case_dir, case_name, is_multiapp):
    """重命名文件并清理不需要的文件"""
    main_file = os.path.join(case_dir, f"main_{case_name}.i")
    if not is_multiapp:
        new_file = os.path.join(case_dir, f"{case_name}.i")
        if os.path.exists(main_file):
            os.rename(main_file, new_file)
        sub_file = os.path.join(case_dir, f"sub_{case_name}.i")
        if os.path.exists(sub_file):
            os.remove(sub_file)

def aligned_length_and_broadcast(parameter_matrix):
    """获取对齐长度，并支持将长度为1的列表广播为共同长度。
    规则：
    - 所有值必须为列表或元组；
    - 取非单例列表的最大长度L；
    - 长度为1的参数广播到L；
    - 若存在两个及以上不同且>1的长度且不一致，则报错。
    - 若全部为单值，则L=1。
    返回：L 以及广播后的副本字典。
    """
    lengths = []
    for key, val in parameter_matrix.items():
        if not isinstance(val, (list, tuple)):
            raise ValueError(f"参数 '{key}' 的值必须是列表或元组")
        lengths.append(len(val))

    multi_lengths = sorted({l for l in lengths if l > 1})
    if len(multi_lengths) > 1:
        raise ValueError(f"存在不一致的参数序列长度：{multi_lengths}")

    L = multi_lengths[0] if multi_lengths else 1
    broadcasted = {}
    for key, val in parameter_matrix.items():
        if len(val) == 1:
            broadcasted[key] = list(itertools.repeat(val[0], L))
        elif len(val) == L:
            broadcasted[key] = list(val)
        else:
            raise ValueError(f"参数 '{key}' 的长度为 {len(val)}，与对齐长度 {L} 不一致")
    return L, broadcasted

def generate_parameter_combinations_aligned(params_dict):
    """按序列“对齐配对”生成参数组合，长度为共同长度（支持单值广播）。"""
    L, broadcasted = aligned_length_and_broadcast(params_dict)
    keys = list(broadcasted.keys())
    combos = []
    for i in range(L):
        combo = {k: broadcasted[k][i] for k in keys}
        combos.append(combo)
    return combos

def generate_study_cases():
    # 校验主程序模板文件
    if not os.path.exists(template_main):
        raise FileNotFoundError("主程序模板文件不存在，请检查路径配置")

    # 判断是否为多程序模式
    is_multiapp = os.path.exists(template_sub) and os.path.abspath(template_main) != os.path.abspath(template_sub)

    # 清理并创建输出目录
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 生成“对齐配对”的参数组合
    all_params = generate_parameter_combinations_aligned(parameter_matrix)

    for idx, params in enumerate(all_params, 1):
        case_name = generate_case_name(params)
        case_dir = os.path.join(output_dir, f"case_{idx:03d}_{case_name}")
        os.makedirs(case_dir, exist_ok=True)

        # 处理主程序文件
        with open(template_main, 'r') as f:
            main_content = generate_header(params) + f.read()

        # 添加checkpoint配置
        main_content = add_checkpoint_to_outputs(main_content)

        # 参数替换
        if is_multiapp:
            main_content = replace_parameters(main_content, params)
        else:
            for param, value in params.items():
                pattern = rf'(\s*){re.escape(param)}\s*=\s*[\d\.eE+-]+(.*?)(\n)'
                replacement = f'\\1{param} = {format_scientific(value)}\\2\\3'
                main_content = re.sub(pattern, replacement, main_content, flags=re.MULTILINE)

        main_output = os.path.join(case_dir, f"main_{case_name}.i")
        with open(main_output, 'w') as f:
            f.write(main_content)

        # 子程序（若为多程序模式）
        if is_multiapp:
            with open(template_sub, 'r') as f:
                sub_content = generate_header(params) + f.read()
            sub_content = replace_parameters(sub_content, params)
            sub_output = os.path.join(case_dir, f"sub_{case_name}.i")
            with open(sub_output, 'w') as f:
                f.write(sub_content)

        # 重命名与清理
        rename_and_cleanup_files(case_dir, case_name, is_multiapp)

        print(f"生成案例 {idx:03d}: {case_name}")
        print(f"  路径: {case_dir}")
        print("  模式: MultiApp" if is_multiapp else "  模式: SingleApp")

if __name__ == '__main__':
    try:
        generate_study_cases()
        print(f"\n所有案例已成功生成至: {os.path.abspath(output_dir)}")
        print(f"总案例数: {len(generate_parameter_combinations_aligned(parameter_matrix))}")
    except Exception as e:
        print(f"\n错误发生: {str(e)}")
        print("故障排查建议:")
        print("1. 检查模板文件路径是否正确")
        print("2. 确认参数名称与模板文件中的变量名完全一致")
        print("3. 验证参数列表长度是否一致（或可广播）")
        print("4. 检查文件系统权限和磁盘空间")


