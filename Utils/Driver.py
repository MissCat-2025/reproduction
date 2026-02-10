#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import json

# ========= 步骤开关：想跑哪一步就设成 True =========
# RUN_STEP1 = True            # 网格生成（单次）
# RUN_STEP1_SERIES = True      # 网格生成（series）


RUN_STEP2 = True            # 运行 parameter_studies
# RUN_STEP2_SERIES = True      # 运行 parameter_studies_series

# RUN_STEP3 = True             # 收敛统计
# RUN_STEP4 = True             # ParaView 单例版
# RUN_STEP4_SERIES = True      # ParaView series 版
# RUN_STEP5 = True             # 时间 + 图片整理

# =========================== Step1 ============================


# Step1 单次参数矩阵和 Checkpoint 配置
template_main_name = "Main.i"
template_sub_name = "Sub.i"
STEP1_PARAM_MATRIX = {
    # "pellet_critical_energy": [3],
    # "PressureFactor": [1e6, 2e6, 3e6, 4e6, 5e6],
    "degradation_factor": [1e-9,1e-8,1e-7,1e-6,1e-5],
}

# Step1 series 参数矩阵
STEP1_SERIES_PARAM_MATRIX = {
    "a2": [3.1748, -0.5, 0.5396842, 1.3868],
    "a3": [0, 0, 0, 0.9106],
    "m": [2.5, 2, 4, 2],
}


STEP1_CHECKPOINT_CONFIG = '''
  [my_checkpoint]
    type = Checkpoint
    time_step_interval = 5
    num_files = 200
    wall_time_interval = 600
  []'''

# =========================== Step2 ============================


MPI_PROCESSES = 15


# =========================== Step3 ============================

# Step3 收敛统计目标子目录（相对于工程目录），为空则默认用 parameter_studies_series
STEP3_STUDIES_SUBDIR = "parameter_studies"

# =========================== Step4 ============================


# ParaView 字段和图像参数
PV_SINGLE_STUDIES_SUBDIR = "parameter_studies"  #要分析的路径 例如 "post_results"

PV_FIELDS_SINGLE = "d:相场变量,hoop_stress:环向应力,radial_stress:径向应力"
PV_FIELDS_SERIES = "d:相场变量,hoop_stress:环向应力,sigma0:断裂强度"
PV_IMAGE_SIZE = "1642,1083"
# 单例 ParaView 的输入/输出根目录（相对于工程目录），为空则默认用 parameter_studies

# =========================== Step5 ============================

# Step5 时间/图片整理目标子目录（相对于工程目录），为空则默认用 parameter_studies_series
STEP5_STUDIES_SUBDIR = PV_SINGLE_STUDIES_SUBDIR
# ParaView / 拼图使用的目标时间（秒）
TARGET_TIMES = [100000, 900000000]
# 时间拼图字体和布局参数
TIME_TITLE_FONT_SIZE = 72
TIME_LABEL_FONT_SIZE = 28
TIME_AXIS_FONT_SIZE = 24
TIME_LABEL_WIDTH = 220
TIME_TITLE_HEIGHT = 80
TIME_AXIS_HEIGHT = 100
TIME_ROW_AXIS_HEIGHT = 100
# ========= 集中参数配置 =========

# Utils 目录：放公共 step 脚本
UTILS_DIR = "/home/yp/projects/reproduction/Utils"

# 工程目录：自动用当前 driver.py 所在目录
PROJECT_BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# MOOSE 相关
MOOSE_APP = "/home/yp/projects/reproduction/reproduction-opt"


STUDIES_DIR_NAME = "parameter_studies"
STUDIES_SERIES_DIR_NAME = "parameter_studies_series"

PV_OUTPUT_DIR_SINGLE = "post_results"
PV_OUTPUT_DIR_SERIES = "parameter_studies_series"


def _step_enabled(name: str, default: bool = False) -> bool:
    return bool(globals().get(name, default))


def run_step(title, cmd, cwd, env):
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)
    result = subprocess.run(cmd, cwd=cwd, env=env)
    if result.returncode != 0:
        print(f"{title} 失败，返回码 = {result.returncode}")
        sys.exit(result.returncode)

def main():
    # 通用环境变量：所有 Utils 脚本共享
    env = os.environ.copy()
    env["PROJECT_BASE_DIR"] = PROJECT_BASE_DIR
    env["MOOSE_APP"] = MOOSE_APP
    env["MPI_PROCESSES"] = str(MPI_PROCESSES)

    # 模板文件名（Main / Sub 可以自定义）
    env["TEMPLATE_MAIN_NAME"] = template_main_name
    env["TEMPLATE_SUB_NAME"] = template_sub_name

    # 统一的目标时间字符串
    env["TARGET_TIMES"] = ",".join(str(t) for t in TARGET_TIMES)
    env["TIME_PICTURE_TARGET_TIMES"] = env["TARGET_TIMES"]

    env["TIME_TITLE_FONT_SIZE"] = str(TIME_TITLE_FONT_SIZE)
    env["TIME_LABEL_FONT_SIZE"] = str(TIME_LABEL_FONT_SIZE)
    env["TIME_AXIS_FONT_SIZE"] = str(TIME_AXIS_FONT_SIZE)
    env["TIME_LABEL_WIDTH"] = str(TIME_LABEL_WIDTH)
    env["TIME_TITLE_HEIGHT"] = str(TIME_TITLE_HEIGHT)
    env["TIME_AXIS_HEIGHT"] = str(TIME_AXIS_HEIGHT)
    env["TIME_ROW_AXIS_HEIGHT"] = str(TIME_ROW_AXIS_HEIGHT)

    # ===== Step 1: 生成参数研究案例 =====
    if _step_enabled("RUN_STEP1"):
        env["STUDIES_SUBDIR"] = STUDIES_DIR_NAME
        env["STEP1_CHECKPOINT_CONFIG"] = STEP1_CHECKPOINT_CONFIG
        env["STEP1_PARAM_MATRIX"] = json.dumps(STEP1_PARAM_MATRIX)
        step1_path = os.path.join(UTILS_DIR, "step1_MeshGenerator.py")
        run_step("Step 1 生成 parameter_studies",
                 [sys.executable, step1_path],
                 cwd=UTILS_DIR, env=env)

    if _step_enabled("RUN_STEP1_SERIES"):
        env["STUDIES_SUBDIR"] = STUDIES_SERIES_DIR_NAME
        env["STEP1_CHECKPOINT_CONFIG"] = STEP1_CHECKPOINT_CONFIG
        env["STEP1_SERIES_PARAM_MATRIX"] = json.dumps(STEP1_SERIES_PARAM_MATRIX)
        step1s_path = os.path.join(UTILS_DIR, "step1_MeshGenerator_series.py")
        run_step("Step 1(series) 生成 parameter_studies_series",
                 [sys.executable, step1s_path],
                 cwd=UTILS_DIR, env=env)

    # ===== Step 2: 运行 MOOSE 案例 =====
    if _step_enabled("RUN_STEP2"):
        env["STUDIES_SUBDIR"] = STUDIES_DIR_NAME
        step2_path = os.path.join(UTILS_DIR, "step2_Run.py")
        run_step("Step 2 运行 parameter_studies",
                 [sys.executable, step2_path],
                 cwd=UTILS_DIR, env=env)

    if _step_enabled("RUN_STEP2_SERIES"):
        env["STUDIES_SUBDIR"] = STUDIES_SERIES_DIR_NAME
        step2s_path = os.path.join(UTILS_DIR, "step2_Run_series.py")
        run_step("Step 2(series) 运行 parameter_studies_series",
                 [sys.executable, step2s_path],
                 cwd=UTILS_DIR, env=env)

    # ===== Step 3: 收敛统计 =====
    if _step_enabled("RUN_STEP3"):
        # 默认统计 parameter_studies_series；如果指定了目录则使用该目录
        studies_subdir = STEP3_STUDIES_SUBDIR if STEP3_STUDIES_SUBDIR else STUDIES_SERIES_DIR_NAME
        env["STUDIES_SUBDIR"] = studies_subdir
        step3_path = os.path.join(UTILS_DIR, "step3_OperationalBriefing.py")
        run_step("Step 3 收敛统计",
                 [sys.executable, step3_path],
                 cwd=UTILS_DIR, env=env)

    # ===== Step 4: ParaView 后处理 =====
    if _step_enabled("RUN_STEP4"):
        # 单例版：如果未指定目录则默认用 parameter_studies；否则使用指定的图片路径
        studies_subdir = PV_SINGLE_STUDIES_SUBDIR if PV_SINGLE_STUDIES_SUBDIR else STUDIES_DIR_NAME
        env["STUDIES_SUBDIR"] = studies_subdir
        env["PV_FIELDS"] = PV_FIELDS_SINGLE
        env["PV_IMAGE_SIZE"] = PV_IMAGE_SIZE
        env["PV_OUTPUT_DIR"] = PV_OUTPUT_DIR_SINGLE
        # 通过 run_paraview.sh 启动 ParaView 专用环境并运行处理脚本
        pv_single = os.path.join(UTILS_DIR, "run_paraview.sh")
        run_step("Step 4 单例 ParaView 截图",
                 ["bash", pv_single],
                 cwd=UTILS_DIR, env=env)

    if _step_enabled("RUN_STEP4_SERIES"):
        env["STUDIES_SUBDIR"] = STUDIES_SERIES_DIR_NAME
        env["PV_SERIES_FIELDS"] = PV_FIELDS_SERIES
        env["PV_IMAGE_SIZE"] = PV_IMAGE_SIZE
        env["PV_OUTPUT_DIR_SERIES"] = PV_OUTPUT_DIR_SERIES
        # run_paraviewSeries.sh 里会调用 step4_paraview_processorSeries.py
        pv_series = os.path.join(UTILS_DIR, "run_paraviewSeries.sh")
        run_step("Step 4(series) ParaView 截图",
                 ["bash", pv_series],
                 cwd=UTILS_DIR, env=env)

    # ===== Step 5: 时间 + 图片整理 =====
    if _step_enabled("RUN_STEP5"):
        # 指定要分析哪个目录：默认 parameter_studies_series，如指定则用 STEP5_STUDIES_SUBDIR
        studies_subdir = STEP5_STUDIES_SUBDIR if STEP5_STUDIES_SUBDIR else STUDIES_SERIES_DIR_NAME
        target_studies_dir = os.path.join(PROJECT_BASE_DIR, studies_subdir)
        env["TARGET_STUDIES_DIR"] = target_studies_dir

        step5_path = os.path.join(UTILS_DIR, "step5_time_picture.py")
        run_step("Step 5 时间与图片整理",
                 [sys.executable, step5_path],
                 cwd=UTILS_DIR, env=env)

    print("\n✅ 全部启用的步骤已执行完毕")

if __name__ == "__main__":
    main()
