#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import json

# ========= 步骤开关：想跑哪一步就设成 True =========
RUN_STEP1 = True            # 网格生成（单次）
RUN_STEP2 = True            # 运行 parameter_studies

# RUN_STEP1_SERIES = True      # 网格生成（series）
# RUN_STEP2_SERIES = True      # 运行 parameter_studies_series

RUN_STEP3 = True             # 收敛统计
RUN_STEP4 = True             # ParaView 单例版
RUN_STEP5 = True             # 时间 + 图片整理
RUN_STEP6 = True             # 全时间域标量导出


# RUN_STEP4_SERIES = True      # ParaView series 版
# =========================== Step1 ============================


# Step1 单次参数矩阵和 Checkpoint 配置
template_main_name = "MainA5.i"
template_sub_name = "Sub.i"
STEP1_PARAM_MATRIX = {
    "pellet_critical_energy": [10],
    # "pellet_critical_fracture_strength": [60e6,70e6,80e6],
    # "CGc": [0.003,0.0035,0.004],
    # "PressureFactor": [0,1e6, 2e6, 3e6, 4e6, 5e6, 6e6],
    # "degradation_factor": [1e-9,1e-8,1e-7,1e-6,1e-5],
    # "PowerTime": [30,60,120,480],
    # "WeibullSeed": [1,2,3,4]
    # "coolant_heat_transfer_coefficient_out": [2500,3000,3500]
    # "Ndt": [50,100,200,300,400],
    # "length_scale_paramete": [7e-5,6e-5,5e-5,4e-5,3e-5]
    # "largestPoreSize0":[20,30,40,50]
    # "LinearPower":[100,85,70,55,40,25,10]
    # "fission_rate":[1.2e+19,2.4e19,3.6e19,4.8e19]
}

# Step1 series 参数矩阵
STEP1_SERIES_PARAM_MATRIX = {
    "a2": [3.1748, -0.5, 0.5396842, 1.3868],
    "a3": [0, 0, 0, 0.9106],
    "m": [2.5, 2, 4, 2],
}

# STEP1_CHECKPOINT_CONFIG = '''
#   [my_checkpoint]
#     type = Checkpoint
#     time_step_interval = 5
#     num_files = 1
#     wall_time_interval = 600
#   []'''
STEP1_CHECKPOINT_CONFIG = '''
 '''

# =========================== Step2 ============================


MPI_PROCESSES = 9


# =========================== Step3 ============================

# Step3 收敛统计目标子目录（相对于工程目录），为空则默认用 parameter_studies_series
STEP3_STUDIES_SUBDIR = "parameter_studies"

# =========================== Step4 ============================


# ParaView 字段和图像参数
PV_SINGLE_STUDIES_SUBDIR = "parameter_studies"  #要分析的路径 例如 "post_results"

# PV_FIELDS_SINGLE = "d:相场变量,hoop_stress:环向应力,radial_stress:径向应力"
PV_FIELDS_SINGLE = "d:相场变量,sigma0:断裂强度,Gc:Gc"

# PV_FIELDS_SERIES = "d:相场变量,hoop_stress:环向应力,sigma0:断裂强度"
PV_IMAGE_SIZE = "1083,1083"
# 单例 ParaView 的输入/输出根目录（相对于工程目录），为空则默认用 parameter_studies

# =========================== Step5 ============================

# Step5 时间/图片整理目标子目录（相对于工程目录），为空则默认用 parameter_studies_series
STEP5_STUDIES_SUBDIR = PV_SINGLE_STUDIES_SUBDIR
# ParaView / 拼图使用的目标时间（秒）
TARGET_TIMES = [50000,200000,2e7]

# =========================== Step7 ============================

# Step7 导出的时刻和网格点场变量数据
DATA_TIMES = [50000, 200000, 2e7]
DATA_FIELDS = "sigma0:断裂强度"

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


PV_FIELDS_SINGLE2 = "pellet_total_strain_energy:总应变能量密度"  # 导出整个时间域的应变能量密度
PV_TIMESERIES_OUTPUT_DIR = PV_OUTPUT_DIR_SINGLE


def _step_enabled(name: str, default: bool = False) -> bool:
    return bool(globals().get(name, default))


def run_step(title, cmd, cwd, env):
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)
    print(f"cwd: {cwd}")
    print("cmd: " + " ".join(cmd))
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

    # 统一的数据导出参数（Step 7）
    if "DATA_TIMES" in globals():
        env["DATA_TARGET_TIMES"] = ",".join(str(t) for t in DATA_TIMES)
    if "DATA_FIELDS" in globals():
        env["DATA_FIELDS"] = DATA_FIELDS
    if "DATA_N" in globals():
        env["DATA_N"] = str(DATA_N)

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
        env["PV_OUTPUT_DIR"] = PV_OUTPUT_DIR_SINGLE

        step5_path = os.path.join(UTILS_DIR, "step5_time_picture.py")
        run_step("Step 5 时间与图片整理",
                 [sys.executable, step5_path],
                 cwd=UTILS_DIR, env=env)

    if _step_enabled("RUN_STEP6"):
        env["STUDIES_SUBDIR"] = STUDIES_DIR_NAME
        env["PV_FIELDS_SINGLE2"] = PV_FIELDS_SINGLE2
        env["PV_TIMESERIES_OUTPUT_DIR"] = PV_TIMESERIES_OUTPUT_DIR
        env["PV_TIMESERIES_MERGED_CSV"] = os.path.join(PROJECT_BASE_DIR, STUDIES_DIR_NAME, PV_TIMESERIES_OUTPUT_DIR, "step6_merged_pellet_total_strain_energy.csv")
        env["PV_PYTHON_SCRIPT"] = os.path.join(UTILS_DIR, "step6_export_timeseries.py")
        pv_single = os.path.join(UTILS_DIR, "run_paraview.sh")
        run_step("Step 6 导出全时间域标量CSV",
                 ["bash", pv_single],
                 cwd=UTILS_DIR, env=env)

    # ===== Step 7: 导出特定时刻全网格点数据 =====
    if _step_enabled("RUN_STEP7"):
        studies_subdir = PV_SINGLE_STUDIES_SUBDIR if PV_SINGLE_STUDIES_SUBDIR else STUDIES_DIR_NAME
        env["STUDIES_SUBDIR"] = studies_subdir
        env["PV_PYTHON_SCRIPT"] = os.path.join(UTILS_DIR, "step7_export_mesh_data.py")
        pv_single = os.path.join(UTILS_DIR, "run_paraview.sh")
        run_step("Step 7 导出特定时刻全网格点数据",
                 ["bash", pv_single],
                 cwd=UTILS_DIR, env=env)

    print("\n✅ 全部启用的步骤已执行完毕")

if __name__ == "__main__":
    main()
