# Utils 脚本说明

本目录用于管理 MOOSE 参数研究的生成、运行与后处理，推荐通过 `Driver.py` 一键控制各步骤。

## 目录结构

- activate_moose.sh：激活 `moose` conda 环境并执行 Python 脚本
- Driver.py：总控入口，按 RUN_STEPx 开关执行各步骤
- extract_pellet_total_strain_energy.py：从 .e 导出 `pellet_total_strain_energy`
- run_paraview.sh：准备 ParaView 环境并执行单例后处理脚本
- run_paraviewSeries.sh：准备 ParaView 环境并执行 series 后处理脚本
- step1_MeshGenerator.py：生成 parameter_studies 案例（参数笛卡尔积）
- step1_MeshGenerator_series.py：生成 parameter_studies_series 案例（参数对齐配对）
- step2_Run.py：批量运行 parameter_studies
- step2_Run_series.py：批量运行 parameter_studies_series
- step3_OperationalBriefing.py：统计收敛与日志信息
- step4_paraview_processor.py：ParaView 单例截图导出
- step4_paraview_processorSeries.py：ParaView 多目标时间截图导出
- step5_time_picture.py：时间/图片整理与组合图生成
- step6_export_timeseries.py：导出全时间域标量 CSV 并汇总

## 快速使用

1. 设置 `Driver.py` 中的步骤开关
2. 在工程目录下运行：

```bash
python /home/yp/projects/reproduction/Utils/Driver.py
```

## 步骤说明

### Step1 生成案例

- 输入模板：`Main.i`、`Sub.i`（可用环境变量覆盖）
- 输出目录：`parameter_studies` 或 `parameter_studies_series`
- 可选：自动插入 `[Outputs]` 的 Checkpoint 配置

### Step2 批量运行

- 自动激活 `moose` 环境
- 支持断点续跑（`*_my_checkpoint_cp`）
- 每个案例写 `run.log`

### Step3 收敛统计

- 解析 `run.log`
- 输出 `convergence_report.csv`
- 从输入文件注释中还原完整参数名

### Step4 ParaView 后处理

- 从 `file_base` 精确定位 `.e`
- 目标时间由 `TARGET_TIMES` 控制
- 输出目录固定在 `case_xxx/post_results/`

### Step5 时间/图片整理

- 汇总 `run.log` 时间信息
- 组合 `post_results/*_images` 输出大图

### Step6 全时间域导出

- 导出 `pellet_total_strain_energy` 等标量
- 输出 `*_step6_timeseries.csv`
- 自动生成跨案例汇总 CSV

## 常用环境变量

由 `Driver.py` 注入，亦可在命令行覆盖：

- PROJECT_BASE_DIR：工程根目录
- STUDIES_SUBDIR：案例目录（parameter_studies 或 parameter_studies_series）
- PV_FIELDS / PV_SERIES_FIELDS：截图字段列表
- PV_IMAGE_SIZE：图片分辨率（如 `1083,1083`）
- PV_OUTPUT_DIR：截图输出目录名（默认 `post_results`）
- PV_FIELDS_SINGLE2：Step6 导出字段
- PV_TIMESERIES_OUTPUT_DIR：Step6 输出目录名
- PV_TIMESERIES_REDUCTION：Step6 场量归约方式（mean/min/max）

## 输出目录约定

- 单案例截图：`case_xxx/post_results/<field>_images/*.png`
- Step6 CSV：`case_xxx/post_results/*_step6_timeseries.csv`
- 汇总 CSV：`parameter_studies/post_results/step6_merged_pellet_total_strain_energy.csv`

