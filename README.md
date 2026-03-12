# PAM4 Jitter Analysis (MATLAB)

本仓库提供一个无需 MATLAB toolbox 的 PAM4 抖动后处理脚本：

- 主脚本：`pam4_jitter_analysis.m`
- 演示脚本：`demo_run.m`

适用于 224G / 112G TX 前仿真波形后处理，输入为两列 CSV：

1. `time`（单位：s）
2. `vdiff`（单位：V）

---

## 快速使用

1. 准备 CSV 文件（两列：`time, vdiff`）。
2. 打开 `demo_run.m`，修改：
   - `csv_file`
   - `fb`
   - `M`（默认 64）
3. 运行 `demo_run.m`。

---

## 输入参数说明

`pam4_jitter_analysis(csv_file, fb, M, use_cru)`

- `csv_file`：输入波形 CSV 路径
- `fb`：波特率（Hz）
- `M`：每 UI 重采样点数（samples/UI），通常 `M=64`
- `use_cru`：是否启用 Golden PLL CRU（`true/false`，默认 `true`）

---

## 输出结果说明

函数返回结构体 `result`，主要字段包括：

- `JRMS_sec / JRMS_ps / JRMS_UI`
- `J3u_sec / J3u_ps / J3u_UI`
- `transition_names`
- `samples_per_class`
- `Tavgi_sec`
- `symbol_means`（V0/V1/V2/V3）
- `thresholds`（th01/th12/th23）
- `resampled`（重采样数据、Ts、M、重采样起点）
- `crossing_not_found`
- `ui_sampling_csv`（每个 UI 的采样相位/电压/符号 CSV）

---

## 15 步流程详解

> 下述编号与脚本内部 Step 标注一致。

### Step 1 读取数据

- 读取 CSV 两列数据（time、vdiff）。
- 清理 `NaN/Inf`。
- 对时间排序，保证单调递增。
- 检查时间是否等间隔：
  - 若不等间隔，先插值到等间隔时间轴。

### Step 1.5 估计首个 crossing（用于重采样起点）

- 根据用户要求，重采样起点改为“**存在 A->B 跳变的第一个过阈值采样点**”。
- 先在原始波形上做 4 电平聚类，得到粗符号 `A/B` 与电平均值。
- 由电平均值计算 `th01/th12/th23`，并按首个 `A->B` 跳变选择对应阈值。
- 使用该 `A->B` 阈值寻找第一个实际过阈值采样点，将其时间设为 `t_start_resample`（不做时间插值）。
- 若未找到，回退到 `t_raw(1)` 并告警。

### Step 2 重采样

- 定义 `UI = 1/fb`。
- 定义 `Ts = UI/M`。
- 构建统一重采样时间轴：`t_uniform = t_start_resample : Ts : t_raw(end)`。
- 使用线性插值得到 `v_uniform`。

### Step 3 UI 分块

- 计算可用 UI 数：`N_UI = floor(length(v_uniform)/M)`。
- 若不能整除，尾部样本丢弃并告警。
- 将波形重排为 `N_UI x M` 的 UI 矩阵。

### Step 4 采样相位选择

- 基线方式是取中心点（`m_center = round(M/2)`）。
- 为解决“最大开口不在中心”的问题，脚本会在 `1..M` 里搜索最优采样相位 `m_opt`：
  - 每个候选相位先做 4 电平聚类
  - 计算电平最小间距与簇内离散度
  - 选择分离度（`min_spacing / mean_within_std`）最大的相位
- 若搜索异常则回退到中心采样。

### Step 5 四电平聚类

- 对 `y_center` 进行 1D 四类聚类（无 toolbox 实现）。
- 聚类中心从低到高排序后映射到符号 `0/1/2/3`。

### Step 6 识别 AAAABB transition（候选）

- 扫描 `symbol`：
  - `symbol[n:n+3] = A`
  - `symbol[n+4:n+5] = B`
  - 且 `A != B`
- 记录候选 transition 的 UI 索引（12 类 A→B）。

### Step 7 计算符号均值

- 计算四电平均值：`V0, V1, V2, V3`。

### Step 8 计算阈值（方法 A）

- `th01 = (V0+V1)/2`
- `th02 = (V0+V2)/2`
- `th03 = (V0+V3)/2`
- `th12 = (V1+V2)/2`
- `th13 = (V1+V3)/2`
- `th23 = (V2+V3)/2`

### Step 9 提取所有 crossing time

- 对每个 `A != B` 的 UI 边界提取 crossing（不提前按 AAAABB 过滤）。
- 在“前一 UI 右半边 + 当前 UI 左半边”拼接窗口中搜索 crossing。
- 通过线性插值计算 crossing time。
- 多候选时，选离边界最近的 crossing。

### Step 10 Golden PLL CRU

- 定义：
  - `fc = fb/13280`
  - `alpha = exp(-2*pi*fc/fb)`
- 构建理想时钟 `t_ideal`。
- 计算 `tie_raw = tcross_abs - t_ideal`。
- 若 `use_cru=true`：PLL 每 UI 更新得到 `tie_LF`，并计算 `tie_HF = tie_raw - tie_LF`。
- 若 `use_cru=false`：跳过 PLL，令 `tie_HF = tie_raw`，后续步骤保持不变。

### Step 11 每类 transition 预处理

- 从 Step 6 候选中按类别提取对应 `tie_HF` 样本。
- 丢弃每类前 500 个样本（PLL 建立段）。
- 若剩余样本 `<100`，给告警。

### Step 12 强制样本数一致

- 统计每类剩余样本 `Ni`。
- 取 `Nmin = min(Ni)`。
- 每类仅保留前 `Nmin` 个样本，保证 12 类完全一致。

### Step 13 每类去均值

- 每类计算均值 `Tavgi`。
- 生成去均值样本：`Δt_i = Si - Tavgi`。

### Step 14 合并

- `Δt_all = [Δt_1; Δt_2; ...; Δt_12]`。

### Step 15 统计

- `JRMS = std(Δt_all)`。
- `t_low = prctile(Δt_all, 0.05%)`。
- `t_high = prctile(Δt_all, 99.95%)`。
- `J3u = t_high - t_low`。
- 输出 sec / ps / UI 三种单位。

---

## 绘图输出

脚本会自动生成：

1. 电平分布图（含聚类可视化）
2. 每类 transition 样本数柱状图
3. `Δt_all` 直方图（标注 0.05% 和 99.95%）

---

## 异常与告警

脚本内包含以下健壮性处理：

- CSV 不存在 / 列数不足
- 时间序列非法（非递增）
- 不等间隔时间自动插值
- 样本过少或 UI 数不足
- reshape 不能整除时尾部丢弃
- crossing 未找到告警
- 每类样本不足告警
- 类别为空时停止并报错

---

## Demo 额外输出

`demo_run.m` 还会：

- 打印符号均值、阈值
- 导出重采样数据 `resampled_waveform.csv`
- 导出每UI采样信息 `ui_sampling_per_ui.csv`（UI索引/采样相位/电压/符号）
- 打印重采样起点 `t_start_resample`

