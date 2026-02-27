# my-codex-test

## `cru_transition_jitter_pipeline.m` 使用说明

本文档说明 `cru_transition_jitter_pipeline.m` 的输入参数、配置项、调用方法、输出字段及含义。

---

## 1. 函数功能概述

`cru_transition_jitter_pipeline` 用于对高速链路波形进行：

1. 基于 UI 的重采样与 4 电平恢复（PAM4）；
2. transition crossing time 提取；
3. 连续时间 Golden PLL 去趋势得到 `e_cru`；
4. 12 类 transition 分类与均衡；
5. 抖动统计（`JRMS`、`J3u`）；
6. EOJ（End-of-Jitter）统计（按 repeat、按类输出）。

---

## 2. 函数原型

```matlab
out = cru_transition_jitter_pipeline(t, v, cfg)
```

- `t`：`Nx1 double`，时间（秒）
- `v`：`Nx1 double`，差分电压（伏特）
- `cfg`：结构体配置

---

## 3. 输入配置 `cfg`

### 3.1 必填字段

- `cfg.UI`
  - 单位间隔（秒）
- `cfg.outputDir`
  - 输出目录（用于保存 PNG/CSV）
- `cfg.context_patterns`
  - 12 类 transition 的判定规则（12 元 cell）

### 3.2 可选字段（含默认值）

- `cfg.M = 32`
  - 每 UI 重采样点数
- `cfg.balance_mode = 'truncate'`
  - 类间样本均衡方式：`'truncate'` 或 `'random'`
- `cfg.rng_seed = 1`
  - `balance_mode='random'` 时的随机种子
- `cfg.edge_dir = 'either'`
  - 过零边沿方向：`'either' | 'rising' | 'falling'`
- `cfg.th_method = 'midpoint'`
  - 阈值方法（当前仅支持 midpoint）
- `cfg.debug_plot = true`
  - 是否输出 CRU/统计调试图
- `cfg.Npat = 8191`
  - pattern 长度（QPRBS13_CEI 用 8191；QPRBS9_CEI 用 511）
- `cfg.eoj_mode = 'max'`
  - EOJ 窗口策略：`'max'`（取所有窗口最大）或 `'first3'`（取首个三 repeat 窗）
- `cfg.eoj_debug_plot = false`
  - 是否输出每类 EOJ 调试图
- `cfg.t0_override = []`
  - 绝对时间锚点（秒），用于 override EOJ repeat 起点
- `cfg.eoj_start_ui_m = []`
  - EOJ 触发起点设为“第 `m` 个 UI 的中点”
  - 当该字段非空时，优先级高于 `t0_override`

### 3.3 `m`（`cfg.eoj_start_ui_m`）取值要求

- `m` 必须是：
  - 标量
  - `double`
  - 有限值
  - 非负整数
- 并且满足：
  - `0 <= m <= N-1`
  - 其中 `N` 为重采样后 UI 数量（函数内部计算）
- 建议：
  - 为获得稳定 EOJ 三 repeat 窗，`m` 不要过大；若剩余长度不足，函数会给 warning

---

## 4. 调用示例

### 4.1 基础调用

```matlab
t = your_time_vector;      % Nx1
d = your_diff_voltage;     % Nx1

cfg = struct();
cfg.UI = 100e-12;
cfg.outputDir = './output_demo';
cfg.context_patterns = {
    [0 1], [1 0], [1 2], [2 1], [2 3], [3 2], ...
    [0 2], [2 0], [1 3], [3 1], [0 3], [3 0] ...
};

out = cru_transition_jitter_pipeline(t, d, cfg);
```

### 4.2 从第 `m` 个 UI 中点作为 EOJ trig 起点

```matlab
cfg.eoj_start_ui_m = 1000;   % 第1000个UI中点作为EOJ起点
out = cru_transition_jitter_pipeline(t, d, cfg);
```

### 4.3 使用绝对时间锚点（次优先级）

```matlab
cfg.t0_override = 12.5e-9;
out = cru_transition_jitter_pipeline(t, d, cfg);
```

> 若同时设置 `eoj_start_ui_m` 与 `t0_override`，将优先使用 `eoj_start_ui_m`。

---

## 5. 输出 `out` 字段说明

### 5.1 基础与CRU相关

- `out.UI`：UI（秒）
- `out.fc`：Golden PLL 截止频率（Hz）
- `out.omega_c`：Golden PLL 角频率（rad/s）
- `out.t_cross_all`：所有有效 transition 的 crossing time（秒）
- `out.e_k_all`：原始 time error（秒）
- `out.e_pll`：PLL 跟踪项（秒）
- `out.e_cru`：高通后残差 jitter（秒）

### 5.2 电平与符号恢复

- `out.sym_code`：每 UI 恢复的符号编码（0~3）
- `out.V_level`：4 个电平估计均值

### 5.3 12 类分类与统计

- `out.Si`：12 类各自 jitter 样本（均衡后）
- `out.Tavgi`：每类均值
- `out.S0i`：每类去均值后样本
- `out.fJ`：12 类 `S0i` 拼接后的总体样本
- `out.JRMS`：`fJ` 的 RMS 抖动
- `out.J3u`：`99.95% - 0.05%` 峰峰值
- `out.counts_before`：均衡前各类样本数
- `out.counts_after`：均衡后各类样本数

### 5.4 事件级输出（EOJ依赖）

- `out.t_cross_event`：事件 crossing 时间（秒）
- `out.ui_index_event`：事件对应全局 UI 索引
- `out.class_id_event`：事件类别（0 或 1~12）
- `out.t0`：EOJ repeat 对齐锚点
- `out.Npat`：pattern 长度
- `out.eoj_start_ui_m`：若启用，记录使用的 `m`

### 5.5 EOJ输出

- `out.EOJ_per_class_UI`：每类 EOJ（UI）
- `out.EOJ_UI`：12 类最大 EOJ（UI）
- `out.EOJ_per_class_s`：每类 EOJ（秒）
- `out.EOJ_s`：总体 EOJ（秒）
- `out.EOJ_num_windows`：每类使用的有效窗口数
- `out.EOJ_notes`：每类备注（如 `max` / `first3` / 无可用窗口）
- `out.EOJ_Tr_by_class`：每类按 repeat 的均值轨迹
- `out.EOJ_repeat_ids_by_class`：每类对应的 repeat 编号
- `out.EOJ_anchor_delta_ui_mean/std/maxabs`：锚点一致性诊断

---

## 6. 文件输出（`cfg.outputDir`）

### 6.1 常规调试图（`cfg.debug_plot=true`）

- `cru_debug.png`
- `histogram_fJ.png`
- `counts_before_after.png`

### 6.2 EOJ输出

- `eoj_per_class.csv`：每类 EOJ 表格
- `eoj_bar.png`：每类 EOJ 柱状图
- `eoj_debug_classi.png`：可选，每类 EOJ 调试图（`cfg.eoj_debug_plot=true`）

---

## 7. 常见注意事项

1. `t`、`v` 必须为列向量，且长度一致。
2. `cfg.context_patterns` 必须是 12 元 cell。
3. 若某一类没有样本，均衡步骤会报错。
4. 若数据中可用 repeat 太少，EOJ 可能出现 `NaN` 或窗口不足告警。
5. 建议优先使用 `eoj_start_ui_m` 来定义“从第 m 个 UI 中点触发”的统计起点。
