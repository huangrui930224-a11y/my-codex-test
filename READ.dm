# CSV抖动统计脚本说明（READ.dm）

本说明文档对应以下 3 个 MATLAB 脚本：

- `jitter_stats_from_csv.m`：计算 `JRMS` 和 `J3u`
- `eoj_stats_from_csv.m`：计算 `EOJ`
- `run_jitter_analysis_from_csv.m`：顶层统一调用接口

---

## 1. 统计方法说明

## 1.1 JRMS / J3u 统计方法（`jitter_stats_from_csv.m`）

输入为 `.csv` 文件中的时间 `t`（秒）和电压 `v`（伏特）列。

核心步骤：

1. 读取 CSV 并自动识别 `t/v` 列（也支持手动指定列名或列索引）。
2. 按 `UI` 和 `M` 做等间隔重采样：
   - `Ts = UI/M`
   - `t_uniform = t(1):Ts:t(end)`
3. 使用全局均值阈值 `thr = mean(v_uniform)` 搜索过零区间。
4. 对每个过零点做线性插值，得到 crossing 时间 `t_cross`。
5. 计算相邻 crossing 间隔偏差：
   - `jitter(k) = (t_cross(k+1)-t_cross(k)) - UI`
6. 计算统计量：
   - `JRMS = std(jitter)`
   - `J3u = prctile(jitter, pr_high) - prctile(jitter, pr_low)`

输出字段：

- `out.JRMS`
- `out.J3u`
- `out.jitter`
- `out.t_cross`
- `out.threshold`

---

## 1.2 EOJ 统计方法（`eoj_stats_from_csv.m`）

输入同样来自 CSV 的 `t/v`。

核心步骤：

1. 读取 CSV 并自动识别 `t/v` 列。
2. 用统一阈值 `thr = mean(v)` 提取 crossing 时间 `t_cross`。
3. 构建 trigger 起点 `t0`：
   - 若给 `cfg.eoj_start_ui_m`：`t0 = t(1) + (m+0.5)*UI`
   - 否则若给 `cfg.t0_override`：`t0 = t0_override`
   - 否则：`t0 = t(1)`
4. 计算 UI 索引、repeat 编号和 repeat 内相对时间：
   - `ui_index = floor((t_cross-t0)/UI)`
   - `repeat_id = floor(ui_index/Npat)+1`
   - `t_rel_UI = (t_cross - t_rep0)/UI`
5. 对每个 repeat 求均值轨迹 `Tr_UI`。
6. 用连续 3 个 repeat 形成窗口，按协议计算：
   - `T1=Tr(r), T2=Tr(r+1), T3=Tr(r+1), T4=Tr(r+2)`
   - `EOJwin = abs((T2-T1) - (T4-T3))`
7. 汇总：
   - `eoj_mode='max'`：取所有窗口最大值
   - `eoj_mode='first3'`：取第一个窗口
   - `EOJ_s = EOJ_UI * UI`

输出字段：

- `out.EOJ_UI`
- `out.EOJ_s`
- `out.EOJ_windows_UI`
- `out.Tr_UI`
- `out.repeat_id`
- `out.t0`
- `out.threshold`

---

## 2. 如何调用脚本

## 2.1 仅计算 JRMS/J3u

```matlab
cfg = struct();
cfg.UI = 100e-12;              % 必填
cfg.M = 32;                    % 可选
cfg.pr_low = 0.05;             % 可选
cfg.pr_high = 99.95;           % 可选
cfg.outputDir = './out_jitter';% 可选
cfg.debug_plot = true;         % 可选

outJ = jitter_stats_from_csv('wave.csv', cfg);
```

---

## 2.2 仅计算 EOJ

```matlab
cfg = struct();
cfg.UI = 100e-12;              % 必填
cfg.Npat = 8191;               % 可选，默认8191
cfg.eoj_mode = 'max';          % 可选：'max' 或 'first3'
cfg.eoj_start_ui_m = 0;        % 可选
% cfg.t0_override = 12e-9;      % 可选
cfg.outputDir = './out_eoj';   % 可选
cfg.debug_plot = true;         % 可选

outE = eoj_stats_from_csv('wave.csv', cfg);
```

---

## 2.3 顶层接口统一调用

```matlab
cfg = struct();
cfg.UI = 100e-12;
cfg.Npat = 8191;

% mode: 'jitter' | 'eoj' | 'both'
out = run_jitter_analysis_from_csv('wave.csv', cfg, 'both');
```

- `mode='jitter'`：返回 jitter 脚本输出
- `mode='eoj'`：返回 eoj 脚本输出
- `mode='both'`：返回结构体 `out.jitter` 与 `out.eoj`

---

## 3. CSV 列识别规则与参数填写

默认自动识别时间/电压列：

- 时间候选：`t/time/time_s/sec/seconds/timestamp`
- 电压候选：`v/voltage/voltage_v/signal/diff/vdiff`

若自动识别失败，可手动指定：

### 3.1 按列名指定

```matlab
cfg.t_col = 'time_s';
cfg.v_col = 'voltage_v';
```

### 3.2 按列索引指定（1-based）

```matlab
cfg.t_col_idx = 1;
cfg.v_col_idx = 2;
```

---

## 4. 必填参数清单（最小可运行）

- `csvPath`：CSV 文件路径
- `cfg.UI`：UI（秒）

> EOJ分析建议补充 `cfg.Npat`（QPRBS13=8191，QPRBS9=511）。
