# my-codex-test

本仓库提供一个可直接运行的 MATLAB EOJ（Even/Odd Jitter）分析流程：

- 主函数：`eoj_pam4_from_csv.m`
- 示例脚本：`demo_eoj_run.m`

输入为两列 CSV：

1. `time`（秒，严格递增）
2. `vdiff`（伏特）

---

## 快速使用

1. 准备 CSV 文件（两列：`time, vdiff`）。
2. 修改 `demo_eoj_run.m` 中的 `csv_file` 路径。
3. 运行：
   ```matlab
   demo_eoj_run
   ```
4. 查看输出：
   - `EOJ_ps`
   - `EOJ_i_ps`
   - `ref_transition`
   - `Tpat_est`
   - `repeat_pair_used`
   - 电平与阈值
   - 采样相位搜索信息

---

## 代码流程（逐步说明）

下面按 `eoj_pam4_from_csv.m` 的主流程说明每一步。

### Step 0) 参数与输入检查

- 通过 `apply_defaults(cfg)` 补齐默认参数：
  - `fb = 112e9`
  - `Npat = 8191`
  - `M = 64`
  - `fc = fb/13280`
  - `discard_first_repeats = 0`
  - `do_plot = true`
  - `verbose = true`
- 通过 `read_csv_two_cols` 读取 CSV，并检查：
  - 列数>=2
  - 时间严格递增
  - NaN 清理后样本足够

### Step 1) 等间隔重采样

- 计算：
  - `UI = 1/fb`
  - `Ts = UI/M`
- 生成统一时间轴：`t_uniform = t(1):Ts:t(end)`
- 线性插值：`v_uniform = interp1(t, v, t_uniform, 'linear')`
- 截断到完整 UI 数：`N_UI * M`

### Step 2) UI 分块 + 最优采样相位搜索

- `Vui = reshape(v_uniform, M, N_UI).'`，每行一个 UI。
- 不再固定中心点；先做相位搜索：
  - 遍历 `mm = 1..M`
  - 取 `x = Vui(:, mm)`
  - 对 `x` 做 `kmeans_1d_no_toolbox(x,4,60)`
  - 计算评分：
    - `min_spacing = min(diff(sort(centers)))`
    - `mean_wstd = mean(within_cluster_std)`
    - `score = min_spacing/(mean_wstd + eps)`
  - 取分数最大相位 `m_opt`
- 若搜索失败，则回退 `m_center = round(M/2)`。
- 最终采样序列：`y = Vui(:, m_used)`。

### Step 3) 四电平符号与电平均值

- 对 `y` 执行 1D kmeans（无 toolbox）分成 4 类。
- 输出 `sym ∈ {0,1,2,3}`。
- 计算四电平均值：`V0..V3`。

### Step 4) 阈值计算

采用 Method A：

- `th01 = (V0+V1)/2`
- `th12 = (V1+V2)/2`
- `th23 = (V2+V3)/2`

### Step 5) transition 定义（12类）

二选一：

1. **手工模式**：使用 `cfg.trans_def(1..12)`。
2. **自动模式**：`infer_transitions_from_symbols` 先在 `offset=0..Npat-1` 上做 pattern 起点搜索（按 AAAABB 覆盖率与命中数选最优 offset），再基于该对齐窗口识别 12 类有向 transition。

### Step 6) crossing 提取

分两部分：

1. **全 crossing 提取**（给 CRU 用）：
   - `collect_all_crossing_events`
   - 在 `th01/th12/th23` 上找所有过阈值 crossing
   - 线性插值求 `tcross_abs`
2. **12类目标 transition crossing 提取**（给 EOJ 用）：
   - 每类按自己的 `thr_type + dir + window` 查找
   - 一个窗口内多 crossing 时取第一个

### Step 7) CRU（golden PLL）校正

- 先把所有 crossing 映射到 UI 索引。
- 计算 raw TIE：`tie_raw = tcross_abs - (t0 + n*UI)`。
- 每 UI 更新一阶 LPF：
  - `alpha = exp(-2*pi*fc/fb)`
  - 有事件用事件值，无事件沿用上一 UI 值（保持带宽正确）
- 对 12 类目标 crossing 做 CRU 校正：
  - `tcross_cru_abs = tcross_abs - tie_LF_ui(n)`

### Step 8) 参考 transition 与 repeat 对选择

- 调用 `choose_ref_transition_and_repeats`。
- 优先默认相邻 repeat 对：`[discard, discard+1]`。
- 若不可用，自动回退到存在数据的相邻 repeat 对。
- 参考周期估计：
  - `T3 = mean(ref @ rep0)`
  - `T4 = mean(ref @ rep1)`
  - `Tpat_est = T4 - T3`

### Step 9) 每类 EOJ_i

- 对每类 `i`：
  - `T1 = mean(i @ rep0)`
  - `T2 = mean(i @ rep1)`
  - `EOJ_i = abs((T2-T1) - Tpat_est)`
- 缺失类记为 `NaN` 并 warning。

### Step 10) 汇总输出

- `EOJ = max(EOJ_i, omitnan)`
- 同时输出秒/ps/UI单位。
- 输出详细结构体字段：
  - `EOJ_i_sec/ps`
  - `T1,T2,T3,T4,Tpat_est`
  - `ref_transition`
  - `repeat_pair_used`
  - `V0..V3`
  - `th01/th12/th23`
  - `symbol_mean_voltage`
  - `crossing_threshold`
  - `phase_search`
  - `all_crossings`

### Step 11) 可视化

`cfg.do_plot=true` 时：

- 中心采样直方图（含电平和阈值）
- `EOJ_i` 柱状图
- 参考 transition 在两次 repeat 的 CRU crossing 分布对比

---

## 常见问题

### 1) 报错：`No valid transition can be used as reference ...`

说明默认参考 repeat 对中没有可用 transition；代码已支持自动回退到其他相邻 repeat 对。建议检查：

- `cfg.discard_first_repeats`
- `cfg.trans_def` 是否准确
- 输入波形覆盖长度是否足够

### 2) crossing 是过零点吗？

不是。当前实现是**过阈值点**（`th01/th12/th23`）并线性插值。

### 3) 如何提高稳定性？

优先提供准确的 `cfg.trans_def`，并根据数据建立过程设置 `discard_first_repeats`。
