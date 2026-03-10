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
  - `ref_repeat_pair_used`
  - `repeat_triplet_used`
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
   - `collect_all_crossing_events_from_symbols`
   - 先按符号跳变边界筛选（仅 `A->B, A!=B`）
   - 再按该跳变对应阈值 (`th01/th12/th23`) 与方向做 crossing 插值
2. **12类目标 transition crossing 提取**（给 EOJ 用）：
   - 先按符号跳变对筛选边界（如 `10`、`13`）
   - 再在该边界两侧局部窗口内（UI n 右半 + UI n+1 左半）按 `thr_type + dir` 做阈值 crossing 插值
   - 同类同 repeat 的所有命中都保留用于均值

### Step 7) 不做 CRU 滤波（按当前需求）

- EOJ 计算直接使用原始 crossing 时间：`tcross_cru_abs = tcross_abs`。
- `tie_LF_ui` 仅保留为全零占位，不参与 EOJ 修正。

### Step 8) 选择3个相邻 repeat 与参考 transition

- 调用 `choose_eoj_repeat_triplet_and_reference`。
- 选择 3 个相邻 repeat：`[repeatA, repeatB, repeatC]`。
- 任选一个参考 transition（手工指定或自动选择），用于周期估计：
  - `T3 = mean(ref @ repeatA)`
  - `T4 = mean(ref @ repeatB)`
  - `Tpat_est = T4 - T3`

### Step 9) 每类 EOJ_i（用另外一对相邻 repeat）

- 对每类 `i`：
  - `T1 = mean(i @ repeatB)`
  - `T2 = mean(i @ repeatC)`
  - `EOJ_i = abs((T2-T1) - (T4-T3))`
- 缺失类记为 `NaN` 并 warning。

### Step 10) 汇总输出

- `EOJ = max(EOJ_i, omitnan)`
- 同时输出秒/ps/UI单位。
- 输出详细结构体字段：
  - `EOJ_i_sec/ps`
  - `T1,T2,T3,T4,Tpat_est`
  - `ref_transition`
  - `repeat_pair_used`
  - `ref_repeat_pair_used`
  - `repeat_triplet_used`
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


### 4) 为什么 `symbol_trace` 看起来都有12类，但仍出现 missing?

`*_symbol_trace.csv` 中的 `is_aaaabb_transition_pos` 是**按 trans_def 窗口的期望位置标注**，不代表该位置一定检测到 crossing。

请以 `*_detected_transition_crossings.csv` 为准查看**实际检测到**的每条 transition crossing（含 repeat、UI、tcross_abs/tcross_cru_abs）。


### 4.1) 3个UI如何判断 AAAABB 上下文

以候选边界 `n->n+1`（`A->B`）为例，本工具采用 **run-length** 判据：

- 左侧 `A` 连续长度 `>=4 UI`
- 右侧 `B` 连续长度 `>=2 UI`

这与 `AAAABB`（边界位于第4个A与第1个B之间）等价。

之所以常说“3个UI上下文”也能判断，是因为边界本身已经贡献了两侧各1个UI：
- 左侧只需再看前面3个UI，就能确认是否达到4个连续A；
- 右侧只需再看后面1个UI，就能确认是否达到2个连续B。

实现上使用 run-length 而不是固定切片，优点是当实际序列为 `AAAAA...BBB...` 时同样可稳定命中，不会因为比 `AAAABB` 更长而漏检。

### 5) 自动 AAAABB 有 missing 如何优化

可增大 `cfg.auto_window_half_width_ui`（默认 1），自动推断时会扩展内部 crossing 搜索窗 `search_begin_ui/search_end_ui = pos(1)±half_window` 以提高命中。

`*_aaaabb_transitions.csv` 中导出的 `begin_ui/end_ui` 定义为 AAAABB 实际跨度：
- `begin_ui`：第一个 A 出现位置（`p0-3`）
- `end_ui`：第二个 B 结束位置（`p0+2`）

并结合 `*_detected_transition_crossings.csv` 检查每类在 `repeat_pair_used` 中是否都检测到事件。

另外会新增 EOJ 结果CSV：
- `*_eoj_summary.csv`：总 EOJ、`Tpat_est`、参考 transition、repeat 对、阈值等汇总
- `*_eoj_per_transition.csv`：12类的 `T1/T2/EOJ_i` 与在 `repeat_pair_used` 中的样本计数
