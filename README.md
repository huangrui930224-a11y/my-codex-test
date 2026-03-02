# 协议版 CSV 抖动统计（JRMS/J3u + EOJ）

本仓库当前提供 3 个对外脚本 + 1 个协议核心：

- `jitter_stats_from_csv.m`（独立调用：JRMS/J3u）
- `eoj_stats_from_csv.m`（独立调用：EOJ）
- `run_jitter_analysis_from_csv.m`（顶层接口）
- `protocol_pipeline_from_csv.m`（内部协议实现）

---

## 1. 与原协议一致的计算流程

核心实现遵循原版流程：

1. UI 重采样得到 `Y(M×N)`
2. UI中心采样 + `kmeans(s,4)` 恢复 `sym_code`
3. 每个 transition 采用专属阈值 `Vth=(V_level(A)+V_level(B))/2`
4. 线性插值提取 `t_cross`
5. `e_k = t_cross - t_ideal`
6. 连续时间 Golden PLL（必须使用 `Delta_t_k = t_cross(k)-t_cross(k-1)`）
7. 12类 context 分类
8. 强制类别均衡（`truncate`/`random`）
9. 计算 `Tavgi/S0i/fJ/JRMS/J3u`
10. 计算 EOJ（repeat-local UI 域，避免秒级大数相减）

---

## 2. 脚本调用

## 2.1 仅 JRMS/J3u

```matlab
cfg = struct();
cfg.UI = 100e-12;
cfg.context_patterns = {
  [0 1],[1 0],[1 2],[2 1],[2 3],[3 2],...
  [0 2],[2 0],[1 3],[3 1],[0 3],[3 0]
};
cfg.outputDir = './out_jitter';

outJ = jitter_stats_from_csv('wave.csv', cfg);
% outJ.JRMS, outJ.J3u
```

## 2.2 仅 EOJ

```matlab
cfg = struct();
cfg.UI = 100e-12;
cfg.context_patterns = {
  [0 1],[1 0],[1 2],[2 1],[2 3],[3 2],...
  [0 2],[2 0],[1 3],[3 1],[0 3],[3 0]
};
cfg.Npat = 8191;         % QPRBS13_CEI
cfg.eoj_mode = 'max';    % 或 'first3'
cfg.outputDir = './out_eoj';

outE = eoj_stats_from_csv('wave.csv', cfg);
% outE.EOJ_UI, outE.EOJ_s
```

## 2.3 顶层统一接口

```matlab
out = run_jitter_analysis_from_csv('wave.csv', cfg, 'both');
% mode = 'jitter' | 'eoj' | 'both'
```

---

## 3. 输入 CSV 自动识别规则

支持从路径读取 `.csv/.txt`，自动识别列：

- 时间列候选：`t/time/time_s/timestamp/sec/seconds`
- 电压列候选：`v/voltage/voltage_v/signal/diff/vdiff`

若自动识别失败，请手动指定：

```matlab
cfg.t_col = 'time_s';
cfg.v_col = 'voltage_v';
% 或
cfg.t_col_idx = 1;
cfg.v_col_idx = 2;
```

---

## 4. 必填参数与常用参数

### 必填

- `cfg.UI`
- `cfg.context_patterns`（12元素 cell）

### 常用可选

- `cfg.outputDir='.'`
- `cfg.M=32`
- `cfg.balance_mode='truncate'|'random'`
- `cfg.rng_seed=1`
- `cfg.edge_dir='either'|'rising'|'falling'`
- `cfg.Npat=8191`（QPRBS13）或 `511`（QPRBS9）
- `cfg.eoj_mode='max'|'first3'`
- `cfg.eoj_start_ui_m`（从第m个UI中点作为EOJ起点）
- `cfg.t0_override`（绝对起点秒值）
- `cfg.debug_plot=true|false`

---

## 5. 输出说明

### `jitter_stats_from_csv`

- `JRMS`、`J3u`
- `fJ`
- `Si/Tavgi/S0i`
- `counts_before/counts_after`

### `eoj_stats_from_csv`

- `EOJ_per_class_UI`、`EOJ_UI`
- `EOJ_per_class_s`、`EOJ_s`
- `EOJ_num_windows`、`EOJ_notes`

> 两个脚本都返回 `full_out`，包含协议核心完整中间量（`e_k_all/e_pll/e_cru/t_cross_all/...`）。
