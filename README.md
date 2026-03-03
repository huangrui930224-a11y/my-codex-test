# my-codex-test

MATLAB QPRBS13/PAM4 分析工具集合。仓库提供从 CSV 波形到符号推断、线性拟合、稳态估计、噪声统计和 SNDR 的模块化 pipeline。

## 总体流程

`qprbs13_pipeline.m` 串联执行：
1. 重采样与周期平均（`qprbs13_resample_module`）
2. 符号反推（`qprbs13_symbol_module`）
3. 符号旋转与矩阵构建（`qprbs13_rotate_x_module`）
4. 线性拟合与误差提取（`qprbs13_linear_fit_module`）
5. 稳态 `v_f` 估计（`qprbs13_steady_state_vf_module`）
6. 电平噪声估计（`qprbs13_level_noise_module`）
7. SNDR 计算（`qprbs13_sndr_module`）

---

## 子模块说明（逐文件）

### `qprbs13_pipeline.m`
- 顶层入口，负责：
  - 统一读取/校验 `cfg`
  - 调用所有子模块
  - 聚合 `out` 结构体
  - 保存 CSV/MAT/PNG 输出
- 输出中包含关键中间量：`Y`, `sym_code`, `x`, `X`, `P/E`, `sigmae`, `sigma_n`, `SNDR` 等。

### `qprbs13_resample_module.m`
- 功能：从 CSV (`t`,`vdiff`) 生成 `y_avg` 和 `Y`。
- 关键步骤：
  - 读取 CSV（支持有表头 `t,vdiff` 和无表头前两列）
  - 数据清洗（有限值过滤、时间排序去重）
  - 起点确定（默认首个 crossing；可切回局部扫描）
  - `numCycles` 周期重采样与平均
- 主要输出：`t0_best`, `y_avg`, `Y`, `Ts`, `Tpat`。

### `qprbs13_symbol_module.m`
- 功能：从 `Y` 中提取每 UI 的代表电压并反推 PAM4 符号。
- 采样模式：
  - `center`
  - `mean`
  - `proxy_opt`（按“类间距最大 / 类内方差最小”代理指标自动选最佳相位）
- 聚类：`kmeans` 或 `gmm`。
- 输出：`sym_code`, `sym_cycle`, `V`（四电平均值）, `ES1`, `ES2`, `RLM`, `x`。

### `qprbs13_rotate_x_module.m`
- 功能：对 `x(n)` 按 `Dp` 旋转得到 `xr`，并构造 `X`。
- 对应公式：Eq.(11-14) / Eq.(11-15)。
- 输出：`xr`, `X`, `Dp`。

### `qprbs13_linear_fit_module.m`
- 功能：线性拟合得到脉冲与误差。
- 关键计算：
  - 构建 `X1 = [X(1:TNp,:); ones(1,N)]`
  - 计算 `P`, `E = P*X1 - Y`
  - 提取 `e(k)`, `sigmae`, `p(k)`, `pmax`
- 对应公式：Eq.(11-16) ~ Eq.(11-18)。

### `qprbs13_steady_state_vf_module.m`
- 功能：用 `NpVf`（默认 20）重算线性拟合脉冲 `p_vf(k)`，估计稳态电压 `v_f`。
- 关键计算：`v_f = sum(p_vf)/M`（Eq.(11-12)）。
- 输出：`vf`, `p_vf`, `NpVf`, `TNpVf`。

### `qprbs13_level_noise_module.m`
- 功能：估计每个 PAM4 电平的噪声 RMS。
- 方法：
  - 只统计连续相同符号长度 `>= minRunLen` 的 run
  - 在固定相位点采样
  - 按 `merge` / `weighted` / `mean` 聚合每个电平 RMS
- 输出：
  - `sigma_levels`（4 电平）
  - `sigma_n = mean(sigma_levels)`
  - `runs_per_level`, `samples_per_level`。

### `qprbs13_sndr_module.m`
- 功能：基于 `pmax`, `sigmae`, `sigma_n` 计算 SNDR。
- 公式：`SNDR = 10*log10(pmax^2/(sigmae^2 + sigma_n^2))`（Eq.(27-5)）。
- 输出：`sndr_db` 及其构成项。

---

## 辅助工具/示例

### `resample_csv_7cycle_avg.m`
- 独立的“重采样 + 7 周期平均”工具，不依赖完整 pipeline。
- 适合先验证 `t0` 对齐和 `Y` 质量。

### `demo_resample_csv_7cycle_avg.m`
- 最小可运行 demo（重采样工具）。

### `DEMO.md`
- 中文调用文档，含完整 pipeline 示例、参数说明、输出文件清单、常见问题。

---

## 快速开始

1. 先看 `DEMO.md`，复制最小示例。
2. 只调三项即可起跑：
   - `cfg.csvPath`
   - `cfg.tStableStart`
   - `cfg.outputDir`
3. 若要提升判决稳定性，可尝试：
   - `cfg.uiSampleMethod = 'proxy_opt'`
   - `cfg.levelRmsMethod = 'weighted'`
