# MATLAB 调用 Demo（用于快速测试）

本文档给出两个可直接复制执行的测试流程：

1. `qprbs13_pipeline.m`（推荐，完整端到端）
2. `resample_csv_7cycle_avg.m`（仅重采样与 7 周期平均）

---

## 0. 前置条件

- 已安装 MATLAB（建议带 Statistics and Machine Learning Toolbox，用于 `kmeans`/`fitgmdist`）。
- 当前仓库目录下已包含：
  - `qprbs13_pipeline.m`
  - `resample_csv_7cycle_avg.m`
- 支持两种 CSV：
  1) 有表头：`t`, `vdiff`（单位 s / V）
  2) 无表头：默认第 1 列为时间，第 2 列为电压。

---

## 1. 最小可运行示例（完整 pipeline）

在 MATLAB Command Window 中执行：

```matlab
cd('/workspace/my-codex-test');

cfg = struct();
cfg.csvPath = '/path/to/your/data.csv';   % 改成你的 CSV 路径
cfg.tStableStart = 1e-6;                  % 当 useFirstCrossing=false 时作为扫描中心
cfg.outputDir = fullfile(pwd, 'demo_out_pipeline');

% 可选参数（不填则走默认）
% cfg.ui = 1/112e9;
% cfg.M = 32;
% cfg.N = 8191;
% cfg.numCycles = 7;
% cfg.dtScan = [];             % 空则使用 Ts/8
% cfg.scanWindowUI = 0.5;
% cfg.useFirstCrossing = true; % true: 起点固定为第一个 crossing time
% cfg.uiSampleMethod = 'center'; % 'center' | 'mean' | 'proxy_opt'
% cfg.clusterMethod = 'kmeans';  % or 'gmm'
% cfg.saveCycles = true;
% cfg.Dp = 4;              % pulse delay for xr/X (Eq.11-14/11-15)
% cfg.Np = 29;              % linear-fit parameter
% cfg.TNp = 2*cfg.Np+1;     % rows from X used in Eq.(11-16)
% cfg.NpVf = 20;             % steady-state refit uses Np=20
% cfg.TNpVf = 2*cfg.NpVf+1;  % rows from X used for vf refit
% cfg.minRunLen = 6;          % >=6 连续相同符号 run 才参与噪声统计
% cfg.fixedSampleMethod = 'center'; % or 'index'
% cfg.fixedSampleIndex = 16;   % 当 method='index' 时生效
% cfg.levelRmsMethod = 'merge';% or 'weighted' or 'mean'

out = qprbs13_pipeline(cfg);

fprintf('Done. t0_best = %.16g s\n', out.t0_best);
fprintf('ES1 = %.6f, ES2 = %.6f, RLM = %.6f\n', out.ES1, out.ES2, out.RLM);
```

### pipeline 输出文件

运行后会在 `cfg.outputDir` 下看到：

- `resample_out.mat`
- `y_avg.csv`
- `symbol_means.csv`
- `s_x_symbols.csv`
- `symbol_voltage_samples.csv`
- `J_scan.png`
- `heatmap_Y.png`
- `s_hist_centers.png`
- `levels_boxplot.png`
- `x_preview.png`
- `xr.csv`
- `X.mat`
- `linear_fit.mat`
- `e_waveform.csv`
- `p_pulse.csv`
- `e_preview.png`
- `p_preview.png`
- `p_pulse_vf.csv`
- `steady_state_vf.mat`
- `p_vf_preview.png`
- `level_noise_rms.csv`
- `level_noise_rms.mat`
- `level_noise_rms.png`
- `sndr.csv`
- `sndr.mat`
- `sndr_summary.png`

---

## 2. 仅验证重采样与 7 周期平均

```matlab
cd('/workspace/my-codex-test');

csvFile = '/path/to/your/data.csv';
tStableStart = 1e-6;

cfg = struct();
cfg.ui = 1/112e9;
cfg.M = 32;
cfg.N = 8191;
cfg.numCycles = 7;
cfg.outputDir = fullfile(pwd, 'demo_out_resample');
% cfg.dtScan = [];          % 空则默认 Ts/4
% cfg.overlayUI = [1, 4096, 8191];

out = resample_csv_7cycle_avg(csvFile, tStableStart, cfg);

fprintf('Done. t0_best = %.16g s\n', out.t0_best);
fprintf('Y size = [%d, %d]\n', size(out.Y,1), size(out.Y,2));
```

### resample 输出文件

- `y_avg.mat`
- `Y.mat`
- `y_avg.csv`
- `J_scan.png`
- `overlay_cycles.png`
- `heatmap_Y.png`

---

## 3. 常见报错与排查

### 报错：CSV 列相关错误
- 若是有表头 CSV，推荐列名为 `t, vdiff`。
- 若无表头 CSV，请确保前两列分别是时间与电压，且数据为数值。

### 报错：`No valid t0 ...` / `Insufficient data coverage ...`
- `tStableStart` 附近数据不够覆盖 7 个周期；
- 尝试：
  - 选更靠中间的 `tStableStart`；
  - 放宽扫描窗 `scanWindowUI`；
  - 使用更长采样数据。

### 报错：聚类后某类计数为 0
- 信号质量或采样点提取方式导致 4 类不可分；
- 尝试：
  - `cfg.uiSampleMethod = 'mean'`；
  - `cfg.clusterMethod = 'gmm'`；
  - 检查输入波形是否确为 PAM4 稳定区间。

---

## 4. 一键复用 demo 模板

如需多次跑不同文件，可以封装成小函数：

```matlab
function out = run_demo(csvPath, tStableStart, outDir)
    cfg = struct();
    cfg.csvPath = csvPath;
    cfg.tStableStart = tStableStart;
    cfg.outputDir = outDir;
    cfg.uiSampleMethod = 'center';
    cfg.clusterMethod = 'kmeans';
    out = qprbs13_pipeline(cfg);
end
```

然后：

```matlab
out = run_demo('/path/to/data.csv', 1e-6, fullfile(pwd,'batch_out_01'));
```


---

## 5. 分模块调用（新增）

现在已拆成多个子模块，顶层 `qprbs13_pipeline` 会自动串联：

- `qprbs13_resample_module(cfg)`：CSV -> 对齐扫描 -> `y_avg` / `Y`
- `qprbs13_symbol_module(Y, cfg)`：`Y` -> `sym_code` / `ES1` / `ES2` / `RLM` / `x`
- `qprbs13_rotate_x_module(x, Dp)`：`x` -> `xr` / `X`
- `qprbs13_linear_fit_module(Y, X, cfg)`：`P` / `E` / `e` / `sigma_e` / `p` / `pmax`
- `qprbs13_steady_state_vf_module(Y, X, M, cfg)`：`Np=20` 重算脉冲并估计 `v_f`
- `qprbs13_level_noise_module(Y, sym_code, cfg)`：4 电平噪声 RMS 与 `sigma_n`
- `qprbs13_sndr_module(pmax, sigma_e, sigma_n)`：计算 `SNDR`

你也可以手动分步调用：

```matlab
cd('/workspace/my-codex-test');

cfg = struct();
cfg.csvPath = '/path/to/your/data.csv';
cfg.tStableStart = 1e-6;
cfg.outputDir = fullfile(pwd, 'demo_out_split');
cfg.ui = 1/112e9;
cfg.M = 32;
cfg.N = 8191;
cfg.numCycles = 7;
cfg.scanWindowUI = 0.5;
cfg.dtScan = (cfg.ui/cfg.M)/8;
cfg.uiSampleMethod = 'center';
cfg.clusterMethod = 'kmeans';
cfg.saveCycles = true;

res = qprbs13_resample_module(cfg);
sym = qprbs13_symbol_module(res.Y, cfg);
rot = qprbs13_rotate_x_module(sym.x, cfg.Dp);
fit1 = qprbs13_linear_fit_module(res.Y, rot.X, cfg);
vf = qprbs13_steady_state_vf_module(res.Y, rot.X, res.M, cfg);
noise = qprbs13_level_noise_module(res.Y, sym.sym_code, cfg);
sndr = qprbs13_sndr_module(fit1.pmax, fit1.sigmae, noise.sigma_n);
```
