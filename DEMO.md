# CRU Transition Jitter Pipeline Demo 文档

本文档提供可直接执行的 Demo 测试流程，用于验证 `cru_transition_jitter_pipeline.m` 在两种典型场景下可运行：

1. **向量输入模式**（直接传 `t`/`v`）
2. **文件路径模式**（直接传 `.csv/.txt`）

---

## 1) 准备环境

在 MATLAB 中将工作目录切换到本仓库路径，确保 `cru_transition_jitter_pipeline.m` 可见。

```matlab
cd('/workspace/my-codex-test')
```

---

## 2) Demo A：向量输入模式（基础冒烟测试）

### 2.1 目标

- 验证主流程可执行：重采样、crossing 提取、CRU、12类统计、EOJ 输出。
- 生成基础输出文件：
  - `cru_debug.png`
  - `histogram_fJ.png`
  - `counts_before_after.png`
  - `eoj_per_class.csv`
  - `eoj_bar.png`

### 2.2 示例代码

```matlab
clear; clc;

% --- 合成数据 ---
UI = 100e-12;
fs = 1e12;                 % 1 ps 采样
Tend = 200e-9;
t = (0:1/fs:Tend).';

% 构造简化 PAM4 风格波形（仅用于demo）
sig = 0.5*square(2*pi*(1/(2*UI))*t);
noise = 0.01*randn(size(t));
v = sig + noise;

% --- 配置 ---
cfg = struct();
cfg.UI = UI;
cfg.outputDir = './demo_output_A';
cfg.context_patterns = {
    [0 1], [1 0], [1 2], [2 1], [2 3], [3 2], ...
    [0 2], [2 0], [1 3], [3 1], [0 3], [3 0] ...
};

cfg.M = 32;
cfg.balance_mode = 'truncate';
cfg.rng_seed = 1;
cfg.edge_dir = 'either';
cfg.th_method = 'midpoint';
cfg.debug_plot = true;
cfg.Npat = 8191;
cfg.eoj_mode = 'max';
cfg.eoj_debug_plot = false;

% 可选：从第 m 个 UI 中点开始 EOJ 统计
cfg.eoj_start_ui_m = 0;

out = cru_transition_jitter_pipeline(t, v, cfg);

disp(out.JRMS);
disp(out.J3u);
disp(out.EOJ_UI);
```

### 2.3 验收点

- `out` 成功返回，且字段如 `JRMS/J3u/EOJ_UI` 非空。
- `demo_output_A` 下存在 PNG 与 CSV 文件。

---

## 3) Demo B：文件路径模式（CSV自动读取测试）

### 3.1 目标

验证文件输入路径模式：

```matlab
out = cru_transition_jitter_pipeline(filePath, [], cfg)
```

并测试自动列识别、按列名/列索引显式指定。

### 3.2 准备 CSV

先将 Demo A 的数据导出为 CSV：

```matlab
Tbl = table(t, v, 'VariableNames', {'time_s','voltage_v'});
writetable(Tbl, 'demo_wave.csv');
```

### 3.3 示例代码：自动列名识别

```matlab
clear cfg;
cfg = struct();
cfg.UI = UI;
cfg.outputDir = './demo_output_B_auto';
cfg.context_patterns = {
    [0 1], [1 0], [1 2], [2 1], [2 3], [3 2], ...
    [0 2], [2 0], [1 3], [3 1], [0 3], [3 0] ...
};

out_auto = cru_transition_jitter_pipeline('demo_wave.csv', [], cfg);
```

### 3.4 示例代码：显式按列名指定

```matlab
cfg.t_col = 'time_s';
cfg.v_col = 'voltage_v';
cfg.outputDir = './demo_output_B_name';
out_name = cru_transition_jitter_pipeline('demo_wave.csv', [], cfg);
```

### 3.5 示例代码：显式按列索引指定

```matlab
cfg = rmfield(cfg, {'t_col','v_col'});
cfg.t_col_idx = 1;
cfg.v_col_idx = 2;
cfg.outputDir = './demo_output_B_idx';
out_idx = cru_transition_jitter_pipeline('demo_wave.csv', [], cfg);
```

### 3.6 验收点

- 三种路径模式都能成功输出 `out_*`。
- 各输出目录均生成 EOJ CSV 与图像。

---

## 4) 常见问题

1. **某些 class 样本为 0 导致报错**
   - 可延长数据长度、检查 `context_patterns`、或确认信号包含足够 transition。

2. **EOJ 窗口不足**
   - 减小 `cfg.eoj_start_ui_m` 或增加采样时长（repeat 数）。

3. **CSV 列无法自动识别**
   - 显式设置 `cfg.t_col/cfg.v_col` 或 `cfg.t_col_idx/cfg.v_col_idx`。

4. **输出图像未生成**
   - 检查 `cfg.outputDir` 写权限，以及 `cfg.debug_plot` 是否开启。

---

## 5) 推荐最小回归清单

每次修改后建议至少执行：

- Demo A（向量输入）
- Demo B（CSV路径输入 + 自动列识别）
- Demo B（CSV路径输入 + 显式列名）

并确认：

- 函数无报错
- `EOJ_UI` 可计算
- 输出文件齐全


---

## 6) 新增独立脚本测试（按统计项拆分）

### 6.1 仅运行 JRMS/J3u

```matlab
cfg = struct();
cfg.UI = 100e-12;
cfg.outputDir = './demo_output_split_jitter';
cfg.context_patterns = {
    [0 1], [1 0], [1 2], [2 1], [2 3], [3 2], ...
    [0 2], [2 0], [1 3], [3 1], [0 3], [3 0] ...
};
out_jitter = jitter_stats_from_csv('demo_wave.csv', cfg);
disp(out_jitter.JRMS);
disp(out_jitter.J3u);
```

### 6.2 仅运行 EOJ

```matlab
cfg.outputDir = './demo_output_split_eoj';
cfg.Npat = 8191;
out_eoj = eoj_stats_from_csv('demo_wave.csv', cfg);
disp(out_eoj.EOJ_UI);
```

### 6.3 顶层接口统一调用

```matlab
cfg.outputDir = './demo_output_split_both';
out_all = run_jitter_analysis_from_csv('demo_wave.csv', cfg, 'both');
```

### 6.4 验收点

- `jitter_stats_from_csv` 可独立输出 `JRMS/J3u`
- `eoj_stats_from_csv` 可独立输出 `EOJ_UI`
- `run_jitter_analysis_from_csv(..., 'both')` 可同时返回两套结果
