# my-codex-test

MATLAB QPRBS13/PAM4 分析工具集合。当前仓库提供从 CSV 波形到符号推断、线性拟合、噪声与 SNDR 计算的模块化 pipeline。

## 文件说明（逐脚本）

### 顶层与文档
- `qprbs13_pipeline.m`
  - 顶层总控入口。
  - 串联调用各子模块：重采样、符号反推、旋转建模、线性拟合、稳态 `v_f`、电平噪声、SNDR。
  - 输出 CSV/MAT/PNG 结果。
- `DEMO.md`
  - 中文调用文档，包含最小运行示例、参数说明、输出文件清单、常见报错。
- `README.md`
  - 本说明文件。

### 重采样相关
- `qprbs13_resample_module.m`
  - 从 CSV 读取波形并执行 `t0` 扫描、7 周期重采样、`y_avg` 与 `Y` 生成。
  - 支持两种 CSV 格式：
    1) 有表头：`t, vdiff`
    2) 无表头：默认第 1 列时间、第 2 列电压。
- `resample_csv_7cycle_avg.m`
  - 独立版重采样工具（不依赖完整 pipeline）。
- `demo_resample_csv_7cycle_avg.m`
  - 重采样工具的最小 demo 脚本。
- `resample_qprbs13_cei_step1.m`
  - 另一个早期/参考重采样函数（保留）。

### 符号与建模分析模块
- `qprbs13_symbol_module.m`
  - 从 `Y` 提取 `s(n)`，聚类得到 `sym_code`，并计算 `ES1/ES2/RLM`、生成 `x(n)`。
- `qprbs13_rotate_x_module.m`
  - 对 `x(n)` 按 `Dp` 旋转得到 `xr`，并构造 `X`（循环结构矩阵）。
- `qprbs13_linear_fit_module.m`
  - 线性拟合求 `P/E`，提取 `e(k)`、`p(k)`，输出 `sigma_e` 与 `pmax`。
- `qprbs13_steady_state_vf_module.m`
  - 用 `NpVf=20` 的重算脉冲估计稳态电压 `v_f`。
- `qprbs13_level_noise_module.m`
  - 在“连续相同符号长度>=minRunLen”区段上按固定相位估计四电平 RMS 噪声与 `sigma_n`。
- `qprbs13_sndr_module.m`
  - 按公式 `SNDR = 10*log10(pmax^2/(sigma_e^2+sigma_n^2))` 计算 SNDR。

## 快速开始

推荐先看 `DEMO.md`，复制其中的最小示例直接运行。

