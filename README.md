# my-codex-test
chatgpt AI assistance

## Top-level auto run (只需要 waveform.csv)

你只需要提供包含 `t` 和 `v` 的 waveform CSV：

```matlab
cfg = struct('M',64, ...
             'save_mat',true, ...
             'output_mat','result.mat', ...
             'save_metrics_csv',true, ...
             'metrics_csv','metrics.csv');
out = qprbs13_cei_run_from_waveform_csv('waveform.csv', cfg);
```

执行流程：
1. 读取 waveform CSV（时间/电压）
2. 重采样
3. 后处理
4. 输出完整 `out` 结构，并可自动保存：
   - `result.mat`（完整结果）
   - `metrics.csv`（关键指标：SNDR、sigma_e、sigma_n、pmax、vf、ES1/ES2/RLM）

## 函数说明

- `resample_qprbs13_cei_step1`：重采样函数
- `qprbs13_cei_tx_postprocess`：后处理主函数（支持直接传 `waveform.csv` 自动执行读取+重采样+后处理）
- `qprbs13_cei_run_from_waveform_csv`：顶层调用函数（单入口 + 自动保存结果）

## 直接调用后处理函数（CSV 自动模式）

```matlab
out = qprbs13_cei_tx_postprocess('waveform.csv', cfg);
```
