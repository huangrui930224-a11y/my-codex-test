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
2. 做重采样（`resample_qprbs13_cei_step1`）
3. 做后处理（`qprbs13_cei_tx_postprocess`）
4. 输出完整 `out` 结构，并可自动保存：
   - `result.mat`（完整结果）
   - `metrics.csv`（关键指标：SNDR、sigma_e、sigma_n、pmax、vf、ES1/ES2/RLM）

## 组件函数

- 重采样函数：`resample_qprbs13_cei_step1`
- 后处理函数：`qprbs13_cei_tx_postprocess`
- CSV入口函数：`qprbs13_cei_tx_postprocess_from_csv`
- 顶层自动调用：`qprbs13_cei_run_from_waveform_csv`
