% demo_run.m
% ------------------------------------------------------------
% 演示如何调用 pam4_jitter_analysis.m
% 请根据实际 CSV 路径与 fb 修改下列参数。

clear; clc;

% ===== 用户配置 =====
csv_file = 'example_waveform.csv';  % 替换为真实 CSV 文件路径
fb = 112e9;                         % 112G PAM4 可用 56e9 或按你的定义填写
M = 64;                             % 每 UI 重采样点数

% ===== 执行分析 =====
result = pam4_jitter_analysis(csv_file, fb, M);

% ===== 结果打印 =====
fprintf('\n[Demo] JRMS = %.6e s, %.3f ps, %.6e UI\n', ...
    result.JRMS_sec, result.JRMS_ps, result.JRMS_UI);
fprintf('[Demo] J3u  = %.6e s, %.3f ps, %.6e UI\n', ...
    result.J3u_sec, result.J3u_ps, result.J3u_UI);

fprintf('[Demo] 每类最终样本数:\n');
for i = 1:numel(result.transition_names)
    fprintf('  %-5s : %d, Tavgi = %.6e s\n', ...
        result.transition_names{i}, result.samples_per_class(i), result.Tavgi_sec(i));
end

% ===== 新增输出示例 =====
fprintf('\n[Demo] 每个符号均值 (V): V0=%.6e, V1=%.6e, V2=%.6e, V3=%.6e\n', ...
    result.symbol_means.V0, result.symbol_means.V1, result.symbol_means.V2, result.symbol_means.V3);
fprintf('[Demo] 跳变阈值 (V): th01=%.6e, th12=%.6e, th23=%.6e\n', ...
    result.thresholds.th01, result.thresholds.th12, result.thresholds.th23);

% 可选：导出重采样数据 CSV（time_s, vdiff_V）
resampled_csv = 'resampled_waveform.csv';
writematrix([result.resampled.time_s, result.resampled.vdiff_V], resampled_csv);
fprintf('[Demo] 已导出重采样数据: %s (N=%d)\n', resampled_csv, numel(result.resampled.vdiff_V));
fprintf('[Demo] 重采样起点 t_start_resample = %.6e s\n', result.resampled.t_start_resample);
