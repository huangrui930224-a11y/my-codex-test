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
