function out = run_jitter_analysis_from_csv(csvPath, cfg, mode)
%RUN_JITTER_ANALYSIS_FROM_CSV 顶层调用接口
% out = run_jitter_analysis_from_csv(csvPath, cfg, mode)
%
% 输入参数说明：
%   csvPath : CSV文件路径（包含时间s与电压v）
%   cfg     : 配置结构体（至少包含 cfg.UI）
%             - 用于JRMS/J3u：cfg.UI
%             - 用于EOJ：cfg.UI, cfg.Npat(可选,默认8191)
%             - 可选列覆盖：cfg.t_col/cfg.v_col 或 cfg.t_col_idx/cfg.v_col_idx
%   mode    : 'jitter' | 'eoj' | 'both'（默认'both'）
%
% 调用行为：
%   - mode='jitter'：仅调用 jitter_stats_from_csv
%   - mode='eoj'   ：仅调用 eoj_stats_from_csv
%   - mode='both'  ：两者都调用并返回 out.jitter/out.eoj

if nargin < 3 || isempty(mode)
    mode = 'both';
end

switch lower(mode)
    case 'jitter'
        out = jitter_stats_from_csv(csvPath, cfg);
    case 'eoj'
        out = eoj_stats_from_csv(csvPath, cfg);
    case 'both'
        out = struct();
        out.jitter = jitter_stats_from_csv(csvPath, cfg);
        out.eoj = eoj_stats_from_csv(csvPath, cfg);
    otherwise
        error('Unsupported mode: %s', mode);
end
end
