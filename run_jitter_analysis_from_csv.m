function out = run_jitter_analysis_from_csv(csvPath, cfg, mode)
%RUN_JITTER_ANALYSIS_FROM_CSV 顶层调用接口（协议版）
% out = run_jitter_analysis_from_csv(csvPath, cfg, mode)
%
% 输入参数：
%   csvPath : .csv/.txt路径（包含时间s、电压v）
%   cfg     : 配置结构体，至少需要：
%             - cfg.UI
%             - cfg.context_patterns (12类transition规则)
%             可选：cfg.outputDir,cfg.M,cfg.balance_mode,cfg.rng_seed,
%                   cfg.edge_dir,cfg.Npat,cfg.eoj_mode,
%                   cfg.eoj_start_ui_m,cfg.t0_override,
%                   cfg.t_col/cfg.v_col 或 cfg.t_col_idx/cfg.v_col_idx
%   mode    : 'jitter' | 'eoj' | 'both'（默认'both'）

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
        error('Unsupported mode: %s. Use ''jitter'', ''eoj'', or ''both''.', mode);
end
end
