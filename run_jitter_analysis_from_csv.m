function out = run_jitter_analysis_from_csv(csvPath, cfg, mode)
%RUN_JITTER_ANALYSIS_FROM_CSV Top-level entry for independent jitter analyses.
%   out = run_jitter_analysis_from_csv(csvPath, cfg, mode)
%
% Inputs:
%   csvPath : string/char, CSV/TXT path containing time and voltage columns.
%   cfg     : config struct used by analysis scripts.
%   mode    : 'jitter' | 'eoj' | 'both' (default: 'both')
%
% Notes:
%   - JRMS/J3u analysis is executed by jitter_stats_from_csv.
%   - EOJ analysis is executed by eoj_stats_from_csv.
%   - CSV columns are auto-detected by the underlying pipeline;
%     user can override via cfg.t_col/cfg.v_col or cfg.t_col_idx/cfg.v_col_idx.

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
