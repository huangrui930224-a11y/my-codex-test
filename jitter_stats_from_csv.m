function out_jitter = jitter_stats_from_csv(csvPath, cfg)
%JITTER_STATS_FROM_CSV Independent callable script for JRMS/J3u analysis.
%   out_jitter = jitter_stats_from_csv(csvPath, cfg)
%
% Required inputs:
%   csvPath : string/char, path to CSV/TXT file.
%   cfg     : struct, must include fields required by
%             cru_transition_jitter_pipeline:
%             - cfg.UI
%             - cfg.outputDir
%             - cfg.context_patterns (12-class rules)
%
% Optional file-column selection (auto-detection by default):
%   cfg.t_col, cfg.v_col           : column names for time/voltage
%   cfg.t_col_idx, cfg.v_col_idx   : column indices for time/voltage
%
% Output:
%   out_jitter.JRMS
%   out_jitter.J3u
%   out_jitter.fJ
%   out_jitter.Tavgi
%   out_jitter.counts_before
%   out_jitter.counts_after
%   out_jitter.full_out            : full pipeline output

validateattributes(csvPath, {'char','string'}, {'nonempty'}, mfilename, 'csvPath', 1);
if ~isstruct(cfg)
    error('cfg must be a struct.');
end

full_out = cru_transition_jitter_pipeline(csvPath, [], cfg);

out_jitter = struct();
out_jitter.JRMS = full_out.JRMS;
out_jitter.J3u = full_out.J3u;
out_jitter.fJ = full_out.fJ;
out_jitter.Tavgi = full_out.Tavgi;
out_jitter.counts_before = full_out.counts_before;
out_jitter.counts_after = full_out.counts_after;
out_jitter.full_out = full_out;
end
