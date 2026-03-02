function out_eoj = eoj_stats_from_csv(csvPath, cfg)
%EOJ_STATS_FROM_CSV Independent callable script for EOJ analysis.
%   out_eoj = eoj_stats_from_csv(csvPath, cfg)
%
% Required inputs:
%   csvPath : string/char, path to CSV/TXT file.
%   cfg     : struct, must include fields required by
%             cru_transition_jitter_pipeline:
%             - cfg.UI
%             - cfg.outputDir
%             - cfg.context_patterns (12-class rules)
%
% EOJ-related optional inputs:
%   cfg.Npat            : 8191 (QPRBS13_CEI) or 511 (QPRBS9_CEI)
%   cfg.eoj_mode        : 'max' | 'first3'
%   cfg.eoj_start_ui_m  : start EOJ from midpoint of m-th UI
%   cfg.t0_override     : absolute trigger anchor time (seconds)
%
% Output:
%   out_eoj.EOJ_per_class_UI
%   out_eoj.EOJ_UI
%   out_eoj.EOJ_per_class_s
%   out_eoj.EOJ_s
%   out_eoj.EOJ_num_windows
%   out_eoj.EOJ_notes
%   out_eoj.full_out

validateattributes(csvPath, {'char','string'}, {'nonempty'}, mfilename, 'csvPath', 1);
if ~isstruct(cfg)
    error('cfg must be a struct.');
end

full_out = cru_transition_jitter_pipeline(csvPath, [], cfg);

out_eoj = struct();
out_eoj.EOJ_per_class_UI = full_out.EOJ_per_class_UI;
out_eoj.EOJ_UI = full_out.EOJ_UI;
out_eoj.EOJ_per_class_s = full_out.EOJ_per_class_s;
out_eoj.EOJ_s = full_out.EOJ_s;
out_eoj.EOJ_num_windows = full_out.EOJ_num_windows;
out_eoj.EOJ_notes = full_out.EOJ_notes;
out_eoj.full_out = full_out;
end
