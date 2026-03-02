function out_eoj = eoj_stats_from_csv(csvPath, cfg)
%EOJ_STATS_FROM_CSV 协议版EOJ统计（独立可调用）
% out_eoj = eoj_stats_from_csv(csvPath, cfg)
%
% 必填参数：
%   csvPath : .csv/.txt路径（含时间s、电压v）
%   cfg.UI
%   cfg.context_patterns (12类规则)
%
% EOJ相关可选参数：
%   cfg.Npat (8191/511)
%   cfg.eoj_mode ('max'|'first3')
%   cfg.eoj_start_ui_m / cfg.t0_override
%   cfg.outputDir
%   cfg.t_col/cfg.v_col 或 cfg.t_col_idx/cfg.v_col_idx

full_out = protocol_pipeline_from_csv(csvPath, cfg);

out_eoj = struct();
out_eoj.EOJ_per_class_UI = full_out.EOJ_per_class_UI;
out_eoj.EOJ_UI = full_out.EOJ_UI;
out_eoj.EOJ_per_class_s = full_out.EOJ_per_class_s;
out_eoj.EOJ_s = full_out.EOJ_s;
out_eoj.EOJ_num_windows = full_out.EOJ_num_windows;
out_eoj.EOJ_notes = full_out.EOJ_notes;
out_eoj.full_out = full_out;
end
