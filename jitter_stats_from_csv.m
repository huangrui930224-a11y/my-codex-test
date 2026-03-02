function out_jitter = jitter_stats_from_csv(csvPath, cfg)
%JITTER_STATS_FROM_CSV 协议版JRMS/J3u统计（独立可调用）
% out_jitter = jitter_stats_from_csv(csvPath, cfg)
%
% 必填参数：
%   csvPath : .csv/.txt路径（含时间s、电压v）
%   cfg.UI
%   cfg.context_patterns (12类规则)
%
% 常用可选参数：
%   cfg.outputDir, cfg.M, cfg.balance_mode, cfg.rng_seed
%   cfg.edge_dir, cfg.Npat, cfg.eoj_mode
%   cfg.t_col/cfg.v_col 或 cfg.t_col_idx/cfg.v_col_idx

full_out = protocol_pipeline_from_csv(csvPath, cfg);

out_jitter = struct();
out_jitter.JRMS = full_out.JRMS;
out_jitter.J3u = full_out.J3u;
out_jitter.fJ = full_out.fJ;
out_jitter.Si = full_out.Si;
out_jitter.Tavgi = full_out.Tavgi;
out_jitter.S0i = full_out.S0i;
out_jitter.counts_before = full_out.counts_before;
out_jitter.counts_after = full_out.counts_after;
out_jitter.full_out = full_out;
end
