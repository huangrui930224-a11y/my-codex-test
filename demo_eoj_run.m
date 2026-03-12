%DEMO_EOJ_RUN Example script for EOJ computation from CSV.
% Update csv_file to your waveform path.

clear; clc;

csv_file = 'waveform.csv'; % <-- replace with your CSV file path

cfg = struct();
cfg.fb = 112e9;
cfg.Npat = 8191;
cfg.M = 64;
cfg.fc = cfg.fb / 13280;
cfg.ref_transition = [];           % auto-select if empty
cfg.discard_first_repeats = 0;
cfg.do_plot = true;
cfg.verbose = true;
cfg.export_csv = true;
cfg.export_dir = '.';
cfg.auto_window_half_width_ui = 2; % widen auto AAAABB window to reduce missing

% Recommended: user-provided transition windows (stable and protocol-driven)
% Example template only (replace begin_ui/end_ui with your known windows).
% Must provide exactly 12 entries when enabled.
use_manual_trans_def = false;
if use_manual_trans_def
    td = repmat(struct('name','','begin_ui',1,'end_ui',1,'thr_type','th12','dir','rise'), 1, 12);
    td(1)  = struct('name','01','begin_ui',100,'end_ui',100,'thr_type','th01','dir','rise');
    td(2)  = struct('name','10','begin_ui',220,'end_ui',220,'thr_type','th01','dir','fall');
    td(3)  = struct('name','12','begin_ui',340,'end_ui',340,'thr_type','th12','dir','rise');
    td(4)  = struct('name','21','begin_ui',460,'end_ui',460,'thr_type','th12','dir','fall');
    td(5)  = struct('name','23','begin_ui',580,'end_ui',580,'thr_type','th23','dir','rise');
    td(6)  = struct('name','32','begin_ui',700,'end_ui',700,'thr_type','th23','dir','fall');
    td(7)  = struct('name','02','begin_ui',820,'end_ui',820,'thr_type','th12','dir','rise');
    td(8)  = struct('name','20','begin_ui',940,'end_ui',940,'thr_type','th12','dir','fall');
    td(9)  = struct('name','13','begin_ui',1060,'end_ui',1060,'thr_type','th23','dir','rise');
    td(10) = struct('name','31','begin_ui',1180,'end_ui',1180,'thr_type','th23','dir','fall');
    td(11) = struct('name','03','begin_ui',1300,'end_ui',1300,'thr_type','th23','dir','rise');
    td(12) = struct('name','30','begin_ui',1420,'end_ui',1420,'thr_type','th23','dir','fall');
    cfg.trans_def = td;
end

out = eoj_pam4_from_csv(csv_file, cfg);

fprintf('\n===== EOJ RESULT =====\n');
fprintf('EOJ      : %.6f ps (%.6e UI)\n', out.EOJ_ps, out.EOJ_UI);
fprintf('EOJ_i(ps):\n');
disp(out.EOJ_i_ps);
fprintf('ref transition index: %d\n', out.ref_transition);
fprintf('Tpat_est : %.6f ns\n', out.Tpat_est * 1e9);
fprintf('repeat pair used: [%d, %d]\n', out.repeat_pair_used(1), out.repeat_pair_used(2));

fprintf('symbol_mean_voltage [V0 V1 V2 V3] (V):\n');
disp(out.symbol_mean_voltage);
fprintf('crossing thresholds [th01 th02 th03 th12 th13 th23] (V):\n');
disp([out.th01, out.th02, out.th03, out.th12, out.th13, out.th23]);

fprintf('sampling phase: m_used=%d, m_opt=%d, m_center=%d, valid=%d\n', ...
    out.phase_search.m_used, out.phase_search.m_opt, out.phase_search.m_center, out.phase_search.valid);
fprintf('phase score (first 16):\n'); disp(out.phase_search.score(1:min(16,end)));

fprintf('Vui size: %d x %d\n', size(out.Vui,1), size(out.Vui,2));
fprintf('symbol phase index used: %d\n', out.symbol_sample.phase_index);
fprintf('symbol sample voltage head (first 10):\n'); disp(out.symbol_sample.voltage(1:min(10,end)).');
fprintf('inferred symbol head (first 32):\n'); disp(out.symbol_inferred(1:min(32,end)).');
fprintf('all crossing threshold times count: th01=%d, th02=%d, th03=%d, th12=%d, th13=%d, th23=%d\n', ...
    numel(out.all_crossing_threshold_times.th01), numel(out.all_crossing_threshold_times.th02), ...
    numel(out.all_crossing_threshold_times.th03), numel(out.all_crossing_threshold_times.th12), ...
    numel(out.all_crossing_threshold_times.th13), numel(out.all_crossing_threshold_times.th23));
fprintf('th01 crossing time head (s):\n'); disp(out.all_crossing_threshold_times.th01(1:min(10,end)).');

fprintf('csv exported files:\n');
disp(out.csv_export_files);
if isfield(out.csv_export_files, 'detected_transition_csv')
    fprintf('detected transition csv: %s\n', out.csv_export_files.detected_transition_csv);
end
