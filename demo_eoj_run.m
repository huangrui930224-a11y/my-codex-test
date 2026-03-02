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
