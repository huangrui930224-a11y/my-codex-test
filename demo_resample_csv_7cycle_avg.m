% Demo script for resample_csv_7cycle_avg
% Update csvFile and tStableStart for your dataset.

clear; clc;

csvFile = 'example.csv';      % CSV must include columns: t, vdiff
tStableStart = 1e-6;          % seconds

cfg = struct();
cfg.ui = 1/112e9;
cfg.M = 32;
cfg.N = 8191;
cfg.numCycles = 7;
cfg.outputDir = fullfile(pwd, 'output_step1');
% cfg.dtScan = [];            % default Ts/4
% cfg.overlayUI = [1, 4096, 8191];

out = resample_csv_7cycle_avg(csvFile, tStableStart, cfg);

fprintf('Done. t0_best = %.16g s\n', out.t0_best);
fprintf('Saved outputs to: %s\n', cfg.outputDir);
fprintf('y_avg length = %d, Y size = [%d, %d]\n', numel(out.y_avg), size(out.Y,1), size(out.Y,2));
