function out = qprbs13_cei_run_from_waveform_csv(wave_csv, cfg)
%QPRBS13_CEI_RUN_FROM_WAVEFORM_CSV Top-level auto runner from waveform CSV.
%   out = qprbs13_cei_run_from_waveform_csv(wave_csv)
%   out = qprbs13_cei_run_from_waveform_csv(wave_csv, cfg)
%
% Workflow:
%   1) Read waveform CSV (t, v)
%   2) Resample
%   3) Post-process
%   4) Save outputs
%
% Notes
%   The read/resample/post-process core is implemented inside
%   qprbs13_cei_tx_postprocess in CSV-auto mode.

    if nargin < 2
        cfg = struct();
    end

    out = qprbs13_cei_tx_postprocess(wave_csv, cfg);

    save_mat = local_get_cfg(cfg, 'save_mat', true);
    if save_mat
        output_mat = local_get_cfg(cfg, 'output_mat', local_default_mat_name(wave_csv));
        save(output_mat, 'out');
    end

    save_metrics_csv = local_get_cfg(cfg, 'save_metrics_csv', true);
    if save_metrics_csv
        metrics_csv = local_get_cfg(cfg, 'metrics_csv', local_default_metrics_name(wave_csv));
        local_write_metrics_csv(metrics_csv, out);
    end
end

function local_write_metrics_csv(csv_path, out)
    metrics = table();
    metrics.SNDR_dB = out.SNDR_dB;
    metrics.sigma_e = out.fit29.sigma_e;
    metrics.sigma_n = out.sigma_n;
    metrics.pmax = out.fit29.pmax;
    metrics.vf = out.fit20.vf;
    metrics.ES1 = out.levels.ES1;
    metrics.ES2 = out.levels.ES2;
    metrics.RLM = out.levels.RLM;
    writetable(metrics, csv_path);
end

function v = local_get_cfg(cfg, name, default_v)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = default_v;
    end
end

function m = local_default_mat_name(wave_csv)
    [p, n, ~] = fileparts(char(wave_csv));
    m = fullfile(p, [n, '_result.mat']);
end

function m = local_default_metrics_name(wave_csv)
    [p, n, ~] = fileparts(char(wave_csv));
    m = fullfile(p, [n, '_metrics.csv']);
end
