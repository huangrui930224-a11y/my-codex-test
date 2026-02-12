function out = qprbs13_cei_run_from_waveform_csv(wave_csv, cfg)
%QPRBS13_CEI_RUN_FROM_WAVEFORM_CSV Top-level auto runner from waveform CSV.
%   out = qprbs13_cei_run_from_waveform_csv(wave_csv)
%   out = qprbs13_cei_run_from_waveform_csv(wave_csv, cfg)
%
% Workflow:
%   1) Read transmitter waveform CSV (t, v)
%   2) Resample waveform (Step-1 function)
%   3) Infer symbols (if no symbols.csv)
%   4) Run post-processing (post-process function)
%   5) Save outputs (.mat + metrics.csv)

    if nargin < 2
        cfg = struct();
    end

    % 1) Read waveform
    [t, v_tx] = local_read_waveform_csv(wave_csv);

    % 2) Resample first (explicit)
    M = local_get_cfg(cfg, 'M', 64);
    [~, ~, info_step1] = resample_qprbs13_cei_step1(t, v_tx, M);

    % 3) Infer one-cycle symbols from waveform center samples
    sym_cycle = local_infer_symbol_cycle_from_step1(info_step1);

    % 4) Post-processing
    out = qprbs13_cei_tx_postprocess(t, v_tx, sym_cycle, cfg);

    % attach top metadata
    out.top = struct('wave_csv', char(wave_csv), ...
                     'M', M, ...
                     't_start', info_step1.t_start, ...
                     't_end', info_step1.t_end);

    % 5) Save files
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

function sym_cycle = local_infer_symbol_cycle_from_step1(info)
    N = info.N;
    y_center = info.y_matrix_full(info.m_center, :).';

    q = quantile(y_center, [0.25, 0.50, 0.75]);
    labels = ones(size(y_center));
    labels(y_center > q(1)) = 2;
    labels(y_center > q(2)) = 3;
    labels(y_center > q(3)) = 4;

    cent = zeros(4,1);
    for i = 1:4
        cent(i) = mean(y_center(labels == i));
    end

    labels2 = zeros(size(labels));
    for n = 1:numel(y_center)
        [~, idx] = min(abs(y_center(n) - cent));
        labels2(n) = idx;
    end

    [~, ord] = sort(cent, 'ascend');
    map = zeros(4,1);
    map(ord(1)) = -1;
    map(ord(2)) = -1/3;
    map(ord(3)) =  1/3;
    map(ord(4)) =  1;
    sym_est = map(labels2);

    sym_mat = reshape(sym_est, N, 20);
    sym_cycle = zeros(N,1);
    levels = [-1, -1/3, 1/3, 1];
    for k = 1:N
        row = sym_mat(k,:);
        cnt = zeros(1,4);
        for j = 1:4
            cnt(j) = sum(abs(row - levels(j)) < 1e-12);
        end
        [~, ix] = max(cnt);
        sym_cycle(k) = levels(ix);
    end
end

function [t, v_tx] = local_read_waveform_csv(csv_path)
    if ~(ischar(csv_path) || isstring(csv_path))
        error('wave_csv must be a file path.');
    end
    if ~isfile(csv_path)
        error('Waveform CSV not found: %s', csv_path);
    end

    T = readtable(csv_path);
    if width(T) < 2
        error('Waveform CSV must contain at least two columns.');
    end

    names = lower(string(T.Properties.VariableNames));
    idx_t = find(names == "t" | names == "time", 1, 'first');
    idx_v = find(names == "v_tx" | names == "v" | names == "voltage", 1, 'first');

    if isempty(idx_t) || isempty(idx_v)
        num_cols = false(1, width(T));
        for k = 1:width(T)
            num_cols(k) = isnumeric(T{:, k});
        end
        num_idx = find(num_cols);
        if numel(num_idx) < 2
            error('Waveform CSV has fewer than two numeric columns.');
        end
        idx_t = num_idx(1);
        idx_v = num_idx(2);
    end

    t = T{:, idx_t};
    v_tx = T{:, idx_v};
    t = t(:);
    v_tx = v_tx(:);
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
