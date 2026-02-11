function out = qprbs13_cei_tx_postprocess_from_csv(wave_csv, sym_src, cfg)
%QPRBS13_CEI_TX_POSTPROCESS_FROM_CSV Read CSV and run full auto processing.
%   out = qprbs13_cei_tx_postprocess_from_csv(wave_csv, sym_src, cfg)
%
% Inputs
%   wave_csv : path to waveform CSV, containing at least two numeric columns:
%              time(s) and v_tx(V).
%              Supported column names (case-insensitive):
%              - time: t, time
%              - voltage: v_tx, v, voltage
%              If names are unavailable, the first two numeric columns are used.
%   sym_src  : one of
%              (a) numeric vector of one-cycle PAM4 symbols (length 8191), or
%              (b) path to CSV file containing one symbol column
%   cfg      : optional struct passed to qprbs13_cei_tx_postprocess
%              extra optional fields:
%                .save_mat      (default false)
%                .output_mat    (default '<wave_csv>_postprocess.mat')
%
% Output
%   out      : same struct returned by qprbs13_cei_tx_postprocess

    if nargin < 3
        cfg = struct();
    end

    [t, v_tx] = read_waveform_csv(wave_csv);
    sym_cycle = read_symbol_source(sym_src);

    out = qprbs13_cei_tx_postprocess(t, v_tx, sym_cycle, cfg);

    save_mat = get_cfg(cfg, 'save_mat', false);
    if save_mat
        output_mat = get_cfg(cfg, 'output_mat', default_mat_name(wave_csv));
        save(output_mat, 'out');
    end
end

function [t, v_tx] = read_waveform_csv(csv_path)
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

    var_names = string(T.Properties.VariableNames);
    var_names_lower = lower(var_names);

    idx_t = find(var_names_lower == "t" | var_names_lower == "time", 1, 'first');
    idx_v = find(var_names_lower == "v_tx" | var_names_lower == "v" | var_names_lower == "voltage", 1, 'first');

    if isempty(idx_t) || isempty(idx_v)
        numeric_cols = false(1, width(T));
        for k = 1:width(T)
            numeric_cols(k) = isnumeric(T{:, k});
        end
        num_idx = find(numeric_cols);
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

function sym_cycle = read_symbol_source(sym_src)
    if isnumeric(sym_src)
        sym_cycle = sym_src(:);
    elseif ischar(sym_src) || isstring(sym_src)
        if ~isfile(sym_src)
            error('Symbol CSV not found: %s', sym_src);
        end
        Ts = readtable(sym_src);
        if width(Ts) < 1
            error('Symbol CSV has no columns.');
        end
        % use first numeric column, or first column if already numeric table
        idx = [];
        for k = 1:width(Ts)
            if isnumeric(Ts{:, k})
                idx = k;
                break;
            end
        end
        if isempty(idx)
            error('Symbol CSV does not contain numeric symbol column.');
        end
        sym_cycle = Ts{:, idx};
        sym_cycle = sym_cycle(:);
    else
        error('sym_src must be numeric vector or CSV path.');
    end

    if numel(sym_cycle) ~= 8191
        error('Symbol sequence length must be 8191 for one QPRBS13-CEI cycle.');
    end

    valid = ismember(sym_cycle, [-1, -1/3, 1/3, 1]);
    if ~all(valid)
        error('Symbols must be in {-1, -1/3, +1/3, +1}.');
    end
end

function v = get_cfg(cfg, name, default_v)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = default_v;
    end
end

function m = default_mat_name(wave_csv)
    [p, n, ~] = fileparts(char(wave_csv));
    m = fullfile(p, [n, '_postprocess.mat']);
end
