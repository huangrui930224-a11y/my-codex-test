function out = qprbs13_cei_tx_postprocess(in1, in2, in3, in4)
%QPRBS13_CEI_TX_POSTPROCESS CEI/QPRBS13 TX post-processing and SNDR.
%   out = qprbs13_cei_tx_postprocess(t, v_tx, sym_cycle, cfg)
%   out = qprbs13_cei_tx_postprocess(t, v_tx, [], cfg)            % auto infer sym
%   out = qprbs13_cei_tx_postprocess(wave_csv, cfg)               % CSV auto mode
%
% CSV auto mode implements: read waveform -> resample -> post-process.
%
% The implementation follows Eq. 11-13~11-22, 16-16~16-19, 27-5.

    % Parse input modes
    if ischar(in1) || isstring(in1)
        % Mode A: qprbs13_cei_tx_postprocess(wave_csv, cfg)
        wave_csv = in1;
        if nargin >= 2 && isstruct(in2)
            cfg = in2;
        else
            cfg = struct();
        end

        [t, v_tx] = local_read_waveform_csv(wave_csv);
        sym_cycle = local_infer_symbol_cycle_from_waveform(t, v_tx, cfg);
        out = run_core(t, v_tx, sym_cycle, cfg);
        out.top = struct('wave_csv', char(wave_csv));
        return;
    end

    % Mode B: qprbs13_cei_tx_postprocess(t, v_tx, sym_cycle, cfg)
    t = in1;
    v_tx = in2;

    if nargin < 3 || isempty(in3)
        sym_cycle = [];
    else
        sym_cycle = in3;
    end

    if nargin < 4 || isempty(in4)
        cfg = struct();
    else
        cfg = in4;
    end

    if isempty(sym_cycle)
        sym_cycle = local_infer_symbol_cycle_from_waveform(t, v_tx, cfg);
    end

    out = run_core(t, v_tx, sym_cycle, cfg);
end

function out = run_core(t, v_tx, sym_cycle, cfg)
    M = get_cfg(cfg, 'M', 64);
    T_Np = get_cfg(cfg, 'T_Np', 29);
    T_Dp = get_cfg(cfg, 'T_Dp', 4);
    T_Dw = get_cfg(cfg, 'T_Dw', 4);
    T_Nw = get_cfg(cfg, 'T_Nw', 20);

    % Step 1: resample waveform and build Y (Eq. 11-13 format)
    [y, Y, info] = resample_qprbs13_cei_step1(t, v_tx, M);
    y_mat_full = info.y_matrix_full;      % M x (20*N)
    m_center = info.m_center;

    N = info.N;
    N1 = info.N1;

    sym_cycle = sym_cycle(:);
    if numel(sym_cycle) ~= N
        error('sym_cycle length must be N=8191.');
    end

    % Build repeated symbol sequence for 20 cycles (aligned with y_mat_full)
    sym_full = repmat(sym_cycle, 20, 1);  % N1 x 1

    % Eq. 16-16~16-19 levels from central sample in each UI
    y_center = y_mat_full(m_center, :).';
    Vn1 = mean(y_center(sym_full == -1));
    Vn13 = mean(y_center(sym_full == -1/3));
    Vp13 = mean(y_center(sym_full ==  1/3));
    Vp1 = mean(y_center(sym_full ==  1));

    Vmid = (Vn1 + Vp1) / 2;
    ES1 = (Vn13 - Vmid) / (Vn1 - Vmid);
    ES2 = (Vp13 - Vmid) / (Vp1 - Vmid);
    RLM = min([3*ES1, 3*ES2, 2 - 3*ES1, 2 - 3*ES2]);

    % For aligned symbols x(n), use -1, -ES1, ES2, +1 (per instruction)
    x = map_symbol_to_x(sym_cycle, ES1, ES2);  % N x 1

    % Step 2: rotate x by pulse delay Dp to xr (Eq. 11-14)
    xr = rotate_left_by_D(x, T_Dp);

    % Step 2: matrix X derived from xr (Eq. 11-15)
    X = build_X_from_xr(xr);

    % Step 3: linear fit with Np=29 (Eq. 11-16, 11-17, 11-18)
    [fit29, P1_29, p29] = linear_fit_block(Y, X, T_Np);
    pmax = max(p29);

    % Step 4: remove transfer function / equalizer solve (Eq. 11-19~11-22)
    p_i = p29((0:T_Np-1)*M + m_center).';
    p_r = rotate_left_by_D(p_i, T_Dw);
    P2 = build_P2_from_pr(p_r);
    P3 = P2(1:T_Nw, :);
    xp = zeros(T_Np, 1);
    xp(T_Dp + 1) = 1;
    w = (P3.' * P3) \ (P3.' * xp);
    q_i = P3 * w;

    % Step 4 (document's next item): repeat linear-fit pulse with Np=20
    [fit20, ~, p20] = linear_fit_block(Y, X, 20);
    vf = sum(p20) / M;  % Eq. 11-12

    % Step 5: sigma_n from fixed point in runs >= 6 identical symbols
    sigma_n = calc_sigma_n_runs(y_center, sym_full, 6);

    % Step 6: SNDR (Eq. 27-5)
    sigma_e = fit29.sigma_e;
    SNDR_dB = 10 * log10((pmax^2) / (sigma_e^2 + sigma_n^2));

    out = struct();
    out.step1 = struct('y', y, 'Y', Y, 'info', info);
    out.levels = struct('Vn1', Vn1, 'Vn13', Vn13, 'Vp13', Vp13, 'Vp1', Vp1, ...
                        'Vmid', Vmid, 'ES1', ES1, 'ES2', ES2, 'RLM', RLM);
    out.fit29 = fit29;
    out.fit29.P1 = P1_29;
    out.fit29.p = p29;
    out.fit29.pmax = pmax;

    out.eq = struct('p_i', p_i, 'p_r', p_r, 'P2', P2, 'P3', P3, ...
                    'xp', xp, 'w', w, 'q_i', q_i);

    out.fit20 = fit20;
    out.fit20.p = p20;
    out.fit20.vf = vf;

    out.sigma_n = sigma_n;
    out.SNDR_dB = SNDR_dB;
end

function sym_cycle = local_infer_symbol_cycle_from_waveform(t, v_tx, cfg)
    M = get_cfg(cfg, 'M', 64);
    [~, ~, info] = resample_qprbs13_cei_step1(t, v_tx, M);

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

function v = get_cfg(cfg, name, default_v)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = default_v;
    end
end

function x = map_symbol_to_x(sym, ES1, ES2)
    x = nan(size(sym));
    x(sym == -1) = -1;
    x(sym == -1/3) = -ES1;
    x(sym ==  1/3) = ES2;
    x(sym ==  1) = 1;
    if any(isnan(x))
        error('sym_cycle contains invalid values; allowed {-1,-1/3,+1/3,+1}.');
    end
end

function xr = rotate_left_by_D(x, D)
    n = numel(x);
    D = mod(D, n);
    xr = [x(D+1:end); x(1:D)];
end

function X = build_X_from_xr(xr)
    N = numel(xr);
    X = zeros(N, N);
    X(1, :) = xr.';
    for r = 2:N
        X(r, :) = circshift(X(r-1, :), [0, 1]);
    end
end

function P2 = build_P2_from_pr(pr)
    T_Np = numel(pr);
    P2 = zeros(T_Np, T_Np);
    P2(1, :) = [pr(1); pr(end:-1:2)].';
    for r = 2:T_Np
        P2(r, :) = circshift(P2(r-1, :), [0, 1]);
    end
end

function [fit, P1, p] = linear_fit_block(Y, X, T_Np)
    N = size(X, 1);
    X1 = [X(1:T_Np, :); ones(1, N)];

    P = Y * X1.' / (X1 * X1.');
    E = P * X1 - Y;

    e = E(:);
    sigma_e = std(e);

    P1 = P(:, 1:T_Np);
    p = P1(:);

    fit = struct('x', [], 'xr', [], 'X', X, 'X1', X1, 'P', P, 'E', E, ...
                 'e', e, 'sigma_e', sigma_e);
end

function sigma_n = calc_sigma_n_runs(y_center, sym_full, min_run)
    levels = [-1; -1/3; 1/3; 1];

    d = [true; diff(sym_full) ~= 0; true];
    idx = find(d);
    starts = idx(1:end-1);
    stops = idx(2:end) - 1;

    rms_each = nan(4, 1);
    for k = 1:4
        lvl = levels(k);
        samples_lvl = [];
        for r = 1:numel(starts)
            s = starts(r);
            e = stops(r);
            if sym_full(s) == lvl && (e - s + 1) >= min_run
                samples_lvl = [samples_lvl; y_center(s:e)]; %#ok<AGROW>
            end
        end
        if isempty(samples_lvl)
            error('No run >= %d found for PAM4 level %g when calculating sigma_n.', ...
                  min_run, lvl);
        end
        mu = mean(samples_lvl);
        rms_each(k) = sqrt(mean((samples_lvl - mu).^2));
    end

    sigma_n = mean(rms_each);
end
