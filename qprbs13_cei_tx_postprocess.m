function out = qprbs13_cei_tx_postprocess(t, v_tx, sym_cycle, cfg)
%QPRBS13_CEI_TX_POSTPROCESS CEI/QPRBS13 TX post-processing and SNDR.
%   out = qprbs13_cei_tx_postprocess(t, v_tx, sym_cycle, cfg)
%
% This implementation follows the equations and symbols shown in the
% provided workflow excerpts (Eq. 11-13~11-22, 16-16~16-19, 27-5).
%
% Inputs
%   t         : sampled time vector (s), strictly increasing
%   v_tx      : TX analog output waveform (V), same length as t
%   sym_cycle : aligned PAM4 symbol sequence for one QPRBS13-CEI cycle,
%               length N=8191, values in {-1, -1/3, +1/3, +1}
%   cfg       : optional struct
%       .M     (default 64)
%       .T_Np  (default 29)   % for Eq. 11-16 linear-fit pulse p(k)
%       .T_Dp  (default 4)
%       .T_Dw  (default 4)
%       .T_Nw  (default 20)
%
% Output fields
%   out.step1.y, out.step1.Y, out.step1.info
%   out.levels.(Vn1, Vn13, Vp13, Vp1, Vmid, ES1, ES2, RLM)
%   out.fit29.(x, xr, X, X1, P, E, e, sigma_e, P1, p, pmax)
%   out.eq.(p_i, p_r, P2, P3, xp, w, q_i)
%   out.fit20.(P, p, vf)
%   out.sigma_n
%   out.SNDR_dB

    if nargin < 4
        cfg = struct();
    end
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
