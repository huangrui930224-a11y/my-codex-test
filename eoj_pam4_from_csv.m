function out = eoj_pam4_from_csv(csv_file, cfg)
%EOJ_PAM4_FROM_CSV Compute EOJ (Even/Odd Jitter) for PAM4 from CSV waveform.
%   out = eoj_pam4_from_csv(csv_file, cfg)
%
% This implementation follows the protocol definition in the task:
% - Trigger/reference clock is from a CRU (golden PLL, first-order LPF,
%   fc = fb/13280, 20 dB/dec).
% - EOJ is computed ONLY from CRU-corrected crossing times, independent of
%   JRMS/J3u.
%
% Why CRU-corrected time axis:
%   EOJ compares repeat-to-repeat pattern timing under recovered-clock
%   reference. Using raw absolute time would include wander already tracked
%   by the CRU and violate the intended measurement reference.
%
% Why PLL must update every UI:
%   Crossings occur sparsely (only at transition UIs). If we update PLL only
%   on crossing events, effective update rate drops and equivalent bandwidth
%   becomes wrong. Therefore we update once per UI, and when no event exists
%   we hold previous estimate (sample-and-hold behavior).

    cfg = apply_defaults(cfg);

    [t, v] = read_csv_two_cols(csv_file);

    UI = 1 / cfg.fb;
    Ts = UI / cfg.M;

    % Step 1: uniform time grid and resampling
    t_uniform = (t(1):Ts:t(end)).';
    v_uniform = interp1(t, v, t_uniform, 'linear');

    N_UI = floor(numel(v_uniform) / cfg.M);
    if N_UI < (cfg.discard_first_repeats + 2) * cfg.Npat
        error('Data does not cover required repeats (need at least repeat0 and repeat1 after discard).');
    end

    v_uniform = v_uniform(1:N_UI * cfg.M);
    t_uniform = t_uniform(1:N_UI * cfg.M);

    % Step 2: UI block + sampling phase optimization
    % Search mm=1..M and choose best phase by separation score:
    % score = min(level spacing) / (mean(within-cluster std) + eps)
    % If no valid phase score, fallback to center sample.
    Vui = reshape(v_uniform, cfg.M, N_UI).';
    m_center = round(cfg.M / 2);
    [m_opt, phase_score, phase_valid] = find_optimal_sampling_phase(Vui, m_center);
    if phase_valid
        m = m_opt;
    else
        m = m_center;
        warning('Phase optimization failed, fallback to center sample m=%d.', m_center);
    end
    y = Vui(:, m);

    % Step 3: infer symbols + 4-level means
    [sym, centers_sorted] = simple_kmeans_1d(y, 4, 100);
    V0 = mean(y(sym == 0));
    V1 = mean(y(sym == 1));
    V2 = mean(y(sym == 2));
    V3 = mean(y(sym == 3));

    % Step 4: thresholds (Method A)
    th01 = (V0 + V1) / 2;
    th12 = (V1 + V2) / 2;
    th23 = (V2 + V3) / 2;

    % Step 5: repeat grouping by nominal Npat*UI
    t0 = t_uniform(1);

    % Build or infer transition definitions (12 classes)
    if isfield(cfg, 'trans_def') && ~isempty(cfg.trans_def)
        trans_def = cfg.trans_def;
    else
        if cfg.verbose
            fprintf('[EOJ] cfg.trans_def not provided, infer transition windows from symbol stream with AAAABB assumption.\n');
        end
        trans_def = infer_transitions_from_symbols(sym, cfg.Npat);
    end

    if numel(trans_def) ~= 12
        error('Need exactly 12 transition classes.');
    end

    % Step 6) Extract ALL crossing events first (for CRU input), then
    % extract 12 target transition crossings for EOJ statistics.
    %
    % This follows the updated requirement:
    %   1) CRU uses all available crossing events.
    %   2) After CRU timeline is obtained, pick AAAABB 12-class transitions.

    use_repeat_start = cfg.discard_first_repeats;
    ui_start_for_use = use_repeat_start * cfg.Npat + 1;

    all_evt = collect_all_crossing_events(t_uniform, v_uniform, cfg.M, th01, th12, th23, ui_start_for_use);

    % Step 6 (target transitions): crossing extraction per class/repeat
    trans = repmat(struct('tcross_abs', [], 'ui_index_global', [], 'repeat_id', [], ...
        'name', '', 'thr_type', '', 'dir', ''), 1, 12);

    max_repeat = floor(N_UI / cfg.Npat) - 1;

    for i = 1:12
        trans(i).name = trans_def(i).name;
        trans(i).thr_type = trans_def(i).thr_type;
        trans(i).dir = trans_def(i).dir;
        th = get_threshold(trans_def(i).thr_type, th01, th12, th23);

        for r = use_repeat_start:max_repeat
            begin_ui = trans_def(i).begin_ui;
            end_ui = trans_def(i).end_ui;
            if begin_ui < 1 || end_ui > cfg.Npat || end_ui < begin_ui
                error('Invalid begin_ui/end_ui in trans_def(%d).', i);
            end

            ui_global_begin = r * cfg.Npat + begin_ui;
            ui_global_end = r * cfg.Npat + end_ui;
            if ui_global_end > N_UI
                continue;
            end

            k0 = (ui_global_begin - 1) * cfg.M + 1;
            k1 = ui_global_end * cfg.M;
            [tcross, kcross, n_cross] = find_first_crossing(t_uniform, v_uniform, k0, k1, th, trans_def(i).dir);

            if n_cross > 1 && cfg.verbose
                fprintf('[EOJ] transition %d (%s), repeat %d: multiple crossings (%d), using first.\n', ...
                    i, trans_def(i).name, r, n_cross);
            end

            if ~isnan(tcross)
                ui_cross = floor((kcross - 1) / cfg.M) + 1; % 1-based UI index
                trans(i).tcross_abs(end+1,1) = tcross; %#ok<AGROW>
                trans(i).ui_index_global(end+1,1) = ui_cross;
                trans(i).repeat_id(end+1,1) = r;
            end
        end
    end

    % Step 7: CRU golden PLL (per-UI update) using ALL crossings.
    fc = cfg.fc;
    alpha = exp(-2 * pi * fc / cfg.fb);

    tie_LF_ui = zeros(N_UI, 1);
    has_event = false(N_UI, 1);
    event_value = nan(N_UI, 1);

    if ~isempty(all_evt.tcross_abs)
        tie_sum = zeros(N_UI, 1);
        tie_cnt = zeros(N_UI, 1);
        for k = 1:numel(all_evt.tcross_abs)
            n_ui = all_evt.ui_index_global(k);
            if n_ui < 1 || n_ui > N_UI
                continue;
            end
            n0 = n_ui - 1;
            t_ideal = t0 + n0 * UI;
            tie_raw = all_evt.tcross_abs(k) - t_ideal;
            tie_sum(n_ui) = tie_sum(n_ui) + tie_raw;
            tie_cnt(n_ui) = tie_cnt(n_ui) + 1;
        end
        idx_evt = find(tie_cnt > 0);
        has_event(idx_evt) = true;
        event_value(idx_evt) = tie_sum(idx_evt) ./ tie_cnt(idx_evt);
    end

    for n = 2:N_UI
        if has_event(n)
            e = event_value(n);
        else
            e = tie_LF_ui(n-1);
        end
        tie_LF_ui(n) = alpha * tie_LF_ui(n-1) + (1 - alpha) * e;
    end

    for i = 1:12
        trans(i).tcross_cru_abs = nan(size(trans(i).tcross_abs));
        for k = 1:numel(trans(i).tcross_abs)
            n_ui = trans(i).ui_index_global(k);
            trans(i).tcross_cru_abs(k) = trans(i).tcross_abs(k) - tie_LF_ui(n_ui);
        end
    end
    % Step 8: reference transition and Tpat_est
    [ref_i, rep0, rep1] = choose_ref_transition_and_repeats(trans, cfg.ref_transition, use_repeat_start, cfg.verbose);

    T3 = mean_or_nan(trans(ref_i).tcross_cru_abs(trans(ref_i).repeat_id == rep0));
    T4 = mean_or_nan(trans(ref_i).tcross_cru_abs(trans(ref_i).repeat_id == rep1));

    if isnan(T3) || isnan(T4)
        error('Reference transition has no valid data in selected adjacent repeats.');
    end
    Tpat_est = T4 - T3;

    % Step 9: EOJ_i
    T1 = nan(1, 12);
    T2 = nan(1, 12);
    EOJ_i = nan(1, 12);
    for i = 1:12
        T1(i) = mean_or_nan(trans(i).tcross_cru_abs(trans(i).repeat_id == rep0));
        T2(i) = mean_or_nan(trans(i).tcross_cru_abs(trans(i).repeat_id == rep1));
        if isnan(T1(i)) || isnan(T2(i))
            warning('Transition %d (%s) missing in repeat0/1. EOJ_i set to NaN.', i, trans(i).name);
            EOJ_i(i) = nan;
        else
            EOJ_i(i) = abs((T2(i) - T1(i)) - Tpat_est);
        end
    end

    % Step 10: output
    out = struct();
    out.EOJ_i_sec = EOJ_i;
    out.EOJ_i_ps = EOJ_i * 1e12;

    if all(isnan(EOJ_i))
        out.EOJ_sec = nan;
        out.EOJ_ps = nan;
        out.EOJ_UI = nan;
        out.valid = false;
    else
        out.EOJ_sec = max(EOJ_i, [], 'omitnan');
        out.EOJ_ps = out.EOJ_sec * 1e12;
        out.EOJ_UI = out.EOJ_sec / UI;
        out.valid = true;
    end

    out.ref_transition = ref_i;
    out.T1 = T1;
    out.T2 = T2;
    out.T3 = T3;
    out.T4 = T4;
    out.Tpat_est = Tpat_est;
    out.repeat_pair_used = [rep0, rep1];

    out.V0 = V0; out.V1 = V1; out.V2 = V2; out.V3 = V3;
    out.symbol_mean_voltage = [V0, V1, V2, V3];
    out.centers_sorted = centers_sorted;
    out.th01 = th01; out.th12 = th12; out.th23 = th23;
    out.crossing_threshold = struct('th01', th01, 'th12', th12, 'th23', th23);

    out.trans_def = trans_def;
    out.trans = trans;
    out.all_crossings = all_evt;
    out.all_crossings_count = numel(all_evt.tcross_abs);
    out.UI = UI;
    out.t0 = t0;
    out.alpha = alpha;
    out.phase_search = struct('m_opt', m_opt, 'm_center', m_center, 'm_used', m, ...
        'valid', phase_valid, 'score', phase_score);

    % Step 11: plots
    if cfg.do_plot
        make_plots(y, [V0 V1 V2 V3], [th01 th12 th23], out, trans, rep0, rep1);
    end

    if cfg.verbose
        fprintf('[EOJ] EOJ = %.6f ps (%.6f UI), ref transition = %d (%s), Tpat_est = %.3f ns\n', ...
            out.EOJ_ps, out.EOJ_UI, ref_i, trans(ref_i).name, out.Tpat_est * 1e9);
        fprintf('[EOJ] Symbol mean voltages [V0 V1 V2 V3] = [%.6g %.6g %.6g %.6g] V\n', ...
            out.symbol_mean_voltage(1), out.symbol_mean_voltage(2), out.symbol_mean_voltage(3), out.symbol_mean_voltage(4));
        fprintf('[EOJ] Crossing thresholds [th01 th12 th23] = [%.6g %.6g %.6g] V\n', ...
            out.th01, out.th12, out.th23);
        fprintf('[EOJ] Sampling phase m_used=%d (m_opt=%d, m_center=%d, phase_valid=%d)\n', ...
            out.phase_search.m_used, out.phase_search.m_opt, out.phase_search.m_center, out.phase_search.valid);
    end
end

function cfg = apply_defaults(cfg)
    if ~isfield(cfg, 'fb') || isempty(cfg.fb), cfg.fb = 112e9; end
    if ~isfield(cfg, 'Npat') || isempty(cfg.Npat), cfg.Npat = 8191; end
    if ~isfield(cfg, 'M') || isempty(cfg.M), cfg.M = 64; end
    if ~isfield(cfg, 'fc') || isempty(cfg.fc), cfg.fc = cfg.fb / 13280; end
    if ~isfield(cfg, 'ref_transition'), cfg.ref_transition = []; end
    if ~isfield(cfg, 'discard_first_repeats') || isempty(cfg.discard_first_repeats)
        cfg.discard_first_repeats = 0;
    end
    if ~isfield(cfg, 'do_plot') || isempty(cfg.do_plot), cfg.do_plot = true; end
    if ~isfield(cfg, 'verbose') || isempty(cfg.verbose), cfg.verbose = true; end
end

function [t, v] = read_csv_two_cols(csv_file)
    x = readmatrix(csv_file);
    if size(x, 2) < 2
        error('CSV must contain at least 2 columns: time, vdiff.');
    end
    t = x(:,1);
    v = x(:,2);
    t = t(:); v = v(:);
    valid = ~(isnan(t) | isnan(v));
    t = t(valid); v = v(valid);
    if numel(t) < 2
        error('Not enough valid samples after NaN filtering.');
    end
    if any(diff(t) <= 0)
        error('time column must be strictly increasing.');
    end
end

function [m_opt, score_vec, valid] = find_optimal_sampling_phase(Vui, m_center)
    [~, M] = size(Vui);
    score_vec = nan(1, M);
    eps0 = 1e-12;

    for mm = 1:M
        x = Vui(:, mm);
        [labels, centers] = kmeans_1d_no_toolbox(x, 4, 60);
        if isempty(labels) || numel(unique(labels)) < 4 || any(isnan(centers))
            continue;
        end

        centers = sort(centers(:).');
        spacing = diff(centers);
        if any(~isfinite(spacing)) || isempty(spacing)
            continue;
        end
        min_spacing = min(spacing);

        wstd = nan(1, 4);
        for kk = 0:3
            xv = x(labels == kk);
            if numel(xv) < 2
                wstd(kk+1) = nan;
            else
                wstd(kk+1) = std(xv);
            end
        end
        mean_wstd = mean(wstd, 'omitnan');
        if ~isfinite(mean_wstd)
            continue;
        end

        score_vec(mm) = min_spacing / (mean_wstd + eps0);
    end

    valid_idx = find(isfinite(score_vec));
    valid = ~isempty(valid_idx);
    if valid
        [~, irel] = max(score_vec(valid_idx));
        m_opt = valid_idx(irel);
    else
        m_opt = m_center;
    end
end

function [labels, centers] = kmeans_1d_no_toolbox(x, K, max_iter)
    [labels, centers] = simple_kmeans_1d(x, K, max_iter);
end

function [labels, centers] = simple_kmeans_1d(x, K, max_iter)
    % Toolbox-free 1D kmeans
    x = x(:);
    q = linspace(0.05, 0.95, K);
    centers = quantile(x, q);
    labels = zeros(size(x));
    for it = 1:max_iter
        old_labels = labels;
        d = abs(x - centers);
        [~, labels] = min(d, [], 2);
        for k = 1:K
            if any(labels == k)
                centers(k) = mean(x(labels == k));
            end
        end
        if isequal(labels, old_labels)
            break;
        end
    end

    [centers_sorted, order] = sort(centers, 'ascend');
    remap = zeros(1, K);
    for k = 1:K
        remap(order(k)) = k - 1; % output labels in 0..3
    end
    labels = remap(labels).';
    labels = labels(:);
    centers = centers_sorted(:).';
end

function trans_def = infer_transitions_from_symbols(sym, Npat)
    % Automatic fallback when cfg.trans_def is missing.
    %
    % IMPORTANT AAAABB classification rule used here:
    % For a transition A->B at boundary n->n+1, we only accept it as a
    % candidate class event when local context is exactly:
    %   sym(n-3:n)   = [A A A A]
    %   sym(n+1:n+2) = [B B]
    % (i.e., AAAABB with the boundary between 4th A and 1st B).
    %
    % This corresponds to "detect AAAABB then classify current transition"
    % and avoids classifying isolated/noisy transitions without context.

    sym = sym(:);
    if numel(sym) < Npat + 2
        error('Not enough UI symbols to infer AAAABB transitions.');
    end
    % Use first full pattern for robust/specific index extraction.
    s = sym(1:Npat);

    % Collect all AAAABB-qualified boundaries in one repeat.
    % boundary UI index n means transition from UI n to n+1.
    class_pos = cell(4, 4);
    for n = 4:(Npat - 2)
        A = s(n);
        B = s(n + 1);
        if A == B
            continue;
        end
        left_ok = all(s(n-3:n) == A);
        right_ok = all(s(n+1:n+2) == B);
        if left_ok && right_ok
            class_pos{A+1, B+1}(end+1) = n; %#ok<AGROW>
        end
    end

    % Build all 12 directed classes explicitly, each class picks one
    % representative window (first AAAABB hit in repeat0).
    trans_def = repmat(struct('name', '', 'begin_ui', 1, 'end_ui', 1, ...
        'thr_type', 'th12', 'dir', 'rise'), 1, 12);

    idx = 0;
    missing = {};
    for a = 0:3
        for b = 0:3
            if a == b
                continue;
            end
            idx = idx + 1;
            nm = sprintf('%d%d', a, b);
            trans_def(idx).name = nm;
            trans_def(idx).dir = ternary(b > a, 'rise', 'fall');
            trans_def(idx).thr_type = choose_thr_type_by_levels(a, b);

            pos = class_pos{a+1, b+1};
            if isempty(pos)
                missing{end+1} = nm; %#ok<AGROW>
                % Fallback to broad repeat window; later crossing stage may
                % still fail and return NaN, but class remains defined.
                trans_def(idx).begin_ui = 1;
                trans_def(idx).end_ui = Npat;
            else
                trans_def(idx).begin_ui = pos(1);
                trans_def(idx).end_ui = pos(1);
            end
        end
    end

    if ~isempty(missing)
        warning('AAAABB auto-infer missing classes in repeat0: %s', strjoin(missing, ', '));
    end
end

function thr_type = choose_thr_type_by_levels(a, b)
    hi = max(a, b);
    lo = min(a, b);
    if lo == 0 && hi == 1
        thr_type = 'th01';
    elseif lo == 1 && hi == 2
        thr_type = 'th12';
    elseif lo == 2 && hi == 3
        thr_type = 'th23';
    else
        % For 2-level jump, pick threshold nearest to midpoint.
        mid = (a + b) / 2;
        if mid < 1
            thr_type = 'th01';
        elseif mid < 2
            thr_type = 'th12';
        else
            thr_type = 'th23';
        end
    end
end


function all_evt = collect_all_crossing_events(t, v, M, th01, th12, th23, ui_start)
    thresholds = [th01, th12, th23];
    tc = [];
    ui = [];

    for k = 1:(numel(v)-1)
        ui_k = floor((k - 1) / M) + 1;
        if ui_k < ui_start
            continue;
        end
        v0 = v(k);
        v1 = v(k+1);
        dv = v1 - v0;
        if dv == 0
            continue;
        end
        for ith = 1:3
            th = thresholds(ith);
            is_cross = (v0 < th && v1 >= th) || (v0 > th && v1 <= th);
            if is_cross
                tcross = t(k) + (th - v0) * (t(k+1) - t(k)) / dv;
                tc(end+1,1) = tcross; %#ok<AGROW>
                ui(end+1,1) = ui_k; %#ok<AGROW>
            end
        end
    end

    all_evt = struct('tcross_abs', tc, 'ui_index_global', ui);
end

function th = get_threshold(thr_type, th01, th12, th23)
    switch lower(thr_type)
        case 'th01'
            th = th01;
        case 'th12'
            th = th12;
        case 'th23'
            th = th23;
        otherwise
            error('Unknown thr_type: %s', thr_type);
    end
end

function [tcross, kcross, n_cross] = find_first_crossing(t, v, k0, k1, th, dir)
    tcross = nan; kcross = nan; n_cross = 0;
    if k1 <= k0 || k1 > numel(v)
        return;
    end
    vv0 = v(k0:k1-1);
    vv1 = v(k0+1:k1);

    switch lower(dir)
        case 'rise'
            idx = find((vv0 < th) & (vv1 >= th));
        case 'fall'
            idx = find((vv0 > th) & (vv1 <= th));
        otherwise
            error('Unknown direction: %s', dir);
    end

    n_cross = numel(idx);
    if isempty(idx)
        return;
    end

    k = k0 + idx(1) - 1;
    dv = v(k+1) - v(k);
    if dv == 0
        return;
    end
    tcross = t(k) + (th - v(k)) * (t(k+1) - t(k)) / dv;
    kcross = k;
end

function [ref_i, rep0, rep1] = choose_ref_transition_and_repeats(trans, ref_cfg, rep_start, verbose)
    % First try protocol-intended pair: repeat0/repeat1 after discard.
    rep0 = rep_start;
    rep1 = rep_start + 1;

    if ~isempty(ref_cfg)
        ref_i = ref_cfg;
        if any(trans(ref_i).repeat_id == rep0) && any(trans(ref_i).repeat_id == rep1)
            return;
        end

        % If user-specified ref missing on (rep0,rep1), keep ref_i but search
        % first adjacent repeat pair where this reference exists.
        [rep0, rep1, ok] = find_adjacent_pair_for_transition(trans(ref_i).repeat_id, rep_start);
        if ok
            if verbose
                warning('ref_transition=%d missing in repeat%d/repeat%d; fallback to repeat%d/repeat%d.', ...
                    ref_i, rep_start, rep_start+1, rep0, rep1);
            end
            return;
        end
        error('Specified ref_transition=%d has no valid adjacent repeat pair.', ref_i);
    end

    % Auto-ref selection: best transition on default pair first.
    best = -1;
    ref_i = 1;
    for i = 1:numel(trans)
        c0 = sum(trans(i).repeat_id == rep0);
        c1 = sum(trans(i).repeat_id == rep1);
        if c0 > 0 && c1 > 0
            score = c0 + c1;
            if score > best
                best = score;
                ref_i = i;
            end
        end
    end
    if best >= 0
        return;
    end

    % Robust fallback: search any adjacent repeat pair and transition.
    best = -1;
    found = false;
    for i = 1:numel(trans)
        [r0, r1, ok, score] = best_adjacent_pair_score(trans(i).repeat_id, rep_start);
        if ok && score > best
            best = score;
            ref_i = i;
            rep0 = r0;
            rep1 = r1;
            found = true;
        end
    end

    if ~found
        error('No valid transition can be used as reference in any adjacent repeat pair.');
    end

    if verbose
        warning('No ref transition in repeat%d/repeat%d. Fallback to repeat%d/repeat%d with transition %d.', ...
            rep_start, rep_start+1, rep0, rep1, ref_i);
    end
end

function [rep0, rep1, ok] = find_adjacent_pair_for_transition(repeat_ids, rep_start)
    rep0 = rep_start;
    rep1 = rep_start + 1;
    ok = false;
    if isempty(repeat_ids)
        return;
    end
    u = unique(repeat_ids(:)).';
    u = u(u >= rep_start);
    for k = 1:(numel(u)-1)
        if u(k+1) == u(k) + 1
            rep0 = u(k);
            rep1 = u(k+1);
            ok = true;
            return;
        end
    end
end

function [rep0, rep1, ok, score] = best_adjacent_pair_score(repeat_ids, rep_start)
    rep0 = rep_start;
    rep1 = rep_start + 1;
    ok = false;
    score = -inf;
    if isempty(repeat_ids)
        return;
    end

    u = unique(repeat_ids(:)).';
    u = u(u >= rep_start);
    for k = 1:(numel(u)-1)
        if u(k+1) == u(k) + 1
            c0 = sum(repeat_ids == u(k));
            c1 = sum(repeat_ids == u(k+1));
            sc = c0 + c1;
            if sc > score
                score = sc;
                rep0 = u(k);
                rep1 = u(k+1);
                ok = true;
            end
        end
    end
end

function v = mean_or_nan(x)
    if isempty(x)
        v = nan;
    else
        v = mean(x);
    end
end

function make_plots(y, V, ths, out, trans, rep0, rep1)
    figure('Name', 'EOJ PAM4 Summary', 'Color', 'w');

    subplot(2,2,1);
    histogram(y, 120);
    hold on;
    xline(V(1), '--', 'V0'); xline(V(2), '--', 'V1');
    xline(V(3), '--', 'V2'); xline(V(4), '--', 'V3');
    xline(ths(1), '-', 'th01'); xline(ths(2), '-', 'th12'); xline(ths(3), '-', 'th23');
    hold off;
    title('Center samples y[n] histogram');
    xlabel('Voltage (V)'); ylabel('Count'); grid on;

    subplot(2,2,2);
    bar(out.EOJ_i_ps);
    title(sprintf('EOJ_i (ps), max EOJ = %.4f ps', out.EOJ_ps));
    xlabel('Transition index'); ylabel('EOJ_i (ps)'); grid on;

    subplot(2,2,[3 4]);
    r = out.ref_transition;
    d0 = trans(r).tcross_cru_abs(trans(r).repeat_id == rep0) * 1e12;
    d1 = trans(r).tcross_cru_abs(trans(r).repeat_id == rep1) * 1e12;
    if isempty(d0), d0 = nan; end
    if isempty(d1), d1 = nan; end
    histogram(d0, 40, 'FaceAlpha', 0.6); hold on;
    histogram(d1, 40, 'FaceAlpha', 0.6);
    hold off;
    legend(sprintf('repeat %d', rep0), sprintf('repeat %d', rep1), 'Location', 'best');
    title(sprintf('Ref transition #%d (%s): tcross_{cru} distribution', r, trans(r).name));
    xlabel('Time (ps)'); ylabel('Count'); grid on;
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
