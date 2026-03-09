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

    % Step 5: time-axis origin anchor uses the first threshold crossing
    % instead of the first sample time (requested behavior).
    t0 = find_first_threshold_crossing_time(t_uniform, v_uniform, th01, th12, th23);

    % Build or infer transition definitions (12 classes)
    if isfield(cfg, 'trans_def') && ~isempty(cfg.trans_def)
        trans_def = cfg.trans_def;
        infer_info = struct('used', false, 'offset_ui', 0, 'coverage', nan, 'hit_count', nan, 'half_window_ui', cfg.auto_window_half_width_ui);
    else
        if cfg.verbose
            fprintf('[EOJ] cfg.trans_def not provided, infer transition windows from symbol stream with AAAABB assumption.\n');
        end
        [trans_def, infer_info] = infer_transitions_from_symbols(sym, cfg.Npat, cfg.auto_window_half_width_ui);
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
            if isfield(trans_def, 'search_begin_ui') && isfield(trans_def, 'search_end_ui')
                begin_ui = trans_def(i).search_begin_ui;
                end_ui = trans_def(i).search_end_ui;
            else
                begin_ui = trans_def(i).begin_ui;
                end_ui = trans_def(i).end_ui;
            end
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
        % No per-UI averaging: if multiple events exist in the same UI,
        % keep the latest event value seen in that UI.
        for k = 1:numel(all_evt.tcross_abs)
            n_ui = all_evt.ui_index_global(k);
            if n_ui < 1 || n_ui > N_UI
                continue;
            end
            n0 = n_ui - 1;
            t_ideal = t0 + n0 * UI;
            tie_raw = all_evt.tcross_abs(k) - t_ideal;
            has_event(n_ui) = true;
            event_value(n_ui) = tie_raw;
        end
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
            warning('Transition %d (%s) missing in repeat%d/repeat%d. EOJ_i set to NaN.', i, trans(i).name, rep0, rep1);
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
    out.infer_info = infer_info;
    out.trans = trans;
    out.all_crossings = all_evt;
    out.all_crossings_count = numel(all_evt.tcross_abs);
    out.UI = UI;
    out.t0 = t0;
    out.alpha = alpha;
    out.phase_search = struct('m_opt', m_opt, 'm_center', m_center, 'm_used', m, ...
        'valid', phase_valid, 'score', phase_score);

    % Additional debug/trace outputs requested
    out.Vui = Vui;
    out.symbol_inferred = sym;
    out.symbol_sample = struct('phase_index', m, 'phase_time_offset_s', (m-1)*(UI/cfg.M), ...
        'voltage', y);

    % All crossing times grouped by threshold
    out.all_crossing_threshold_times = struct( ...
        'th01', all_evt.tcross_abs(all_evt.th_id == 1), ...
        'th12', all_evt.tcross_abs(all_evt.th_id == 2), ...
        'th23', all_evt.tcross_abs(all_evt.th_id == 3));

    % Optional CSV exports for symbol trace / AAAABB markers / crossings
    if cfg.export_csv
        out.csv_export_files = export_eoj_csv_outputs(out, cfg, csv_file, N_UI);
    else
        out.csv_export_files = struct();
    end

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
        if out.infer_info.used
            fprintf('[EOJ] AAAABB infer offset=%d UI, coverage=%d/12, hits=%d, half_window=%d UI\n', ...
                out.infer_info.offset_ui, out.infer_info.coverage, out.infer_info.hit_count, out.infer_info.half_window_ui);
        end
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
    if ~isfield(cfg, 'export_csv') || isempty(cfg.export_csv), cfg.export_csv = false; end
    if ~isfield(cfg, 'export_dir') || isempty(cfg.export_dir), cfg.export_dir = '.'; end
    if ~isfield(cfg, 'auto_window_half_width_ui') || isempty(cfg.auto_window_half_width_ui)
        cfg.auto_window_half_width_ui = 1;
    end
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

function [trans_def, infer_info] = infer_transitions_from_symbols(sym, Npat, half_w)
    % Automatic fallback when cfg.trans_def is missing.
    %
    % IMPORTANT AAAABB classification rule used here:
    % For a transition A->B at boundary n->n+1, we only accept it as a
    % candidate class event when:
    %   - A-side consecutive run-length at boundary n is >= 4 UI
    %   - B-side consecutive run-length at boundary n+1 is >= 2 UI
    % which is equivalent to requiring the local pattern AAAABB with the
    % boundary between the 4th A and 1st B.
    %
    % Why this can be judged with "3 UI context":
    % the boundary UI itself already provides 1 UI on each side. Therefore,
    % checking AAAABB only needs:
    %   - 3 additional UI on the left of boundary (to reach 4A total)
    %   - 1 additional UI on the right of boundary (to reach 2B total)
    % implemented below by run-lengths, which is robust to longer runs too.
    %
    % Optimization (Plan B): search pattern start offset first, then run
    % AAAABB detection on the best-aligned pattern window.

    sym = sym(:);
    if numel(sym) < Npat + 2
        error('Not enough UI symbols to infer AAAABB transitions.');
    end
    half_w = max(0, round(half_w));

    max_off = min(Npat - 1, numel(sym) - Npat);
    best_cov = -1;
    best_hits = -1;
    best_off = 0;
    best_class_pos = cell(4,4);

    for off = 0:max_off
        s = sym(off + 1 : off + Npat);
        class_pos = cell(4,4);
        hit_count = 0;

        left_run = ones(1, Npat);
        for k = 2:Npat
            if s(k) == s(k-1)
                left_run(k) = left_run(k-1) + 1;
            end
        end
        right_run = ones(1, Npat);
        for k = Npat-1:-1:1
            if s(k) == s(k+1)
                right_run(k) = right_run(k+1) + 1;
            end
        end

        for n = 4:(Npat - 2)
            A = s(n);
            B = s(n + 1);
            if A == B
                continue;
            end
            left_ok = left_run(n) >= 4;
            right_ok = right_run(n+1) >= 2;
            if left_ok && right_ok
                class_pos{A+1, B+1}(end+1) = n; %#ok<AGROW>
                hit_count = hit_count + 1;
            end
        end

        cov = 0;
        for a = 0:3
            for b = 0:3
                if a == b, continue; end
                if ~isempty(class_pos{a+1,b+1})
                    cov = cov + 1;
                end
            end
        end

        if (cov > best_cov) || (cov == best_cov && hit_count > best_hits)
            best_cov = cov;
            best_hits = hit_count;
            best_off = off;
            best_class_pos = class_pos;
        end
    end

    class_pos = best_class_pos;
    infer_info = struct('used', true, 'offset_ui', best_off, ...
        'coverage', best_cov, 'hit_count', best_hits, 'half_window_ui', half_w);

% Build all 12 directed classes explicitly, each class picks one
    % representative window (first AAAABB hit in repeat0).
    trans_def = repmat(struct('name', '', 'begin_ui', 1, 'end_ui', 1, ...
        'search_begin_ui', 1, 'search_end_ui', 1, ...
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
                trans_def(idx).search_begin_ui = 1;
                trans_def(idx).search_end_ui = Npat;
            else
                p0 = pos(1); % boundary between 4th A and 1st B
                trans_def(idx).begin_ui = max(1, p0 - 3);
                trans_def(idx).end_ui = min(Npat, p0 + 2);
                trans_def(idx).search_begin_ui = max(1, p0 - half_w);
                trans_def(idx).search_end_ui = min(Npat, p0 + half_w);
            end
        end
    end

    if ~isempty(missing)
        warning('AAAABB auto-infer missing classes (offset=%d): %s', infer_info.offset_ui, strjoin(missing, ', '));
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



function t0 = find_first_threshold_crossing_time(t, v, th01, th12, th23)
    thresholds = [th01, th12, th23];
    t0 = t(1);

    for k = 1:(numel(v)-1)
        v0 = v(k);
        v1 = v(k+1);
        dv = v1 - v0;
        if dv == 0
            continue;
        end

        tcands = [];
        for ith = 1:3
            th = thresholds(ith);
            is_cross = (v0 < th && v1 >= th) || (v0 > th && v1 <= th);
            if is_cross
                tc = t(k) + (th - v0) * (t(k+1) - t(k)) / dv;
                tcands(end+1) = tc; %#ok<AGROW>
            end
        end

        if ~isempty(tcands)
            t0 = min(tcands);
            return;
        end
    end
end

function all_evt = collect_all_crossing_events(t, v, M, th01, th12, th23, ui_start)
    thresholds = [th01, th12, th23];
    tc = [];
    ui = [];
    th_id = [];
    th_val = [];

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
                th_id(end+1,1) = ith; %#ok<AGROW>
                th_val(end+1,1) = th; %#ok<AGROW>
            end
        end
    end

    all_evt = struct('tcross_abs', tc, 'ui_index_global', ui, ...
        'th_id', th_id, 'th_value', th_val);
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


function files = export_eoj_csv_outputs(out, cfg, csv_file, N_UI)
    if ~exist(cfg.export_dir, 'dir')
        mkdir(cfg.export_dir);
    end

    [~, stem, ~] = fileparts(csv_file);
    if isempty(stem)
        stem = 'eoj_output';
    end

    % 1) Symbol trace CSV: inferred symbols + sampling phase/voltage + AAAABB marker
    ui_global = (1:N_UI).';
    repeat_id = floor((ui_global - 1) / cfg.Npat);
    ui_in_repeat = mod(ui_global - 1, cfg.Npat) + 1;

    % expected position markers from transition definitions
    mark = false(N_UI,1);
    label = repmat({''}, N_UI, 1);
    for i = 1:numel(out.trans_def)
        b = out.trans_def(i).begin_ui;
        e = out.trans_def(i).end_ui;
        mask = (ui_in_repeat >= b) & (ui_in_repeat <= e);
        mark(mask) = true;
        for k = find(mask).'
            if isempty(label{k})
                label{k} = out.trans_def(i).name;
            else
                label{k} = [label{k}, '|', out.trans_def(i).name]; %#ok<AGROW>
            end
        end
    end

    % actual detected crossing markers from extracted transition events
    det_mark = false(N_UI,1);
    det_label = repmat({''}, N_UI, 1);
    for i = 1:numel(out.trans)
        ui_evt = out.trans(i).ui_index_global(:);
        for jj = 1:numel(ui_evt)
            ug = ui_evt(jj);
            if ug < 1 || ug > N_UI
                continue;
            end
            det_mark(ug) = true;
            if isempty(det_label{ug})
                det_label{ug} = out.trans(i).name;
            else
                det_label{ug} = [det_label{ug}, '|', out.trans(i).name]; %#ok<AGROW>
            end
        end
    end

    T_sym = table(ui_global, repeat_id, ui_in_repeat, out.symbol_inferred(:), ...
        repmat(out.symbol_sample.phase_index, N_UI, 1), ...
        repmat(out.symbol_sample.phase_time_offset_s, N_UI, 1), ...
        out.symbol_sample.voltage(:), mark, string(label), det_mark, string(det_label), ...
        'VariableNames', {'ui_global','repeat_id','ui_in_repeat','symbol', ...
        'sample_phase_index','sample_phase_offset_s','sample_voltage_V', ...
        'is_aaaabb_transition_pos','aaaabb_transition_name', ...
        'is_detected_transition_pos','detected_transition_name'});

    sym_csv = fullfile(cfg.export_dir, [stem, '_symbol_trace.csv']);
    writetable(T_sym, sym_csv);

    % 2) AAAABB transition definition CSV
    N = numel(out.trans_def);
    idx_i = (1:N).';
    name = strings(N,1);
    begin_ui = zeros(N,1);
    end_ui = zeros(N,1);
    thr = strings(N,1);
    dir = strings(N,1);
    for i = 1:N
        name(i) = string(out.trans_def(i).name);
        begin_ui(i) = out.trans_def(i).begin_ui;
        end_ui(i) = out.trans_def(i).end_ui;
        thr(i) = string(out.trans_def(i).thr_type);
        dir(i) = string(out.trans_def(i).dir);
    end
    T_def = table(idx_i, name, begin_ui, end_ui, thr, dir, ...
        'VariableNames', {'transition_index','name','begin_ui','end_ui','thr_type','dir'});
    def_csv = fullfile(cfg.export_dir, [stem, '_aaaabb_transitions.csv']);
    writetable(T_def, def_csv);

    % 3) Actual detected transition crossings CSV (per event)
    rows = 0;
    for i = 1:numel(out.trans)
        rows = rows + numel(out.trans(i).tcross_abs);
    end
    tr_idx = zeros(rows,1);
    tr_name = strings(rows,1);
    rep_id = zeros(rows,1);
    ui_evt = zeros(rows,1);
    tc_abs = nan(rows,1);
    tc_cru = nan(rows,1);
    ridx = 0;
    for i = 1:numel(out.trans)
        n = numel(out.trans(i).tcross_abs);
        if n == 0, continue; end
        rr = (ridx+1):(ridx+n);
        tr_idx(rr) = i;
        tr_name(rr) = string(out.trans(i).name);
        rep_id(rr) = out.trans(i).repeat_id(:);
        ui_evt(rr) = out.trans(i).ui_index_global(:);
        tc_abs(rr) = out.trans(i).tcross_abs(:);
        tc_cru(rr) = out.trans(i).tcross_cru_abs(:);
        ridx = ridx + n;
    end
    T_det = table(tr_idx, tr_name, rep_id, ui_evt, tc_abs, tc_cru, ...
        'VariableNames', {'transition_index','transition_name','repeat_id', ...
        'ui_index_global','tcross_abs_s','tcross_cru_abs_s'});
    det_csv = fullfile(cfg.export_dir, [stem, '_detected_transition_crossings.csv']);
    writetable(T_det, det_csv);

    % 4) All crossing threshold times CSV
    th_name = strings(numel(out.all_crossings.th_id),1);
    for k = 1:numel(th_name)
        if out.all_crossings.th_id(k) == 1
            th_name(k) = "th01";
        elseif out.all_crossings.th_id(k) == 2
            th_name(k) = "th12";
        else
            th_name(k) = "th23";
        end
    end
    T_cross = table(out.all_crossings.tcross_abs(:), out.all_crossings.ui_index_global(:), ...
        out.all_crossings.th_id(:), th_name, out.all_crossings.th_value(:), ...
        'VariableNames', {'tcross_abs_s','ui_index_global','threshold_id','threshold_name','threshold_value_V'});
    cross_csv = fullfile(cfg.export_dir, [stem, '_all_crossings.csv']);
    writetable(T_cross, cross_csv);

    files = struct('symbol_trace_csv', sym_csv, ...
        'aaaabb_transition_csv', def_csv, ...
        'detected_transition_csv', det_csv, ...
        'crossing_csv', cross_csv);
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
