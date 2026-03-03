function result = pam4_jitter_analysis(csv_file, fb, M)
% PAM4 抖动后处理主函数（224G/112G TX 前仿真）
% ------------------------------------------------------------
% 输入:
%   csv_file : CSV 路径，包含两列 [time(s), vdiff(V)]
%   fb       : 符号率/波特率 (Hz)
%   M        : 每 UI 重采样点数，默认 64
%
% 输出:
%   result   : 结构体，包含 JRMS/J3u/每类样本数/每类均值 等结果
%
% 说明:
%   - 不依赖 toolbox（手写 1D kmeans 与分位数函数）
%   - 流程严格按用户给定的 15 个步骤实现

    if nargin < 3 || isempty(M)
        M = 64;
    end

    %% ==============================
    % Step 1 读取数据
    % ===============================
    if ~isfile(csv_file)
        error('CSV 文件不存在: %s', csv_file);
    end

    raw = readmatrix(csv_file);
    if size(raw, 2) < 2
        error('CSV 至少需要两列: time, vdiff');
    end

    t_raw = raw(:, 1);
    v_raw = raw(:, 2);

    % 清理 NaN/Inf
    valid = isfinite(t_raw) & isfinite(v_raw);
    t_raw = t_raw(valid);
    v_raw = v_raw(valid);

    if numel(t_raw) < 10
        error('有效样本过少，无法分析。');
    end

    % 时间排序（防止乱序）
    [t_raw, idx_sort] = sort(t_raw);
    v_raw = v_raw(idx_sort);

    dt = diff(t_raw);
    dt_med = median(dt);
    if dt_med <= 0
        error('时间序列非法（非递增）。');
    end

    tol = 1e-3; % 0.1% 容差
    is_uniform = max(abs(dt - dt_med)) <= tol * dt_med;

    if ~is_uniform
        fprintf('[Step1] 检测到时间不等间距，进行等间距插值。\n');
        t_eq = (t_raw(1):dt_med:t_raw(end)).';
        v_eq = interp1(t_raw, v_raw, t_eq, 'linear', 'extrap');
        t_raw = t_eq;
        v_raw = v_eq;
    else
        fprintf('[Step1] 时间等间距，直接使用原始数据。\n');
    end

    %% ==============================
    % Step 2 重采样
    % ===============================
    UI = 1 / fb;
    Ts = UI / M;

    t_uniform = (t_raw(1):Ts:t_raw(end)).';
    if numel(t_uniform) < M * 20
        warning('重采样后样本偏少，统计可信度可能不足。');
    end
    v_uniform = interp1(t_raw, v_raw, t_uniform, 'linear', 'extrap');

    %% ==============================
    % Step 3 UI 分块
    % ===============================
    N_total = numel(v_uniform);
    N_UI = floor(N_total / M);

    if N_UI < 20
        error('可用 UI 数量过少（N_UI=%d），无法稳定统计。', N_UI);
    end

    rem_samp = mod(N_total, M);
    if rem_samp ~= 0
        warning('重采样点数无法整除 M，尾部丢弃 %d 个样本。', rem_samp);
    end

    v_trim = v_uniform(1:N_UI * M);
    t_trim = t_uniform(1:N_UI * M);

    Vui = reshape(v_trim, M, N_UI).'; %#ok<NASGU>
    Tui = reshape(t_trim, M, N_UI).';

    %% ==============================
    % Step 4 中心采样
    % ===============================
    m = round(M / 2);
    y = reshape(v_trim, M, N_UI).';
    y_center = y(:, m);

    %% ==============================
    % Step 5 4 电平聚类
    % ===============================
    [cluster_idx, centers] = kmeans_1d_no_toolbox(y_center, 4, 100);

    % 对电平中心排序并映射到 symbol {0,1,2,3}
    [cent_sorted, ord] = sort(centers(:), 'ascend');
    map_old2new = zeros(4, 1);
    for ii = 1:4
        map_old2new(ord(ii)) = ii - 1;
    end

    symbol = zeros(size(cluster_idx));
    for ii = 1:numel(cluster_idx)
        symbol(ii) = map_old2new(cluster_idx(ii));
    end

    %% ==============================
    % Step 6 识别 AAAABB transition（先记录候选UI索引，统计放在CRU之后）
    % ===============================
    transition_names = generate_transition_names(); % 固定 12 类
    class_indices_candidate = cell(12, 1); % 每类先存候选 UI index

    for n = 1:(N_UI - 5)
        A = symbol(n);
        % AAAA
        if all(symbol(n:n+3) == A)
            B = symbol(n + 4);
            if (symbol(n + 5) == B) && (A ~= B)
                cls = transition_class_id(A, B);
                if cls > 0
                    % crossing 发生在 A->B 切换处，即第 n+4 个 UI 的开始边界附近
                    ui_idx = n + 4;
                    class_indices_candidate{cls}(end+1, 1) = ui_idx; %#ok<AGROW>
                end
            end
        end
    end

    %% ==============================
    % Step 7 计算电平均值
    % ===============================
    V = zeros(4, 1);
    for lv = 0:3
        yy = y_center(symbol == lv);
        if isempty(yy)
            error('电平 %d 没有样本，无法计算阈值。', lv);
        end
        V(lv + 1) = mean(yy);
    end
    V0 = V(1); V1 = V(2); V2 = V(3); V3 = V(4); %#ok<NASGU>

    %% ==============================
    % Step 8 计算阈值（方法A）
    % ===============================
    th01 = (V(1) + V(2)) / 2;
    th12 = (V(2) + V(3)) / 2;
    th23 = (V(3) + V(4)) / 2;

    %% ==============================
    % Step 9 先提取所有 crossing time（不按 AAAABB 分类）
    % ===============================
    tcross_abs = nan(N_UI, 1);       % 每UI边界（对应UI索引）的 crossing 时间
    cross_A = nan(N_UI, 1);          % crossing 前符号
    cross_B = nan(N_UI, 1);          % crossing 后符号
    crossing_not_found = 0;

    for ui_idx = 2:N_UI
        A = symbol(ui_idx - 1);
        B = symbol(ui_idx);

        % 仅在发生电平变化时提取 crossing
        if A == B
            continue;
        end

        th = choose_threshold(A, B, th01, th12, th23);

        % crossing 常发生在 UI 边界附近，使用跨边界窗口搜索
        t_boundary = Tui(ui_idx, 1);
        vwin = [y(ui_idx - 1, :), y(ui_idx, :)];
        twin = [Tui(ui_idx - 1, :), Tui(ui_idx, :)];

        [tc, ok] = find_crossing_near_boundary(twin, vwin, th, A, B, t_boundary);
        if ok
            tcross_abs(ui_idx) = tc;
            cross_A(ui_idx) = A;
            cross_B(ui_idx) = B;
        else
            crossing_not_found = crossing_not_found + 1;
        end
    end

    if crossing_not_found > 0
        warning('有 %d 个电平变化 UI 未找到 crossing。', crossing_not_found);
    end

    %% ==============================
    % Step 10 Golden PLL CRU
    % ===============================
    fc = fb / 13280;
    alpha = exp(-2 * pi * fc / fb);

    t_ideal = ((1:N_UI).' - 1) * UI + t_trim(1);

    tie_raw = nan(N_UI, 1);
    has_cross = ~isnan(tcross_abs);
    tie_raw(has_cross) = tcross_abs(has_cross) - t_ideal(has_cross);

    tie_LF = zeros(N_UI, 1);
    tie_HF = nan(N_UI, 1);

    for n = 2:N_UI
        if has_cross(n)
            e = tie_raw(n);
        else
            e = tie_LF(n - 1);
        end
        tie_LF(n) = alpha * tie_LF(n - 1) + (1 - alpha) * e;
    end

    valid_raw = has_cross;
    tie_HF(valid_raw) = tie_raw(valid_raw) - tie_LF(valid_raw);

    %% ==============================
    % Step 11 每类 transition 预处理
    % ===============================
    class_samples = cell(12, 1);
    class_mean_pre = nan(12, 1);
    drop_n = 500;

    for cls = 1:12
        idx_cand = class_indices_candidate{cls};
        idx_cand = idx_cand(~isnan(tie_HF(idx_cand))); % 先经过CRU并且有crossing的样本
        % 再次校验 A->B 与类别一致（避免误配）
        if ~isempty(idx_cand)
            ab = parse_transition_name(transition_names{cls});
            keep = (cross_A(idx_cand) == ab(1)) & (cross_B(idx_cand) == ab(2));
            idx_all = idx_cand(keep);
        else
            idx_all = [];
        end
        s = tie_HF(idx_all);

        if numel(s) > drop_n
            s = s(drop_n + 1:end);
        else
            s = [];
        end

        if numel(s) < 100
            warning('类别 %s 丢弃前500后剩余样本仅 %d (<100)。', transition_names{cls}, numel(s));
        end

        class_samples{cls} = s;
        if ~isempty(s)
            class_mean_pre(cls) = mean(s);
        end
    end

    %% ==============================
    % Step 12 强制样本数一致
    % ===============================
    Ni = cellfun(@numel, class_samples);
    if any(Ni == 0)
        error('至少有一个 transition 类别无可用样本，无法执行等样本截断。');
    end

    Nmin = min(Ni);
    for cls = 1:12
        class_samples{cls} = class_samples{cls}(1:Nmin);
    end
    Ni_final = cellfun(@numel, class_samples);

    %% ==============================
    % Step 13 每类去均值
    % ===============================
    Tavgi = zeros(12, 1);
    dt_class = cell(12, 1);

    for cls = 1:12
        s = class_samples{cls};
        Tavgi(cls) = mean(s);
        dt_class{cls} = s - Tavgi(cls);
    end

    %% ==============================
    % Step 14 合并
    % ===============================
    dt_all = vertcat(dt_class{:});

    %% ==============================
    % Step 15 统计
    % ===============================
    JRMS = std(dt_all);
    t_low = prctile_no_toolbox(dt_all, 0.05);
    t_high = prctile_no_toolbox(dt_all, 99.95);
    J3u = t_high - t_low;

    JRMS_ps = JRMS * 1e12;
    J3u_ps = J3u * 1e12;

    JRMS_UI = JRMS / UI;
    J3u_UI = J3u / UI;

    %% 绘图
    make_plots(y_center, symbol, transition_names, Ni_final, dt_all, t_low, t_high);

    %% 输出结果
    result = struct();
    result.fb = fb;
    result.UI = UI;
    result.fc = fc;
    result.alpha = alpha;
    result.level_centers = cent_sorted;
    result.symbol_means = struct('V0', V(1), 'V1', V(2), 'V2', V(3), 'V3', V(4));
    result.thresholds = struct('th01', th01, 'th12', th12, 'th23', th23);
    % 新增输出：重采样数据（time, vdiff）
    result.resampled = struct('time_s', t_uniform, 'vdiff_V', v_uniform, ...
                              'time_trim_s', t_trim, 'vdiff_trim_V', v_trim, ...
                              'M', M, 'Ts', Ts);
    result.transition_names = transition_names;
    result.samples_per_class = Ni_final;
    result.Tavgi_sec = Tavgi;
    result.JRMS_sec = JRMS;
    result.JRMS_ps = JRMS_ps;
    result.JRMS_UI = JRMS_UI;
    result.J3u_sec = J3u;
    result.J3u_ps = J3u_ps;
    result.J3u_UI = J3u_UI;
    result.crossing_not_found = crossing_not_found;

    fprintf('\n=========== PAM4 Jitter Analysis Result ===========\n');
    fprintf('JRMS = %.6e s | %.3f ps | %.6e UI\n', JRMS, JRMS_ps, JRMS_UI);
    fprintf('J3u  = %.6e s | %.3f ps | %.6e UI\n', J3u, J3u_ps, J3u_UI);
    fprintf('符号均值: V0=%.6e, V1=%.6e, V2=%.6e, V3=%.6e (V)\n', V(1), V(2), V(3), V(4));
    fprintf('阈值: th01=%.6e, th12=%.6e, th23=%.6e (V)\n', th01, th12, th23);
    fprintf('重采样数据点数: full=%d, trim=%d, Ts=%.6e s, M=%d\n', ...
        numel(v_uniform), numel(v_trim), Ts, M);
    fprintf('---------------------------------------------------\n');
    fprintf('每类 transition 最终样本数 (已强制一致 Nmin=%d):\n', Nmin);
    for cls = 1:12
        fprintf('  %-5s : %d\n', transition_names{cls}, Ni_final(cls));
    end
    fprintf('每类 Tavgi (sec):\n');
    for cls = 1:12
        fprintf('  %-5s : %.6e\n', transition_names{cls}, Tavgi(cls));
    end
    fprintf('===================================================\n\n');
end

function names = generate_transition_names()
    names = {'0->1','0->2','0->3','1->0','1->2','1->3', ...
             '2->0','2->1','2->3','3->0','3->1','3->2'};
end

function cls = transition_class_id(A, B)
    table = [0 1;0 2;0 3;1 0;1 2;1 3;2 0;2 1;2 3;3 0;3 1;3 2];
    cls = 0;
    for i = 1:size(table,1)
        if table(i,1) == A && table(i,2) == B
            cls = i;
            return;
        end
    end
end

function th = choose_threshold(A, B, th01, th12, th23)
    lo = min(A, B);
    hi = max(A, B);

    if lo == 0 && hi == 1
        th = th01;
    elseif lo == 1 && hi == 2
        th = th12;
    elseif lo == 2 && hi == 3
        th = th23;
    elseif lo == 0 && hi == 2
        th = th12;
    elseif lo == 1 && hi == 3
        th = th23;
    elseif lo == 0 && hi == 3
        th = th12;
    else
        error('未知的 A/B 组合: %d -> %d', A, B);
    end
end

function [tc, ok] = find_crossing_near_boundary(t, v, th, A, B, t_boundary)
    tc = nan;
    ok = false;

    % 方向约束的候选 crossing（上升沿 / 下降沿）
    if B > A
        cand = find(v(1:end-1) < th & v(2:end) >= th);
    else
        cand = find(v(1:end-1) > th & v(2:end) <= th);
    end

    % 若方向约束下没有，回退到任意过阈值
    if isempty(cand)
        cand = find((v(1:end-1)-th).*(v(2:end)-th) <= 0);
    end

    if isempty(cand)
        return;
    end

    % 对所有候选做线性插值，选择离边界最近的 crossing
    tc_all = nan(numel(cand), 1);
    for ii = 1:numel(cand)
        k = cand(ii);
        dv = v(k+1) - v(k);
        if abs(dv) < eps
            tc_all(ii) = t(k);
        else
            tc_all(ii) = t(k) + (th - v(k)) * (t(k+1) - t(k)) / dv;
        end
    end

    [~, id_best] = min(abs(tc_all - t_boundary));
    tc = tc_all(id_best);
    ok = true;
end

function [idx, centers] = kmeans_1d_no_toolbox(x, k, max_iter)
    x = x(:);
    n = numel(x);
    if n < k
        error('kmeans 样本不足。');
    end

    xs = sort(x);
    centers = zeros(k,1);
    for i = 1:k
        pos = round((i - 0.5) * n / k);
        pos = max(1, min(n, pos));
        centers(i) = xs(pos);
    end

    idx = ones(n,1);
    for it = 1:max_iter
        dist = abs(x - centers.');
        [~, new_idx] = min(dist, [], 2);

        new_centers = centers;
        for j = 1:k
            members = x(new_idx == j);
            if isempty(members)
                new_centers(j) = xs(randi(n));
            else
                new_centers(j) = mean(members);
            end
        end

        if all(new_idx == idx)
            centers = new_centers;
            break;
        end

        idx = new_idx;
        centers = new_centers;
    end
end

function q = prctile_no_toolbox(x, p)
    x = sort(x(:));
    n = numel(x);
    if n == 0
        q = nan;
        return;
    end
    if p <= 0
        q = x(1);
        return;
    elseif p >= 100
        q = x(end);
        return;
    end

    pos = (p/100) * (n - 1) + 1;
    lo = floor(pos);
    hi = ceil(pos);
    w = pos - lo;

    if lo == hi
        q = x(lo);
    else
        q = x(lo) * (1 - w) + x(hi) * w;
    end
end

function ab = parse_transition_name(name)
    % 将 'A->B' 解析为 [A B]（兼容无 toolbox 环境）
    tok = regexp(name, '^(\d)->(\d)$', 'tokens', 'once');
    if isempty(tok)
        error('非法 transition 名称: %s', name);
    end
    ab = [str2double(tok{1}), str2double(tok{2})];
end

function make_plots(y_center, symbol, transition_names, Ni, dt_all, t_low, t_high)
    figure('Name', 'PAM4 电平分布', 'Color', 'w');
    subplot(2,1,1);
    histogram(y_center, 120);
    grid on;
    xlabel('中心采样电压 (V)');
    ylabel('计数');
    title('中心采样电压分布');

    subplot(2,1,2);
    scatter(1:numel(y_center), y_center, 4, symbol, 'filled');
    grid on;
    xlabel('UI 索引');
    ylabel('中心采样电压 (V)');
    title('中心采样与聚类结果');

    figure('Name', '每类样本数量', 'Color', 'w');
    bar(Ni);
    grid on;
    set(gca, 'XTick', 1:12, 'XTickLabel', transition_names, 'XTickLabelRotation', 45);
    ylabel('样本数');
    title('12 类 transition 最终样本数（已对齐）');

    figure('Name', '\Delta t_all 直方图', 'Color', 'w');
    histogram(dt_all, 150);
    hold on;
    yl = ylim;
    plot([t_low t_low], yl, 'r--', 'LineWidth', 1.5);
    plot([t_high t_high], yl, 'm--', 'LineWidth', 1.5);
    grid on;
    xlabel('\Delta t (s)');
    ylabel('计数');
    title('\Delta t_{all} 直方图（标出 0.05% 与 99.95%）');
    legend('hist', '0.05%', '99.95%', 'Location', 'best');
end
