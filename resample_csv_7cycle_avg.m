function out = resample_csv_7cycle_avg(csvFile, tStableStart, cfg)
%RESAMPLE_CSV_7CYCLE_AVG Resample CSV(t, vdiff) and average 7 pattern cycles.
%   out = resample_csv_7cycle_avg(csvFile, tStableStart)
%   out = resample_csv_7cycle_avg(csvFile, tStableStart, cfg)
%
% Inputs
%   csvFile      : CSV path with fixed column names: t, vdiff (s, V)
%   tStableStart : user-provided stable start time (s), may be off UI boundary
%   cfg          : optional struct fields
%       .ui         (default 1/112e9)
%       .M          (default 32)
%       .N          (default 8191)
%       .numCycles  (default 7)
%       .outputDir  (default pwd)
%       .dtScan     (default Ts/4)
%       .overlayUI  (default [1, round(N/2), N])
%
% Outputs
%   out struct fields
%       .t0_best, .y_avg, .Y, .ui, .M, .N, .Ts, .Tpat
%       .y_cycles (MN x numCycles), .J_scan, .t0_candidates
%
% Files written to cfg.outputDir
%   y_avg.mat, Y.mat, y_avg.csv
%   J_scan.png, overlay_cycles.png, heatmap_Y.png

    if nargin < 3
        cfg = struct();
    end

    validateattributes(csvFile, {'char', 'string'}, {'scalartext'}, mfilename, 'csvFile');
    validateattributes(tStableStart, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'tStableStart');

    % ---- Defaults ----
    ui = get_cfg(cfg, 'ui', 1/112e9);
    M = get_cfg(cfg, 'M', 32);
    N = get_cfg(cfg, 'N', 8191);
    numCycles = get_cfg(cfg, 'numCycles', 7);
    outputDir = get_cfg(cfg, 'outputDir', pwd);

    validateattributes(ui, {'numeric'}, {'real', 'finite', 'positive', 'scalar'});
    validateattributes(M, {'numeric'}, {'real', 'finite', 'positive', 'integer', 'scalar'});
    validateattributes(N, {'numeric'}, {'real', 'finite', 'positive', 'integer', 'scalar'});
    validateattributes(numCycles, {'numeric'}, {'real', 'finite', 'positive', 'integer', 'scalar'});

    M = double(M);
    N = double(N);
    numCycles = double(numCycles);
    ui = double(ui);
    tStableStart = double(tStableStart);

    Ts = ui / M;
    Tpat = N * ui;
    MN = M * N;

    dtScan = get_cfg(cfg, 'dtScan', Ts / 4);
    overlayUI = get_cfg(cfg, 'overlayUI', [1, round(N/2), N]);

    validateattributes(dtScan, {'numeric'}, {'real', 'finite', 'positive', 'scalar'});

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % ---- Read CSV ----
    T = readtable(csvFile);
    reqVars = {'t', 'vdiff'};
    for i = 1:numel(reqVars)
        if ~ismember(reqVars{i}, T.Properties.VariableNames)
            error('CSV must contain column "%s".', reqVars{i});
        end
    end

    t = double(T.t(:));
    vdiff = double(T.vdiff(:));

    if numel(t) ~= numel(vdiff) || numel(t) < 2
        error('Columns t and vdiff must have the same length and contain at least 2 samples.');
    end

    finiteMask = isfinite(t) & isfinite(vdiff);
    t = t(finiteMask);
    vdiff = vdiff(finiteMask);

    if numel(t) < 2
        error('Not enough finite samples after removing NaN/Inf.');
    end

    % ---- If t is not strictly increasing: sort + unique ----
    [tSorted, idxSort] = sort(t, 'ascend');
    vSorted = vdiff(idxSort);
    [tUniq, idxUniq] = unique(tSorted, 'stable');
    vUniq = vSorted(idxUniq);

    t = double(tUniq);
    vdiff = double(vUniq);

    if numel(t) < 2
        error('Time vector has fewer than 2 unique samples after sort/unique.');
    end

    % ---- Local fine alignment scan ----
    t0Min = tStableStart - 0.5 * ui;
    t0Max = tStableStart + 0.5 * ui;
    t0Candidates = (t0Min:dtScan:t0Max).';
    if isempty(t0Candidates) || t0Candidates(end) < t0Max
        t0Candidates = [t0Candidates; t0Max];
    end

    J = inf(size(t0Candidates));

    for i = 1:numel(t0Candidates)
        t0 = t0Candidates(i);
        tEndNeed = t0 + (numCycles-1) * Tpat + (MN-1) * Ts;

        if t0 < t(1) || tEndNeed > t(end)
            continue;
        end

        yCycles = zeros(MN, numCycles, 'double');
        valid = true;
        for c = 0:(numCycles-1)
            tgrid = t0 + c * Tpat + (0:MN-1).' * Ts;
            yc = interp1(t, vdiff, tgrid, 'linear');
            if any(isnan(yc))
                valid = false;
                break;
            end
            yCycles(:, c+1) = yc;
        end

        if ~valid
            continue;
        end

        yAvgCand = mean(yCycles, 2);
        diffMat = yCycles - yAvgCand;
        J(i) = sum(diffMat(:).^2);
    end

    [Jmin, iBest] = min(J);
    if isinf(Jmin)
        error(['No valid t0 found in scan window. Ensure data covers scan candidates and ', ...
               'required cycles around tStableStart.']);
    end

    t0_best = t0Candidates(iBest);

    % ---- Coverage check for chosen t0_best ----
    tNeedStart = t0_best;
    tNeedEnd = t0_best + (numCycles-1) * Tpat + (MN-1) * Ts;
    if tNeedStart < t(1) || tNeedEnd > t(end)
        error(['Insufficient data coverage. Required t range: [%.16g, %.16g] s; ', ...
               'available range: [%.16g, %.16g] s.'], tNeedStart, tNeedEnd, t(1), t(end));
    end

    % ---- Final resampling at t0_best ----
    y_cycles = zeros(MN, numCycles, 'double');
    for c = 0:(numCycles-1)
        tgrid = t0_best + c * Tpat + (0:MN-1).' * Ts;
        yc = interp1(t, vdiff, tgrid, 'linear');
        if any(isnan(yc))
            error('Interpolation produced NaN at cycle %d. Check time coverage.', c);
        end
        y_cycles(:, c+1) = yc;
    end

    y_avg = mean(y_cycles, 2);
    Y = reshape(y_avg, M, N);

    out = struct();
    out.t0_best = t0_best;
    out.y_avg = y_avg;
    out.Y = Y;
    out.ui = ui;
    out.M = M;
    out.N = N;
    out.Ts = Ts;
    out.Tpat = Tpat;
    out.y_cycles = y_cycles;
    out.J_scan = J;
    out.t0_candidates = t0Candidates;

    % ---- Save MAT files ----
    save(fullfile(outputDir, 'y_avg.mat'), 'y_avg');
    save(fullfile(outputDir, 'Y.mat'), 'Y');

    % ---- Save CSV ----
    k = (1:MN).';
    ui_index = floor((k-1) / M) + 1;
    sample_in_ui = mod(k-1, M) + 1;
    y = y_avg;
    Tout = table(k, ui_index, sample_in_ui, y);
    writetable(Tout, fullfile(outputDir, 'y_avg.csv'));

    % ---- Plot 1: J scan ----
    f1 = figure('Visible', 'off', 'Color', 'w');
    plot(t0Candidates, J, '-', 'LineWidth', 1.2);
    hold on;
    plot(t0_best, Jmin, 'ro', 'MarkerFaceColor', 'r');
    hold off;
    grid on;
    xlabel('t0 (s)');
    ylabel('J(t0)');
    title('Alignment scan: J vs t0');
    saveas(f1, fullfile(outputDir, 'J_scan.png'));
    close(f1);

    % ---- Plot 2: overlay cycles for selected UI columns ----
    overlayUI = unique(max(1, min(N, round(double(overlayUI(:).')))));
    if isempty(overlayUI)
        overlayUI = [1, round(N/2), N];
    end

    nShow = numel(overlayUI);
    f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 420*nShow, 360]);
    tiledlayout(1, nShow, 'Padding', 'compact', 'TileSpacing', 'compact');

    s = (1:M).';
    for ii = 1:nShow
        uiCol = overlayUI(ii);
        idx = (uiCol-1)*M + (1:M);
        nexttile;
        hold on;
        for c = 1:numCycles
            plot(s, y_cycles(idx, c), '-', 'LineWidth', 0.9);
        end
        plot(s, y_avg(idx), 'k-', 'LineWidth', 2.0);
        hold off;
        grid on;
        xlabel('sample\_in\_ui');
        ylabel('v');
        title(sprintf('UI column %d', uiCol));
        if ii == nShow
            legend([arrayfun(@(x) sprintf('cycle %d', x), 1:numCycles, 'UniformOutput', false), {'avg'}], ...
                   'Location', 'best');
        end
    end
    saveas(f2, fullfile(outputDir, 'overlay_cycles.png'));
    close(f2);

    % ---- Plot 3: heatmap of Y ----
    f3 = figure('Visible', 'off', 'Color', 'w');
    imagesc(1:N, 1:M, Y);
    axis xy;
    xlabel('ui\_index (1..N)');
    ylabel('sample\_in\_ui (1..M)');
    title('Heatmap of Y = reshape(y\_avg, M, N)');
    colormap(turbo);
    colorbar;
    saveas(f3, fullfile(outputDir, 'heatmap_Y.png'));
    close(f3);
end

function v = get_cfg(cfg, name, defaultValue)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = defaultValue;
    end
end
