function res = qprbs13_resample_module(cfg)
%QPRBS13_RESAMPLE_MODULE
% Resampling module:
% CSV(t,vdiff) -> scan t0 -> resample numCycles -> y_avg -> Y.

    validate_cfg(cfg);

    [t, vdiff] = read_csv_data(cfg.csvPath);
    [t0Candidates, JScan, t0Best, const] = scan_t0(t, vdiff, cfg);
    [yCycles, yAvg, Y] = resample_with_t0(t, vdiff, t0Best, const, cfg.numCycles);

    res = struct();
    res.t0_best = t0Best;
    res.t0_candidates = t0Candidates;
    res.J_scan = JScan;
    res.ui = const.ui;
    res.M = const.M;
    res.N = const.N;
    res.Ts = const.Ts;
    res.Tpat = const.Tpat;
    res.numCycles = cfg.numCycles;
    res.y_avg = yAvg;
    res.Y = Y;
    if isfield(cfg, 'saveCycles') && logical(cfg.saveCycles)
        res.y_cycles = yCycles;
    end
end

function validate_cfg(cfg)
    req = {'csvPath','tStableStart','ui','M','N','numCycles','dtScan','scanWindowUI'};
    for i = 1:numel(req)
        if ~isfield(cfg, req{i}) || isempty(cfg.(req{i}))
            error('qprbs13_resample_module requires cfg.%s.', req{i});
        end
    end
end

function [t, vdiff] = read_csv_data(csvPath)
    T = readtable(csvPath, 'ReadVariableNames', true);

    % Support both styles:
    % 1) Header present with columns: t, vdiff
    % 2) No header, first column = time, second column = voltage
    varNames = lower(T.Properties.VariableNames);
    idxT = find(strcmp(varNames, 't'), 1);
    idxV = find(strcmp(varNames, 'vdiff'), 1);
    hasNamedCols = ~isempty(idxT) && ~isempty(idxV);

    if hasNamedCols
        t = double(T{:, idxT});
        vdiff = double(T{:, idxV});
    else
        raw = table2array(T);
        if size(raw, 2) < 2
            error('CSV must contain at least 2 columns: time and voltage.');
        end
        warning('CSV without recognized headers (t, vdiff). Using first two columns as [t, vdiff].');
        t = double(raw(:,1));
        vdiff = double(raw(:,2));
    end

    t = t(:);
    vdiff = vdiff(:);

    if numel(t) ~= numel(vdiff) || numel(t) < 2
        error('t and vdiff must have same length and at least 2 rows.');
    end

    finiteMask = isfinite(t) & isfinite(vdiff);
    t = t(finiteMask);
    vdiff = vdiff(finiteMask);

    if numel(t) < 2
        error('Not enough finite samples after removing NaN/Inf.');
    end

    if ~all(diff(t) > 0)
        warning('t is not strictly increasing. Applying sort + unique by t.');
        [tSort, idxSort] = sort(t, 'ascend');
        vSort = vdiff(idxSort);
        [t, idxUniq] = unique(tSort, 'stable');
        vdiff = vSort(idxUniq);
    end

    if numel(t) < 2
        error('Time vector has <2 unique points after sort+unique.');
    end
end

function [t0Candidates, J, t0Best, const] = scan_t0(t, vdiff, cfg)
    const = struct();
    const.ui = double(cfg.ui);
    const.M = double(cfg.M);
    const.N = double(cfg.N);
    const.Ts = const.ui / const.M;
    const.Tpat = const.N * const.ui;
    const.MN = const.M * const.N;

    useFirstCrossing = true;
    if isfield(cfg, 'useFirstCrossing') && ~isempty(cfg.useFirstCrossing)
        useFirstCrossing = logical(cfg.useFirstCrossing);
    end

    if useFirstCrossing
        t0Cross = find_first_crossing_time(t, vdiff, cfg);
        t0Candidates = double(t0Cross);
        J = 0;
    else
        t0Min = double(cfg.tStableStart - cfg.scanWindowUI * const.ui);
        t0Max = double(cfg.tStableStart + cfg.scanWindowUI * const.ui);

        t0Candidates = (t0Min:cfg.dtScan:t0Max).';
        if isempty(t0Candidates) || t0Candidates(end) < t0Max
            t0Candidates = [t0Candidates; t0Max]; %#ok<AGROW>
        end

        J = inf(size(t0Candidates), 'double');
    end
    for i = 1:numel(t0Candidates)
        t0 = t0Candidates(i);
        tEnd = t0 + (cfg.numCycles - 1) * const.Tpat + (const.MN - 1) * const.Ts;
        if t0 < t(1) || tEnd > t(end)
            continue;
        end

        yCycles = zeros(const.MN, cfg.numCycles, 'double');
        valid = true;
        for c = 0:(cfg.numCycles-1)
            tgrid = t0 + c * const.Tpat + (0:const.MN-1).' * const.Ts;
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

        yAvgTmp = mean(yCycles, 2);
        d = yCycles - yAvgTmp;
        if ~useFirstCrossing
            J(i) = sum(d(:).^2);
        else
            J(i) = 0;
        end
    end

    [jMin, idxMin] = min(J);
    if isinf(jMin)
        reqStart = t0Candidates(1);
        reqEnd = t0Candidates(end) + (cfg.numCycles - 1) * const.Tpat + (const.MN - 1) * const.Ts;
        error(['No valid t0 in scan window. Need data to cover about [%.16g, %.16g] s. ', ...
               'Available [%.16g, %.16g] s.'], reqStart, reqEnd, t(1), t(end));
    end

    t0Best = t0Candidates(idxMin);

    reqStart = t0Best;
    reqEnd = t0Best + (cfg.numCycles - 1) * const.Tpat + (const.MN - 1) * const.Ts;
    if reqStart < t(1) || reqEnd > t(end)
        error(['Insufficient data coverage. Required [%.16g, %.16g] s, ', ...
               'available [%.16g, %.16g] s.'], reqStart, reqEnd, t(1), t(end));
    end
end


function tCross = find_first_crossing_time(t, v, cfg)
% Find the first transition crossing using adaptive threshold (not fixed zero).
% Steps:
%   1) Quantize waveform into 4 amplitude levels by kmeans.
%   2) Find first adjacent samples whose level code changes.
%   3) Use threshold = midpoint between the two involved level centers.
%   4) Interpolate crossing time where waveform crosses that threshold.
    t = double(t(:));
    v = double(v(:));

    if numel(t) ~= numel(v) || numel(t) < 2
        error('t and v must have same length and at least 2 points.');
    end

    % Optional override for fixed threshold crossing.
    useFixedThreshold = isfield(cfg, 'crossingThreshold') && ~isempty(cfg.crossingThreshold);
    if useFixedThreshold
        thr = double(cfg.crossingThreshold);
        validateattributes(thr, {'numeric'}, {'real','finite','scalar'});
        idx = find((v(1:end-1) - thr) .* (v(2:end) - thr) <= 0, 1, 'first');
        if isempty(idx)
            error('No crossing found for cfg.crossingThreshold=%g.', thr);
        end
        tCross = interpolate_crossing(t(idx), t(idx+1), v(idx), v(idx+1), thr);
        return;
    end

    try
        [idxRaw, centersRaw] = kmeans(v, 4, 'Replicates', 5, 'MaxIter', 200);
    catch
        % Fallback for environments lacking Statistics Toolbox.
        idxRaw = discretize(v, quantile(v, [0 0.25 0.5 0.75 1]));
        centersRaw = accumarray(max(1,min(4,idxRaw(~isnan(idxRaw)))), v(~isnan(idxRaw)), [4,1], @mean, mean(v));
        idxRaw(isnan(idxRaw)) = 1;
    end

    [centersSorted, ord] = sort(double(centersRaw(:)), 'ascend');
    mapOldToNew = zeros(numel(ord), 1);
    for k = 1:numel(ord)
        mapOldToNew(ord(k)) = k;
    end
    idxLevel = mapOldToNew(idxRaw);

    idx = find(idxLevel(1:end-1) ~= idxLevel(2:end), 1, 'first');
    if isempty(idx)
        % Last fallback: use midpoint between global min/max and find first crossing.
        thr = (min(v) + max(v)) / 2;
        idx = find((v(1:end-1) - thr) .* (v(2:end) - thr) <= 0, 1, 'first');
        if isempty(idx)
            error('No symbol transition crossing found in waveform.');
        end
    else
        c1 = centersSorted(idxLevel(idx));
        c2 = centersSorted(idxLevel(idx+1));
        thr = (c1 + c2) / 2;
    end

    tCross = interpolate_crossing(t(idx), t(idx+1), v(idx), v(idx+1), thr);
end

function tc = interpolate_crossing(t1, t2, v1, v2, thr)
    if v1 == thr
        tc = t1;
        return;
    end
    if v2 == thr
        tc = t2;
        return;
    end
    if v2 == v1
        tc = (t1 + t2) / 2;
        return;
    end
    tc = t1 + (thr - v1) * (t2 - t1) / (v2 - v1);
end

function [yCycles, yAvg, Y] = resample_with_t0(t, vdiff, t0Best, const, numCycles)
    yCycles = zeros(const.MN, numCycles, 'double');
    for c = 0:(numCycles-1)
        tgrid = t0Best + c * const.Tpat + (0:const.MN-1).' * const.Ts;
        yc = interp1(t, vdiff, tgrid, 'linear');
        if any(isnan(yc))
            error('Interpolation returned NaN at cycle %d.', c);
        end
        yCycles(:, c+1) = yc;
    end

    yAvg = mean(yCycles, 2);
    % Eq(11-13): Y(:,n) is one UI column with M samples.
    Y = reshape(yAvg, const.M, const.N);
end
