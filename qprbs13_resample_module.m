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
        t0Cross = find_first_crossing_time(t, vdiff);
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


function tCross = find_first_crossing_time(t, v)
    t = double(t(:));
    v = double(v(:));

    iZero = find(v == 0, 1, 'first');
    if ~isempty(iZero)
        tCross = t(iZero);
        return;
    end

    idx = find(v(1:end-1) .* v(2:end) < 0, 1, 'first');
    if isempty(idx)
        error('No crossing found in waveform. Cannot set start point to first crossing time.');
    end

    t1 = t(idx); t2 = t(idx+1);
    v1 = v(idx); v2 = v(idx+1);
    tCross = t1 + (0 - v1) * (t2 - t1) / (v2 - v1);
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
