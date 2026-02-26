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
    T = readtable(csvPath);
    requiredCols = {'t', 'vdiff'};
    for i = 1:numel(requiredCols)
        if ~ismember(requiredCols{i}, T.Properties.VariableNames)
            error('CSV must contain column "%s".', requiredCols{i});
        end
    end

    t = double(T.t(:));
    vdiff = double(T.vdiff(:));

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

    t0Min = double(cfg.tStableStart - cfg.scanWindowUI * const.ui);
    t0Max = double(cfg.tStableStart + cfg.scanWindowUI * const.ui);

    t0Candidates = (t0Min:cfg.dtScan:t0Max).';
    if isempty(t0Candidates) || t0Candidates(end) < t0Max
        t0Candidates = [t0Candidates; t0Max]; %#ok<AGROW>
    end

    J = inf(size(t0Candidates), 'double');
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
        J(i) = sum(d(:).^2);
    end

    [jMin, idxMin] = min(J);
    if isinf(jMin)
        reqStart = t0Min;
        reqEnd = t0Max + (cfg.numCycles - 1) * const.Tpat + (const.MN - 1) * const.Ts;
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
