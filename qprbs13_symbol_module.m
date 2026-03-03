function sym = qprbs13_symbol_module(Y, cfg)
%QPRBS13_SYMBOL_MODULE
% Symbol inference module:
% Y(MxN) -> s(n) -> 4-class clustering -> symbol means -> ES1/ES2/RLM -> x(n).

    validate_cfg(cfg);

    [s, symCode, symCycle, levels, clusterInfo, counts] = infer_symbols_from_Y(Y, cfg);
    [V, ES1, ES2, RLM, x] = compute_es_and_x(s, symCode);

    sym = struct();
    sym.s = s;
    sym.sym_code = symCode;
    sym.sym_cycle = symCycle;
    sym.levels = levels;
    sym.cluster = clusterInfo;
    sym.V = V;
    sym.ES1 = ES1;
    sym.ES2 = ES2;
    sym.RLM = RLM;
    sym.counts = counts;
    sym.x = x;
end

function validate_cfg(cfg)
    req = {'uiSampleMethod','clusterMethod'};
    for i = 1:numel(req)
        if ~isfield(cfg, req{i}) || isempty(cfg.(req{i}))
            error('qprbs13_symbol_module requires cfg.%s.', req{i});
        end
    end
end

function [s, symCode, symCycle, levels, clusterInfo, counts] = infer_symbols_from_Y(Y, cfg)
    M = size(Y, 1);

    uiMethod = lower(char(cfg.uiSampleMethod));
    switch uiMethod
        case 'center'
            m0 = round(M/2);
            s = double(Y(m0, :).');
        case 'mean'
            m0 = NaN;
            s = double(mean(Y, 1).');
        case 'proxy_opt'
            [m0, s] = select_phase_by_proxy(Y, cfg);
        otherwise
            error('Unsupported cfg.uiSampleMethod. Use ''center'', ''mean'', or ''proxy_opt''.');
    end

    [idx, centers] = cluster_4(s, cfg);

    [centSorted, ord] = sort(double(centers), 'ascend');
    mapOldClassToCode = zeros(4,1);
    for k = 1:4
        mapOldClassToCode(ord(k)) = k - 1;
    end

    symCode = mapOldClassToCode(idx);

    levels = double([-1; -1/3; 1/3; 1]);
    symCycle = levels(symCode + 1);

    counts = [sum(symCode==0); sum(symCode==1); sum(symCode==2); sum(symCode==3)];
    if any(counts == 0)
        error('At least one symbol class is empty after clustering. counts = [%d %d %d %d].', counts(1), counts(2), counts(3), counts(4));
    end

    clusterInfo = struct();
    clusterInfo.method = lower(char(cfg.clusterMethod));
    clusterInfo.centers_sorted = centSorted;
    clusterInfo.idx_raw = idx;
    clusterInfo.ui_sample_method = uiMethod;
    clusterInfo.selected_phase_index = m0;
end

function [idx, centers] = cluster_4(s, cfg)
    switch lower(char(cfg.clusterMethod))
        case 'kmeans'
            [idx, centers] = kmeans(s, 4, 'Replicates', 10, 'MaxIter', 200);
            centers = centers(:);
        case 'gmm'
            gm = fitgmdist(s, 4, ...
                'RegularizationValue', 1e-6, ...
                'Replicates', 5, ...
                'Options', statset('MaxIter', 500));
            idx = cluster(gm, s);
            centers = gm.mu(:);
        otherwise
            error('Unsupported cfg.clusterMethod. Use ''kmeans'' or ''gmm''.');
    end
end

function [mBest, sBest] = select_phase_by_proxy(Y, cfg)
% Proxy-optimal phase (method 3): maximize class-separation / within-class spread.
    M = size(Y, 1);
    bestScore = -inf;
    mBest = 1;
    sBest = double(Y(1,:).');

    for m = 1:M
        s = double(Y(m,:).');
        [idx, centers] = cluster_4(s, cfg);
        c = sort(double(centers(:)), 'ascend');
        sep = min(diff(c));

        w = zeros(4,1);
        for k = 1:4
            vals = s(idx == k);
            if numel(vals) < 2
                w(k) = 0;
            else
                w(k) = std(vals, 0);
            end
        end
        within = mean(w);
        score = sep / (within + eps);

        if score > bestScore
            bestScore = score;
            mBest = m;
            sBest = s;
        end
    end
end

function [V, ES1, ES2, RLM, x] = compute_es_and_x(s, symCode)
    V_m1 = mean(s(symCode == 0));
    V_m1_3 = mean(s(symCode == 1));
    V_p1_3 = mean(s(symCode == 2));
    V_p1 = mean(s(symCode == 3));

    Vmid = (V_m1 + V_p1) / 2;

    % Eq(16-16 ~ 16-19)
    den1 = (V_m1 - Vmid);
    den2 = (V_p1 - Vmid);
    if abs(den1) < eps || abs(den2) < eps
        error('Degenerate voltage levels: denominator near zero while computing ES1/ES2.');
    end

    ES1 = (V_m1_3 - Vmid) / den1;
    ES2 = (V_p1_3 - Vmid) / den2;
    RLM = min([3*ES1, 3*ES2, 2 - 3*ES1, 2 - 3*ES2]);

    x = zeros(size(symCode), 'double');
    x(symCode == 0) = -1;
    x(symCode == 1) = -ES1;
    x(symCode == 2) = ES2;
    x(symCode == 3) = 1;

    V = struct('m1', V_m1, 'm1_3', V_m1_3, 'p1_3', V_p1_3, 'p1', V_p1, 'mid', Vmid);
end
