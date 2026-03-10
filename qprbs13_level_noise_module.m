function noiseOut = qprbs13_level_noise_module(Y, symCode, cfg)
%QPRBS13_LEVEL_NOISE_MODULE
% Estimate noise RMS per PAM4 level using runs with >= minRunLen identical symbols.
%
% For each level (0..3):
%   1) Find runs of consecutive identical symbols with length >= minRunLen.
%   2) At fixed UI phase sample index m0, take one voltage sample per UI in each run.
%   3) Compute RMS deviation relative to mean:
%      - 'merge'   : merge all eligible run samples of this level, then RMS-deviation once
%      - 'weighted': RMS-deviation per run, then weighted average by run length
%      - 'mean'    : RMS-deviation per run, then arithmetic mean
% Then sigma_n = mean([sigma_L0..sigma_L3]).
%
% Inputs
%   Y       : MxN matrix
%   symCode : Nx1 symbol code in {0,1,2,3}
%   cfg optional:
%       .minRunLen         (default 6)
%       .fixedSampleMethod ('center'|'index', default 'center')
%       .fixedSampleIndex  (default round(M/2), used when method='index')
%       .levelRmsMethod    ('merge'|'weighted'|'mean', default 'merge')
%
% Outputs
%   noiseOut struct:
%       .m0, .minRunLen, .method
%       .sigma_levels (4x1 for levels 0..3)
%       .sigma_n
%       .samples_per_level (4x1)
%       .runs_per_level (4x1)

    Y = double(Y);
    symCode = double(symCode(:));

    [M, N] = size(Y);
    if numel(symCode) ~= N
        error('symCode length (%d) must equal Y columns N (%d).', numel(symCode), N);
    end
    if any(~ismember(symCode, [0,1,2,3]))
        error('symCode must contain only 0/1/2/3.');
    end

    minRunLen = double(get_cfg(cfg, 'minRunLen', 6));
    fixedSampleMethod = lower(char(get_cfg(cfg, 'fixedSampleMethod', 'center')));
    levelRmsMethod = lower(char(get_cfg(cfg, 'levelRmsMethod', 'merge')));

    validateattributes(minRunLen, {'numeric'}, {'real','finite','integer','>=',1,'scalar'});

    switch fixedSampleMethod
        case 'center'
            m0 = round(M/2);
        case 'index'
            m0 = double(get_cfg(cfg, 'fixedSampleIndex', round(M/2)));
            validateattributes(m0, {'numeric'}, {'real','finite','integer','>=',1,'<=',M,'scalar'});
        otherwise
            error('cfg.fixedSampleMethod must be ''center'' or ''index''.');
    end

    if ~ismember(levelRmsMethod, {'merge','weighted','mean'})
        error('cfg.levelRmsMethod must be ''merge'', ''weighted'', or ''mean''.');
    end

    sigmaLevels = nan(4,1);
    samplesPerLevel = zeros(4,1);
    runsPerLevel = zeros(4,1);

    for lev = 0:3
        runList = find_runs(symCode == lev, minRunLen);
        runsPerLevel(lev+1) = size(runList,1);

        if isempty(runList)
            error('Level %d has no runs with length >= %d.', lev, minRunLen);
        end

        runSigmas = zeros(size(runList,1),1);
        runLens = zeros(size(runList,1),1);
        merged = [];

        for r = 1:size(runList,1)
            s = runList(r,1);
            e = runList(r,2);
            idx = s:e;
            vals = Y(m0, idx).';

            runLens(r) = numel(vals);
            samplesPerLevel(lev+1) = samplesPerLevel(lev+1) + runLens(r);

            if runLens(r) < 2
                runSigmas(r) = 0;
            else
                % RMS deviation as in sigma = sqrt(1/N * sum((Vi - Vmean)^2))
                runSigmas(r) = std(vals, 1);
            end
            merged = [merged; vals]; %#ok<AGROW>
        end

        switch levelRmsMethod
            case 'merge'
                if numel(merged) < 2
                    sigmaLevels(lev+1) = 0;
                else
                    sigmaLevels(lev+1) = std(merged, 1);
                end
            case 'weighted'
                sigmaLevels(lev+1) = sum(runSigmas .* runLens) / sum(runLens);
            case 'mean'
                sigmaLevels(lev+1) = mean(runSigmas);
        end
    end

    sigma_n = mean(sigmaLevels);

    noiseOut = struct();
    noiseOut.m0 = m0;
    noiseOut.minRunLen = minRunLen;
    noiseOut.method = levelRmsMethod;
    noiseOut.sigma_levels = sigmaLevels;
    noiseOut.sigma_n = sigma_n;
    noiseOut.samples_per_level = samplesPerLevel;
    noiseOut.runs_per_level = runsPerLevel;
end

function runs = find_runs(mask, minLen)
    mask = logical(mask(:));
    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    lens = ends - starts + 1;
    keep = lens >= minLen;
    runs = [starts(keep), ends(keep)];
end

function v = get_cfg(cfg, name, defaultValue)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = defaultValue;
    end
end
