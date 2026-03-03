function out = qprbs13_pipeline(cfg)
%QPRBS13_PIPELINE Top-level orchestration.
% Split architecture:
%   1) qprbs13_resample_module.m : CSV -> t0 scan -> y_avg/Y
%   2) qprbs13_symbol_module.m   : Y -> symbols/ES/RLM/x

    cfg = normalize_cfg(cfg);

    res = qprbs13_resample_module(cfg);
    sym = qprbs13_symbol_module(res.Y, cfg);
    rot = qprbs13_rotate_x_module(sym.x, cfg.Dp);
    fitOut = qprbs13_linear_fit_module(res.Y, rot.X, cfg);
    vfOut = qprbs13_steady_state_vf_module(res.Y, rot.X, res.M, cfg);
    noiseOut = qprbs13_level_noise_module(res.Y, sym.sym_code, cfg);
    sndrOut = qprbs13_sndr_module(fitOut.pmax, fitOut.sigmae, noiseOut.sigma_n);

    out = struct();
    out.t0_best = res.t0_best;
    out.t0_candidates = res.t0_candidates;
    out.J_scan = res.J_scan;
    out.ui = res.ui;
    out.M = res.M;
    out.N = res.N;
    out.Ts = res.Ts;
    out.Tpat = res.Tpat;
    out.numCycles = res.numCycles;
    out.y_avg = res.y_avg;
    out.Y = res.Y;
    if isfield(res, 'y_cycles')
        out.y_cycles = res.y_cycles;
    end

    out.s = sym.s;
    out.sym_code = sym.sym_code;
    out.sym_cycle = sym.sym_cycle;
    out.levels = sym.levels;
    out.cluster = sym.cluster;
    out.selected_phase_index = sym.cluster.selected_phase_index;
    out.V = sym.V;
    out.ES1 = sym.ES1;
    out.ES2 = sym.ES2;
    out.RLM = sym.RLM;
    out.counts = sym.counts;
    out.x = sym.x;
    out.xr = rot.xr;
    out.X = rot.X;
    out.Dp = rot.Dp;

    out.Np = fitOut.Np;
    out.TNp = fitOut.TNp;
    out.X1 = fitOut.X1;
    out.P = fitOut.P;
    out.E = fitOut.E;
    out.e = fitOut.e;
    out.sigmae = fitOut.sigmae;
    out.P1 = fitOut.P1;
    out.p = fitOut.p;
    out.pmax = fitOut.pmax;

    out.NpVf = vfOut.NpVf;
    out.TNpVf = vfOut.TNpVf;
    out.p_vf = vfOut.p_vf;
    out.vf = vfOut.vf;

    out.noise = noiseOut;
    out.sigma_L = noiseOut.sigma_levels;
    out.sigma_n = noiseOut.sigma_n;

    out.SNDR = sndrOut.sndr_db;
    out.sndr = sndrOut;

    write_outputs(out, cfg.outputDir);
    plot_outputs(out, cfg.outputDir);
end

function cfg = normalize_cfg(cfg)
    if nargin < 1 || isempty(cfg)
        error('cfg is required.');
    end

    requiredFields = {'csvPath', 'tStableStart', 'outputDir'};
    for i = 1:numel(requiredFields)
        f = requiredFields{i};
        if ~isfield(cfg, f) || isempty(cfg.(f))
            error('cfg.%s is required.', f);
        end
    end

    cfg.csvPath = char(cfg.csvPath);
    cfg.outputDir = char(cfg.outputDir);
    cfg.tStableStart = double(cfg.tStableStart);

    cfg.ui = double(get_cfg(cfg, 'ui', 1/112e9));
    cfg.M = double(get_cfg(cfg, 'M', 32));
    cfg.N = double(get_cfg(cfg, 'N', 8191));
    cfg.numCycles = double(get_cfg(cfg, 'numCycles', 7));
    cfg.scanWindowUI = double(get_cfg(cfg, 'scanWindowUI', 0.5));
    cfg.uiSampleMethod = lower(char(get_cfg(cfg, 'uiSampleMethod', 'center')));
    cfg.clusterMethod = lower(char(get_cfg(cfg, 'clusterMethod', 'kmeans')));
    cfg.saveCycles = logical(get_cfg(cfg, 'saveCycles', false));
    cfg.Dp = double(get_cfg(cfg, 'Dp', 4));
    cfg.Np = double(get_cfg(cfg, 'Np', 29));
    cfg.TNp = double(get_cfg(cfg, 'TNp', 2*cfg.Np + 1));
    cfg.NpVf = double(get_cfg(cfg, 'NpVf', 20));
    cfg.TNpVf = double(get_cfg(cfg, 'TNpVf', 2*cfg.NpVf + 1));
    cfg.minRunLen = double(get_cfg(cfg, 'minRunLen', 6));
    cfg.fixedSampleMethod = lower(char(get_cfg(cfg, 'fixedSampleMethod', 'center')));
    cfg.fixedSampleIndex = double(get_cfg(cfg, 'fixedSampleIndex', round(cfg.M/2)));
    cfg.levelRmsMethod = lower(char(get_cfg(cfg, 'levelRmsMethod', 'merge')));

    validateattributes(cfg.ui, {'numeric'}, {'real','finite','positive','scalar'});
    validateattributes(cfg.M, {'numeric'}, {'real','finite','positive','integer','scalar'});
    validateattributes(cfg.N, {'numeric'}, {'real','finite','positive','integer','scalar'});
    validateattributes(cfg.numCycles, {'numeric'}, {'real','finite','positive','integer','scalar'});
    validateattributes(cfg.tStableStart, {'numeric'}, {'real','finite','scalar'});
    validateattributes(cfg.scanWindowUI, {'numeric'}, {'real','finite','positive','scalar'});
    validateattributes(cfg.Dp, {'numeric'}, {'real','finite','integer','nonnegative','scalar'});
    validateattributes(cfg.Np, {'numeric'}, {'real','finite','integer','positive','scalar'});
    validateattributes(cfg.TNp, {'numeric'}, {'real','finite','integer','positive','scalar'});
    validateattributes(cfg.NpVf, {'numeric'}, {'real','finite','integer','positive','scalar'});
    validateattributes(cfg.TNpVf, {'numeric'}, {'real','finite','integer','positive','scalar'});
    validateattributes(cfg.minRunLen, {'numeric'}, {'real','finite','integer','>=',1,'scalar'});
    validateattributes(cfg.fixedSampleIndex, {'numeric'}, {'real','finite','integer','>=',1,'<=',cfg.M,'scalar'});

    if ~(strcmp(cfg.uiSampleMethod, 'center') || strcmp(cfg.uiSampleMethod, 'mean') || strcmp(cfg.uiSampleMethod, 'proxy_opt'))
        error('cfg.uiSampleMethod must be ''center'', ''mean'' or ''proxy_opt''.');
    end
    if ~(strcmp(cfg.clusterMethod, 'kmeans') || strcmp(cfg.clusterMethod, 'gmm'))
        error('cfg.clusterMethod must be ''kmeans'' or ''gmm''.');
    end
    if cfg.TNp > cfg.N
        error('cfg.TNp must be <= cfg.N.');
    end
    if cfg.TNpVf > cfg.N
        error('cfg.TNpVf must be <= cfg.N.');
    end
    if ~(strcmp(cfg.fixedSampleMethod, 'center') || strcmp(cfg.fixedSampleMethod, 'index'))
        error('cfg.fixedSampleMethod must be ''center'' or ''index''.');
    end
    if ~(strcmp(cfg.levelRmsMethod, 'merge') || strcmp(cfg.levelRmsMethod, 'weighted') || strcmp(cfg.levelRmsMethod, 'mean'))
        error('cfg.levelRmsMethod must be ''merge'', ''weighted'' or ''mean''.');
    end

    Ts = cfg.ui / cfg.M;
    cfg.dtScan = double(get_cfg(cfg, 'dtScan', Ts / 8));
    validateattributes(cfg.dtScan, {'numeric'}, {'real','finite','positive','scalar'});

    if ~exist(cfg.outputDir, 'dir')
        mkdir(cfg.outputDir);
    end
end

function write_outputs(out, outputDir)
    save(fullfile(outputDir, 'resample_out.mat'), 'out');

    MN = out.M * out.N;
    k = (1:MN).';
    ui_index = floor((k - 1) / out.M) + 1;
    sample_in_ui = mod(k - 1, out.M) + 1;
    y = out.y_avg;
    T_y = table(k, ui_index, sample_in_ui, y);
    writetable(T_y, fullfile(outputDir, 'y_avg.csv'));

    symbol_code = (0:3).';
    symbol_value = out.levels;
    count = out.counts;
    mean_voltage = [out.V.m1; out.V.m1_3; out.V.p1_3; out.V.p1];
    T_sym = table(symbol_code, symbol_value, count, mean_voltage);
    writetable(T_sym, fullfile(outputDir, 'symbol_means.csv'));

    n = (1:out.N).';
    sym_code = out.sym_code;
    sym_cycle = out.sym_cycle;
    s = out.s;
    x = out.x;
    T_sx = table(n, sym_code, sym_cycle, s, x);
    writetable(T_sx, fullfile(outputDir, 's_x_symbols.csv'));

    xr = out.xr;
    T_xr = table(n, xr);
    writetable(T_xr, fullfile(outputDir, 'xr.csv'));

    X = out.X; %#ok<NASGU>
    save(fullfile(outputDir, 'X.mat'), 'X');

    linear_fit = struct('Np', out.Np, 'TNp', out.TNp, 'sigmae', out.sigmae, 'pmax', out.pmax, ...
                        'X1', out.X1, 'P', out.P, 'E', out.E, 'P1', out.P1); %#ok<NASGU>
    save(fullfile(outputDir, 'linear_fit.mat'), 'linear_fit');

    k_e = (1:numel(out.e)).';
    e = out.e;
    T_e = table(k_e, e);
    writetable(T_e, fullfile(outputDir, 'e_waveform.csv'));

    k_p = (1:numel(out.p)).';
    p = out.p;
    T_p = table(k_p, p);
    writetable(T_p, fullfile(outputDir, 'p_pulse.csv'));

    k_p_vf = (1:numel(out.p_vf)).';
    p_vf = out.p_vf;
    T_p_vf = table(k_p_vf, p_vf);
    writetable(T_p_vf, fullfile(outputDir, 'p_pulse_vf.csv'));

    steady_state = struct('vf', out.vf, 'NpVf', out.NpVf, 'TNpVf', out.TNpVf); %#ok<NASGU>
    save(fullfile(outputDir, 'steady_state_vf.mat'), 'steady_state');

    level = (0:3).';
    sigma = out.sigma_L;
    runs = out.noise.runs_per_level;
    samples = out.noise.samples_per_level;
    T_noise = table(level, sigma, runs, samples);
    writetable(T_noise, fullfile(outputDir, 'level_noise_rms.csv'));

    noise_summary = struct('sigma_L', out.sigma_L, 'sigma_n', out.sigma_n, ...
                           'm0', out.noise.m0, 'method', out.noise.method, ...
                           'minRunLen', out.noise.minRunLen); %#ok<NASGU>
    save(fullfile(outputDir, 'level_noise_rms.mat'), 'noise_summary');

    sndr_db = out.SNDR;
    pmax = out.sndr.pmax;
    sigmae = out.sndr.sigmae;
    sigman = out.sndr.sigman;
    T_sndr = table(sndr_db, pmax, sigmae, sigman);
    writetable(T_sndr, fullfile(outputDir, 'sndr.csv'));

    sndr = out.sndr; %#ok<NASGU>
    save(fullfile(outputDir, 'sndr.mat'), 'sndr');
end

function plot_outputs(out, outputDir)
    f1 = figure('Visible','off','Color','w');
    plot(out.t0_candidates, out.J_scan, 'b-', 'LineWidth', 1.2);
    hold on;
    [jMin, iMin] = min(out.J_scan);
    plot(out.t0_candidates(iMin), jMin, 'ro', 'MarkerFaceColor', 'r');
    hold off;
    grid on;
    xlabel('t0 (s)'); ylabel('J(t0)'); title('J scan');
    saveas(f1, fullfile(outputDir, 'J_scan.png'));
    close(f1);

    f2 = figure('Visible','off','Color','w');
    imagesc(1:out.N, 1:out.M, out.Y);
    axis xy; colorbar; colormap(turbo);
    xlabel('UI index n'); ylabel('sample in UI m');
    title('heatmap Y (Eq.11-13 reshape)');
    saveas(f2, fullfile(outputDir, 'heatmap_Y.png'));
    close(f2);

    f3 = figure('Visible','off','Color','w');
    histogram(out.s, max(40, round(sqrt(numel(out.s)))), 'Normalization', 'pdf');
    hold on;
    yl = ylim;
    for i = 1:4
        xline(out.cluster.centers_sorted(i), 'r--', 'LineWidth', 1.4);
    end
    ylim(yl);
    hold off;
    grid on;
    xlabel('s(n)'); ylabel('PDF');
    title('Histogram of s with 4 centers');
    saveas(f3, fullfile(outputDir, 's_hist_centers.png'));
    close(f3);

    f4 = figure('Visible','off','Color','w');
    group = categorical(out.sym_code, 0:3, {'code0','code1','code2','code3'});
    boxplot(out.s, group);
    hold on;
    yline(out.V.m1, 'k-', 'V_{-1}');
    yline(out.V.m1_3, 'b-', 'V_{-1/3}');
    yline(out.V.p1_3, 'm-', 'V_{+1/3}');
    yline(out.V.p1, 'r-', 'V_{+1}');
    hold off;
    grid on;
    ylabel('s(n)');
    title('s grouped by sym\_code with symbol means');
    saveas(f4, fullfile(outputDir, 'levels_boxplot.png'));
    close(f4);

    f5 = figure('Visible','off','Color','w');
    nShow = min(500, numel(out.x));
    stem(1:nShow, out.x(1:nShow), 'filled', 'MarkerSize', 3);
    grid on;
    xlabel('n'); ylabel('x(n)');
    title(sprintf('x preview (first %d)', nShow));
    saveas(f5, fullfile(outputDir, 'x_preview.png'));
    close(f5);

    f6 = figure('Visible','off','Color','w');
    nShowE = min(2000, numel(out.e));
    plot(1:nShowE, out.e(1:nShowE), 'LineWidth', 1.0);
    grid on;
    xlabel('k'); ylabel('e(k)');
    title(sprintf('e(k) preview, sigma_e = %.4g', out.sigmae));
    saveas(f6, fullfile(outputDir, 'e_preview.png'));
    close(f6);

    f7 = figure('Visible','off','Color','w');
    nShowP = min(2000, numel(out.p));
    plot(1:nShowP, out.p(1:nShowP), 'LineWidth', 1.0);
    hold on; yline(out.pmax, 'r--', 'pmax'); hold off;
    grid on;
    xlabel('k'); ylabel('p(k)');
    title('Linear fit pulse response p(k) preview');
    saveas(f7, fullfile(outputDir, 'p_preview.png'));
    close(f7);

    f8 = figure('Visible','off','Color','w');
    nShowPv = min(2000, numel(out.p_vf));
    plot(1:nShowPv, out.p_vf(1:nShowPv), 'LineWidth', 1.0);
    grid on;
    xlabel('k'); ylabel('p_{vf}(k)');
    title(sprintf('p_{vf}(k) preview, v_f = %.6g', out.vf));
    saveas(f8, fullfile(outputDir, 'p_vf_preview.png'));
    close(f8);

    f9 = figure('Visible','off','Color','w');
    bar(0:3, out.sigma_L(:));
    grid on;
    xlabel('PAM4 level code'); ylabel('RMS');
    title(sprintf('Level noise RMS, sigma_n = %.6g', out.sigma_n));
    saveas(f9, fullfile(outputDir, 'level_noise_rms.png'));
    close(f9);

    f10 = figure('Visible','off','Color','w');
    axis off;
    txt = sprintf('SNDR = %.4f dB\npmax = %.4g\nsigma_e = %.4g\nsigma_n = %.4g', out.SNDR, out.sndr.pmax, out.sndr.sigmae, out.sndr.sigman);
    text(0.1, 0.6, txt, 'FontSize', 12);
    title('SNDR summary');
    saveas(f10, fullfile(outputDir, 'sndr_summary.png'));
    close(f10);
end

function v = get_cfg(cfg, name, defaultValue)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = defaultValue;
    end
end
