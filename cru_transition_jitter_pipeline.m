function out = cru_transition_jitter_pipeline(t, v, cfg)
%CRU_TRANSITION_JITTER_PIPELINE Transition jitter pipeline with continuous-time Golden PLL.
%   out = cru_transition_jitter_pipeline(t, v, cfg)
%
% Inputs
%   t   : Nx1 double, time in seconds
%   v   : Nx1 double, differential voltage signal
%   cfg : configuration struct
%       Required:
%         cfg.UI
%         cfg.outputDir
%         cfg.context_patterns
%       Optional:
%         cfg.M = 32
%         cfg.balance_mode = 'truncate' | 'random'
%         cfg.rng_seed = 1
%         cfg.edge_dir = 'either' | 'rising' | 'falling'
%         cfg.th_method = 'midpoint'
%         cfg.debug_plot = true
%
% Output fields
%   out.UI, out.fc, out.omega_c, out.t_cross_all, out.e_k_all,
%   out.e_pll, out.e_cru, out.sym_code, out.V_level,
%   out.Si, out.Tavgi, out.S0i, out.fJ, out.JRMS, out.J3u,
%   out.counts_before, out.counts_after
%
% Notes
%   This implementation uses strict continuous-time first-order Golden PLL:
%       e_pll(k) = exp(-omega_c*Delta_t_k)*e_pll(k-1) + ...
%                  (1-exp(-omega_c*Delta_t_k))*e_k_all(k)
%   where Delta_t_k is the actual transition interval between crossings.

% ---------- Input validation ----------
validateattributes(t, {'double'}, {'column', 'real', 'finite', 'nonempty'}, mfilename, 't', 1);
validateattributes(v, {'double'}, {'column', 'real', 'finite', 'numel', numel(t)}, mfilename, 'v', 2);
if ~isstruct(cfg)
    error('cfg must be a struct.');
end
required_fields = {'UI', 'outputDir', 'context_patterns'};
for kf = 1:numel(required_fields)
    if ~isfield(cfg, required_fields{kf})
        error('cfg.%s is required.', required_fields{kf});
    end
end

cfg = apply_defaults(cfg);
validateattributes(cfg.UI, {'double'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'cfg.UI');

if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

% ---------- Step A: UI segmentation and level restoration ----------
UI = double(cfg.UI);
M  = double(cfg.M);
Ts = UI / M;

t_uniform = (t(1):Ts:t(end)).';
v_uniform = interp1(t, v, t_uniform, 'linear');

N = floor(numel(t_uniform) / M);
if N < 2
    error('Not enough samples after resampling to form at least 2 UIs.');
end

Y = reshape(v_uniform(1:N*M), M, N);
m0 = round(M / 2);
s = Y(m0, :).';

[idx, centers] = kmeans(s, 4, 'Replicates', 8, 'MaxIter', 1000);
[sorted_centers, order] = sort(centers(:), 'ascend');
% Map cluster IDs into symbol codes {0,1,2,3} in ascending voltage order.
sym_code = zeros(N, 1, 'double');
for i = 1:4
    sym_code(idx == order(i)) = i - 1;
end

V_level = zeros(4,1,'double');
for i = 0:3
    members = s(sym_code == i);
    if isempty(members)
        error('No samples mapped to symbol level %d.', i);
    end
    V_level(i+1) = mean(members, 'double');
end

% ---------- Step B: crossing extraction with transition-specific threshold ----------
t_cross_all = zeros(N-1, 1, 'double');
e_k_all     = zeros(N-1, 1, 'double');
trans_n     = zeros(N-1, 1, 'double');
trans_A     = zeros(N-1, 1, 'double');
trans_B     = zeros(N-1, 1, 'double');

k_valid = 0;
for n = 2:N
    A = sym_code(n-1);
    B = sym_code(n);
    if A == B
        continue;
    end

    if ~edge_direction_match(A, B, cfg.edge_dir)
        continue;
    end

    if ~strcmpi(cfg.th_method, 'midpoint')
        error('Unsupported cfg.th_method: %s. Only ''midpoint'' is implemented.', cfg.th_method);
    end

    Vth = (V_level(A+1) + V_level(B+1)) / 2.0;

    x = Y(:, n);
    t_ui0 = t_uniform((n-1)*M + 1);
    t_local = t_ui0 + (0:M-1).' * Ts;

    cross_idx = find((x(1:end-1) - Vth) .* (x(2:end) - Vth) <= 0, 1, 'first');
    if isempty(cross_idx)
        continue;
    end

    v1 = x(cross_idx);
    v2 = x(cross_idx + 1);
    t1 = t_local(cross_idx);
    t2 = t_local(cross_idx + 1);

    if v2 == v1
        t_cross = (t1 + t2) / 2.0;
    else
        alpha = (Vth - v1) / (v2 - v1);
        t_cross = t1 + alpha * (t2 - t1);
    end

    t_ideal = (n - 1) * UI;
    e_k = t_cross - t_ideal;

    k_valid = k_valid + 1;
    t_cross_all(k_valid) = t_cross;
    e_k_all(k_valid) = e_k;
    trans_n(k_valid) = n;
    trans_A(k_valid) = A;
    trans_B(k_valid) = B;
end

if k_valid < 2
    error('Not enough valid crossings to run CRU (need >=2, got %d).', k_valid);
end

t_cross_all = t_cross_all(1:k_valid);
e_k_all = e_k_all(1:k_valid);
trans_n = trans_n(1:k_valid);
trans_A = trans_A(1:k_valid);
trans_B = trans_B(1:k_valid);

% ---------- Step C: optimized continuous-time Golden PLL (strict Delta_t_k) ----------
fb = 1.0 / UI;
fc = fb / 13280.0;
omega_c = 2.0 * pi * fc;

K = numel(e_k_all);
e_pll = zeros(K,1,'double');
e_pll(1) = e_k_all(1);

for k = 2:K
    Delta_t_k = t_cross_all(k) - t_cross_all(k-1);
    a_k = exp(-omega_c * Delta_t_k);
    e_pll(k) = a_k * e_pll(k-1) + (1.0 - a_k) * e_k_all(k);
end

e_cru = e_k_all - e_pll;

% ---------- Step D: 12-class transition classification ----------
class_id = zeros(K,1,'double');
for k = 1:K
    n = trans_n(k);
    context.prev2 = get_sym(sym_code, n-2);
    context.prev1 = get_sym(sym_code, n-1);
    context.curr  = get_sym(sym_code, n);
    context.next1 = get_sym(sym_code, n+1);

    cid = identify_class(context, cfg.context_patterns);
    if cid >= 1 && cid <= 12
        class_id(k) = cid;
    end
end

valid_cls = class_id >= 1 & class_id <= 12;
class_id = class_id(valid_cls);
e_cru_cls = e_cru(valid_cls);

Si = cell(12,1);
for i = 1:12
    Si{i} = e_cru_cls(class_id == i);
end

% ---------- Step E: force balanced sample counts ----------
counts_before = zeros(12,1,'double');
for i = 1:12
    counts_before(i) = numel(Si{i});
end

if any(counts_before == 0)
    error('At least one class has zero samples before balancing. counts_before=%s', mat2str(counts_before.'));
end

Nmin = min(counts_before);
if strcmpi(cfg.balance_mode, 'truncate')
    for i = 1:12
        Si{i} = Si{i}(1:Nmin);
    end
elseif strcmpi(cfg.balance_mode, 'random')
    rng(cfg.rng_seed);
    for i = 1:12
        idx_sel = randperm(numel(Si{i}), Nmin);
        Si{i} = Si{i}(idx_sel);
    end
else
    error('Unsupported cfg.balance_mode: %s', cfg.balance_mode);
end

counts_after = Nmin * ones(12,1,'double');

% ---------- Step F: statistics ----------
Tavgi = zeros(12,1,'double');
S0i = cell(12,1);
for i = 1:12
    Tavgi(i) = mean(Si{i}, 'double');
    S0i{i} = Si{i} - Tavgi(i);
end

fJ = vertcat(S0i{:});
JRMS = std(fJ, 0, 'all');
J3u = prctile(fJ, 99.95) - prctile(fJ, 0.05);

% ---------- Package outputs ----------
out = struct();
out.UI = UI;
out.fc = fc;
out.omega_c = omega_c;
out.t_cross_all = t_cross_all;
out.e_k_all = e_k_all;
out.e_pll = e_pll;
out.e_cru = e_cru;
out.sym_code = sym_code;
out.V_level = V_level;
out.Si = Si;
out.Tavgi = Tavgi;
out.S0i = S0i;
out.fJ = fJ;
out.JRMS = JRMS;
out.J3u = J3u;
out.counts_before = counts_before;
out.counts_after = counts_after;

% ---------- Debug plots ----------
if cfg.debug_plot
    make_debug_plots(cfg.outputDir, e_k_all, e_pll, e_cru, fJ, counts_before, counts_after);
end

end

% ============================== Local functions ==============================

function cfg = apply_defaults(cfg)
if ~isfield(cfg, 'M'), cfg.M = 32; end
if ~isfield(cfg, 'balance_mode'), cfg.balance_mode = 'truncate'; end
if ~isfield(cfg, 'rng_seed'), cfg.rng_seed = 1; end
if ~isfield(cfg, 'edge_dir'), cfg.edge_dir = 'either'; end
if ~isfield(cfg, 'th_method'), cfg.th_method = 'midpoint'; end
if ~isfield(cfg, 'debug_plot'), cfg.debug_plot = true; end

cfg.M = double(cfg.M);
cfg.rng_seed = double(cfg.rng_seed);
end

function tf = edge_direction_match(A, B, edge_dir)
switch lower(edge_dir)
    case 'either'
        tf = true;
    case 'rising'
        tf = B > A;
    case 'falling'
        tf = B < A;
    otherwise
        error('Unsupported cfg.edge_dir: %s', edge_dir);
end
end

function s = get_sym(sym_code, idx)
if idx < 1 || idx > numel(sym_code)
    s = NaN;
else
    s = sym_code(idx);
end
end

function cid = identify_class(context, patterns)
% Identify class id in [1..12] using cfg.context_patterns.
% Supported pattern forms per class:
%   - function_handle: fn(context) -> logical
%   - numeric row/vector with wildcard NaN/-1:
%       [A B]           where A=prev1, B=curr
%       [P A B]         where P=prev2, A=prev1, B=curr
%       [P A B N]       where N=next1
%   - numeric matrix: each row is an alternative rule
%   - struct with fields among {prev2,prev1,curr,next1}; missing fields ignored

cid = 0;
if ~iscell(patterns) || numel(patterns) ~= 12
    error('cfg.context_patterns must be a 12-element cell array.');
end

for i = 1:12
    p = patterns{i};
    if match_pattern(context, p)
        cid = i;
        return;
    end
end
end

function tf = match_pattern(context, p)
if isa(p, 'function_handle')
    tf = logical(p(context));
    return;
end

if isnumeric(p)
    if isvector(p)
        tf = match_numeric_rule(context, p(:).');
        return;
    end
    tf = false;
    for r = 1:size(p,1)
        if match_numeric_rule(context, p(r,:))
            tf = true;
            return;
        end
    end
    return;
end

if isstruct(p)
    tf = true;
    fns = fieldnames(p);
    for k = 1:numel(fns)
        fn = fns{k};
        if ~isfield(context, fn)
            tf = false;
            return;
        end
        if ~is_wildcard(p.(fn)) && context.(fn) ~= p.(fn)
            tf = false;
            return;
        end
    end
    return;
end

error('Unsupported context pattern type: %s', class(p));
end

function tf = match_numeric_rule(context, row)
vals = row(~isnan(row)); %#ok<NASGU>
L = numel(row);
switch L
    case 2
        rule = [context.prev1, context.curr];
    case 3
        rule = [context.prev2, context.prev1, context.curr];
    case 4
        rule = [context.prev2, context.prev1, context.curr, context.next1];
    otherwise
        error('Numeric pattern rows must be length 2, 3, or 4.');
end

tf = true;
for i = 1:L
    if ~is_wildcard(row(i)) && rule(i) ~= row(i)
        tf = false;
        return;
    end
end
end

function tf = is_wildcard(x)
tf = isnan(x) || isequal(x, -1);
end

function make_debug_plots(outputDir, e_k_all, e_pll, e_cru, fJ, counts_before, counts_after)
fig1 = figure('Visible', 'off');
plot(e_k_all, 'LineWidth', 1); hold on;
plot(e_pll, 'LineWidth', 1);
plot(e_cru, 'LineWidth', 1);
grid on;
xlabel('Transition index k');
ylabel('Time error (s)');
legend({'e\_k\_all','e\_pll','e\_cru'}, 'Location', 'best');
title('CRU debug: raw, PLL track, and high-pass residual');
saveas(fig1, fullfile(outputDir, 'cru_debug.png'));
close(fig1);

fig2 = figure('Visible', 'off');
histogram(fJ, 100);
grid on;
xlabel('Jitter (s)');
ylabel('Count');
title('Histogram of fJ');
saveas(fig2, fullfile(outputDir, 'histogram_fJ.png'));
close(fig2);

fig3 = figure('Visible', 'off');
b = bar([(1:12).', counts_before, counts_after], 'grouped'); %#ok<NASGU>
grid on;
xlabel('Class id');
ylabel('Sample count');
legend({'before', 'after'}, 'Location', 'best');
title('Class counts before/after balancing');
saveas(fig3, fullfile(outputDir, 'counts_before_after.png'));
close(fig3);
end

%% ============================== Demo example ==============================
%{
% Demo usage (replace with your own t/v data):
%
% t = (0:1e-12:20e-9).';
% UI = 100e-12;
% v = 0.4 * square(2*pi*(1/UI/2)*t);   % Example only
%
% cfg = struct();
% cfg.UI = UI;
% cfg.outputDir = './output_demo';
% cfg.context_patterns = {
%     [0 1], [1 0], [1 2], [2 1], [2 3], [3 2], ...
%     [0 2], [2 0], [1 3], [3 1], [0 3], [3 0] ...
% };
% cfg.M = 32;
% cfg.balance_mode = 'truncate'; % or 'random'
% cfg.rng_seed = 1;
% cfg.edge_dir = 'either';
% cfg.th_method = 'midpoint';
% cfg.debug_plot = true;
%
% out = cru_transition_jitter_pipeline(t, v, cfg);
% disp(out.JRMS);
% disp(out.J3u);
%}
