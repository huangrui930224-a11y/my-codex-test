function out = cru_transition_jitter_pipeline(t, v, cfg)
%CRU_TRANSITION_JITTER_PIPELINE CRU + 12-class transition jitter + EOJ pipeline.
%   out = cru_transition_jitter_pipeline(t, v, cfg)
%
% Required cfg fields:
%   cfg.UI, cfg.outputDir, cfg.context_patterns
% Optional cfg fields:
%   cfg.M=32, cfg.balance_mode='truncate'|'random', cfg.rng_seed=1,
%   cfg.edge_dir='either'|'rising'|'falling', cfg.th_method='midpoint',
%   cfg.debug_plot=true, cfg.Npat=8191 (or 511), cfg.eoj_mode='max'|'first3',
%   cfg.eoj_debug_plot=false, cfg.t0_override=[]

% ---------- Input validation ----------
validateattributes(t, {'double'}, {'column','real','finite','nonempty'}, mfilename, 't', 1);
validateattributes(v, {'double'}, {'column','real','finite','numel',numel(t)}, mfilename, 'v', 2);
if ~isstruct(cfg), error('cfg must be a struct.'); end
required_fields = {'UI','outputDir','context_patterns'};
for i = 1:numel(required_fields)
    if ~isfield(cfg, required_fields{i}), error('cfg.%s is required.', required_fields{i}); end
end

cfg = apply_defaults(cfg);
validateattributes(cfg.UI, {'double'}, {'scalar','real','finite','positive'}, mfilename, 'cfg.UI');
if ~exist(cfg.outputDir, 'dir'), mkdir(cfg.outputDir); end

% ---------- Step A: UI segmentation and level restoration ----------
UI = double(cfg.UI);
M  = double(cfg.M);
Ts = UI / M;

t_uniform = (t(1):Ts:t(end)).';
v_uniform = interp1(t, v, t_uniform, 'linear');
N = floor(numel(t_uniform)/M);
if N < 2, error('Not enough samples after resampling.'); end

Y = reshape(v_uniform(1:N*M), M, N);
m0 = round(M/2);
s = Y(m0,:).';

[idx, centers] = kmeans(s, 4, 'Replicates', 8, 'MaxIter', 1000);
[~, order] = sort(centers(:), 'ascend');
sym_code = zeros(N,1,'double');
for i = 1:4
    sym_code(idx == order(i)) = i - 1;
end

V_level = zeros(4,1,'double');
for i = 0:3
    members = s(sym_code == i);
    if isempty(members), error('No samples in symbol level %d.', i); end
    V_level(i+1) = mean(members, 'double');
end

% ---------- Step B: crossing extraction (transition-specific threshold) ----------
t_cross_all = zeros(N-1,1,'double');
e_k_all     = zeros(N-1,1,'double');
trans_n     = zeros(N-1,1,'double');

k_valid = 0;
for n = 2:N
    A = sym_code(n-1); B = sym_code(n);
    if A == B, continue; end
    if ~edge_direction_match(A, B, cfg.edge_dir), continue; end
    if ~strcmpi(cfg.th_method, 'midpoint')
        error('Unsupported cfg.th_method: %s. Only ''midpoint'' is implemented.', cfg.th_method);
    end

    Vth = (V_level(A+1) + V_level(B+1)) / 2.0;
    x = Y(:,n);
    t_ui0 = t_uniform((n-1)*M + 1);
    t_local = t_ui0 + (0:M-1).' * Ts;

    cross_idx = find((x(1:end-1)-Vth).*(x(2:end)-Vth) <= 0, 1, 'first');
    if isempty(cross_idx), continue; end

    v1 = x(cross_idx); v2 = x(cross_idx+1);
    t1 = t_local(cross_idx); t2 = t_local(cross_idx+1);
    if v2 == v1
        t_cross = (t1+t2)/2.0;
    else
        alpha = (Vth-v1)/(v2-v1);
        t_cross = t1 + alpha*(t2-t1);
    end

    t_ideal = (n-1) * UI;
    e_k = t_cross - t_ideal;

    k_valid = k_valid + 1;
    t_cross_all(k_valid) = t_cross;
    e_k_all(k_valid) = e_k;
    trans_n(k_valid) = n;
end

if k_valid < 2, error('Not enough valid crossings to run CRU.'); end
t_cross_all = t_cross_all(1:k_valid);
e_k_all = e_k_all(1:k_valid);
trans_n = trans_n(1:k_valid);

% ---------- Step C: strict continuous-time Golden PLL using Delta_t_k ----------
fb = 1.0/UI;
fc = fb/13280.0;
omega_c = 2*pi*fc;

K = numel(e_k_all);
e_pll = zeros(K,1,'double');
e_pll(1) = e_k_all(1);
for k = 2:K
    Delta_t_k = t_cross_all(k) - t_cross_all(k-1);
    a_k = exp(-omega_c * Delta_t_k);
    e_pll(k) = a_k*e_pll(k-1) + (1-a_k)*e_k_all(k);
end
e_cru = e_k_all - e_pll;

% ---------- Step D: classify 12 transition classes ----------
class_id_event = zeros(K,1,'double');
for k = 1:K
    n = trans_n(k);
    context.prev2 = get_sym(sym_code, n-2);
    context.prev1 = get_sym(sym_code, n-1);
    context.curr  = get_sym(sym_code, n);
    context.next1 = get_sym(sym_code, n+1);
    cid = identify_class(context, cfg.context_patterns);
    if cid >= 1 && cid <= 12
        class_id_event(k) = cid;
    end
end

valid_cls = class_id_event >= 1 & class_id_event <= 12;
Si = cell(12,1);
for i = 1:12
    Si{i} = e_cru(class_id_event == i);
end

% ---------- Step E: balance sample counts ----------
counts_before = zeros(12,1,'double');
for i = 1:12, counts_before(i) = numel(Si{i}); end
if any(counts_before == 0)
    error('At least one class has zero samples before balancing. counts_before=%s', mat2str(counts_before.'));
end

Nmin = min(counts_before);
if strcmpi(cfg.balance_mode, 'truncate')
    for i = 1:12, Si{i} = Si{i}(1:Nmin); end
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

% ---------- Output packing ----------
out = struct();
out.UI = UI;
out.Npat = double(cfg.Npat);
% UI#0 start anchor. If cfg.t0_override is provided, it explicitly defines
% the first trigger/repeat anchor (practice-B controlled alignment).
if ~isempty(cfg.t0_override)
    out.t0 = double(cfg.t0_override);
else
    out.t0 = double(t_uniform(1));
end
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

% Event-level outputs required by EOJ
out.t_cross_event = t_cross_all;          % Kx1 crossing time (s)
out.ui_index_event = trans_n - 1;         % Kx1 global UI boundary index
out.class_id_event = class_id_event;      % Kx1, 0 for not in class 1..12

% ---------- EOJ sub-module ----------
out = compute_eoj(out, cfg);

% ---------- Debug plots ----------
if cfg.debug_plot
    make_debug_plots(cfg.outputDir, e_k_all, e_pll, e_cru, fJ, counts_before, counts_after);
end

end

% ============================== Local functions ==============================

function cfg = apply_defaults(cfg)
if ~isfield(cfg,'M'), cfg.M = 32; end
if ~isfield(cfg,'balance_mode'), cfg.balance_mode = 'truncate'; end
if ~isfield(cfg,'rng_seed'), cfg.rng_seed = 1; end
if ~isfield(cfg,'edge_dir'), cfg.edge_dir = 'either'; end
if ~isfield(cfg,'th_method'), cfg.th_method = 'midpoint'; end
if ~isfield(cfg,'debug_plot'), cfg.debug_plot = true; end
if ~isfield(cfg,'Npat'), cfg.Npat = 8191; end
if ~isfield(cfg,'eoj_mode'), cfg.eoj_mode = 'max'; end
if ~isfield(cfg,'eoj_debug_plot'), cfg.eoj_debug_plot = false; end
if ~isfield(cfg,'t0_override'), cfg.t0_override = []; end
cfg.M = double(cfg.M);
cfg.rng_seed = double(cfg.rng_seed);
cfg.Npat = double(cfg.Npat);
if ~isempty(cfg.t0_override)
    cfg.t0_override = double(cfg.t0_override);
    validateattributes(cfg.t0_override, {'double'}, {'scalar','real','finite'}, mfilename, 'cfg.t0_override');
end
end

function tf = edge_direction_match(A, B, edge_dir)
switch lower(edge_dir)
    case 'either',  tf = true;
    case 'rising',  tf = B > A;
    case 'falling', tf = B < A;
    otherwise, error('Unsupported cfg.edge_dir: %s', edge_dir);
end
end

function s = get_sym(sym_code, idx)
if idx < 1 || idx > numel(sym_code), s = NaN; else, s = sym_code(idx); end
end

function cid = identify_class(context, patterns)
cid = 0;
if ~iscell(patterns) || numel(patterns) ~= 12
    error('cfg.context_patterns must be a 12-element cell array.');
end
for i = 1:12
    p = patterns{i};
    if match_pattern(context, p), cid = i; return; end
end
end

function tf = match_pattern(context, p)
if isa(p,'function_handle')
    tf = logical(p(context)); return;
end
if isnumeric(p)
    if isvector(p), tf = match_numeric_rule(context, p(:).'); return; end
    tf = false;
    for r = 1:size(p,1)
        if match_numeric_rule(context, p(r,:)), tf = true; return; end
    end
    return;
end
if isstruct(p)
    tf = true; fns = fieldnames(p);
    for k = 1:numel(fns)
        fn = fns{k};
        if ~isfield(context, fn), tf = false; return; end
        if ~is_wildcard(p.(fn)) && context.(fn) ~= p.(fn), tf = false; return; end
    end
    return;
end
error('Unsupported context pattern type: %s', class(p));
end

function tf = match_numeric_rule(context, row)
L = numel(row);
switch L
    case 2, rule = [context.prev1, context.curr];
    case 3, rule = [context.prev2, context.prev1, context.curr];
    case 4, rule = [context.prev2, context.prev1, context.curr, context.next1];
    otherwise, error('Numeric pattern rows must be length 2, 3, or 4.');
end
tf = true;
for i = 1:L
    if ~is_wildcard(row(i)) && rule(i) ~= row(i), tf = false; return; end
end
end

function tf = is_wildcard(x)
tf = isnan(x) || isequal(x,-1);
end

function out = compute_eoj(out, cfg)
%COMPUTE_EOJ EOJ statistics using repeat-local UI domain to avoid big-number subtraction.
% T1~T4 protocol mapping for class i and repeat window (r,r+1,r+2):
%   T1 = Tr(i,r),   T2 = Tr(i,r+1)   => P2 = T2 - T1
%   T3 = Tr(i,r+1), T4 = Tr(i,r+2)   => P3 = T4 - T3
%   EOJ_i_win = abs(P2 - P3), EOJ_i = max over windows (or first window for first3)

UI = double(out.UI);
if isfield(out, 'Npat')
    Npat = double(out.Npat);
elseif isfield(cfg, 'Npat')
    Npat = double(cfg.Npat);
else
    error('Npat not found in out/cfg.');
end
if Npat <= 0, error('Npat must be positive.'); end

% Required event-level fields
for f = {'t_cross_event','ui_index_event','class_id_event'}
    if ~isfield(out, f{1}), error('out.%s is required for EOJ.', f{1}); end
end

t_cross = double(out.t_cross_event(:));
ui_idx = double(out.ui_index_event(:));
cid = double(out.class_id_event(:));
valid = cid >= 1 & cid <= 12;

t_cross = t_cross(valid);
ui_idx = ui_idx(valid);
cid = cid(valid);

repeat_id = floor(ui_idx ./ Npat) + 1;

% Repeat start time anchor.
if isfield(out,'t0')
    t0 = double(out.t0);
else
    % Fallback: align by earliest repeat using event-time/UI-index relation.
    % This keeps a consistent t0 for all repeats.
    [~, imin] = min(t_cross);
    t0 = t_cross(imin) - ui_idx(imin)*UI;
end

t_rep0 = t0 + (repeat_id - 1) * Npat * UI;
% IMPORTANT: all EOJ math is performed in UI units (not large second values).
t_rel_UI = (t_cross - t_rep0) ./ UI;

EOJ_per_class_UI = nan(12,1);
EOJ_per_class_s  = nan(12,1);
num_windows_used = zeros(12,1,'double');
notes = strings(12,1);
Tr_by_class = cell(12,1);
rep_by_class = cell(12,1);

for i = 1:12
    m = (cid == i);
    rep_i = repeat_id(m);
    trel_i = t_rel_UI(m);
    if isempty(rep_i)
        notes(i) = "no events";
        continue;
    end

    [rep_unique, ~, g] = unique(rep_i);
    Tr = accumarray(g, trel_i, [], @mean);

    % Require consecutive repeats: (r,r+1,r+2)
    eoj_wins = [];
    for k = 1:(numel(rep_unique)-2)
        r1 = rep_unique(k); r2 = rep_unique(k+1); r3 = rep_unique(k+2);
        if (r2 == r1+1) && (r3 == r2+1)
            % Protocol-consistent definitions:
            % T1=Tr(r), T2=Tr(r+1), T3=Tr(r+1), T4=Tr(r+2)
            T1 = Tr(k);
            T2 = Tr(k+1);
            T3 = Tr(k+1);
            T4 = Tr(k+2);
            P2_UI = T2 - T1;
            P3_UI = T4 - T3;
            eoj_wins(end+1,1) = abs(P2_UI - P3_UI); %#ok<AGROW>
        end
    end

    Tr_by_class{i} = Tr;
    rep_by_class{i} = rep_unique;

    if isempty(eoj_wins)
        notes(i) = "no consecutive 3-repeat window";
        continue;
    end

    if strcmpi(cfg.eoj_mode, 'first3')
        EOJ_per_class_UI(i) = eoj_wins(1);
        notes(i) = "first3";
    else
        EOJ_per_class_UI(i) = max(eoj_wins);
        notes(i) = "max";
    end
    EOJ_per_class_s(i) = EOJ_per_class_UI(i) * UI;
    num_windows_used(i) = numel(eoj_wins);

    if cfg.eoj_debug_plot
        fig = figure('Visible','off');
        subplot(2,1,1);
        plot(rep_unique, Tr, '-o'); grid on;
        xlabel('repeat id'); ylabel('Tr(i,r) in UI');
        title(sprintf('Class %d repeat-mean time (UI domain)', i));
        subplot(2,1,2);
        plot(eoj_wins, '-o'); grid on;
        xlabel('window index'); ylabel('|P2-P3| (UI)');
        title(sprintf('Class %d EOJ windows', i));
        saveas(fig, fullfile(cfg.outputDir, sprintf('eoj_debug_class%d.png', i)));
        close(fig);
    end
end

out.EOJ_per_class_UI = EOJ_per_class_UI;
out.EOJ_UI = max(EOJ_per_class_UI, [], 'omitnan');
out.EOJ_per_class_s = EOJ_per_class_s;
out.EOJ_s = out.EOJ_UI * UI;
out.EOJ_num_windows = num_windows_used;
out.EOJ_notes = notes;
out.EOJ_Tr_by_class = Tr_by_class;
out.EOJ_repeat_ids_by_class = rep_by_class;

% Anchor consistency diagnostic (practice-B correctness check):
% delta_ui = (t_cross - (t0 + ui_idx*UI))/UI should stay bounded and centered.
delta_ui = (t_cross - (t0 + ui_idx .* UI)) ./ UI;
out.EOJ_anchor_delta_ui_mean = mean(delta_ui, 'omitnan');
out.EOJ_anchor_delta_ui_std = std(delta_ui, 0, 'omitnan');
out.EOJ_anchor_delta_ui_maxabs = max(abs(delta_ui), [], 'omitnan');
if out.EOJ_anchor_delta_ui_maxabs > 0.5
    warning('EOJ anchor check: max abs delta_ui = %.4f UI (>0.5 UI). Check cfg.t0_override alignment.', ...
        out.EOJ_anchor_delta_ui_maxabs);
end

% Save CSV
class_id = (1:12).';
T = table(class_id, EOJ_per_class_UI, EOJ_per_class_s, num_windows_used, notes, ...
    'VariableNames', {'class_id','EOJ_UI','EOJ_s','num_windows_used','note'});
writetable(T, fullfile(cfg.outputDir, 'eoj_per_class.csv'));

% Save bar chart
figb = figure('Visible','off');
bar(class_id, EOJ_per_class_UI);
grid on; xlabel('class id'); ylabel('EOJ (UI)');
title('EOJ per class (UI domain)');
saveas(figb, fullfile(cfg.outputDir, 'eoj_bar.png'));
close(figb);
end

function make_debug_plots(outputDir, e_k_all, e_pll, e_cru, fJ, counts_before, counts_after)
fig1 = figure('Visible','off');
plot(e_k_all, 'LineWidth', 1); hold on;
plot(e_pll, 'LineWidth', 1);
plot(e_cru, 'LineWidth', 1);
grid on;
xlabel('Transition index k'); ylabel('Time error (s)');
legend({'e\_k\_all','e\_pll','e\_cru'}, 'Location', 'best');
title('CRU debug: raw, PLL track, and high-pass residual');
saveas(fig1, fullfile(outputDir, 'cru_debug.png')); close(fig1);

fig2 = figure('Visible','off');
histogram(fJ, 100);
grid on; xlabel('Jitter (s)'); ylabel('Count');
title('Histogram of fJ');
saveas(fig2, fullfile(outputDir, 'histogram_fJ.png')); close(fig2);

fig3 = figure('Visible','off');
bar([(1:12).', counts_before, counts_after], 'grouped');
grid on; xlabel('Class id'); ylabel('Sample count');
legend({'before','after'}, 'Location', 'best');
title('Class counts before/after balancing');
saveas(fig3, fullfile(outputDir, 'counts_before_after.png')); close(fig3);
end

%% ============================== Demo example ==============================
%{
% Demo A: normal pipeline usage (replace t/v with measured waveform)
% t = (0:1e-12:200e-9).';
% UI = 100e-12;
% v = 0.4 * square(2*pi*(1/UI/2)*t);
% cfg = struct();
% cfg.UI = UI;
% cfg.outputDir = './output_demo';
% cfg.context_patterns = {
%   [0 1],[1 0],[1 2],[2 1],[2 3],[3 2],[0 2],[2 0],[1 3],[3 1],[0 3],[3 0]
% };
% cfg.M = 32; cfg.balance_mode='truncate'; cfg.rng_seed=1;
% cfg.edge_dir='either'; cfg.th_method='midpoint'; cfg.debug_plot=true;
% cfg.Npat = 8191; cfg.eoj_mode='max'; cfg.eoj_debug_plot=false;
% out = cru_transition_jitter_pipeline(t, v, cfg);
% disp(out.EOJ_UI); disp(out.EOJ_s);
%
% Demo B: synthetic EOJ-chain check directly via compute_eoj input fields.
% out2 = struct(); out2.UI = UI; out2.Npat = 8191; out2.t0 = 0;
% K = 240; out2.ui_index_event = (0:K-1).';
% rep = floor(out2.ui_index_event/out2.Npat)+1;
% out2.class_id_event = mod((1:K).',12)+1;
% base = out2.ui_index_event*UI;
% drift = 0.01*UI*rep;         % slow drift
% wiggle = 0.001*UI*sin(2*pi*rep/5); % repeat-dependent perturbation -> EOJ
% out2.t_cross_event = base + drift + wiggle;
% cfg2 = struct('outputDir','./output_demo','eoj_mode','max','eoj_debug_plot',false,'Npat',8191);
% out2 = compute_eoj(out2, cfg2);
% disp(out2.EOJ_per_class_UI);
%}
