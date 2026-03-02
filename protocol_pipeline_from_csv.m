function out = protocol_pipeline_from_csv(csvPath, cfg)
%PROTOCOL_PIPELINE_FROM_CSV Protocol-aligned CRU transition jitter + EOJ core.
% This helper implements the original protocol-style flow:
% 1) UI resampling and PAM4 level recovery
% 2) per-transition threshold crossing extraction
% 3) continuous-time Golden PLL using Delta_t_k
% 4) 12-class classification + balancing
% 5) JRMS/J3u from balanced fJ
% 6) EOJ per class from repeat-local UI domain

cfg = apply_defaults(cfg);
[t, v] = read_tv_from_csv(csvPath, cfg);
validateattributes(t, {'double'}, {'column','real','finite','nonempty'});
validateattributes(v, {'double'}, {'column','real','finite','numel',numel(t)});

UI = double(cfg.UI);
M = double(cfg.M);
Ts = UI / M;

t_uniform = (t(1):Ts:t(end)).';
v_uniform = interp1(t, v, t_uniform, 'linear');
N = floor(numel(t_uniform)/M);
if N < 2, error('Not enough samples after resampling.'); end

Y = reshape(v_uniform(1:N*M), M, N);
m0 = round(M/2);
s = Y(m0,:).';

[idx, c] = kmeans(s, 4, 'Replicates', 8, 'MaxIter', 1000);
[~, ord] = sort(c(:), 'ascend');
sym_code = zeros(N,1,'double');
for i = 1:4
    sym_code(idx == ord(i)) = i - 1;
end

V_level = zeros(4,1,'double');
for i = 0:3
    members = s(sym_code == i);
    if isempty(members), error('No samples for symbol level %d.', i); end
    V_level(i+1) = mean(members, 'double');
end

% crossing extraction
maxK = N-1;
t_cross_all = zeros(maxK,1,'double');
e_k_all = zeros(maxK,1,'double');
trans_n = zeros(maxK,1,'double');
k_valid = 0;
for n = 2:N
    A = sym_code(n-1); B = sym_code(n);
    if A == B, continue; end
    if ~edge_match(A, B, cfg.edge_dir), continue; end
    Vth = (V_level(A+1) + V_level(B+1)) / 2.0;

    x = Y(:,n);
    t_ui0 = t_uniform((n-1)*M + 1);
    t_local = t_ui0 + (0:M-1).' * Ts;
    i0 = find((x(1:end-1)-Vth).*(x(2:end)-Vth) <= 0, 1, 'first');
    if isempty(i0), continue; end

    v1 = x(i0); v2 = x(i0+1);
    t1 = t_local(i0); t2 = t_local(i0+1);
    if v2 == v1
        t_cross = (t1+t2)/2;
    else
        alpha = (Vth-v1)/(v2-v1);
        t_cross = t1 + alpha*(t2-t1);
    end

    k_valid = k_valid + 1;
    t_cross_all(k_valid) = t_cross;
    e_k_all(k_valid) = t_cross - (n-1)*UI;
    trans_n(k_valid) = n;
end
if k_valid < 2, error('Not enough valid crossings for CRU.'); end

t_cross_all = t_cross_all(1:k_valid);
e_k_all = e_k_all(1:k_valid);
trans_n = trans_n(1:k_valid);

% continuous-time golden PLL with Delta_t
fb = 1/UI;
fc = fb/13280;
omega_c = 2*pi*fc;
K = numel(e_k_all);
e_pll = zeros(K,1,'double');
e_pll(1) = e_k_all(1);
for k = 2:K
    Delta_t_k = t_cross_all(k) - t_cross_all(k-1);
    a_k = exp(-omega_c*Delta_t_k);
    e_pll(k) = a_k*e_pll(k-1) + (1-a_k)*e_k_all(k);
end
e_cru = e_k_all - e_pll;

% classify 12 classes
class_id_event = zeros(K,1,'double');
for k = 1:K
    n = trans_n(k);
    ctx.prev2 = get_sym(sym_code, n-2);
    ctx.prev1 = get_sym(sym_code, n-1);
    ctx.curr = get_sym(sym_code, n);
    ctx.next1 = get_sym(sym_code, n+1);
    cid = identify_class(ctx, cfg.context_patterns);
    if cid >= 1 && cid <= 12
        class_id_event(k) = cid;
    end
end

Si = cell(12,1);
counts_before = zeros(12,1,'double');
for i = 1:12
    Si{i} = e_cru(class_id_event == i);
    counts_before(i) = numel(Si{i});
end
if any(counts_before == 0)
    error('At least one class has zero samples before balancing.');
end
Nmin = min(counts_before);
if strcmpi(cfg.balance_mode, 'truncate')
    for i = 1:12, Si{i} = Si{i}(1:Nmin); end
elseif strcmpi(cfg.balance_mode, 'random')
    rng(cfg.rng_seed);
    for i = 1:12
        p = randperm(numel(Si{i}), Nmin);
        Si{i} = Si{i}(p);
    end
else
    error('Unsupported cfg.balance_mode: %s', cfg.balance_mode);
end
counts_after = Nmin * ones(12,1,'double');

Tavgi = zeros(12,1,'double');
S0i = cell(12,1);
for i = 1:12
    Tavgi(i) = mean(Si{i}, 'double');
    S0i{i} = Si{i} - Tavgi(i);
end
fJ = vertcat(S0i{:});
JRMS = std(fJ, 0, 'all');
J3u = prctile(fJ, 99.95) - prctile(fJ, 0.05);

% EOJ from event stream (class-specific)
out = struct();
out.UI = UI;
out.Npat = double(cfg.Npat);
if ~isempty(cfg.eoj_start_ui_m)
    out.t0 = double(t_uniform(1) + (cfg.eoj_start_ui_m + 0.5)*UI);
elseif ~isempty(cfg.t0_override)
    out.t0 = double(cfg.t0_override);
else
    out.t0 = double(t_uniform(1));
end
out.t_cross_event = t_cross_all;
out.ui_index_event = trans_n - 1;
out.class_id_event = class_id_event;

out = compute_eoj_core(out, cfg);

% package jitter outputs
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

if cfg.debug_plot
    if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
    fig = figure('Visible','off');
    plot(e_k_all); hold on; plot(e_pll); plot(e_cru); grid on;
    legend({'e_k','e_pll','e_cru'}, 'Location', 'best');
    saveas(fig, fullfile(cfg.outputDir, 'cru_debug.png')); close(fig);

    fig = figure('Visible','off'); histogram(fJ,100); grid on;
    saveas(fig, fullfile(cfg.outputDir, 'histogram_fJ.png')); close(fig);

    fig = figure('Visible','off');
    bar([(1:12).', counts_before, counts_after], 'grouped'); grid on;
    saveas(fig, fullfile(cfg.outputDir, 'counts_before_after.png')); close(fig);
end
end

function out = compute_eoj_core(out, cfg)
UI = double(out.UI);
Npat = double(out.Npat);
cid = out.class_id_event(:);
t_cross = out.t_cross_event(:);
ui_idx = out.ui_index_event(:);
valid = cid >= 1 & cid <= 12;
cid = cid(valid); t_cross = t_cross(valid); ui_idx = ui_idx(valid);

repeat_id = floor(ui_idx ./ Npat) + 1;
t0 = double(out.t0);
t_rep0 = t0 + (repeat_id - 1) * Npat * UI;
t_rel_UI = (t_cross - t_rep0) ./ UI;

EOJ_per_class_UI = nan(12,1); EOJ_per_class_s = nan(12,1);
num_w = zeros(12,1,'double'); notes = strings(12,1);
for i = 1:12
    m = (cid == i);
    rep_i = repeat_id(m); tr_i = t_rel_UI(m);
    if isempty(rep_i), notes(i) = "no events"; continue; end
    [rep_u,~,g] = unique(rep_i);
    Tr = accumarray(g, tr_i, [], @mean);
    wins = [];
    for k = 1:(numel(rep_u)-2)
        r1=rep_u(k); r2=rep_u(k+1); r3=rep_u(k+2);
        if r2==r1+1 && r3==r2+1
            T1=Tr(k); T2=Tr(k+1); T3=Tr(k+1); T4=Tr(k+2);
            wins(end+1,1)=abs((T2-T1) - (T4-T3)); %#ok<AGROW>
        end
    end
    if isempty(wins), notes(i) = "no 3-repeat window"; continue; end
    if strcmpi(cfg.eoj_mode, 'first3')
        EOJ_per_class_UI(i) = wins(1); notes(i) = "first3";
    else
        EOJ_per_class_UI(i) = max(wins); notes(i) = "max";
    end
    EOJ_per_class_s(i) = EOJ_per_class_UI(i) * UI;
    num_w(i) = numel(wins);
end

out.EOJ_per_class_UI = EOJ_per_class_UI;
out.EOJ_UI = max(EOJ_per_class_UI, [], 'omitnan');
out.EOJ_per_class_s = EOJ_per_class_s;
out.EOJ_s = out.EOJ_UI * UI;
out.EOJ_num_windows = num_w;
out.EOJ_notes = notes;

if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
T = table((1:12).', EOJ_per_class_UI, EOJ_per_class_s, num_w, notes, ...
    'VariableNames', {'class_id','EOJ_UI','EOJ_s','num_windows_used','note'});
writetable(T, fullfile(cfg.outputDir, 'eoj_per_class.csv'));
fig = figure('Visible','off'); bar((1:12).', EOJ_per_class_UI); grid on;
saveas(fig, fullfile(cfg.outputDir, 'eoj_bar.png')); close(fig);
end

function [t, v] = read_tv_from_csv(csvPath, cfg)
if ~(ischar(csvPath) || (isstring(csvPath) && isscalar(csvPath)))
    error('csvPath must be char/string.');
end
csvPath = char(csvPath);
if ~exist(csvPath, 'file'), error('File not found: %s', csvPath); end

if endsWith(lower(csvPath), '.csv') || endsWith(lower(csvPath), '.txt')
    try
        T = readtable(csvPath);
        [t, v] = extract_from_table(T, cfg);
    catch
        M = readmatrix(csvPath);
        [t, v] = extract_from_matrix(M, cfg);
    end
else
    error('Unsupported file extension. Use .csv or .txt');
end

t = double(t(:)); v = double(v(:));
if numel(t) ~= numel(v), error('Loaded t and v lengths differ.'); end
end

function [t, v] = extract_from_table(T, cfg)
vars = string(T.Properties.VariableNames);
varsL = lower(vars);
if isfield(cfg,'t_col') && isfield(cfg,'v_col') && ~isempty(cfg.t_col) && ~isempty(cfg.v_col)
    t = T.(cfg.t_col); v = T.(cfg.v_col); return;
end
tCand = ["t","time","time_s","timestamp","sec","seconds"];
vCand = ["v","voltage","voltage_v","signal","diff","vdiff"];
it = find(ismember(varsL,tCand),1,'first');
iv = find(ismember(varsL,vCand),1,'first');
if ~isempty(it) && ~isempty(iv)
    t = T.(vars(it)); v = T.(vars(iv)); return;
end
M = table2array(T); [t, v] = extract_from_matrix(M, cfg);
end

function [t, v] = extract_from_matrix(M, cfg)
if size(M,2) < 2, error('CSV must contain at least 2 columns for t and v.'); end
if isfield(cfg,'t_col_idx') && isfield(cfg,'v_col_idx') && ~isempty(cfg.t_col_idx) && ~isempty(cfg.v_col_idx)
    t = M(:, cfg.t_col_idx); v = M(:, cfg.v_col_idx);
else
    t = M(:,1); v = M(:,2);
end
end

function cfg = apply_defaults(cfg)
if ~isstruct(cfg), error('cfg must be a struct.'); end
req = {'UI','context_patterns'};
for i=1:numel(req)
    if ~isfield(cfg, req{i}), error('cfg.%s is required.', req{i}); end
end
if ~isfield(cfg,'outputDir'), cfg.outputDir = '.'; end
if ~isfield(cfg,'M'), cfg.M = 32; end
if ~isfield(cfg,'balance_mode'), cfg.balance_mode = 'truncate'; end
if ~isfield(cfg,'rng_seed'), cfg.rng_seed = 1; end
if ~isfield(cfg,'edge_dir'), cfg.edge_dir = 'either'; end
if ~isfield(cfg,'Npat'), cfg.Npat = 8191; end
if ~isfield(cfg,'eoj_mode'), cfg.eoj_mode = 'max'; end
if ~isfield(cfg,'eoj_start_ui_m'), cfg.eoj_start_ui_m = []; end
if ~isfield(cfg,'t0_override'), cfg.t0_override = []; end
if ~isfield(cfg,'debug_plot'), cfg.debug_plot = false; end
if ~isfield(cfg,'t_col'), cfg.t_col = ''; end
if ~isfield(cfg,'v_col'), cfg.v_col = ''; end
if ~isfield(cfg,'t_col_idx'), cfg.t_col_idx = []; end
if ~isfield(cfg,'v_col_idx'), cfg.v_col_idx = []; end
end

function tf = edge_match(A,B,edge_dir)
switch lower(edge_dir)
    case 'either', tf = true;
    case 'rising', tf = B > A;
    case 'falling', tf = B < A;
    otherwise, error('Unsupported cfg.edge_dir: %s', edge_dir);
end
end

function s = get_sym(sym_code, idx)
if idx < 1 || idx > numel(sym_code), s = NaN; else, s = sym_code(idx); end
end

function cid = identify_class(context, patterns)
if ~iscell(patterns) || numel(patterns) ~= 12
    error('cfg.context_patterns must be a 12-element cell array.');
end
cid = 0;
for i = 1:12
    if match_pattern(context, patterns{i}), cid = i; return; end
end
end

function tf = match_pattern(context, p)
if isa(p, 'function_handle')
    tf = logical(p(context)); return;
end
if isnumeric(p)
    if isvector(p), tf = match_numeric(context, p(:).'); return; end
    tf = false;
    for r = 1:size(p,1)
        if match_numeric(context, p(r,:)), tf = true; return; end
    end
    return;
end
if isstruct(p)
    tf = true;
    fns = fieldnames(p);
    for k = 1:numel(fns)
        fn = fns{k};
        if ~isfield(context, fn), tf = false; return; end
        if ~is_wild(p.(fn)) && context.(fn) ~= p.(fn), tf = false; return; end
    end
    return;
end
error('Unsupported pattern type: %s', class(p));
end

function tf = match_numeric(context, row)
L = numel(row);
switch L
    case 2, rule = [context.prev1, context.curr];
    case 3, rule = [context.prev2, context.prev1, context.curr];
    case 4, rule = [context.prev2, context.prev1, context.curr, context.next1];
    otherwise, error('Numeric pattern length must be 2/3/4.');
end
tf = true;
for i=1:L
    if ~is_wild(row(i)) && rule(i) ~= row(i), tf = false; return; end
end
end

function tf = is_wild(x)
tf = isnan(x) || isequal(x,-1);
end
