function out = eoj_stats_from_csv(csvPath, cfg)
%EOJ_STATS_FROM_CSV  EOJ统计（独立脚本）
% out = eoj_stats_from_csv(csvPath, cfg)
%
% 输入参数：
%   csvPath (char/string) : .csv 文件路径，文件包含时间(s)和电压(v)
%   cfg.UI    (double, 必填) : UI（秒）
%   cfg.Npat  (double, 可选, 默认8191) : repeat长度（QPRBS13=8191, QPRBS9=511）
%   cfg.eoj_mode (char, 可选, 'max'|'first3', 默认'max')
%   cfg.eoj_start_ui_m (double, 可选) : 从第m个UI中点作为trigger起点
%   cfg.t0_override (double, 可选) : 绝对trigger起点（秒）
%   cfg.t_col/cfg.v_col 或 cfg.t_col_idx/cfg.v_col_idx : 指定列
%   cfg.debug_plot, cfg.outputDir : 可选调试输出
%
% 输出：
%   out.EOJ_UI, out.EOJ_s
%   out.EOJ_windows_UI
%   out.Tr_UI, out.repeat_id
%   out.t0

if ~isstruct(cfg), error('cfg must be a struct.'); end
if ~isfield(cfg,'UI'), error('cfg.UI is required.'); end
if ~isfield(cfg,'Npat'), cfg.Npat = 8191; end
if ~isfield(cfg,'eoj_mode'), cfg.eoj_mode = 'max'; end
if ~isfield(cfg,'debug_plot'), cfg.debug_plot = false; end
if ~isfield(cfg,'outputDir'), cfg.outputDir = '.'; end
if ~isfield(cfg,'eoj_start_ui_m'), cfg.eoj_start_ui_m = []; end
if ~isfield(cfg,'t0_override'), cfg.t0_override = []; end

[t, v] = read_tv_from_csv(csvPath, cfg);
UI = double(cfg.UI);
Npat = double(cfg.Npat);

% crossing抽取（统一门限）
thr = mean(v, 'omitnan');
idx = find((v(1:end-1)-thr).*(v(2:end)-thr) <= 0);
if numel(idx) < 4
    error('Not enough crossings to compute EOJ.');
end

t_cross = zeros(numel(idx),1);
for k = 1:numel(idx)
    i = idx(k);
    t1=t(i); t2=t(i+1); v1=v(i); v2=v(i+1);
    if v2==v1
        t_cross(k) = 0.5*(t1+t2);
    else
        a = (thr-v1)/(v2-v1);
        t_cross(k) = t1 + a*(t2-t1);
    end
end

% trigger起点
if ~isempty(cfg.eoj_start_ui_m)
    m = double(cfg.eoj_start_ui_m);
    validateattributes(m, {'double'}, {'scalar','integer','nonnegative','finite'});
    t0 = t(1) + (m + 0.5) * UI;
elseif ~isempty(cfg.t0_override)
    t0 = double(cfg.t0_override);
else
    t0 = t(1);
end

ui_index = floor((t_cross - t0) ./ UI);
repeat_id = floor(ui_index ./ Npat) + 1;
t_rep0 = t0 + (repeat_id - 1) * Npat * UI;
t_rel_UI = (t_cross - t_rep0) ./ UI;

[rep_u, ~, g] = unique(repeat_id);
Tr_UI = accumarray(g, t_rel_UI, [], @mean);

if numel(rep_u) < 3
    error('Need at least 3 repeats to compute EOJ windows.');
end

eoj_w = [];
for k = 1:(numel(rep_u)-2)
    r1=rep_u(k); r2=rep_u(k+1); r3=rep_u(k+2);
    if r2==r1+1 && r3==r2+1
        T1=Tr_UI(k); T2=Tr_UI(k+1); T3=Tr_UI(k+1); T4=Tr_UI(k+2);
        P2=T2-T1; P3=T4-T3;
        eoj_w(end+1,1)=abs(P2-P3); %#ok<AGROW>
    end
end
if isempty(eoj_w)
    error('No consecutive 3-repeat windows for EOJ.');
end

if strcmpi(cfg.eoj_mode, 'first3')
    EOJ_UI = eoj_w(1);
else
    EOJ_UI = max(eoj_w);
end
EOJ_s = EOJ_UI * UI;

out = struct();
out.EOJ_UI = EOJ_UI;
out.EOJ_s = EOJ_s;
out.EOJ_windows_UI = eoj_w;
out.Tr_UI = Tr_UI;
out.repeat_id = rep_u;
out.t0 = t0;
out.threshold = thr;

if cfg.debug_plot
    if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
    f = figure('Visible','off');
    subplot(2,1,1); plot(rep_u, Tr_UI, '-o'); grid on; title('Tr per repeat (UI)');
    subplot(2,1,2); plot(eoj_w, '-o'); grid on; title('EOJ windows (UI)');
    saveas(f, fullfile(cfg.outputDir, 'eoj_stats_debug.png'));
    close(f);
end
end

function [t, v] = read_tv_from_csv(csvPath, cfg)
if ~(ischar(csvPath) || (isstring(csvPath) && isscalar(csvPath)))
    error('csvPath must be char/string.');
end
csvPath = char(csvPath);
if ~exist(csvPath, 'file'), error('File not found: %s', csvPath); end

opts = detectImportOptions(csvPath, 'NumHeaderLines', 0);
T = readtable(csvPath, opts);
vars = string(T.Properties.VariableNames);
varsL = lower(vars);

if isfield(cfg,'t_col') && isfield(cfg,'v_col') && ~isempty(cfg.t_col) && ~isempty(cfg.v_col)
    t = double(T.(cfg.t_col)); v = double(T.(cfg.v_col));
    t=t(:); v=v(:); return;
end

tCand = ["t","time","time_s","sec","seconds","timestamp"];
vCand = ["v","voltage","voltage_v","signal","diff","vdiff"];
it = find(ismember(varsL, tCand), 1, 'first');
iv = find(ismember(varsL, vCand), 1, 'first');
if ~isempty(it) && ~isempty(iv)
    t = double(T.(vars(it))); v = double(T.(vars(iv)));
    t=t(:); v=v(:); return;
end

M = table2array(T);
if size(M,2) < 2, error('CSV must contain at least 2 columns.'); end
if isfield(cfg,'t_col_idx') && isfield(cfg,'v_col_idx') && ~isempty(cfg.t_col_idx) && ~isempty(cfg.v_col_idx)
    t = double(M(:, cfg.t_col_idx));
    v = double(M(:, cfg.v_col_idx));
else
    t = double(M(:,1));
    v = double(M(:,2));
end
t=t(:); v=v(:);
end
