function out = jitter_stats_from_csv(csvPath, cfg)
%JITTER_STATS_FROM_CSV  JRMS/J3u统计（独立脚本）
% out = jitter_stats_from_csv(csvPath, cfg)
%
% 输入参数：
%   csvPath (char/string) : .csv 文件路径，文件需包含时间(s)和电压(v)列。
%   cfg.UI    (double, 必填) : UI（秒）
%   cfg.M     (double, 可选, 默认32) : 每UI重采样点数
%   cfg.pr_low/high (double, 可选, 默认0.05/99.95) : J3u分位点
%   cfg.debug_plot (logical, 可选, 默认false)
%   cfg.outputDir (char/string, 可选, 默认'.')
%   cfg.t_col/cfg.v_col (可选) : 指定列名
%   cfg.t_col_idx/cfg.v_col_idx (可选) : 指定列索引(1-based)
%
% 输出：
%   out.JRMS   : RMS抖动（秒）
%   out.J3u    : 峰峰值抖动（秒）= prctile(high)-prctile(low)
%   out.jitter : 逐transition时间误差（秒）
%   out.t_cross: crossing time（秒）

if ~isstruct(cfg), error('cfg must be a struct.'); end
if ~isfield(cfg,'UI'), error('cfg.UI is required.'); end
if ~isfield(cfg,'M'), cfg.M = 32; end
if ~isfield(cfg,'pr_low'), cfg.pr_low = 0.05; end
if ~isfield(cfg,'pr_high'), cfg.pr_high = 99.95; end
if ~isfield(cfg,'debug_plot'), cfg.debug_plot = false; end
if ~isfield(cfg,'outputDir'), cfg.outputDir = '.'; end

[t, v] = read_tv_from_csv(csvPath, cfg);
validateattributes(t, {'double'}, {'column','real','finite','nonempty'});
validateattributes(v, {'double'}, {'column','real','finite','numel',numel(t)});

UI = double(cfg.UI);
M  = double(cfg.M);
Ts = UI / M;

t_uniform = (t(1):Ts:t(end)).';
v_uniform = interp1(t, v, t_uniform, 'linear', 'extrap');

thr = mean(v_uniform, 'omitnan');
idx = find((v_uniform(1:end-1)-thr).*(v_uniform(2:end)-thr) <= 0);
if numel(idx) < 3
    error('Not enough crossings found from waveform.');
end

t_cross = zeros(numel(idx),1);
for k = 1:numel(idx)
    i = idx(k);
    t1 = t_uniform(i); t2 = t_uniform(i+1);
    v1 = v_uniform(i); v2 = v_uniform(i+1);
    if v2 == v1
        t_cross(k) = 0.5*(t1+t2);
    else
        a = (thr-v1)/(v2-v1);
        t_cross(k) = t1 + a*(t2-t1);
    end
end

% 一阶差分得到transition间时间误差
jitter = diff(t_cross) - UI;
JRMS = std(jitter, 0, 'omitnan');
J3u = prctile(jitter, cfg.pr_high) - prctile(jitter, cfg.pr_low);

out = struct();
out.JRMS = JRMS;
out.J3u = J3u;
out.jitter = jitter;
out.t_cross = t_cross;
out.threshold = thr;

if cfg.debug_plot
    if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
    f = figure('Visible','off');
    subplot(2,1,1); plot(t_uniform, v_uniform); hold on; yline(thr,'r--'); grid on; title('Waveform and threshold');
    subplot(2,1,2); histogram(jitter, 100); grid on; title('Jitter histogram');
    saveas(f, fullfile(cfg.outputDir, 'jitter_stats_debug.png'));
    close(f);
end
end

function [t, v] = read_tv_from_csv(csvPath, cfg)
if ~(ischar(csvPath) || (isstring(csvPath) && isscalar(csvPath)))
    error('csvPath must be char/string.');
end
csvPath = char(csvPath);
if ~exist(csvPath, 'file')
    error('File not found: %s', csvPath);
end

opts = detectImportOptions(csvPath, 'NumHeaderLines', 0);
T = readtable(csvPath, opts);
vars = string(T.Properties.VariableNames);
varsL = lower(vars);

if isfield(cfg,'t_col') && isfield(cfg,'v_col') && ~isempty(cfg.t_col) && ~isempty(cfg.v_col)
    t = double(T.(cfg.t_col));
    v = double(T.(cfg.v_col));
    t = t(:); v = v(:); return;
end

tCand = ["t","time","time_s","sec","seconds","timestamp"];
vCand = ["v","voltage","voltage_v","signal","diff","vdiff"];
it = find(ismember(varsL, tCand), 1, 'first');
iv = find(ismember(varsL, vCand), 1, 'first');
if ~isempty(it) && ~isempty(iv)
    t = double(T.(vars(it))); v = double(T.(vars(iv)));
    t = t(:); v = v(:); return;
end

M = table2array(T);
if size(M,2) < 2
    error('CSV must contain at least 2 columns.');
end
if isfield(cfg,'t_col_idx') && isfield(cfg,'v_col_idx') && ~isempty(cfg.t_col_idx) && ~isempty(cfg.v_col_idx)
    t = double(M(:, cfg.t_col_idx));
    v = double(M(:, cfg.v_col_idx));
else
    t = double(M(:,1));
    v = double(M(:,2));
end

t = t(:); v = v(:);
end
