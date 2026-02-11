function [y, Y, info] = resample_qprbs13_cei_step1(t, v_tx, M)
%RESAMPLE_QPRBS13_CEI_STEP1 Step-1 resampling per CEI workflow.
%   [y, Y, info] = resample_qprbs13_cei_step1(t, v_tx, M)
%
% Inputs
%   t    : time vector (s)
%   v_tx : TX analog output waveform (V)
%   M    : samples per UI (default 64)
%
% Outputs
%   y    : resampled waveform sequence, column vector, length M*(20*N)
%   Y    : M x N matrix (UI-blocked, averaged across 20 analyzed cycles)
%   info : struct with constants and intermediate indexing metadata
%          includes y_matrix_full (M x 20*N) and m_center index
%
% Global parameters (must-use)
%   Data rate    : 224G PAM4
%   fb           : 112e9 (symbol rate)
%   UI           : 1/fb
%   Sequence     : QPRBS13-CEI
%   N            : 8191 symbols / cycle
%   N1           : 20*N analyzed symbols (20 cycles)

    if nargin < 3 || isempty(M)
        M = 64;
    end

    % ---- Fixed constants from specification ----
    fb = 112e9;
    UI = 1 / fb;
    N = 8191;
    N1 = 20 * N;

    % ---- Basic input checks ----
    t = t(:);
    v_tx = v_tx(:);

    if numel(t) ~= numel(v_tx)
        error('t and v_tx must have the same length.');
    end
    if numel(t) < 2
        error('t and v_tx must contain at least 2 samples.');
    end
    if any(diff(t) <= 0)
        error('t must be strictly increasing.');
    end
    if ~(isscalar(M) && M == floor(M) && M > 0)
        error('M must be a positive integer.');
    end

    % ---- Step 1: choose a valid 20*N-UI window from captured waveform ----
    % The waveform may be cropped and t(1) may not be 0 or a cycle boundary.
    total_ui = floor((t(end) - t(1)) / UI);
    if total_ui < N1
        error('Input waveform is too short. Need at least 20*N UIs of data.');
    end

    % Use a centered N1-UI window to reduce edge sensitivity.
    start_ui = floor((total_ui - N1) / 2);
    t_start = t(1) + start_ui * UI;
    t_end = t_start + N1 * UI;

    % ---- Generate M uniformly spaced sample offsets inside each UI ----
    % (Left-closed, right-open per UI to avoid duplicated boundaries.)
    ui_offsets = (0:M-1).' * (UI / M);  % M x 1

    % ---- Build sampling grid for 20*N UIs ----
    ui_starts = t_start + (0:N1-1) * UI;      % 1 x N1
    t_sample = ui_offsets + ui_starts;         % M x N1 (implicit expansion)

    % ---- Ensure interpolation range is covered ----
    if t_sample(1,1) < t(1) || t_sample(end,end) > t(end)
        error('Input waveform does not cover required 20*N UI sampling window.');
    end

    % ---- Resample by interpolation ----
    y_mat = interp1(t, v_tx, t_sample, 'linear');  % M x N1

    % Sequence y(k): concatenate UI blocks column-by-column.
    y = y_mat(:);

    % Build Y with dimension M x N (UI-blocked), using 20-cycle average.
    y_3d = reshape(y_mat, M, N, 20);
    Y = mean(y_3d, 3);

    % ---- Metadata ----
    info = struct();
    info.data_rate = 224e9;
    info.fb = fb;
    info.UI = UI;
    info.N = N;
    info.N1 = N1;
    info.M = M;
    info.t_start = t_start;
    info.t_end = t_end;
    info.ui_offsets = ui_offsets;
    info.sample_time_grid = t_sample;
    info.y_matrix_full = y_mat;
    info.m_center = round(M / 2);
end
