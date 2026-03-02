function fitOut = qprbs13_linear_fit_module(Y, xOrX, cfg)
%QPRBS13_LINEAR_FIT_MODULE
% Linear fit on resampled data using Eq.(11-16)~(11-18).
%
% Inputs
%   Y : MxN matrix from y(k) reshape (Eq.11-13)
%   xOrX : either xr as Nx1 vector, or a pre-built X matrix
%   cfg.Np  : pulse span parameter (default 29)
%   cfg.TNp : number of selected rows from X (default 2*Np+1)
%
% Outputs
%   fitOut.P      : M x (TNp+1) coefficient matrix (Eq.11-16)
%   fitOut.E      : M x N error matrix (Eq.11-17)
%   fitOut.e      : error waveform read column-wise from E
%   fitOut.sigmae : std(e)
%   fitOut.P1     : M x TNp matrix (Eq.11-18)
%   fitOut.p      : pulse response read column-wise from P1
%   fitOut.pmax   : max(p)
%   fitOut.X1     : (TNp+1) x N fitting matrix [X(1:TNp,:); ones(1,N)]

    Y = double(Y);
    xOrX = double(xOrX);

    if ndims(Y) ~= 2
        error('Y must be a 2-D matrix.');
    end

    [~, N] = size(Y);
    Np = double(get_cfg(cfg, 'Np', 29));
    TNp = double(get_cfg(cfg, 'TNp', 2*Np + 1));

    validateattributes(Np, {'numeric'}, {'real','finite','integer','positive','scalar'});
    validateattributes(TNp, {'numeric'}, {'real','finite','integer','positive','scalar'});

    if TNp > N
        error('TNp=%d exceeds N=%d. Reduce TNp or increase N.', TNp, N);
    end

    % Build X1 using first TNp rows of X plus a row of ones.
    X1 = build_x1(xOrX, N, TNp);

    % Eq.(11-16): P = Y*X1^T*(X1*X1^T)^(-1)
    G = X1 * X1.';
    if rcond(G) < 1e-12
        warning('X1*X1'' is ill-conditioned (rcond=%g). Using pinv for stability.', rcond(G));
        PinvTerm = pinv(G);
    else
        PinvTerm = inv(G);
    end
    P = Y * X1.' * PinvTerm;

    % Eq.(11-17): E = P*X1 - Y
    E = P * X1 - Y;
    e = E(:);            % read column-wise
    sigmae = std(e, 0);  % sample std

    % Eq.(11-18): P1 is first TNp columns of P
    P1 = P(:, 1:TNp);
    p = P1(:);           % read column-wise
    pmax = max(p);

    fitOut = struct();
    fitOut.Np = Np;
    fitOut.TNp = TNp;
    fitOut.X1 = X1;
    fitOut.P = P;
    fitOut.E = E;
    fitOut.e = e;
    fitOut.sigmae = sigmae;
    fitOut.P1 = P1;
    fitOut.p = p;
    fitOut.pmax = pmax;
end

function v = get_cfg(cfg, name, defaultValue)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = defaultValue;
    end
end

function X1 = build_x1(xOrX, N, TNp)
    if isvector(xOrX)
        xr = xOrX(:);
        if numel(xr) ~= N
            error('Dimension mismatch: Y is MxN but xr must have N elements.');
        end
        ii = (0:TNp-1).';
        jj = (0:N-1);
        idx = mod(jj - ii, N) + 1;
        Xsel = xr(idx);
    else
        if ndims(xOrX) ~= 2
            error('xOrX must be a vector xr or a 2-D matrix X.');
        end
        [rowsX, colsX] = size(xOrX);
        if colsX ~= N
            error('Dimension mismatch: Y is MxN but X must have N columns.');
        end
        if rowsX < TNp
            error('X has %d rows, but TNp=%d rows are required.', rowsX, TNp);
        end
        Xsel = xOrX(1:TNp, :);
    end

    X1 = [Xsel; ones(1, N)];
end
