function rot = qprbs13_rotate_x_module(x, Dp)
%QPRBS13_ROTATE_X_MODULE Rotate aligned symbols and build circulant-like matrix.
% Eq(11-14): xr = rotate x by pulse delay Dp.
% Eq(11-15): X is an N-by-N matrix built from xr with cyclic right shifts per row.

    x = double(x(:));
    N = numel(x);
    validateattributes(Dp, {'numeric'}, {'real','finite','integer','nonnegative','scalar'});

    if N == 0
        error('x must be non-empty.');
    end

    d = mod(double(Dp), N);

    % Eq(11-14): xr = [x(N-d+1:N), x(1:N-d)]'  (for d=0, xr=x)
    xr = circshift(x, d);

    % Eq(11-15):
    % row1 = xr(1..N)
    % row2 = xr(N), xr(1), ..., xr(N-1)
    % ...
    ii = (0:N-1).';
    jj = (0:N-1);
    idx = mod(jj - ii, N) + 1;
    X = xr(idx);

    rot = struct();
    rot.xr = xr;
    rot.X = X;
    rot.Dp = d;
end
