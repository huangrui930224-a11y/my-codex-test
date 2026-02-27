function vfOut = qprbs13_steady_state_vf_module(Y, xr, M, cfg)
%QPRBS13_STEADY_STATE_VF_MODULE
% Recompute linear-fit pulse p(k) with Np=20, then compute steady-state voltage v_f.
% Eq.(11-12): v_f = (1/M) * sum_k p(k)
%
% Inputs
%   Y   : MxN resampled matrix
%   xr  : rotated symbol vector used to build fitting rows
%   M   : samples per UI
%   cfg : optional fields
%       .NpVf   (default 20)
%       .TNpVf  (default 2*NpVf+1)
%
% Outputs
%   vfOut.vf          : steady-state voltage estimate
%   vfOut.NpVf        : Np used for this re-fit
%   vfOut.TNpVf       : TNp used for this re-fit
%   vfOut.p_vf        : pulse response waveform (column-wise from P1)
%   vfOut.fit_linear  : full output struct from qprbs13_linear_fit_module

    validateattributes(M, {'numeric'}, {'real','finite','positive','integer','scalar'});
    M = double(M);

    NpVf = double(get_cfg(cfg, 'NpVf', 20));
    TNpVf = double(get_cfg(cfg, 'TNpVf', 2*NpVf + 1));

    validateattributes(NpVf, {'numeric'}, {'real','finite','integer','positive','scalar'});
    validateattributes(TNpVf, {'numeric'}, {'real','finite','integer','positive','scalar'});

    cfgVf = struct();
    cfgVf.Np = NpVf;
    cfgVf.TNp = TNpVf;
    fitVf = qprbs13_linear_fit_module(Y, xr, cfgVf);

    p_vf = fitVf.p;
    vf = sum(p_vf) / M;

    vfOut = struct();
    vfOut.vf = vf;
    vfOut.NpVf = NpVf;
    vfOut.TNpVf = TNpVf;
    vfOut.p_vf = p_vf;
    vfOut.fit_linear = fitVf;
end

function v = get_cfg(cfg, name, defaultValue)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        v = cfg.(name);
    else
        v = defaultValue;
    end
end
