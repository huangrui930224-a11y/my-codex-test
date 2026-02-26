function sndrOut = qprbs13_sndr_module(pmax, sigmae, sigman)
%QPRBS13_SNDR_MODULE Compute SNDR in dB.
% Eq.(27-5):
%   SNDR = 10*log10( pmax^2 / (sigma_e^2 + sigma_n^2) )

    pmax = double(pmax);
    sigmae = double(sigmae);
    sigman = double(sigman);

    validateattributes(pmax, {'numeric'}, {'real','finite','scalar'});
    validateattributes(sigmae, {'numeric'}, {'real','finite','nonnegative','scalar'});
    validateattributes(sigman, {'numeric'}, {'real','finite','nonnegative','scalar'});

    signalPower = pmax.^2;
    noiseDistPower = sigmae.^2 + sigman.^2;

    if noiseDistPower <= 0
        error('sigma_e^2 + sigma_n^2 must be > 0 to compute SNDR.');
    end

    sndr_db = 10 * log10(signalPower / noiseDistPower);

    sndrOut = struct();
    sndrOut.pmax = pmax;
    sndrOut.sigmae = sigmae;
    sndrOut.sigman = sigman;
    sndrOut.sndr_db = sndr_db;
end
