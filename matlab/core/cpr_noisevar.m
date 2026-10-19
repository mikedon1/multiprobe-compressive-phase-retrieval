function wvar = cpr_noisevar(c, z, x)
%CPR_NOISEVAR  Variance of the complex noise added before detection, for SNR c.snrdB.
%   wvar = cpr_noisevar(c, z, x)   z: noiseless combined signals (one per row), x: the pattern
% cpr_settings().snrModel selects the definition:
%   'detector' (default): a fixed detector noise. The SNR is the reference ratio of the average
%       power of ONE probe at ONE tone with the full transmit power, mean(|x|^2), to the noise
%       variance. The K tones of a reading share the transmit power (amplitude 1/sqrt(K) each);
%       since the sensing rows use unit-amplitude tones, this is applied as a noise variance K
%       times larger. Combining P probes with random relative phases then raises the average
%       signal power of a reading about P times, as it would in hardware.
%   'normalized': the noise is set relative to the average power of the rows actually measured,
%       mean(|z|^2), so every configuration has the same per-row SNR (earlier versions).
S = cpr_settings();
model = 'detector'; if isfield(S, 'snrModel'), model = S.snrModel; end
switch model
    case 'detector'
        K = 1; if isfield(c, 'K'), K = c.K; end
        wvar = 10^(-c.snrdB/10) * K * mean(abs(x(:)).^2);
    case 'normalized'
        wvar = 10^(-c.snrdB/10) * mean(abs(z(:)).^2);
    otherwise
        error('cpr_noisevar: unknown snrModel %s', model);
end
end
