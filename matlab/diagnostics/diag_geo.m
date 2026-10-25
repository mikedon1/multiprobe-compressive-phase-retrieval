function diag_geo()
%DIAG_GEO (2026-10-21) Sampling geometry vs operator. Section VII pipeline (fwc_recover), ideal FAR-field point probes with the
%CORRECT convention (rows Ri * diag(j^(l+TE)), so the probes see the dipole's true field), dipole coefficients, base setting
%(16 freq, K=4, m/n 0.3, L=12, 50 dB, solverSNRdB 50), seed 1.
%  G1 Section VI geometry: full great circle, theta_s = 0:11.25:348.75 (32 positions x 2 pol = 64 ports; theta_s > 180 deg
%     is (360 - theta_s, phi = 180) with the across (theta-hat) component reversed, as vsh_scan), 16 rotations psi = 0:11.25:168.75
%  G2 half arc pole to pole: theta = linspace(0.01, 179.99, 16), 32 ports, 32 rotations over 360 (corrected T2)
%(T1 corrected, arc 15-165 deg, 32 rotations: -4.3 dB, diag_t1c.)
pf = fopen(fullfile(tempdir, 'diag_geo_progress.txt'), 'w');
say = @(varargin) fprintf(pf, '%s  %s\n', datestr(now, 'HH:MM:SS'), sprintf(varargin{:}));
try
    cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper'); maxNumCompThreads(2);
    S = load(fullfile(getenv('CPR_RESULTS'), 'fwc_validation.mat')); V = S.V; modes = V.A.modes; nf = 16; L = 12;
    fac = (1i .^ (modes(:, 1) + (modes(:, 3) == 1))).'; a = V.A.a(:, 1:nf);
    % G1
    ts = (0:31).' * 11.25; far = ts > 180; th = ts; th(far) = 360 - ts(far); ph = zeros(size(ts)); ph(far) = 180;
    th = min(max(th, 0.01), 179.99) * pi/180; ph = ph * pi/180; sg = ones(size(ts)); sg(far) = -1;
    [Yt, Yp] = vsh_eval(L, th, ph); R1 = [sg .* Yt; Yp] .* fac;
    geos = {'G1 Sec VI geometry (full circles, 16 rot / 180 deg)', R1, (0:15) * 11.25};
    th2 = linspace(0.01, 179.99, 16).' * pi/180; [Yt2, Yp2] = vsh_eval(L, th2, zeros(size(th2)));
    geos(end+1, :) = {'G2 half arc pole to pole (32 rot / 360 deg)', [Yt2; Yp2] .* fac, (0:31) * 11.25};
    for i = 1:size(geos, 1)
        R = repmat(geos{i, 2}, 1, 1, nf); psi = geos{i, 3}; np = size(R, 1); nrot = numel(psi);
        E = exp(-1i * modes(:, 2).' .* (psi(:) * pi/180)); Xm = zeros(np, nrot, nf);
        for k = 1:nf, for l = 1:nrot, Xm(:, l, k) = (R(:, :, k) .* E(l, :)) * a(:, k); end, end
        D = struct('Xring', Xm, 'R', R, 'modes', modes, 'psi', psi * pi/180, 'f', V.freq(1:nf).', 'aTrue', a);
        c = cpr_cfg('arch', 'IC', 'snrdB', 50, 'nd', nf, 'K', 4, 'ratio', 0.3, 'vshL', L, 'maxTry', 10, 'solverSNRdB', 50);
        r = fwc_recover(c, 1, D);
        say('%-52s per-frequency %6.1f dB, tries %2d, nres %6.1f dB, m %d, %5.0f s', geos{i, 1}, r.nmse.perfreqScale, r.info.numTry, r.info.nres_dB, r.m, r.time);
    end
    say('DIAG GEO DONE');
catch ME
    say('ERROR: %s', ME.message);
    for k = 1:min(6, numel(ME.stack)), say('   in %s line %d', ME.stack(k).name, ME.stack(k).line); end
end
fclose(pf);
end
