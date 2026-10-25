function diag_geo3()
%DIAG_GEO2 (2026-10-21) Does a FULL ring help WITHOUT polar coverage? (diag_geo: full circle incl. poles -36.5 dB; half arc
%pole to pole -30.1; half arc 15-165 -4.3.) Same setup as diag_geo (Sec VII pipeline, ideal far-field probes, correct
%convention, dipole, base, seed 1).
%  G3 full ring with the built arc's 15-deg polar gaps: theta = 15:10:165 at phi = 0 AND at phi = 180 (the arc and its
%     mirror, 32 positions x 2 pol = 64 ports), 16 rotations over 180 deg. This is what the existing 32-rotation HFSS
%     data can emulate (the arc at psi + 180 is the far half of a full ring at psi).
pf = fopen(fullfile(tempdir, 'diag_geo3_progress.txt'), 'w');
say = @(varargin) fprintf(pf, '%s  %s\n', datestr(now, 'HH:MM:SS'), sprintf(varargin{:}));
try
    cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper'); maxNumCompThreads(2);
    S = load(fullfile(getenv('CPR_RESULTS'), 'fwc_validation.mat')); V = S.V; modes = V.A.modes; nf = 16; L = 12;
    fac = (1i .^ (modes(:, 1) + (modes(:, 3) == 1))).'; a = V.A.a(:, 1:nf);
    th = linspace(5, 175, 16).' * pi/180; ph = zeros(16, 1);
    [Yt, Yp] = vsh_eval(L, th, ph); Rg = [Yt; Yp] .* fac;
    R = repmat(Rg, 1, 1, nf); psi = (0:31) * 11.25; np = size(R, 1); nrot = numel(psi);
    E = exp(-1i * modes(:, 2).' .* (psi(:) * pi/180)); Xm = zeros(np, nrot, nf);
    for k = 1:nf, for l = 1:nrot, Xm(:, l, k) = (R(:, :, k) .* E(l, :)) * a(:, k); end, end
    D = struct('Xring', Xm, 'R', R, 'modes', modes, 'psi', psi * pi/180, 'f', V.freq(1:nf).', 'aTrue', a);
    c = cpr_cfg('arch', 'IC', 'snrdB', 50, 'nd', nf, 'K', 4, 'ratio', 0.3, 'vshL', L, 'maxTry', 10, 'solverSNRdB', 50);
    r = fwc_recover(c, 1, D);
    say('G4 half arc 5-175 deg (32 rot / 360 deg)  per-frequency %6.1f dB, tries %2d, nres %6.1f dB, m %d, %5.0f s', ...
        r.nmse.perfreqScale, r.info.numTry, r.info.nres_dB, r.m, r.time);
    say('DIAG GEO2 DONE');
catch ME
    say('ERROR: %s', ME.message);
    for k = 1:min(6, numel(ME.stack)), say('   in %s line %d', ME.stack(k).name, ME.stack(k).line); end
end
fclose(pf);
end
