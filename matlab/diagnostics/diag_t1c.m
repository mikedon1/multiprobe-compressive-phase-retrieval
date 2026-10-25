function diag_t1c()
%DIAG_T1C (2026-10-21) CORRECTED ideal-probe test. R maps OUTGOING-wave coefficients a; the far field is c = a .* j^(l+1) (TE)
%or a .* j^l (TM) (fwc_aut_modes). diag_cfg T1-T3 and diag_seeds('ideal') used the far-field matrix Ri on a directly, i.e.
%a different pattern. Here the ideal operator is Ri * diag(j^(l + TE)), so the ideal probes see the dipole's true field.
%(1) per-port check: R vs the corrected ideal rows after a per-port complex gain; (2) T1 corrected, seeds 1-2.
pf = fopen(fullfile(tempdir, 'diag_t1c_progress.txt'), 'w');
say = @(varargin) fprintf(pf, '%s  %s\n', datestr(now, 'HH:MM:SS'), sprintf(varargin{:}));
try
    cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper'); maxNumCompThreads(2);
    S = load(fullfile(getenv('CPR_RESULTS'), 'fwc_validation.mat')); V = S.V; modes = V.A.modes; nf = 16; L = 12;
    th = (15:10:165).' * pi/180; [Yt, Yp] = vsh_eval(L, th, zeros(size(th)));
    fac = (1i .^ (modes(:, 1) + (modes(:, 3) == 1))).';
    for sgn = [1 -1]
        Ri = [Yt; Yp] .* (fac .^ sgn);
        for f = [1 8 16]
            R = V.R(:, :, f); rr = zeros(32, 1);
            for p = 1:32
                cr = abs(Ri * R(p, :)') ./ (vecnorm(Ri, 2, 2) * norm(R(p, :))); [~, q] = max(cr);
                k = (Ri(q, :) * R(p, :)') / (Ri(q, :) * Ri(q, :)'); rr(p) = 20*log10(norm(R(p, :) - k' * Ri(q, :)) / norm(R(p, :)));
            end
            say('factor j^(%+d(l+TE)) f%-2d: R vs ideal after per-port gain, residual median %6.1f dB, worst %6.1f dB', sgn, f, median(rr), max(rr));
        end
    end
    % pick the convention fwc_aut_modes states (a = c ./ j^(...)  ->  c = a .* j^(...)), sgn = +1
    Ri = repmat([Yt; Yp] .* fac, 1, 1, nf);
    a = V.A.a(:, 1:nf); E = exp(-1i * modes(:, 2).' .* (V.psi(:) * pi/180)); Xm = zeros(32, 32, nf);
    for k = 1:nf, for l = 1:32, Xm(:, l, k) = (Ri(:, :, k) .* E(l, :)) * a(:, k); end, end
    D = struct('Xring', Xm, 'R', Ri, 'modes', modes, 'psi', V.psi(:).' * pi/180, 'f', V.freq(1:nf).', 'aTrue', a);
    for seed = 1:2
        c = cpr_cfg('arch', 'IC', 'snrdB', 50, 'nd', nf, 'K', 4, 'ratio', 0.3, 'vshL', L, 'maxTry', 10, 'solverSNRdB', 50);
        r = fwc_recover(c, seed, D);
        say('T1 corrected (ideal probes, true dipole field), seed %d: per-frequency %6.1f dB, tries %2d, nres %6.1f dB, %5.0f s', ...
            seed, r.nmse.perfreqScale, r.info.numTry, r.info.nres_dB, r.time);
    end
    say('DIAG T1C DONE');
catch ME
    say('ERROR: %s', ME.message);
    for k = 1:min(6, numel(ME.stack)), say('   in %s line %d', ME.stack(k).name, ME.stack(k).line); end
end
fclose(pf);
end
