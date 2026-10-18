function r = cpr_sphere_run(r1, seed, varargin)
%CPR_SPHERE_RUN  Phaseless compressive phase retrieval (IC architecture by default) on the CRLH near field
%on a sphere of radius r1 (grid and frequencies from cpr_settings.m): the compressive side of the
%two-scan comparison (Section VII). Extra name/value pairs override cpr_cfg. Saves
%<results>/cpr_sphere_r<r1 cm>_s<seed>_snr<dB>_pol<n>[_<arch>_m<ratio>].mat.
outDir = cpr_paths().results;
maxNumCompThreads(cpr_settings().threads);
c = cpr_cfg('meas', 'phaseless', 'arch', 'IC', varargin{:});
D = crlh_sphere(r1, c.nw, c.nd, c.f0idx);
X = D.X; if c.npol == 2, X = cat(3, D.Xt, D.Xp); end      % dual polarization: arc theta and phi components
r = cpr_run_pattern(c, seed, X, D.f);
r.r1 = r1;
fprintf('CPR IC on sphere r1 = %.2f m, seed %d, npol %d: per-lambda %.1f dB (vector %.1f dB), global %.1f dB (vector %.1f dB), %.0f s, m = %d\n', ...
    r1, seed, c.npol, r.score.perlambda, r.scoreVec.perlambda, r.score.global, r.scoreVec.global, r.time, r.m);
% IC keeps the original name (Fig. 10's file); other architectures/ratios get a suffix
tag = ''; if ~strcmp(c.arch, 'IC'), tag = sprintf('_%s_m%02d', c.arch, round(100*c.ratio)); end
if strcmp(c.weights, 'gauss'), tag = [tag '_gauss']; end         % never overwrite the switch results
save(fullfile(outDir, sprintf('cpr_sphere_r%02d_s%d_snr%d_pol%d%s.mat', round(100*r1), seed, round(c.snrdB), c.npol, tag)), 'r', '-v7.3');
end
