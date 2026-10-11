function twoscan_run(antenna, r1, r2, snr, tidx, varargin)
%IFT_SPHERE_FLOOR  Two-sphere IFT with the same restart rule as CPR: keep starting (quasi-
%Newton, random starts) until the data residual reaches the noise floor, up to MaxStarts.
%Per-frequency results are saved after each frequency (resumable) to
%<results>/ift_floor_<antenna>_r<r1>_<r2>_snr<snr>.mat.
%   twoscan_run('crlh', 0.15, 0.30, 50, 1:16)
% StopRes_dB default -52: at 50 dB SNR correct solutions reach ~-53 dB and wrong ones stall at
% -38 to -48 dB (60-start run, the paper). Noise: rng(1000 + t) per frequency.
p = inputParser; p.addParameter('MaxStarts', cpr_settings().nearfield.twoScanStarts); p.addParameter('StopRes_dB', -52);
p.addParameter('r0', 0.08); p.addParameter('extra', 5); p.parse(varargin{:}); o = p.Results;
outDir = cpr_paths().results;
maxNumCompThreads(cpr_settings().threads);
f = fullfile(outDir, sprintf('ift_floor_%s_r%02d_%02d_snr%d.mat', antenna, round(100*r1), round(100*r2), snr));
if exist(f, 'file'), S = load(f); R = S.R; else, R = struct([]); end
nd = max(tidx);
nw = cpr_settings().nw;
D1 = crlh_sphere(r1, nw, nd, 1); D2 = crlh_sphere(r2, nw, nd, 1);
s2 = D2.scale / D1.scale; S = D1.S;
for t = tidx
    if ~isempty(R) && any([R.t] == t), continue; end
    rng(1000 + t);
    [a1t, a1p] = noisy(D1.Et(:,:,t), D1.Ep(:,:,t), snr);
    [a2t, a2p] = noisy(s2*D2.Et(:,:,t), s2*D2.Ep(:,:,t), snr);
    Nc = min(S.N, ceil(D1.k(t)*o.r0) + o.extra);
    [Et, Ep, info] = twoscan_solve(a1t, a1p, a2t, a2p, S, D1.k(t), r1, r2, 'Nc', Nc, 'Method', 'bfgs', ...
        'MaxIter', 2000, 'Restarts', o.MaxStarts, 'StopRes_dB', o.StopRes_dB, 'Seed', 2000 + t);
    a = [Et(:); Ep(:)]; b = reshape([D1.Et(:,:,t), D1.Ep(:,:,t)], [], 1);
    a2 = a * exp(1j*angle(a'*b));
    E = swe_cprh(Et, Ep, deg2rad(D1.phi)); Xh = reshape(E(D1.map), nw, nw/2);
    X = D1.X(:,:,t); Xh = Xh * exp(1j*angle(Xh(:)'*X(:)));
    e = struct('t', t, 'f', D1.f(t), 'Nc', Nc, 'res_dB', info.res_dB, 'starts', info.starts, ...
        'bestStart', info.restart, 'time', info.time, 'fieldErr', 20*log10(norm(a2 - b)/norm(b)), ...
        'cprhErr', 20*log10(norm(Xh(:) - X(:))/norm(X(:))), 'Xh', Xh);
    if isempty(R), R = e; else, R(end+1) = e; end %#ok<AGROW>
    save(f, 'R');
    fprintf('%s f %.2f GHz Nc %d: residual %.1f dB after %d starts (best %d), E_cprh error %.1f dB, field %.1f dB, %.0fs\n', ...
        antenna, e.f/1e9, Nc, e.res_dB, e.starts, e.bestStart, e.cprhErr, e.fieldErr, e.time);
end
end

function [at, ap] = noisy(Et, Ep, snr)
v = 10^(-snr/10) * mean(abs([Et(:); Ep(:)]).^2);
at = abs(Et + sqrt(v/2)*(randn(size(Et)) + 1j*randn(size(Et))));
ap = abs(Ep + sqrt(v/2)*(randn(size(Ep)) + 1j*randn(size(Ep))));
end
