function twoscan_continue(antenna, r1, r2, snr, seq, varargin)
%TWOSCAN_CONTINUE  Frequency continuation for the two-scan technique (Section VIII): at frequencies where twoscan_run did
%not reach its target after its random starts, the FIRST start is the converged coefficient vector of the neighbouring
%frequency, and later starts are random as before.
%   twoscan_continue('crlh', 0.15, 0.30, 50, [16 15 14 13])   the paper's two chains (after twoscan_run ... 1:16)
%   twoscan_continue('crlh', 0.15, 0.30, 50, [11 10 9 8 7])
% seq(1) is a frequency index that reached the target in twoscan_run; it is solved again with its original settings
% (same noise, seed and starts, so the same solution) to obtain its coefficients. Each later index starts from the
% solution of the previous one. Same data, noise realization (rng(1000 + t)), Nc rule, target and start budget as
% twoscan_run. The entries of seq(2:end) replace those in <results>/ift_floor_<antenna>_r<r1>_<r2>_snr<snr>.mat, which
% must not have reached the target; the file before the first continuation is kept as ..._randomstarts.mat. Every
% entry gets the field 'init' ('random' or 'from <f> GHz').
p = inputParser; p.addParameter('MaxStarts', cpr_settings().nearfield.twoScanStarts); p.addParameter('StopRes_dB', -52);
p.addParameter('r0', 0.08); p.addParameter('extra', 5); p.parse(varargin{:}); o = p.Results;
maxNumCompThreads(cpr_settings().threads);
f = fullfile(cpr_paths().results, sprintf('ift_floor_%s_r%02d_%02d_snr%d.mat', antenna, round(100*r1), round(100*r2), snr));
bak = strrep(f, '.mat', '_randomstarts.mat');
if ~exist(bak, 'file'), copyfile(f, bak); end
S0 = load(f); R = S0.R; tBuild = S0.tBuild;
if ~isfield(R, 'init'), [R.init] = deal('random'); end
nw = cpr_settings().nw; nd = max([R.t]);
D1 = crlh_sphere(r1, nw, nd, 1); D2 = crlh_sphere(r2, nw, nd, 1); s2 = D2.scale / D1.scale; S = D1.S;
cPrev = [];
for i = 1:numel(seq)
    t = seq(i);
    rng(1000 + t);
    [a1t, a1p] = noisy(D1.Et(:,:,t), D1.Ep(:,:,t), snr);
    [a2t, a2p] = noisy(s2*D2.Et(:,:,t), s2*D2.Ep(:,:,t), snr);
    Nc = min(S.N, ceil(D1.k(t)*o.r0) + o.extra);
    args = {'Nc', Nc, 'Method', 'bfgs', 'MaxIter', 2000, 'Restarts', o.MaxStarts, 'StopRes_dB', o.StopRes_dB, 'Seed', 2000 + t};
    if i > 1, args = [args, {'InitCoef', cPrev}]; init = sprintf('from %.2f GHz', D1.f(seq(i-1))/1e9); else, init = 'random'; end
    [Et, Ep, info] = twoscan_solve(a1t, a1p, a2t, a2p, S, D1.k(t), r1, r2, args{:});
    cPrev = info.c;
    if i == 1
        if info.res_dB >= o.StopRes_dB, error('twoscan_continue: the anchor t = %d does not reach the target', t); end
        continue
    end
    a = [Et(:); Ep(:)]; b = reshape([D1.Et(:,:,t), D1.Ep(:,:,t)], [], 1); a2 = a * exp(1j*angle(a'*b));
    E = swe_cprh(Et, Ep, deg2rad(D1.phi)); Xh = reshape(E(D1.map), nw, nw/2);
    X = D1.X(:,:,t); Xh = Xh * exp(1j*angle(Xh(:)'*X(:)));
    e = struct('t', t, 'f', D1.f(t), 'Nc', Nc, 'res_dB', info.res_dB, 'starts', info.starts, 'bestStart', info.restart, ...
        'time', info.time, 'fieldErr', 20*log10(norm(a2 - b)/norm(b)), 'cprhErr', 20*log10(norm(Xh(:) - X(:))/norm(X(:))), ...
        'Xh', Xh, 'init', init);
    j = find([R.t] == t);
    if isempty(j), error('twoscan_continue: t = %d is not in %s (run twoscan_run first)', t, f); end
    if R(j).res_dB < o.StopRes_dB && strcmp(R(j).init, 'random')
        error('twoscan_continue: t = %d already reached the target with random starts', t);
    end
    R(j) = orderfields(e, R(j));
    save(f, 'R', 'tBuild');
    fprintf('%s f %.2f GHz Nc %d (%s): residual %.1f dB after %d starts (best %d), E_cprh error %.1f dB, %.0fs\n', ...
        antenna, e.f/1e9, Nc, init, e.res_dB, e.starts, e.bestStart, e.cprhErr, e.time);
end
end

function [at, ap] = noisy(Et, Ep, snr)
v = 10^(-snr/10) * mean(abs([Et(:); Ep(:)]).^2);
at = abs(Et + sqrt(v/2)*(randn(size(Et)) + 1j*randn(size(Et))));
ap = abs(Ep + sqrt(v/2)*(randn(size(Ep)) + 1j*randn(size(Ep))));
end
