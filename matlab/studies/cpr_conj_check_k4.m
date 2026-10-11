function cpr_conj_check_k4()
%CPR_CONJ_CHECK_K4  Conjugate check for multi-tone readings (K = 4, the base value).
%cpr_conj_check.m covers K = 1 only. With K > 1 there is no exact real-times-diagonal
%factorization in general, so three conjugate-related candidates are tested and the best-fitting
%one (lowest amplitude misfit, Eq. (25) of the paper) is reported:
%  c1  x*                       plain conjugate
%  c2  D^-1 (D x)*              D fitted to the K = 4 operator itself (row + column phases,
%                               spanning tree, as in cpr_conj_check)
%  c3  D1^-1 (D1 x)*            D1 fitted to the K = 1 operator of the SAME hardware (same seed,
%                               same cables): the alternative that exists when readings use one tone
%Saves <results>/conj_check_k4.mat.
outDir = cpr_paths().results;
maxNumCompThreads(4);
c0 = cpr_cfg('meas', 'phaseless', 'arch', 'IS');
[X, f] = pattern_crlh(c0.nw, c0.nd, c0.f0idx); x = X(:);
cases = {'IS', {}; 'RS', {'rnRS', 1/32}; 'RS', {'rnRS', 1/16}; 'RS', {'rnRS', 1/8}; ...
         'IC', {'nlegs', 0}; 'IC', {'nlegs', 8}; 'RC', {'nlegs', 0}; 'RC', {'nlegs', 8}};
seeds = 1:5;
Rs = struct('arch', {}, 'ov', {}, 'm1', {}, 'm2', {}, 'm3', {}, 'relImag', {});
fprintf('%-4s %-10s  misfit (dB), min over seeds: x*   D(K=4)   D(K=1)   | max relImag of K=4 fit\n', '', '');
for i = 1:size(cases, 1)
    m1 = zeros(size(seeds)); m2 = m1; m3 = m1; ri = m1;
    for s = seeds
        c4 = cpr_cfg('meas', 'phaseless', 'arch', cases{i,1}, 'K', 4, cases{i,2}{:});
        S4 = make_readings(c4, f, s);
        c1 = cpr_cfg('meas', 'phaseless', 'arch', cases{i,1}, 'K', 1, cases{i,2}{:});
        S1 = make_readings(c1, f, s);
        y = abs(S4.Phi*x);
        mf = @(xs) 20*log10(norm(abs(S4.Phi*xs) - y) / norm(y));
        m1(s) = mf(conj(x));
        [d4, p4] = coldir(S4.Phi); m2(s) = mf(conj(d4) .* conj(d4 .* x));
        [ii, jj, v] = find(S4.Phi); b = conj(p4(ii)) .* v .* conj(d4(jj)); ri(s) = max(abs(imag(b))) / max(abs(b));
        d1 = coldir(S1.Phi); m3(s) = mf(conj(d1) .* conj(d1 .* x));
    end
    Rs(end+1) = struct('arch', cases{i,1}, 'ov', {cases{i,2}}, 'm1', m1, 'm2', m2, 'm3', m3, 'relImag', ri); %#ok<AGROW>
    best = min([m1; m2; m3], [], 1);
    fprintf('%-4s %-10s  %6.1f  %6.1f  %6.1f   | %.2f    best candidate per seed: %s dB\n', cases{i,1}, ov2str(cases{i,2}), ...
        min(m1), min(m2), min(m3), max(ri), mat2str(round(best, 1)));
end
save(fullfile(outDir, 'conj_check_k4.mat'), 'Rs');
end

function [d, p] = coldir(Phi)
[m, n] = size(Phi);
[i, j, v] = find(Phi);
G = graph(i, m + j, [], m + n);
cc = conncomp(G); [~, firstNode] = unique(cc, 'first');
G = addedge(G, repmat(m + n + 1, numel(firstNode), 1), firstNode(:));
E = bfsearch(G, m + n + 1, 'edgetonew');
ph = zeros(m + n + 1, 1);
A = sparse(i, j, angle(v), m, n);
for e = 1:size(E, 1)
    u = E(e,1); w = E(e,2);
    if u == m + n + 1, ph(w) = 0; continue; end
    if u <= m, ph(w) = A(u, w - m) - ph(u); else, ph(w) = A(w, u - m) - ph(u); end
end
p = exp(1i*ph(1:m)); d = exp(1i*ph(m+1:m+n));
end

function s = ov2str(ov)
s = ''; for i = 1:2:numel(ov), s = [s sprintf('%s=%g ', ov{i}, ov{i+1})]; end %#ok<AGROW>
end
