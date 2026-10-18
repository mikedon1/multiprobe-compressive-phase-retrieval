function cpr_conj_check()
%CPR_CONJ_CHECK  Conjugate ambiguity of the sensing operators (50 cm cables, on/off switches, K = 1),
%Section III-C.
%  Factorization: D = the phase of the first nonzero entry in each column of Phi. If Phi = B D
%  with B real, then x# = D^-1 (D x)* gives the same |Phi x|. Reported per operator:
%    relImag = max|imag(B)| / max|B| (0: Phi factors exactly), and
%    misfit  = || |Phi x#| - |Phi x| || / || Phi x || (dB; amplitude-domain).
%  Sparsity of x# vs. x (IS, 20 cable realizations): coefficients holding 99% of the energy in
%  the joint basis (2D DFT x delay-compensated DCT), after the K = 1 per-wavelength neighbor
%  alignment, minimized over the delay tau (the solver re-chooses both).
%Saves <results>/conj_check_v2.mat.
outDir = cpr_paths().results;
maxNumCompThreads(2);
c0 = cpr_cfg('meas', 'phaseless', 'arch', 'IS');
[X, f] = pattern_crlh(c0.nw, c0.nd, c0.f0idx); x = X(:);
cases = {'IS', {}; 'RS', {'rnRS', 1/32}; 'RS', {'rnRS', 1/16}; 'RS', {'rnRS', 1/8}; ...
         'IC', {'nlegs', 0}; 'IC', {'nlegs', 8}; 'RC', {'nlegs', 0}; 'RC', {'nlegs', 8}};
seeds = 1:10;
Rs = struct('arch', {}, 'ov', {}, 'relImag', {}, 'misfit', {});
fprintf('%-4s %-12s  max relImag   misfit (dB) per seed\n', '', '');
for i = 1:size(cases, 1)
    ri = zeros(size(seeds)); mf = ri;
    for s = seeds
        c = cpr_cfg('meas', 'phaseless', 'arch', cases{i,1}, 'K', 1, cases{i,2}{:});
        S = make_readings(c, f, s);
        [ri(s), mf(s)] = misfit(S.Phi, x);
    end
    Rs(end+1) = struct('arch', cases{i,1}, 'ov', {cases{i,2}}, 'relImag', ri, 'misfit', mf); %#ok<AGROW>
    fprintf('%-4s %-12s  %.1e   %s  (range %.1f to %.1f)\n', cases{i,1}, ov2str(cases{i,2}), max(ri), ...
        sprintf('%.1f ', mf), min(mf), max(mf));
end
% --- sparsity of the conjugate alternative (IS) ---
nr = 20; kx = k99(X, f); kc = zeros(1, nr);
for s = 1:nr
    c = cpr_cfg('meas', 'phaseless', 'arch', 'IS', 'K', 1);
    S = make_readings(c, f, s);
    d = coldir(S.Phi);
    xs = conj(d) .* conj(d .* x);
    kc(s) = k99(reshape(xs, size(X)), f);
end
fprintf('\nIS 99%% coefficients: x %d, x# mean %.0f (range %d to %d) over %d cable realizations\n', ...
    kx, mean(kc), min(kc), max(kc), nr);
save(fullfile(outDir, 'conj_check_v2.mat'), 'Rs', 'kx', 'kc');
end

function [d, p] = coldir(Phi)
% Phaseless K = 1 rows carry a known random phase e^{j psi} (make_readings), so the factorization
% to test is Phi = P B D (P, D diagonal phases, B real). Phase(Phi_ij) = pi_i + delta_j is solved
% exactly on a spanning tree of the bipartite row-column graph (one free phase per component);
% every other entry then tests the factorization.
[m, n] = size(Phi);
[i, j, v] = find(Phi);
G = graph(i, m + j, [], m + n);
cc = conncomp(G); [~, firstNode] = unique(cc, 'first');
G = addedge(G, repmat(m + n + 1, numel(firstNode), 1), firstNode(:));   % virtual root
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

function [ri, mf] = misfit(Phi, x)
[d, p] = coldir(Phi);
[i, j, v] = find(Phi);
b = conj(p(i)) .* v .* conj(d(j)); ri = max(abs(imag(b))) / max(abs(b));
xs = conj(d) .* conj(d .* x);          % D^-1 (D x)*, |d| = 1
y = abs(Phi*x);
mf = 20*log10(norm(abs(Phi*xs) - y) / norm(y));
end

function K = k99(X, f)
[nw, nh, nd] = size(X);
Y = X;                                  % K = 1 per-wavelength gauge: align each slice to its neighbor
for k = 2:nd
    a = Y(:,:,k-1); b = Y(:,:,k); Y(:,:,k) = b * exp(-1i*angle(a(:)'*b(:)));
end
X2 = reshape(fft2(Y) / sqrt(nw*nh), nw*nh, nd);
K = inf;
for tau = (-20:0.25:20)*1e-9
    S2 = X2 * conj(cpr_basis('dct', f, tau));
    e = sort(abs(S2(:)).^2, 'descend');
    K = min(K, find(cumsum(e) >= 0.99*sum(e), 1));
end
end

function s = ov2str(ov)
s = ''; for i = 1:2:numel(ov), s = [s sprintf('%s=%g ', ov{i}, ov{i+1})]; end %#ok<AGROW>
end
