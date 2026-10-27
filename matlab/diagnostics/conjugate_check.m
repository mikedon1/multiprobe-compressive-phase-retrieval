function conjugate_check()
%CONJUGATE_CHECK  Ambiguity check at the paper's defaults (supplement, conjugate paragraph): dual polarization, 32 x 16 grid,
%16 wavelengths, VSH basis L = 12. For each alternative pattern that phaseless data cannot tell apart
%from the truth in general, report (a) how far it lies outside the span of the basis (dB, energy outside
%/ total) and (b) how well it fits the measured amplitudes (misfit || |Phi x_alt| - |Phi x| || / ||Phi x||,
%dB; -inf = fits exactly). An alternative is excluded if (a) or (b) is large.
%  Alternatives: independent phase per rotation (phi column; the IC/RC block ambiguity), the plain
%  conjugate conj(x), and the cable-phase conjugate x# = D^-1 (D x)* (D from the factorization of Phi).
%Prints a table and saves <results>/conjugate_check.mat.
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper');
maxNumCompThreads(1); warning('off', 'all');
c0 = cpr_cfg('arch', 'IS');
[X, f] = pattern_crlh_vsh(c0.nw, c0.nd, c0.f0idx, 16, c0.npol);     % nw x nh x (npol*nd)
[nw, nh, ns] = size(X); nd = ns / c0.npol; x = X(:);
Pb = vsh_basis(nw, nh, c0.vshL, c0.npol);                             % (npol*n2) x p, orthonormal
stack = @(v) reshape(permute(reshape(v, nw*nh, nd, c0.npol), [1 3 2]), [], nd);   % -> (npol*n2) x nd
outside = @(v) 20*log10(norm(stack(v) - Pb*(Pb'*stack(v)), 'fro') / norm(v));
fprintf('truth outside the VSH span: %.1f dB\n', outside(x));
fprintf('plain conjugate conj(x) outside the span: %.1f dB\n', outside(conj(x)));
rng(1); Rb = zeros(1, 3);
for s = 1:3
    ph = exp(1i*2*pi*rand(1, nh));                                    % one phase per rotation (phi column), all pols and wavelengths
    xb = X .* reshape(ph, 1, nh, 1); Rb(s) = outside(xb(:));
end
fprintf('independent phase per rotation outside the span: %s dB\n', mat2str(round(Rb, 1)));
archs = {'IS', {}; 'RS', {}; 'IC', {}; 'RC', {}};
T = struct('arch', {}, 'K', {}, 'outside', {}, 'misfit', {}, 'misfitRot', {});
fprintf('\n%-4s %-3s  x# outside span (dB)   x# misfit (dB)   rotation-phase misfit (dB)\n', '', 'K');
for a = 1:size(archs, 1)
    for K = [1 4]
        o = zeros(1, 3); mf = o; mr = o;
        for s = 1:3
            c = cpr_cfg('arch', archs{a, 1}, 'K', K);
            c2 = c; [c2.nw, c2.nh, ~] = size(X); c2.nd = nd;
            S = make_readings(c2, f, s);
            z = S.Phi * x; d = coldir(S.Phi);
            xs = conj(d) .* conj(d .* x);
            o(s) = outside(xs);
            mf(s) = 20*log10(norm(abs(S.Phi*xs) - abs(z)) / norm(z));
            ph = exp(1i*2*pi*rand(1, nh)); xb = X .* reshape(ph, 1, nh, 1);
            mr(s) = 20*log10(norm(abs(S.Phi*xb(:)) - abs(z)) / norm(z));
        end
        T(end+1) = struct('arch', archs{a, 1}, 'K', K, 'outside', o, 'misfit', mf, 'misfitRot', mr); %#ok<AGROW>
        fprintf('%-4s %-3d  %-22s %-16s %s\n', archs{a, 1}, K, mat2str(round(o, 1)), mat2str(round(mf, 1)), mat2str(round(mr, 1)));
    end
end
save(fullfile(cpr_paths().results, 'conjugate_check.mat'), 'T', 'Rb');
fprintf('CONJUGATE CHECK DONE\n');
end

function d = coldir(Phi)
% unit-modulus column phases D such that Phi D^-1 is as real as possible (spanning-forest propagation)
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
d = exp(1i*ph(m+1:m+n));
end
