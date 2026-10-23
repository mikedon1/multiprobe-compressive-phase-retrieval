function [S, Y, Icur, info] = mom_wire_sparams(wires, freq, opts)
%MOM_WIRE_SPARAMS  Thin-wire method of moments (Pocklington EFIE, reduced kernel, Galerkin with triangular basis functions) for straight
%wires with a delta-gap port at the centre of every wire. Time convention e^{j omega t}, as in HFSS.
%   [S, Y, Icur, info] = mom_wire_sparams(wires, freq)
%   wires   struct array: A, B (3x1 end points, m), a (radius, m), M (number of segments, even)
%   freq    Hz (vector)
%   opts    .Z0 (50), .nq (quadrature points per segment, 8), .nqNear (points for self and adjacent segment pairs, 48)
% S, Y    nw x nw x nf port matrices (ports are the centre nodes of the wires, 1 V delta gap; Y = admittance, S for Z0)
% Icur    N x nw x nf current coefficients (triangular basis, interior nodes of each wire) for each port excited in turn, all other ports
%         shorted (voltage 0); info gives the node bookkeeping (info.portNode, info.wireOfNode, info.nodePos).
% Basis f_j (1 at node j, 0 at the neighbours); no current at the wire ends. Z = j w mu (t.t') int int f f' G + 1/(j w eps) int int f' f'' G with
% G = exp(-j k R)/(4 pi R), R = sqrt(|r - r'|^2 + a^2). Not suitable for wires closer than about 2 radii or for electrically thick wires.
if nargin < 3, opts = struct(); end
Z0 = getf(opts, 'Z0', 50); nq = getf(opts, 'nq', 8); nqNear = getf(opts, 'nqNear', 48);
c0 = 299792458; mu0 = 4e-7*pi; eps0 = 1/(mu0*c0^2);
nw = numel(wires); M = [wires.M]; if any(mod(M, 2)), error('mom_wire_sparams: M must be even'); end
segoff = [0 cumsum(M(1:end-1))]; Sg = sum(M);
st = zeros(Sg, 3); tdir = zeros(Sg, 3); dl = zeros(Sg, 1); wof = zeros(Sg, 1); arad = zeros(Sg, 1);
for w = 1:nw
    A = wires(w).A(:).'; B = wires(w).B(:).'; L = norm(B - A); t = (B - A)/L; D = L/M(w);
    ii = segoff(w) + (1:M(w)); st(ii, :) = A + ((0:M(w)-1)' * D) * t; tdir(ii, :) = repmat(t, M(w), 1); dl(ii) = D; wof(ii) = w; arad(ii) = wires(w).a;
end
% nodes: interior nodes j = 1..M-1 of each wire lie between segment j (left, f rising) and segment j+1 (right, f falling)
nodeW = []; pr = []; pd = [];
for w = 1:nw, j = (1:M(w)-1)'; nodeW = [nodeW; repmat(w, M(w)-1, 1)]; pr = [pr; segoff(w) + j]; pd = [pd; segoff(w) + j + 1]; end %#ok<AGROW>
N = numel(nodeW); portNode = zeros(1, nw);
for w = 1:nw, portNode(w) = find(nodeW == w, 1) + M(w)/2 - 1; end
[xg, wg] = gauleg(nq); ug = (xg + 1)/2; wgt = wg/2;                           % on [0, 1]
[xh, wh] = gauleg(nqNear); uh = (xh + 1)/2; wht = wh/2;
Pq = zeros(Sg, nq, 3); for d = 1:3, Pq(:, :, d) = st(:, d) + dl .* tdir(:, d) * ug.'; end
Wg = dl .* wgt.'; Wr = Wg .* ug.'; Wd = Wg .* (1 - ug).';                       % Sg x nq weights (length included)
T = tdir(pr, :) * tdir(pr, :).';                                                % t.t' for the nodes (same wire direction for both segments)
nf = numel(freq); S = zeros(nw, nw, nf); Y = S; Icur = zeros(N, nw, nf);
retZ = getf(opts, 'returnZ', false); Zall = [];                                % opts.returnZ: info.Z = N x N x nf impedance matrices
E = zeros(N, nw); for w = 1:nw, E(portNode(w), w) = 1; end
a2 = arad.^2;
for ifq = 1:nf
    k = 2*pi*freq(ifq)/c0; om = 2*pi*freq(ifq);
    Mrr = zeros(Sg); Mrd = Mrr; Mdr = Mrr; Mdd = Mrr; Gi = Mrr;
    ch = 40;
    for p0 = 1:ch:Sg
        pp = p0:min(p0 + ch - 1, Sg); cp = numel(pp);
        R2 = zeros(cp, nq, Sg, nq);
        for d = 1:3
            R2 = R2 + (reshape(Pq(pp, :, d), [cp nq 1 1]) - reshape(Pq(:, :, d), [1 1 Sg nq])).^2;
        end
        R = sqrt(R2 + reshape(a2, [1 1 Sg 1]));  G = exp(-1j*k*R) ./ (4*pi*R);
        for yy = 1:3
            switch yy, case 1, Wq = Wr; case 2, Wq = Wd; otherwise, Wq = Wg; end
            tmp = sum(G .* reshape(Wq, [1 1 Sg nq]), 4);                         % cp x nq x Sg
            for xx = 1:3
                switch xx, case 1, Wp = Wr(pp, :); case 2, Wp = Wd(pp, :); otherwise, Wp = Wg(pp, :); end
                v = reshape(sum(tmp .* Wp, 2), [cp Sg]);
                if yy == 3 && xx == 3, Gi(pp, :) = v;
                elseif yy == 1 && xx == 1, Mrr(pp, :) = v; elseif yy == 2 && xx == 1, Mrd(pp, :) = v;
                elseif yy == 1 && xx == 2, Mdr(pp, :) = v; elseif yy == 2 && xx == 2, Mdd(pp, :) = v; end
            end
        end
    end
    % self and adjacent segment pairs of the same wire: finer quadrature
    for p = 1:Sg
        for q = max(1, p-1):min(Sg, p+1)
            if wof(q) ~= wof(p), continue, end
            [Mrr(p, q), Mrd(p, q), Mdr(p, q), Mdd(p, q), Gi(p, q)] = nearpair(st(p, :), tdir(p, :), dl(p), st(q, :), tdir(q, :), dl(q), arad(p), k, uh, wht);
        end
    end
    Z = 1j*om*mu0 * T .* (Mrr(pr, pr) + Mrd(pr, pd) + Mdr(pd, pr) + Mdd(pd, pd)) ...
        + (1/(1j*om*eps0)) * ((1./dl(pr)) .* Gi(pr, pr) .* (1./dl(pr)).' - (1./dl(pr)) .* Gi(pr, pd) .* (1./dl(pd)).' ...
                              - (1./dl(pd)) .* Gi(pd, pr) .* (1./dl(pr)).' + (1./dl(pd)) .* Gi(pd, pd) .* (1./dl(pd)).');
    if retZ, Zall(:, :, ifq) = Z; end %#ok<AGROW>
    Ic = Z \ E; Yp = E.' * Ic;
    Icur(:, :, ifq) = Ic; Y(:, :, ifq) = Yp;
    S(:, :, ifq) = (eye(nw) + Z0*Yp) \ (eye(nw) - Z0*Yp);
end
nodePos = zeros(N, 3);
for n = 1:N, nodePos(n, :) = st(pd(n), :); end
segNodeL = zeros(Sg, 1); segNodeR = zeros(Sg, 1);                           % node (index into the unknowns) at the start / end of each segment, 0 at a wire end
for w = 1:nw
    nodeFirst = find(nodeW == w, 1);
    for i = 1:M(w)
        if i >= 2, segNodeL(segoff(w) + i) = nodeFirst + i - 2; end
        if i <= M(w) - 1, segNodeR(segoff(w) + i) = nodeFirst + i - 1; end
    end
end
info = struct('portNode', portNode, 'wireOfNode', nodeW, 'nodePos', nodePos, 'nodeDir', tdir(pr, :), 'N', N, 'segments', Sg, 'Z', Zall, ...
    'segStart', st, 'segDir', tdir, 'segLen', dl, 'segWire', wof, 'segNodeL', segNodeL, 'segNodeR', segNodeR);
end

function [rr, rd, dr, dd, g] = nearpair(sp, tp, Dp, sq, tq, Dq, a, k, u, w)
% integrals over a pair of nearby segments with nq points per segment (rising r = u, falling d = 1 - u)
P = sp + (u * Dp) * tp; Q = sq + (u * Dq) * tq; nq = numel(u);
R = sqrt(sum((reshape(P, [nq 1 3]) - reshape(Q, [1 nq 3])).^2, 3) + a^2); G = exp(-1j*k*R) ./ (4*pi*R);
wp = Dp * w; wq = Dq * w;
rr = (wp .* u).' * G * (wq .* u); rd = (wp .* u).' * G * (wq .* (1 - u)); dr = (wp .* (1 - u)).' * G * (wq .* u);
dd = (wp .* (1 - u)).' * G * (wq .* (1 - u)); g = wp.' * G * wq;
end

function [x, w] = gauleg(n)
i = 1:n-1; b = i ./ sqrt(4*i.^2 - 1); [V, D] = eig(diag(b, 1) + diag(b, -1)); [x, ix] = sort(diag(D)); w = 2 * V(1, ix).'.^2; x = x(:); w = w(:);
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
