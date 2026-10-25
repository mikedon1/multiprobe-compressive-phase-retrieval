function Tm = mom_make_ring(Vm, opts)
%MOM_MAKE_RING  True full-ring data of Section VII from the thin-wire method of moments (64 probes + AUT, 16 rotations over 180 deg).
%   Tm = mom_make_ring(Vm, opts)
% Vm    a mom_make_validation result (arch): its AUT modes Vm.A, port map Vm.map, frequencies Vm.freq and arch responses Vm.R are reused
% opts  .Mprobe (8), .Maut (16), .Z0 (50), .L (12), .psiDeg (0:11.25:168.75), .verbose (true)
% Tm.Xring  64 x nrot x nf  AUT -> probe transmissions from the exact solution of the AUT and all 64 probes (others loaded with Z0);
%           rows 1..32 near half in the order of Vm.map, rows 33..64 far half in the same order (the layout of fwc_ring_true)
% Tm.psi    rotation angles (rad)
% Tm.Rring  64 x nmode x nf  embedded responses of the ports of the FULL ring alone (reciprocity, other probes loaded), same row order:
%           the calibration of the real ring, including the coupling between the halves
% Tm.Rarch  the same rows from the arch response and the 180-deg symmetry, [R; R.*(-1)^m] (the calibration used with the HFSS data)
% Tm.A, Tm.freq, Tm.map as in Vm; Tm.nmse_db / Tm.nmse_arch_db: model error of the truth (one factor per frequency) with Rring / Rarch.
if nargin < 2, opts = struct(); end
Mp = gf(opts, 'Mprobe', 8); Ma = gf(opts, 'Maut', 16); Z0 = gf(opts, 'Z0', 50); L = gf(opts, 'L', 12);
psiDeg = gf(opts, 'psiDeg', (0:15) * 11.25); verbose = gf(opts, 'verbose', true);
f = Vm.freq(:).'; nf = numel(f); np = height(Vm.map); nrot = numel(psiDeg); A = Vm.A; nmode = size(A.modes, 1);
c0 = 299792458;
wopt = struct('ring', 'full', 'Mprobe', Mp, 'Maut', Ma);
% ---- 1. transmissions from the exact solution, rotation by rotation ----
[~, lab0] = mom_ring_wires(0, setfield(wopt, 'aut', true)); %#ok<SFLD>
ord = rowOrder(lab0, Vm.map, np);
Xring = zeros(2*np, nrot, nf);
for ir = 1:nrot
    t0 = tic; wires = mom_ring_wires(psiDeg(ir), setfield(wopt, 'aut', true)); %#ok<SFLD>
    S = mom_wire_sparams(wires, f, struct('Z0', Z0));
    Xring(:, ir, :) = reshape(S(ord, 1, :), 2*np, 1, nf);
    if verbose, fprintf('MoM ring rotation %2d of %d (psi %.2f deg): %.0f s\n', ir, nrot, psiDeg(ir), toc(t0)); end
end
% ---- 2. embedded responses of the 64 ports of the ring alone ----
[wR, labR] = mom_ring_wires(0, setfield(wopt, 'aut', false)); %#ok<SFLD>
ordR = rowOrder(labR, Vm.map, np);
[~, ~, ~, iR2] = mom_wire_sparams(wR, f, struct('returnZ', true));
[xq, wq] = gauleg(6); uq = (xq + 1)/2; wqt = wq/2;
N = iR2.N; Er = zeros(N, numel(wR)); for w = 1:numel(wR), Er(iR2.portNode(w), w) = 1; end
nSeg = iR2.segments; Pq = zeros(nSeg, numel(uq), 3);
for d = 1:3, Pq(:, :, d) = iR2.segStart(:, d) + (iR2.segLen .* iR2.segDir(:, d)) * uq.'; end
wseg = iR2.segLen .* wqt.';
pts = reshape(Pq, [], 3); pp = [pts(:, 3), pts(:, 1), pts(:, 2)]; rr = sqrt(sum(pp.^2, 2)); tpp = acos(pp(:, 3)./rr); ppp = atan2(pp(:, 2), pp(:, 1));
thh = [cos(tpp).*cos(ppp), cos(tpp).*sin(ppp), -sin(tpp)]; phh = [-sin(ppp), cos(ppp), zeros(size(ppp))];
td = repmat(reshape(iR2.segDir, [nSeg 1 3]), 1, numel(uq), 1); td = reshape(td, [], 3); tdp = [td(:, 3), td(:, 1), td(:, 2)];
tT = sum(tdp .* thh, 2); tP = sum(tdp .* phh, 2);
Rring = zeros(2*np, nmode, nf);
for q = 1:nf
    k = 2*pi*f(q)/c0;
    Z = iR2.Z(:, :, q); Ir = (Z + Z0*(Er*Er.')) \ Er;
    [Et, Ep] = sph_modes(L, k, rr, tpp, ppp, 'out');
    Em = Et .* tT + Ep .* tP;
    for w = 1:numel(wR)
        Iw = Ir(:, w); Isg = zeros(nSeg, numel(uq));
        for s = 1:nSeg
            iL = iR2.segNodeL(s); iRt = iR2.segNodeR(s); IL = 0; IRt = 0;
            if iL > 0, IL = Iw(iL); end
            if iRt > 0, IRt = Iw(iRt); end
            Isg(s, :) = IL*(1 - uq.') + IRt*uq.';
        end
        wv = reshape(Isg .* wseg, [], 1);
        ipos = find(ordR == w, 1); if isempty(ipos), continue, end
        Rring(ipos, :, q) = wv.' * Em;
    end
end
sgn = reshape((-1) .^ A.modes(:, 2), 1, []);
Rarch = [Vm.R; Vm.R .* sgn];
Tm = struct('Xring', Xring, 'psi', psiDeg * pi/180, 'Rring', Rring, 'Rarch', Rarch, 'A', A, 'freq', f, 'map', Vm.map);
% ---- self check: model of the paper against the exact solution, both calibrations ----
Tm.nmse_db = modelErr(Tm, Rring); Tm.nmse_arch_db = modelErr(Tm, Rarch);
if verbose
    fprintf('MoM ring: model error of the truth, ring calibration %.1f to %.1f dB, arch+mirror calibration %.1f to %.1f dB (one factor per frequency)\n', ...
        min(Tm.nmse_db), max(Tm.nmse_db), min(Tm.nmse_arch_db), max(Tm.nmse_arch_db));
end
end

function e = modelErr(Tm, R)
nf = numel(Tm.freq); [np2, nrot, ~] = size(Tm.Xring); modes = Tm.A.modes; nm = size(modes, 1); e = zeros(1, nf);
for f = 1:nf
    Gf = zeros(np2*nrot, nm);
    for l = 1:nrot, Gf((l-1)*np2 + (1:np2), :) = R(:, :, f) .* exp(-1i * modes(:, 2).' * Tm.psi(l)); end
    xm = Gf * Tm.A.a(:, f); xd = reshape(Tm.Xring(:, :, f), [], 1); k = (xm' * xd) / (xm' * xm);
    e(f) = 20*log10(norm(k*xm - xd) / norm(xd));
end
end

function ord = rowOrder(lab, map, np)
% wire index of every row: near half in the order of map, then far half in the same order
ord = zeros(1, 2*np);
for h = 1:2
    hs = 'near'; if h == 2, hs = 'far'; end
    for r = 1:np
        j = find(strcmp({lab.kind}, 'probe') & strcmp({lab.half}, hs) & strcmp({lab.pol}, map.pol(r)) & [lab.idx] == map.ring_index(r));
        if numel(j) ~= 1, error('mom_make_ring: probe %s %s %d not unique', hs, map.pol(r), map.ring_index(r)); end
        ord((h-1)*np + r) = j;
    end
end
end

function [x, w] = gauleg(n)
i = 1:n-1; b = i ./ sqrt(4*i.^2 - 1); [V, D] = eig(diag(b, 1) + diag(b, -1)); [x, ix] = sort(diag(D)); w = 2 * V(1, ix).'.^2; x = x(:); w = w(:);
end

function v = gf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
