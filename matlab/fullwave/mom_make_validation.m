function V = mom_make_validation(Vh, opts)
%MOM_MAKE_VALIDATION  The full-wave data of Section VII from the thin-wire method of moments instead of HFSS: the same structure as fwc_validate_hfss.
%   V = mom_make_validation(Vh, opts)
% Vh      an fwc_validate_hfss result (only its port map V.map, rotation angles V.psi (deg) and frequencies V.freq are used, so the ports, rotations and
%         frequencies are those of the paper)
% opts    .Mprobe (8) and .Maut (16) segments per probe and for the AUT, .Z0 (50), .L (12), .Lreg unused, .rotations (default all of Vh.psi), .verbose
% V.Smeas   nport x nrot x nf: transmissions from the AUT to the probes, from the exact full solution (AUT + all probes, others loaded with Z0)
% V.R       nport x nmode x nf: embedded responses of the probes to the outgoing modes. For each probe, the reaction of its array-alone current (the other probes
%           loaded) with the field of the mode along the wire (Lorentz reciprocity; tangential field components, the probes are tangent to the sphere)
% V.A       AUT modes (fwc_aut_modes) of the far field of the AUT-alone current, on the grid of the HFSS export
% Same constants per frequency as the HFSS route (they are absorbed by the one factor per frequency of the recovery and its score).
if nargin < 2, opts = struct(); end
Mp = gf(opts, 'Mprobe', 8); Ma = gf(opts, 'Maut', 16); Z0 = gf(opts, 'Z0', 50); L = gf(opts, 'L', 12); verbose = gf(opts, 'verbose', true);
f = Vh.freq(:).'; nf = numel(f); psi = Vh.psi(:).'; rot = gf(opts, 'rotations', 1:numel(psi)); nrot = numel(rot); np = height(Vh.map);
c0 = 299792458; mu0 = 4e-7*pi;
wopt = struct('ring', 'arch', 'Mprobe', Mp, 'Maut', Ma);
% ---- 1. transmissions from the exact full solution, rotation by rotation ----
[w0, lab0] = mom_ring_wires(0, setfield(wopt, 'aut', true)); %#ok<SFLD>
ord = zeros(1, np); for r = 1:np, ord(r) = find(strcmp({lab0.pol}, Vh.map.pol(r)) & [lab0.idx] == Vh.map.ring_index(r) & strcmp({lab0.kind}, 'probe')); end
Smeas = zeros(np, nrot, nf);
for ir = 1:nrot
    t0 = tic; wires = mom_ring_wires(psi(rot(ir)), setfield(wopt, 'aut', true)); %#ok<SFLD>
    S = mom_wire_sparams(wires, f, struct('Z0', Z0));
    Smeas(:, ir, :) = reshape(S(ord, 1, :), np, 1, nf);
    if verbose, fprintf('MoM rotation %2d of %d (psi %.2f deg): %.0f s\n', ir, nrot, psi(rot(ir)), toc(t0)); end
end
% ---- 2. AUT modes from the far field of the AUT-alone current ----
wA = mom_ring_wires(0, struct('ring', 'none', 'aut', true, 'Maut', Ma));
[~, ~, IA, iA] = mom_wire_sparams(wA, f);
Pa = load(fullfile(cpr_paths().data, 'hfss_fullarc', 'hfss_aut.mat'));                      % grid of the HFSS far-field export
th = Pa.theta(:) * pi/180; ph = Pa.phi(:).' * pi/180; [TH, PH] = ndgrid(th, ph);
rh = [sin(TH(:)).*cos(PH(:)), sin(TH(:)).*sin(PH(:)), cos(TH(:))]; tth = [cos(TH(:)).*cos(PH(:)), cos(TH(:)).*sin(PH(:)), -sin(TH(:))]; pph = [-sin(PH(:)), cos(PH(:)), zeros(numel(TH), 1)];
[xq, wq] = gauleg(6); uq = (xq + 1)/2; wqt = wq/2;
Tf.x = Pa.theta(:); Tf.phi = Pa.phi(:); Tf.freq = f(:); Tf.c.rETheta = zeros(numel(th), numel(ph), nf); Tf.c.rEPhi = Tf.c.rETheta;
for q = 1:nf
    k = 2*pi*f(q)/c0; om = 2*pi*f(q); Cur = zeros(numel(TH), 3);
    [Iseg, P, tq, wt] = seg_current(iA, IA(:, 1, q), uq, wqt);                            % values along the AUT segments
    for s = 1:size(P, 1)
        ph0 = exp(1j*k*(rh * P(s, :).')); Cur = Cur + (Iseg(s) * wt(s)) * ph0 .* tq(s, :);
    end
    Cperp = Cur - sum(Cur .* rh, 2) .* rh;                                              % transverse part of the radiation integral
    pref = -1j*om*mu0/(4*pi);
    Tf.c.rETheta(:, :, q) = reshape(pref * sum(Cperp .* tth, 2), numel(th), numel(ph)); Tf.c.rEPhi(:, :, q) = reshape(pref * sum(Cperp .* pph, 2), numel(th), numel(ph));
end
A = fwc_aut_modes(Tf, L);
% ---- 3. embedded responses of the probes (array alone) ----
wR = mom_ring_wires(0, setfield(wopt, 'aut', false)); [~, labR] = mom_ring_wires(0, setfield(wopt, 'aut', false)); %#ok<SFLD>
ordR = zeros(1, np); for r = 1:np, ordR(r) = find(strcmp({labR.pol}, Vh.map.pol(r)) & [labR.idx] == Vh.map.ring_index(r)); end
[~, ~, ~, iR] = mom_wire_sparams(wR, f(1), struct('returnZ', false));                    % geometry bookkeeping only
nmode = size(A.modes, 1); R = zeros(np, nmode, nf);
[Spp, ~, ~, iR2] = mom_wire_sparams(wR, f, struct('returnZ', true)); %#ok<ASGLU>
N = iR2.N; Er = zeros(N, numel(wR)); for w = 1:numel(wR), Er(iR2.portNode(w), w) = 1; end
nSeg = iR2.segments; Pq = zeros(nSeg, numel(uq), 3);
for d = 1:3, Pq(:, :, d) = iR2.segStart(:, d) + (iR2.segLen .* iR2.segDir(:, d)) * uq.'; end
wseg = iR2.segLen .* wqt.';                                                             % nSeg x nq weights (length included)
% points in the rotation-axis frame: (x', y', z') = (z, x, y) of the HFSS frame
pts = reshape(Pq, [], 3); pp = [pts(:, 3), pts(:, 1), pts(:, 2)]; rr = sqrt(sum(pp.^2, 2)); tpp = acos(pp(:, 3)./rr); ppp = atan2(pp(:, 2), pp(:, 1));
thh = [cos(tpp).*cos(ppp), cos(tpp).*sin(ppp), -sin(tpp)]; phh = [-sin(ppp), cos(ppp), zeros(size(ppp))];
td = repmat(reshape(iR2.segDir, [nSeg 1 3]), 1, numel(uq), 1); td = reshape(td, [], 3); tdp = [td(:, 3), td(:, 1), td(:, 2)];
tT = sum(tdp .* thh, 2); tP = sum(tdp .* phh, 2);
for q = 1:nf
    k = 2*pi*f(q)/c0;
    Z = iR2.Z(:, :, q); Ir = (Z + Z0*(Er*Er.')) \ Er;                                    % probe w driven (1 V), the others loaded with Z0
    [Et, Ep] = sph_modes(L, k, rr, tpp, ppp, 'out');                                    % tangential fields of the outgoing modes at the wire points
    Em = Et .* tT + Ep .* tP;                                                           % npts x nmode: component along the wire
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
        R(ipos, :, q) = wv.' * Em;
    end
end
V.Smeas = Smeas; V.R = R; V.A = A; V.psi = psi(rot); V.freq = f; V.map = Vh.map;
% ---- self check: the model of the paper against the exact solution ----
Sp = fwc_predict(A.a, A.modes, R, psi(rot)*pi/180, 1); [nm, kap, ~] = fwc_fit_scale(Sp, Smeas); %#ok<ASGLU>
V.nmse_db = nm; V.kappa = kap; V.Spred = Sp;
if verbose, fprintf('MoM validation: AUT fit residual %.1f to %.1f dB; model vs exact solution, one factor per frequency: %.1f to %.1f dB\n', 20*log10(min(A.relres)), 20*log10(max(A.relres)), min(nm(:)), max(nm(:))); end
end

function [Iseg, P, tq, wt] = seg_current(info, I, uq, wqt)
% current at the quadrature points of every segment (triangular basis), positions, directions, weights (length included)
nS = info.segments; nq = numel(uq); Iseg = zeros(nS*nq, 1); P = zeros(nS*nq, 3); tq = P; wt = zeros(nS*nq, 1); k = 0;
for s = 1:nS
    iL = info.segNodeL(s); iR = info.segNodeR(s); IL = 0; IR = 0; if iL > 0, IL = I(iL); end, if iR > 0, IR = I(iR); end
    for i = 1:nq
        k = k + 1; Iseg(k) = IL*(1 - uq(i)) + IR*uq(i); P(k, :) = info.segStart(s, :) + uq(i) * info.segLen(s) * info.segDir(s, :); tq(k, :) = info.segDir(s, :); wt(k) = info.segLen(s) * wqt(i);
    end
end
end

function [x, w] = gauleg(n)
i = 1:n-1; b = i ./ sqrt(4*i.^2 - 1); [V, D] = eig(diag(b, 1) + diag(b, -1)); [x, ix] = sort(diag(D)); w = 2 * V(1, ix).'.^2; x = x(:); w = w(:);
end

function v = gf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
