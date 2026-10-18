function D = crlh_sphere(r, nw, nd, f0idx, tol_dB)
%CRLH_SPHERE  The CRLH antenna's field on a sphere of radius r (m), from its openEMS far field
%through the vector spherical-wave expansion (two-scan technique, Section VII).
%
%   D = crlh_sphere(r, nw, nd, f0idx)      r = Inf gives the far field (paper pattern)
%
% D.X   : nw x nw/2 x nd E_cprh on the paper's grid (theta 0-360 along great circles x phi
%         0-180, same mapping and 90-degree phi shift as pattern_crlh.m),
%         normalized so max|X| = 1 over all frequencies at THIS radius (D.scale holds the factor)
% D.Et, D.Ep : (nw/2+1) x nw x nd tangential components on the physical tensor grid
%         D.theta = 0:dt:180 deg, D.phi = 0:dt:360-dt deg (dt = 360/nw), same scale as D.X
% D.c   : spherical-wave coefficients at this radius (one column per frequency), D.S the
%         expansion on the tensor grid (N = D.S.N), D.Neff per frequency, D.f, D.k
% Degrees whose far-field energy is below tol_dB (default -100) re the strongest are dropped
% (swe_trunc.m) before propagating inward.
if nargin < 5, tol_dB = -100; end
F = crlh_data();
f = F.freq(f0idx + (0:nd-1)); f = f(:).'; k = 2*pi*f/299792458;
% fit on the openEMS grid (theta 0:2:180, phi 0:2:358)
thF = F.theta(:); phF = F.phi(1:end-1); phF = phF(:);
NF = 60;
persistent SFc
if isempty(SFc), SFc = swe_setup(thF, phF, NF); end
SF = SFc;
key = @(Sx, n, ty, m) sub2ind([NF, 2, 2*NF + 1], n, ty, m + NF + 1);
L = zeros(NF, 2, 2*NF + 1); L(key(SF, SF.n, SF.type, SF.m(SF.mi).')) = 1:numel(SF.n);
% synthesis grid
dt = 360/nw; th = deg2rad(0:dt:180).'; ph = deg2rad(0:dt:360-dt).';
N = min(NF, floor((numel(ph) - 1)/2)); N = min(N, numel(th) - 1);
S = swe_setup(th, ph, N);
D.Et = zeros(numel(th), numel(ph), nd); D.Ep = D.Et; D.c = zeros(numel(S.n), nd); D.Neff = zeros(1, nd);
for t = 1:nd
    i = f0idx + t - 1;
    cF = swe_fit(SF, F.E_theta{i}(:, 1:end-1), F.E_phi{i}(:, 1:end-1));
    [cF, D.Neff(t)] = swe_trunc(SF, cF, tol_dB);
    if D.Neff(t) > N
        warning('crlh_sphere: f = %.2f GHz needs N = %d > %d supported by the nw = %d grid', f(t)/1e9, D.Neff(t), N, nw);
    end
    % map the degree-<=N part of the coefficients onto the synthesis expansion S
    c = cF(L(key(S, S.n, S.type, S.m(S.mi).')));
    c = c .* swe_radial(S, k(t), r) ./ swe_radial(S, k(t), Inf);
    D.c(:, t) = c;
    [D.Et(:,:,t), D.Ep(:,:,t)] = swe_synth(S, c);
end
% paper grid: great circles, theta_g 0..360-dt along phi_g 0..180-dt, 90-degree phi shift
nh = nw/2; tg = (0:nw-1)*dt; pg = (0:nh-1)*dt;
D.X = zeros(nw, nh, nd);
[TG, PG] = ndgrid(tg, pg);
up = TG > 180;
TH = TG; TH(up) = 360 - TG(up);
PB = PG; PB(up) = PG(up) + 180;
PH = mod(PB - 90, 360);
it = round(TH/dt) + 1; ip = round(PH/dt) + 1;
lin = sub2ind([numel(th) numel(ph)], it, ip);
for t = 1:nd
    E = swe_cprh(D.Et(:,:,t), D.Ep(:,:,t), ph.');
    D.X(:,:,t) = E(lin);
end
% Dual polarization (Section VII): the components a dual-polarized probe on the arc measures,
% along the arc tangent (theta port) and the arc normal (phi port). On the far half of each
% great circle (theta_g > 180) both point opposite to the physical theta-hat / phi-hat, so the
% arc components are the physical ones times -1 there (continuous along each great circle).
D.arcSign = ones(size(TG)); D.arcSign(up) = -1;
D.Xt = zeros(nw, nh, nd); D.Xp = D.Xt;
for t = 1:nd
    [D.Xt(:,:,t), D.Xp(:,:,t)] = sphere_to_arc(D.Et(:,:,t), D.Ep(:,:,t), lin, D.arcSign);
end
D.scale = max(abs(D.X(:)));
D.X = D.X / D.scale; D.Et = D.Et / D.scale; D.Ep = D.Ep / D.scale; D.c = D.c / D.scale;
D.Xt = D.Xt / D.scale; D.Xp = D.Xp / D.scale;
D.S = S; D.f = f; D.k = k; D.r = r; D.theta = rad2deg(th).'; D.phi = rad2deg(ph).'; D.map = lin;
end
