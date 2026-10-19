function D = array_data(nx, nz, dcm, steerDeg)
%ARRAY_DATA  Far field of a planar array of short dipoles, the electrically larger second example antenna
%(2026-10-19), in the same format as crlh_data (theta 0-pi, phi 0-2pi on a 2-degree grid, the 32 CRLH frequencies).
%   D = array_data()                    5 x 5 array, 6 cm spacing, beam steered 20 deg
%   D = array_data(nx, nz, dcm, steerDeg)
% Geometry (source frame): nx x nz y-directed Hertzian dipoles on a grid in the x-z plane, spacing dcm, centred
% on the origin. Excitation: a cosine-on-pedestal taper, fixed phase-shifter phases that steer the beam
% steerDeg off broadside toward +x at the band centre (so the beam squints with frequency), and small fixed random
% amplitude/phase errors (seed 1) that break the symmetries. Far field (e^{-jkr} convention, common factors
% dropped): E = sum_n w_n exp(1j k r_hat . r_n) [ (y_hat . theta_hat) theta_hat + (y_hat . phi_hat) phi_hat ].
if nargin < 1, nx = 5; end
if nargin < 2, nz = 5; end
if nargin < 3, dcm = 6; end
if nargin < 4, steerDeg = 20; end
persistent cache key
k0 = sprintf('%d_%d_%g_%g', nx, nz, dcm, steerDeg);
if ~isempty(cache) && strcmp(key, k0), D = cache; return, end
C = crlh_data();                                         % grid and frequencies of the first example
th = C.theta(:); ph = C.phi(:).';
[TH, PH] = ndgrid(th, ph);
ux = sin(TH).*cos(PH); uy = sin(TH).*sin(PH); uz = cos(TH);
eth_y = cos(TH).*sin(PH); eph_y = cos(PH);               % y_hat . theta_hat, y_hat . phi_hat
[ix, iz] = ndgrid(((1:nx) - (nx+1)/2) * dcm/100, ((1:nz) - (nz+1)/2) * dcm/100);
rs = RandStream('mt19937ar', 'Seed', 1);
taper = (0.5 + 0.5*cos(pi*ix/(nx*dcm/100))) .* (0.5 + 0.5*cos(pi*iz/(nz*dcm/100)));
err = (1 + 0.05*randn(rs, nx, nz)) .* exp(1i*(5*pi/180)*randn(rs, nx, nz));
fc = mean(C.freq(1:16)); kc = 2*pi*fc/299792458;
w = taper .* err .* exp(-1i*kc*ix*sind(steerDeg));      % fixed phase-shifter phases (set at the band centre)
D.theta = C.theta; D.phi = C.phi; D.freq = C.freq;
nf = numel(C.freq); D.E_theta = cell(1, nf); D.E_phi = cell(1, nf);
for t = 1:nf
    k = 2*pi*C.freq(t)/299792458;
    AF = zeros(size(TH));
    for n = 1:numel(ix), AF = AF + w(n) * exp(1i*k*(ux*ix(n) + uz*iz(n))); end
    D.E_theta{t} = AF .* eth_y; D.E_phi{t} = AF .* eph_y;
end
D.r0 = max(hypot(ix(:), iz(:))) + 0.005;                  % minimum-sphere radius (m), element half-size 0.5 cm
D.info = struct('nx', nx, 'nz', nz, 'dcm', dcm, 'steerDeg', steerDeg);
cache = D; key = k0;
end
