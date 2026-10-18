function [X, f] = pattern_crlh(nw, nd, f0idx)
%PATTERN_CRLH  CRLH leaky-wave antenna test pattern on the paper's scan grid.
%
%   [X, f] = pattern_crlh(nw, nd, f0idx)
%
% X is nw x nw/2 x nd complex (scan theta 0-360 deg, scan phi 0-180 deg), normalized so that
% max|X| = 1, and f (1 x nd, Hz) are its frequencies, starting at frequency index f0idx
% (50 MHz steps from 1.90 GHz). The data are the circular component
% E_cprh = (E_theta + 1i*E_phi).*exp(1i*phi)/sqrt(2) of the openEMS far field in
% data/crlh_far_field.mat. A scan angle theta > 180 deg is the direction
% (360 - theta, phi + 180 deg), so each great circle through the poles is covered once.
if nargin < 3, f0idx = 1; end
D = crlh_data();
nh = nw/2;
[Th, Ph] = meshgrid(D.theta(1:end-1)*180/pi, D.phi*180/pi);   % 181 x 90 source grid (deg)
[Xq, Yq] = meshgrid(0:180/nh:180 - 180/nh, 0:360/nw:360 - 360/nw);
[nth, nph] = size(D.E_theta{1});                                  % 91 x 181
X = zeros(nw, nh, nd);
for t = 1:nd
    i = f0idx + t - 1;
    Eu = circshift(D.E_cprh{i}(:, 1:end-1), 45, 2);               % 90-degree phi shift of the 180 unique
    E = [Eu, Eu(:, 1)];                                           %   columns (0-358 deg); 360 deg re-appended
    G = zeros(nph, nth - 1);                                      % 181 x 90: scan theta 0-360
    G(1:nth, 1:nth-1) = E(1:nth, 1:nth-1);
    G(nth+1:nph, 1:nth-1) = flipud(E(1:nth-1, nth:nph-1));
    X(:, :, t) = interp2(Th, Ph, G, Xq, Yq);
end
mX = abs(X); mX = mX / max(mX(:));                          % magnitude normalized, phase kept (as in the
X = cos(angle(X)).*mX + sin(angle(X)).*mX*1i;               % earlier pattern generator, so results match it exactly)
X = X / max(abs(X(:)));
f = D.freq(f0idx + (0:nd-1)); f = f(:).';
end
