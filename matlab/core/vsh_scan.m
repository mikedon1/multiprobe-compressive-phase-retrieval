function [th, phs, sg] = vsh_scan(nw, nh)
%VSH_SCAN  True directions of the scan grid of pattern_crlh.m (nw x nh, column-major).
%   [th, phs, sg] = vsh_scan(nw, nh)
% Scan angle (theta_s, phi_s): for theta_s <= 180 deg the direction is (theta_s, phi_s), otherwise
% (360 - theta_s, phi_s + 180 deg) with the across-circle unit vector reversed (sg = -1). phs is
% the azimuth in the source (openEMS) frame, which pattern_crlh rotates by 90 deg.
if nargin < 2, nh = nw/2; end
[ts, ps] = ndgrid((0:nw-1)*2*pi/nw, (0:nh-1)*2*pi/nw);
ts = ts(:); ps = ps(:);
far = ts > pi + 1e-9;
th = ts; th(far) = 2*pi - ts(far);
p = ps; p(far) = ps(far) + pi;
phs = mod(p - pi/2, 2*pi);
sg = ones(size(ts)); sg(far) = -1;
end
