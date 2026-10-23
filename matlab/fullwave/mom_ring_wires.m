function [wires, lab] = mom_ring_wires(psiDeg, opts)
%MOM_RING_WIRES  Thin-wire geometry of the Section VII model: AUT (a wire dipole tilted 10 deg from the rotation axis y, offset (2, 1, 2) cm at psi = 0,
%rotated about y by psi) and the probe dipoles of the arch or the full ring (theta-hat dipoles at R = 100 cm, phi-hat dipoles along x at R = 101 cm,
%theta = 15:10:165 deg; the far half of the full ring is the near half turned 180 deg about y).
%   [wires, lab] = mom_ring_wires(psiDeg, opts)
% opts: .ring 'arch' (32 probes) | 'full' (64) | 'none'; .aut true (default); .Mprobe (6), .Maut (14), .a (0.082 cm, in m), .Lprobe (4.91 cm), .Laut (13.625 cm)
% lab: struct array with fields kind ('aut'|'probe'), half ('near'|'far'|''), pol ('theta'|'phi'|''), idx (0..15); wire 1 is the AUT when present.
if nargin < 2, opts = struct(); end
ring = gf(opts, 'ring', 'arch'); useAut = gf(opts, 'aut', true); Mp = gf(opts, 'Mprobe', 6); Ma = gf(opts, 'Maut', 14); a = gf(opts, 'a', 0.082e-2);
Lp = gf(opts, 'Lprobe', 4.91e-2); La = gf(opts, 'Laut', 13.625e-2);
wires = struct('A', {}, 'B', {}, 'a', {}, 'M', {}); lab = struct('kind', {}, 'half', {}, 'pol', {}, 'idx', {});
if useAut
    C0 = [2 1 2] * 1e-2; u0 = [0 cosd(10) sind(10)]; ps = psiDeg*pi/180;
    ry = @(p) [p(1)*cos(ps) + p(3)*sin(ps), p(2), -p(1)*sin(ps) + p(3)*cos(ps)];
    C = ry(C0); u = ry(u0);
    wires(end+1) = struct('A', (C - La/2*u).', 'B', (C + La/2*u).', 'a', a, 'M', Ma); lab(end+1) = struct('kind', 'aut', 'half', '', 'pol', '', 'idx', 0);
end
halves = {}; if strcmp(ring, 'arch'), halves = {'near'}; elseif strcmp(ring, 'full'), halves = {'near', 'far'}; end
for h = 1:numel(halves)
    sg = 1; if strcmp(halves{h}, 'far'), sg = -1; end
    for i = 0:15
        th = (15 + 10*i)*pi/180; ph = sg*th;                                 % position angle in the y-z plane, from +y toward +z
        for pol = {'theta', 'phi'}
            if strcmp(pol{1}, 'theta'), R = 1.00; ax = [0 -sin(ph) cos(ph)]; else, R = 1.01; ax = [1 0 0]; end
            c = [0 R*cos(ph) R*sin(ph)];
            wires(end+1) = struct('A', (c - Lp/2*ax).', 'B', (c + Lp/2*ax).', 'a', a, 'M', Mp); %#ok<AGROW>
            lab(end+1) = struct('kind', 'probe', 'half', halves{h}, 'pol', pol{1}, 'idx', i); %#ok<AGROW>
        end
    end
end
end

function v = gf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
