function [X, g] = dipole_ic_pattern(L_lambda, alpha_deg, varargin)
%DIPOLE_IC_PATTERN  Idealized IC measurement "true pattern" for the redesign.
%
%   [X, g] = DIPOLE_IC_PATTERN(L_lambda, alpha_deg, ...)
%
%   Complex probe responses X (nTheta x nPhi) for a centre-fed thin dipole
%   AUT measured by the Section VI geometry (also the geometry of the HFSS model):
%
%     rotation axis +y; probes on an arc in the y-z plane, theta from +y;
%     probe i at  p_i = R (0, cos th_i, sin th_i),  axis a_i = (0, -sin th_i, cos th_i)
%     (tangent to the arc = perpendicular to its sweep circle).
%     AUT at the origin, axis d = (0, cos alpha, sin alpha): along +y tilted by
%     alpha toward +z. Rotation L turns the AUT by psi_L = (L-1)*dpsi,
%     right-handed about +y.
%
%   X(i,L) = a_i' . E_AUT(rotated by psi_L)(p_i), evaluated in the AUT frame
%   as b . E(q) with q = Ry(-psi_L) p_i, b = Ry(-psi_L) a_i. Probe
%   polarization is included, not assumed.
%
%   Methods ('Method'):
%     'analytic' (default)  exact far field of a sinusoidal-current dipole,
%         E ~ [cos(kL/2 cos t) - cos(kL/2)] / sin^2 t * (cos t r - d) * e^{-jkr}/r
%     'integral'  numerical radiation integral of the same current with exact
%         source-to-probe distances (-jw A_perp). Slower; used to validate
%         the analytic form at the finite probe radius.
%
%   Options:
%     'NTheta' 16, 'NPhi' 32, 'ThetaDeg' [20 160], 'PhiSpanDeg' 360
%         (dpsi = PhiSpanDeg/NPhi), 'Rcm' 100, 'LambdaCm' 10.9,
%     'ProbePsiDeg' 0  rotate the whole probe arc about +y by this azimuth
%         before measuring -- used to place a LEG probe (psi = +k*dpsi)
%         and check leg(L) == anchor(L-k) (validate_dipole_pattern.m).
%     'OffsetCm' [0 0 0]  AUT displacement from the rotation centre.
%     'NSeg' 400  segments for 'integral'.
%
%   g returns the geometry (thetas, psis, probe positions/axes, dipole axis)
%   and the covered solid angle.

    p = inputParser;
    p.addParameter('NTheta', 16);
    p.addParameter('NPhi', 32);
    p.addParameter('ThetaDeg', [20 160]);
    p.addParameter('PhiSpanDeg', 360);
    p.addParameter('Rcm', 100);
    p.addParameter('LambdaCm', 10.9);
    p.addParameter('Method', 'analytic');
    p.addParameter('ProbePsiDeg', 0);
    p.addParameter('OffsetCm', [0 0 0]);
    p.addParameter('NSeg', 400);
    p.parse(varargin{:});
    o = p.Results;

    k   = 2*pi / o.LambdaCm;
    Lcm = L_lambda * o.LambdaCm;
    th  = linspace(o.ThetaDeg(1), o.ThetaDeg(2), o.NTheta).';
    dpsi = o.PhiSpanDeg / o.NPhi;
    psi = (0:o.NPhi-1) * dpsi;

    ca = cosd(alpha_deg); sa = sind(alpha_deg);
    d  = [0; ca; sa];

    % probe arc, optionally rotated by ProbePsiDeg about +y
    P0 = o.Rcm * [zeros(1, o.NTheta); cosd(th).'; sind(th).'];   % 3 x nTheta
    A0 = [zeros(1, o.NTheta); -sind(th).'; cosd(th).'];
    P0 = roty(o.ProbePsiDeg) * P0;
    A0 = roty(o.ProbePsiDeg) * A0;

    X = complex(zeros(o.NTheta, o.NPhi));
    off = o.OffsetCm(:);
    for L = 1:o.NPhi
        Rm = roty(-psi(L));
        Q = Rm * P0 - off;          % probe positions relative to the AUT, AUT frame
        B = Rm * A0;
        switch lower(o.Method)
            case 'analytic'
                X(:, L) = analytic_response(Q, B, d, k, Lcm);
            case 'integral'
                X(:, L) = integral_response(Q, B, d, k, Lcm, o.NSeg);
            otherwise
                error('dipole_ic_pattern:method', 'Unknown Method %s', o.Method);
        end
    end

    g = struct();
    g.theta = th; g.psi = psi; g.dpsi = dpsi;
    g.P = P0; g.A = A0; g.d = d;
    g.L_lambda = L_lambda; g.alpha = alpha_deg;
    g.Rcm = o.Rcm; g.LambdaCm = o.LambdaCm;
    g.pathRadius = o.Rcm * sind(th);
    g.solidAngleFrac = (cosd(o.ThetaDeg(1)) - cosd(o.ThetaDeg(2))) / 2 * (o.PhiSpanDeg / 360);
end


function R = roty(deg)
    c = cosd(deg); s = sind(deg);
    R = [c 0 s; 0 1 0; -s 0 c];
end


function s = analytic_response(Q, B, d, k, Lcm)
    r  = sqrt(sum(Q.^2, 1));
    rh = Q ./ r;
    u  = d.' * rh;                     % cos(theta')
    sin2 = 1 - u.^2;
    a = k * Lcm / 2;
    num = cos(a * u) - cos(a);
    F = zeros(size(u));
    ok = sin2 > 1e-10;
    F(ok) = num(ok) ./ sin2(ok);
    % on-axis limit (u -> +-1): num/sin2 -> a*sin(a)/2
    F(~ok) = a * sin(a) / 2;
    % E direction (cos t * r_hat - d), times spherical-wave factor
    Evec = rh .* u - d;
    Evec = Evec .* (F .* exp(-1i * k * r) ./ r);
    s = sum(B .* Evec, 1).';
end


function s = integral_response(Q, B, d, k, Lcm, nseg)
    ds = Lcm / nseg;
    z = ((1:nseg) - 0.5) * ds - Lcm/2;          % segment centres along d
    I = sin(k * (Lcm/2 - abs(z)));
    np = size(Q, 2);
    s = complex(zeros(np, 1));
    for j = 1:np
        q = Q(:, j);
        src = d * z;                            % 3 x nseg
        Rv = q - src;
        Rn = sqrt(sum(Rv.^2, 1));
        Aj = sum(I .* exp(-1i * k * Rn) ./ Rn) * ds * d;   % vector potential (const dropped)
        rh = q / norm(q);
        Aperp = Aj - (rh.' * Aj) * rh;
        E = -1i * Aperp;                        % omega*mu/(4 pi) dropped
        s(j) = B(:, j).' * E;
    end
end
