function [tp, pp, Et, Ep, n] = fwc_hfss2rot(th, ph, E)
%FWC_HFSS2ROT  HFSS-frame directions and Cartesian vectors -> the AUT rotation-axis frame.
%   [tp, pp, Et, Ep] = fwc_hfss2rot(th, ph, E)
% th, ph: HFSS spherical angles (rad, measured from HFSS +z / +x). E: N x 3 Cartesian vector in HFSS
% axes. The AUT rotates about HFSS +y, so the mode frame is x' = z, y' = x, z' = y (right-handed);
% a rotation by psi about HFSS +y is a rotation by psi about z'.
% Returns the primed angles (tp from z', pp from x') and the tangential components of E along the
% primed theta-hat and phi-hat. n: N x 3 unit directions in the primed frame.
th = th(:); ph = ph(:);
n = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
n = n(:, [3 1 2]); E = E(:, [3 1 2]);
tp = acos(max(-1, min(1, n(:, 3))));
pp = atan2(n(:, 2), n(:, 1));
pole = sin(tp) < 1e-9; pp(pole) = 0;
tt = [cos(tp).*cos(pp), cos(tp).*sin(pp), -sin(tp)];
pv = [-sin(pp), cos(pp), zeros(size(pp))];
Et = sum(E .* tt, 2); Ep = sum(E .* pv, 2);
end
