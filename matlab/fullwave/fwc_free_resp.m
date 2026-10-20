function Rf = fwc_free_resp(V)
%FWC_FREE_RESP  Free-space probe model for Section VII (round 7, 2026-10-20): every port is an ideal point probe at its
%nominal position, reading the field of each outgoing AUT mode along its polarization, with no coupling or scattering.
%   Rf = fwc_free_resp(V)        V from fwc_validate_hfss; Rf has the size and port/mode order of V.R
% Geometry from step84_build_dualpol_fullarc.py (gated positions): probe i of the ring at theta_i = 15 + 10 i deg in
% the HFSS x = 0 plane, theta-hat port at R = 100 cm, x-directed phi-hat port at R = 101 cm. In the mode frame of
% fwc_hfss2rot (x' = z, y' = x, z' = y; the AUT turns about z' = HFSS y) this is (r, theta' = theta_i, phi' = 0),
% with the theta port along theta-hat' and the phi port along phi-hat'. All theta ports are copies of one anchor and
% all phi ports copies of a cross seed made from it, and both map their positive axis to the negative unit vector,
% so the polarity is common to all ports and is absorbed by the per-frequency factor.
% Rf(p, q, f) = u_p . E_q^out(r_p) with E_q^out the outgoing mode q of sph_modes ('out'), the same modes and order
% as fwc_probe_resp. The scale per frequency differs from the reciprocity responses by one constant per frequency.
c0 = 299792458;
M = V.map; np = height(M); nm = size(V.R, 2); nf = numel(V.freq);
L = max(V.A.modes(:, 1));
th = (15 + 10 * M.ring_index(:)) * pi/180;
isPhi = strcmp(string(M.pol(:)), "phi");
r = 1.00 + 0.01 * isPhi;                          % m
Rf = zeros(np, nm, nf);
for f = 1:nf
    k = 2*pi*V.freq(f) / c0;
    for p = 1:np
        [Et, Ep, ~, ~, modes] = sph_modes(L, k, r(p), th(p), 0, 'out');
        if f == 1 && p == 1 && ~isequal(modes(:, 1:3), V.A.modes(:, 1:3))
            error('fwc_free_resp: mode order differs from the AUT modes');
        end
        if isPhi(p), Rf(p, :, f) = Ep; else, Rf(p, :, f) = Et; end
    end
end
end
