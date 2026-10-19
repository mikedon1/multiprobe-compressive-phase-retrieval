function A = cpr_operator(c, Phi, nw, nh, ns, M)
%CPR_OPERATOR  Sensing operator H = Phi * Psi for the angular basis c.angBasis.
%   'dft' : 2D DFT over the scan grid (SparseDftLinTrans)
%   'vsh' : far-field vector spherical harmonics up to degree c.vshL (VshLinTrans); with
%           c.npol = 2 both polarizations share the mode coefficients
ang = 'dft'; if isfield(c, 'angBasis'), ang = c.angBasis; end
npol = 1; if isfield(c, 'npol'), npol = c.npol; end
switch ang
    case 'dft', A = SparseDftLinTrans(Phi, nw, nh, ns, M, c.sqMode);
    case 'vsh', A = VshLinTrans(Phi, nw, nh, ns, M, c.vshL, npol);
    otherwise, error('cpr_operator: unknown angBasis %s', ang);
end
end
