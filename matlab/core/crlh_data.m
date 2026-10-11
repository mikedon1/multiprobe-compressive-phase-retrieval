function D = crlh_data()
%CRLH_DATA  Load the CRLH far-field data once (cached in memory).
%   D.theta (91x1, rad, 0-pi), D.phi (181x1, rad, 0-2pi), D.freq (32x1, Hz),
%   D.E_theta, D.E_phi, D.E_cprh : 1 x 32 cells of 91 x 181 complex fields.
persistent Dc
if isempty(Dc)
    S = load(fullfile(cpr_paths().data, 'crlh_far_field.mat'));
    c = S.crlh;
    P = repmat(c.phi(:).', numel(c.theta), 1);                   % phi grid (varies along columns)
    c.E_cprh = cellfun(@(et, ep) (et + 1i*ep).*exp(1i*P)/sqrt(2), c.E_theta, c.E_phi, 'UniformOutput', false);
    Dc = c;
end
D = Dc;
end
