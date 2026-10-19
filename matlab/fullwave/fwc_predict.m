function S = fwc_predict(a, modes, R, psi, sense)
%FWC_PREDICT  Predicted AUT->probe transmission from AUT mode coefficients and probe responses.
%   S = fwc_predict(a, modes, R, psi, sense)
% a: nmode x nf AUT mode coefficients (fwc_aut_modes .a, AUT at rotation 0). modes: [l m type].
% R: nport x nmode x nf embedded probe responses (fwc_probe_resp, one row per port).
% psi: rotation angles (rad). sense: +1 if the AUT rotation by psi multiplies mode (l, m) by
% exp(-1i*m*psi) (a right-handed rotation about z'), -1 for the opposite sense.
% S: nport x numel(psi) x nf, up to one complex constant per frequency.
m = modes(:, 2);
[np, ~, nf] = size(R);
S = zeros(np, numel(psi), nf);
for f = 1:nf
    S(:, :, f) = R(:, :, f) * (a(:, f) .* exp(-1i * sense * m * psi(:).'));
end
end
