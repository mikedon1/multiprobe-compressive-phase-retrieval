function [nmse_db, kappa, perPort] = fwc_fit_scale(Spred, Smeas)
%FWC_FIT_SCALE  Compare predicted and measured transmissions up to one complex constant per frequency.
%   [nmse_db, kappa, perPort] = fwc_fit_scale(Spred, Smeas)
% Spred, Smeas: nport x nrot x nf. kappa(f) = least-squares scale of Spred onto Smeas.
% nmse_db(f) = 10 log10 ||kappa Spred - Smeas||^2 / ||Smeas||^2 over all ports and rotations;
% perPort (nport x nf, dB) the same per port (normalized by that port's own energy).
nf = size(Smeas, 3); kappa = zeros(1, nf); nmse_db = zeros(1, nf);
perPort = zeros(size(Smeas, 1), nf);
for f = 1:nf
    p = Spred(:, :, f); y = Smeas(:, :, f);
    kappa(f) = (p(:)' * y(:)) / (p(:)' * p(:));
    e = kappa(f) * p - y;
    nmse_db(f) = 10*log10(sum(abs(e(:)).^2) / sum(abs(y(:)).^2));
    perPort(:, f) = 10*log10(sum(abs(e).^2, 2) ./ sum(abs(y).^2, 2));
end
end
