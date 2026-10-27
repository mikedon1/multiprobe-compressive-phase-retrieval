function fwc_joint_delay()
%FWC_JOINT_DELAY  One-dimensional joint delay optimum (review item 5.1, Eq. R6) for the ten HFSS ring-calibration trials: the delay tau and one
%complex constant that minimize the NMSE of all coefficients together, against the score that takes the delay from the fitted factors.
cpr_setup(); maxNumCompThreads(1);
V = load(fullfile(cpr_paths().data, 'fullwave', 'fwc_validation.mat')); V = V.V; a = V.A.a; f = V.freq(:).';
T = fwc_load_trials('^fwc_trials_true_hfss_ring\.mat$'); dfr = f - f(1); e0 = zeros(1, numel(T)); e1 = e0; tau = e0; tauFit = e0;
sc = @(ah) 10*log10(norm((ah(:)'*a(:))/(ah(:)'*ah(:))*ah(:) - a(:))^2 / norm(a(:))^2);
for k = 1:numel(T)
    ah = T(k).aHat; kf = T(k).perfreqFactor(:).'; ph = unwrap(angle(kf)); p = polyfit(dfr, ph, 1);
    e0(k) = T(k).nmse.globalDelay; tauFit(k) = p(1)/(2*pi) * 1e9;
    obj = @(tn) sc(ah .* exp(1i*2*pi*dfr*tn*1e-9)); g = -20:0.05:20; v = arrayfun(obj, g); [~, i] = min(v);
    tn = fminsearch(obj, g(i)); tau(k) = tn; e1(k) = obj(tn);
end
fprintf('fitted-factor delay score (as in the paper): mean %.2f dB (%s)\n', mean(e0), mat2str(round(e0,2)));
fprintf('joint delay optimum:                          mean %.2f dB (%s), tau %s ns (factor-fit tau %s ns)\n', mean(e1), mat2str(round(e1,2)), mat2str(round(tau,2)), mat2str(round(tauFit,2)));
end
