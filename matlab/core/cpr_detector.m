function [y, on] = cpr_detector(c, S, x, z, w, f)
%CPR_DETECTOR  Phaseless data with a detector of finite response (2026-10-25, study pr_detbw). Returns the amplitudes y that
%replace abs(z + w) when c.detBwGHz > 0 (lowpass) or c.detIntPs > 0 (boxcar); on = true when the response was applied.
%   [y, on] = cpr_detector(c, S, x, z, w, f)     S = make_readings output, x the pattern vector, z = S.Phi*x, w the noise
% An ideal square-law detector reports |sum_k u_k|^2 for the tone components u_k of a reading row. A linear detector response
% h with transfer H weights each beat note of tones k and k' (frequency difference D = f_k - f_k') by H(D):
%   q_r = sum_{k,k'} H(f_k - f_k') u_k conj(u_k')      (the sample-time phases are already inside u_k, make_readings)
% H(-D) = conj(H(D)), so q is real. Models: 'lowpass' H = 1/(1 + iD/B), B = detBwGHz (a first-order response, whose group delay
% the solver does not know); 'boxcar' H = sinc(D T), T = detIntPs (centered window; an uncentered one adds a known delay).
% The noise enters as for the ideal detector, q + (|z+w|^2 - |z|^2), unfiltered (it is small at the SNRs of the studies), and
% y = sqrt(max(q, 0)). The solver keeps the ideal model |A x| (a mismatch, as for timing jitter and cable errors).
on = false; y = abs(z + w);
bw = 0; if isfield(c, 'detBwGHz'), bw = c.detBwGHz; end
ti = 0; if isfield(c, 'detIntPs'), ti = c.detIntPs; end
mdl = 'lowpass'; if isfield(c, 'detModel'), mdl = c.detModel; end
if strcmp(mdl, 'lowpass') && bw <= 0, return, end
if strcmp(mdl, 'boxcar') && ti <= 0, return, end
if ~isfield(S.ent, 'vTrue'), error('cpr_detector: S.ent.vTrue missing (make_readings of 2026-10-25 or later)'); end
e = S.ent; nd = numel(f); f = f(:);
v = e.vTrue .* x(e.j);\nu = accumarray([e.i, e.tone], real(v), [S.m, nd]) + 1i * accumarray([e.i, e.tone], imag(v), [S.m, nd]);   % m x nd tone components
D = f - f.';                                                                % D(k, k') = f_k - f_k'
switch mdl
    case 'lowpass', H = 1 ./ (1 + 1i * D / (bw * 1e9));
    case 'boxcar',  H = sinc(D * ti * 1e-12);                               % sinc(x) = sin(pi x)/(pi x)
    otherwise, error('cpr_detector: unknown detModel %s', mdl);
end
q = real(sum((u * H) .* conj(u), 2));                                       % u H u^H per row
q = q + (abs(z + w).^2 - abs(z).^2);
y = sqrt(max(q, 0)); on = true;
end