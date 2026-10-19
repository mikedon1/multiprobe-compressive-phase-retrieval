function D = fwc_compare_mesh(fileA, fileB)
%FWC_COMPARE_MESH  How much do the AUT->probe transmissions change between two solves of the same
%geometry (e.g. loose vs tight convergence)?
%   D = fwc_compare_mesh(fileA, fileB)
% For S(p, 1), p = 2..N, per frequency: magnitude change (dB), phase change (deg), and the relative
% complex change 20 log10 |SB - SA| / |SA| (dB), raw and after removing one complex constant per frequency
% common to all probes (what a per-frequency calibration constant absorbs). Co-pol ports 2..17 and
% cross-pol ports 18..33 are reported separately (dual ring).
[SA, f] = read_touchstone(fileA); SB = read_touchstone(fileB);
a = squeeze(SA(2:end, 1, :)); b = squeeze(SB(2:end, 1, :));
D.freq = f(:).';
D.dmag = 20*log10(abs(b) ./ abs(a)); D.dph = angle(b ./ a) * 180/pi;
D.rel = 20*log10(abs(b - a) ./ abs(a));
k = sum(conj(a) .* b, 1) ./ sum(abs(a).^2, 1);            % one constant per frequency
D.relConst = 20*log10(abs(b - k .* a) ./ abs(a));
groups = {'co (ports 2-17)', 1:16; 'cross (ports 18-33)', 17:32};
for g = 1:2
    r = groups{g, 2};
    fprintf('%-20s |dmag| max %.2f dB, rms %.2f dB; |dphase| max %.2f deg; rel change median %.1f dB (worst %.1f); after per-freq constant median %.1f dB (worst %.1f)\n', ...
        groups{g, 1}, max(abs(D.dmag(r, :)), [], 'all'), rms(D.dmag(r, :), 'all'), max(abs(D.dph(r, :)), [], 'all'), ...
        median(D.rel(r, :), 'all'), max(D.rel(r, :), [], 'all'), median(D.relConst(r, :), 'all'), max(D.relConst(r, :), [], 'all'));
end
fprintf('co ports, per frequency, rel change after constant (median over ports): %s dB\n', mat2str(round(median(D.relConst(1:16, :), 1)*10)/10));
end
