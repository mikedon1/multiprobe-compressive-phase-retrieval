function Tr = fwc_build_ring(outFile)
%FWC_BUILD_RING  True-ring data of Section VII from the packed HFSS exports: the AUT-to-probe transmissions of the 65-port model (16 rotations)
%and the 64-port embedded-response calibration of the ring alone (from its near fields), with the model check of Section VII-A.
%   Tr = fwc_build_ring            reads data/hfss_fullring (fwc_pack_ring) and data/fullwave/fwc_validation.mat; writes
%                                  data/fullwave/fwc_ring_hfss.mat (double), which run_fullwave and fwc_ring_data use directly
% Tr.Xring, Tr.psi, Tr.freq as in fwc_ring_true; Tr.Rring 64 x nmode x nf (ring calibration); Tr.check lists the model error of the truth
% (one complex factor per frequency) for the ring calibration, the opposite rotation sense, free-field probes and arch + symmetry.
if nargin < 1 || isempty(outFile), outFile = fullfile(cpr_paths().data, 'fullwave', 'fwc_ring_hfss.mat'); end
S0 = load(fullfile(cpr_paths().data, 'fullwave', 'fwc_validation.mat')); V = S0.V;
Tr = fwc_ring_true(V);
[Rring, info] = fwc_ring_calibration(V, fullfile(cpr_paths().data, 'hfss_fullring'));
Tr.Rring = Rring; Tr.info = info;
modes = V.A.modes; sgn = reshape((-1) .^ modes(:, 2), 1, []); nf = size(Tr.Xring, 3); nr = size(Tr.Xring, 2); np2 = size(Tr.Xring, 1); nm = size(modes, 1);
Rfree = fwc_free_resp(V); Rfree = [Rfree; Rfree .* sgn]; Rarch = [V.R; V.R .* sgn];
cases = {'ring calibration', Rring, -1; 'opposite rotation sense', Rring, +1; 'free-field probes', Rfree, -1; 'arch + symmetry', Rarch, -1};
Tr.check = strings(size(cases, 1), 1);
for c = 1:size(cases, 1)
    R = cases{c, 2}; sg = cases{c, 3}; e = zeros(1, nf); en = zeros(1, nf);
    for f = 1:nf
        Gf = zeros(np2*nr, nm);
        for l = 1:nr, Gf((l-1)*np2 + (1:np2), :) = R(:, :, f) .* exp(sg*1i*modes(:, 2).' * Tr.psi(l)); end
        xm = Gf * V.A.a(:, f); xd = reshape(Tr.Xring(:, :, f), [], 1); k = (xm' * xd) / (xm' * xm); e(f) = norm(k*xm - xd)^2; en(f) = norm(xd)^2;
    end
    Tr.check(c) = sprintf('%-26s: per frequency %.1f to %.1f dB, overall %.2f dB', cases{c, 1}, min(10*log10(e./en)), max(10*log10(e./en)), 10*log10(sum(e)/sum(en)));
    disp(Tr.check(c));
end
save(outFile, 'Tr', '-v7');
end
