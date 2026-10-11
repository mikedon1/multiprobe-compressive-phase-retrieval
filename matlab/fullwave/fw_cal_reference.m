function fw_cal_reference()
%FW_CAL_REFERENCE  Leg calibration factors from the reference-dipole HFSS solve (Section VI).
%The reference dipole lies along the rotation axis, so its E_theta does not depend on the
%rotation and the anchor ring probe (theta = 90 deg) and the three legs see the same field.
%Each leg's factor at each frequency is G = leg / anchor. No AUT data are used. Also reports the
%residual mismatch between each leg and its anchor entry in the AUT data, without and with G.
%Saves <results>/fw_cal_reference.mat with G (nlegs x nf).
outDir = cpr_paths().results;
[Xp, info] = load_hfss_mf(fullfile(cpr_paths().data, 'hfss', 'sweep_reference_rot00_psi000.00.s20p'));
nTh = 16; anchor = 16; nl = size(Xp, 1) - nTh;
G = complex(zeros(nl, numel(info.freq)));
for k = 1:nl
    G(k, :) = squeeze(Xp(nTh + k, 1, :)).' ./ squeeze(Xp(anchor, 1, :)).';
end
save(fullfile(outDir, 'fw_cal_reference.mat'), 'G');
fprintf('reference |G| by leg (rows) and frequency:\n'); disp(round(abs(G), 3));
fprintf('reference phase(G) (deg):\n'); disp(round(angle(G)*180/pi, 1));
L = load(fullfile(outDir, 'fw_data_hfss.mat')); D = L.D;
a = squeeze(D.Xring(D.anchor, :, :));
for k = 1:nl                                             % ||l_q - kappa_q x_a|| / ||l_q|| (Section VI)
    l = squeeze(D.Xleg(k, :, :)); ref = circshift(a, -D.legShift*k, 1);
    fprintf('leg %d: mismatch with the anchor entry %.1f%% without and %.2f%% with the calibration factor\n', k, ...
        100*norm(l - ref, 'fro')/norm(l, 'fro'), 100*norm(l - G(k, :).*ref, 'fro')/norm(l, 'fro'));
end
end
