function summary_fullwave()
%SUMMARY_FULLWAVE  Section VII numbers from the trials of run_fullwave (<results>/fwc_trials_true_*.mat): per-frequency NMSE of the refined
%estimates and of PR-GAMP alone, the refinement gain and rule, constant+delay and global-phase scores, accepted residuals and restarts for the
%ten HFSS trials, then the means of the other sets (m/n 0.6, perturbed calibration, moment-method ring).
cpr_setup(); maxNumCompThreads(1);
T = fwc_load_trials('^fwc_trials_true_hfss_ring\.mat$');
ref = arrayfun(@(t) t.nmse.perfreqScale, T); raw = arrayfun(@(t) t.nmse.rawPerfreqScale, T);
gd = arrayfun(@(t) t.nmse.globalDelay, T); gp = arrayfun(@(t) t.nmse.globalPhase, T); mo = arrayfun(@(t) t.nmse.model, T);
nres = arrayfun(@(t) t.nres_dB, T); tries = arrayfun(@(t) t.numTry, T); lin = @(x) 10*log10(mean(10.^(x/10)));
fprintf('seeds %s\n', mat2str([T.seed]));
fprintf('refined   per-frequency: %s | mean %.2f (linear %.2f), best %.2f, worst %.2f\n', mat2str(round(ref,1)), mean(ref), lin(ref), min(ref), max(ref));
fprintf('PR-GAMP   per-frequency: %s | mean %.2f (linear %.2f), best %.2f, worst %.2f\n', mat2str(round(raw,1)), mean(raw), lin(raw), min(raw), max(raw));
fprintf('gain refined-raw per seed: %s | mean %.2f dB, linear-mean gain %.2f dB, worst seed change %.2f dB (rule: mean >= 1 dB and no seed worse by more than 1 dB)\n', mat2str(round(raw-ref,2)), mean(raw-ref), lin(raw)-lin(ref), min(raw-ref));
fprintf('constant+delay: mean %.2f, best %.2f, worst %.2f | one global phase: mean %.2f | model error of the truth: mean %.2f\n', mean(gd), min(gd), max(gd), mean(gp), mean(mo));
fprintf('accepted residual of PR-GAMP (dB): %s | restarts: %s\n', mat2str(round(nres,1)), mat2str(tries));
sets = {'^fwc_trials_true_hfss_ring_r60\.mat$', 'HFSS ring, m/n 0.6'; '^fwc_trials_true_az0\.5_rel-Inf_hfss_ring_pa05\.mat$', 'HFSS ring, azimuth 0.5 deg'; ...
    '^fwc_trials_true_az1_rel-Inf_hfss_ring_pa1\.mat$', 'HFSS ring, azimuth 1 deg'; '^fwc_trials_true_az0_rel-40_hfss_ring_pr40\.mat$', 'HFSS ring, entries -40 dB'; ...
    '^fwc_trials_true_az0_rel-30_hfss_ring_pr30\.mat$', 'HFSS ring, entries -30 dB'; '^fwc_trials_true_hfss_arch\.mat$', 'HFSS ring, half-ring calibration'; '^fwc_trials_true_hfss_arch_r60\.mat$', 'HFSS ring, half-ring calibration, m/n 0.6'; '^fwc_trials_true_momring_ring_r30\.mat$', 'MoM ring, ring calibration, m/n 0.3'; ...
    '^fwc_trials_true_momring_arch_r30\.mat$', 'MoM ring, arch calibration, m/n 0.3'; '^fwc_trials_true_momring_ring_r60\.mat$', 'MoM ring, ring calibration, m/n 0.6'; ...
    '^fwc_trials_true_momring_arch_r60\.mat$', 'MoM ring, arch calibration, m/n 0.6'};
for i = 1:size(sets, 1)
    try, U = fwc_load_trials(sets{i, 1}); catch, fprintf('%-38s (not run)\n', sets{i, 2}); continue, end
    v = arrayfun(@(t) t.nmse.perfreqScale, U); r = arrayfun(@(t) t.nmse.rawPerfreqScale, U);
    fprintf('%-38s seeds %-12s refined mean %7.2f | PR-GAMP alone %7.2f | constant+delay %7.2f\n', sets{i, 2}, mat2str([U.seed]), mean(v), mean(r), mean(arrayfun(@(t) t.nmse.globalDelay, U)));
end
end
