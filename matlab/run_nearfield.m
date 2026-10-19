function run_nearfield()
%RUN_NEARFIELD  Section VII: two-scan technique on two spheres (15 and 30 cm) vs. compressive phase
%retrieval on the 15 cm sphere (IC and IS, both polarizations), then the figures.
%Settings: cpr_settings.m (pattern grid, frequencies, SNR, maximum two-scan starts, trials).
cpr_setup();
if isempty(which('fminunc')) || ~license('test', 'Optimization_Toolbox')
    warning('run_nearfield:toolbox', ['No Optimization Toolbox: the two-scan benchmark uses the built-in BFGS ' ...
        '(twoscan_solve ''Optimizer'' ''own'') instead of fminunc.']);
end
S = cpr_settings(); snr = S.nearfield.snr;
twoscan_run('crlh', 0.15, 0.30, snr, 1:S.nd);
for seed = 1:S.trials.nearfield                    % both at the base ratio m/n 0.3 (user, 2026-10-19: IS 0.6 dropped)
    cpr_sphere_run(0.15, seed, 'npol', 2, 'snrdB', snr);                                    % IC, m/n 0.3
    cpr_sphere_run(0.15, seed, 'npol', 2, 'snrdB', snr, 'arch', 'IS', 'ratio', 0.3);
    cpr_sphere_run(0.15, seed, 'npol', 2, 'snrdB', snr, 'arch', 'RS', 'ratio', 0.3);   % all four (user, 2026-10-19)
    cpr_sphere_run(0.15, seed, 'npol', 2, 'snrdB', snr, 'arch', 'RC', 'ratio', 0.3);
end
cpr_fig_sec5();
end