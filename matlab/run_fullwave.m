function run_fullwave()
%RUN_FULLWAVE  Section VII: ten trials of the emulated full ring with embedded-probe calibration (HFSS data), the four
%perturbed-calibration runs (three trials each), then the cut figures. Settings: cpr_settings.m.
%Data: data/fullwave/fwc_validation.mat (processed HFSS data: free-space AUT coefficients, embedded probe responses,
%solved ring data). fwc_validate_hfss rebuilds it from the raw HFSS exports (data/hfss_fullarc or CPR_HFSS).
%The earlier leg-probe model (cpr_fw_trials, data/hfss) is kept in matlab/fullwave but is not used by the paper.
cpr_setup();
S = cpr_settings();
fwc_run_trials(1:S.trials.fullwave);                         % -> <results>/fwc_trials.mat (main and constant+delay scores)
fwc_run_trials(1:3, 'PerturbAzDeg', 0.5);
fwc_run_trials(1:3, 'PerturbAzDeg', 1);
fwc_run_trials(1:3, 'PerturbRelDb', -40);
fwc_run_trials(1:3, 'PerturbRelDb', -30);
cpr_fwc_cuts();                                              % worst trial, azimuth and elevation cuts at 2.725 GHz
end
