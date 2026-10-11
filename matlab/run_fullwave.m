function run_fullwave()
%RUN_FULLWAVE  Section VI: HFSS data checks, reference-dipole leg calibration, ten calibrated trials
%(and the same trials without calibration), then the cut figures. Settings: cpr_settings.m.
%Data: data/hfss (32 rotations of the 20-port HFSS model and the reference-dipole solve).
cpr_setup();
S = cpr_settings();
cpr_fw_diag();                                   % reads the Touchstone files, writes fw_data_hfss.mat
fw_cal_reference();                              % leg calibration factors G (reference-dipole solve)
cal = fullfile(cpr_paths().results, 'fw_cal_reference.mat');
cpr_fw_trials(1:S.trials.fullwave, 'calop', 'CalFile', cal);
cpr_fw_trials(1:S.trials.fullwave, 'uncal');
Rt = load(fullfile(cpr_paths().results, 'fw_trials_calop.mat')); R = Rt.R;
sc = [R.score]; [~, i] = max([sc.global]);               % the worst calibrated trial is the one plotted
Dl = load(fullfile(cpr_paths().results, 'fw_data_hfss.mat'));
nf = size(R(i).Xh, 3);
cpr_fw_cuts('calop', R(i).seed, mean(Dl.D.f(1:nf))/1e9);
end