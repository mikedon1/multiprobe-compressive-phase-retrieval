function run_fullwave()
%RUN_FULLWAVE  Section VII: the TRUE full ring (65-port HFSS model, 16 rotations over 180 deg) with the ring-alone embedded-probe calibration
%and the refinement of the estimates (fwc_polish), then the moment-method ring as the accurate-data check, then the cut figures.
%   - ten trials at m/n 0.3 (S.trials.fullwave), three at m/n 0.6, and four perturbed-calibration runs of three trials each
%     (azimuth 0.5 and 1 deg, entries -40 and -30 dB); results fwc_trials_true_hfss_ring*.mat in <results>
%   - the half-ring calibration (arch response + 180-deg symmetry) on the same HFSS data, ten trials at m/n 0.3 and three at 0.6 (supplement)
%   - the same on the moment-method ring (fwc_trials_true_momring_*): ring calibration and arch + symmetry calibration, m/n 0.3 and 0.6
%   - cpr_fwc_cuts: worst trial, azimuth and elevation cuts at 2.725 GHz
%Data (data/fullwave, built from data/hfss_fullring and data/hfss_fullarc; see README): fwc_ring_data('hfss'), fwc_ring_data('mom').
%Number of trials and the extra sets: cpr_settings (paper 10 and all sets; quick 1 trial of the HFSS ring and of the MoM ring, no extras). summary_fullwave prints the numbers of the paper from the saved results.
cpr_setup();
S = cpr_settings(); n = S.trials.fullwave; n3 = min(3, n);
Tr = fwc_ring_data('hfss');
base = {'Ring', 'true', 'RingData', Tr, 'RingCal', 'ring', 'Polish', true};
fwc_run_trials(1:n, base{:}, 'OutTag', '_hfss_ring');
if S.fullwaveExtras
    fwc_run_trials(1:n3, base{:}, 'ratio', 0.6, 'OutTag', '_hfss_ring_r60');
    fwc_run_trials(1:n, 'Ring', 'true', 'RingData', Tr, 'RingCal', 'arch', 'Polish', true, 'OutTag', '_hfss_arch');        % half-ring calibration (supplement)
    fwc_run_trials(1:n3, 'Ring', 'true', 'RingData', Tr, 'RingCal', 'arch', 'Polish', true, 'ratio', 0.6, 'OutTag', '_hfss_arch_r60');
    fwc_run_trials(1:n3, base{:}, 'PerturbAzDeg', 0.5, 'OutTag', '_hfss_ring_pa05');
    fwc_run_trials(1:n3, base{:}, 'PerturbAzDeg', 1, 'OutTag', '_hfss_ring_pa1');
    fwc_run_trials(1:n3, base{:}, 'PerturbRelDb', -40, 'OutTag', '_hfss_ring_pr40');
    fwc_run_trials(1:n3, base{:}, 'PerturbRelDb', -30, 'OutTag', '_hfss_ring_pr30');
end
Vm = fwc_ring_data('momval'); Tm = fwc_ring_data('mom');
cals = {'ring'}; if S.fullwaveExtras, cals = {'ring', 'arch'}; end
for cal = cals
    mom = {'Validation', Vm, 'Ring', 'true', 'RingData', Tm, 'RingCal', cal{1}, 'Polish', true};
    fwc_run_trials(1:n, mom{:}, 'OutTag', ['_momring_' cal{1} '_r30']);
    if S.fullwaveExtras, fwc_run_trials(1:n3, mom{:}, 'ratio', 0.6, 'OutTag', ['_momring_' cal{1} '_r60']); end
end
cpr_fwc_cuts();
summary_fullwave();
end
