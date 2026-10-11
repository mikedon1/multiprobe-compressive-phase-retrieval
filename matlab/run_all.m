function run_all(varargin)
%RUN_ALL  Run every simulation behind the paper and its supplementary material, then draw the figures.
%
%   run_all                              everything with the 'paper' settings (takes days on a laptop)
%   run_all('Preset', 'quick')           everything at reduced size (about 1.5 hours; results go to
%                                        results_quick/, and the numbers do NOT match the paper)
%   run_all('Skip', {'nearfield', 'fullwave', 'checks'})   leave parts out
%
% All settings (pattern size, restarts, trials, sweep grids, whether figures are also written as
% PDF and PNG) are at the top of cpr_settings.m. Everything is checkpointed, so an interrupted run
% resumes where it stopped. A single study can be run with run_study('pr_ratio') and read back with
% cpr_load('pr_ratio').
p = inputParser; p.addParameter('Skip', {}); p.addParameter('Preset', '');
p.parse(varargin{:});
skip = p.Results.Skip;
cpr_setup();
if ~isempty(p.Results.Preset), cpr_settings(p.Results.Preset); end

% ---- Simulation studies (core/ and studies/): one study = one sweep, definitions in cpr_jobs.m ----
% pr_* are phaseless (Section IV and supplement S4-S5), cs_* complex measurements (supplement S1-S3).
%   *_ratio     sampling ratio m/n, all four architectures
%   *_P, *_K    active probes, tones per reading
%   *_rn        installed-probe fraction of RS and RC
%   *_nd        wavelengths, joint vs. one at a time
%   *_snr       SNR, all four architectures
%   *_legs      leg probes of IC and RC
%   *_size      angular resolution
%   *_examples  IS runs that keep the recovered pattern, used by the pattern-cut figures
for meas = {'pr', 'cs'}
    for kind = {'ratio', 'P', 'K', 'rn', 'nd', 'snr', 'legs', 'size', 'examples'}
        run_study([meas{1} '_' kind{1}]);
    end
end
run_study('pr_gauss');                % on/off switches vs Gaussian weights (Section IV-A)
run_study('pr_calerr');               % sensitivity to cable phase errors (Section IV-C)
cpr_robustness_report();              % prints both tables

% ---- Figures (figures/); saved to <results>/figs/paper when S.figures.print is true ----
cpr_paper_figs('IV');                 % Section IV: pr2_nd, pr2_K, pr2_ratio_ISRS, pr2_ratio, pr2_legs, pr2_snr
                                      %   and supplement S4-S5: pr2_P, pr2_size
cpr_paper_figs('supplement');         % supplement S1-S3: cs2_ratio, cs2_P, cs2_K, cs2_rn, cs2_nd, cs2_snr, cs2_size, ...
cpr_cuts(1, 'pr');                    % Fig. 4 pattern cut, elevation (pr2_cut)
cpr_cuts(1, 'pr', 'phi');             % supplement pattern cut, azimuth (pr2_cut_phi)
cpr_cuts(1, 'cs');                    % supplement pattern cut, complex measurements (cs2_cut)
cpr_ic_bad();                         % IC without leg probes: phase recovered per cut only (IC_bad)

% ---- Section V: comparison with the two-scan technique (nearfield/) -------------------------------
% run_nearfield: spherical-wave expansion, near-field data, two-scan technique and the compressive
% runs on the same data; writes ift_freq (Fig. 8) and ift_accel (supplement S8, cpr_fig_sec5).
if ~ismember('nearfield', skip), run_nearfield(); end

% ---- Section VI: full-wave example (fullwave/) -----------------------------------------------------
% run_fullwave: reads the HFSS Touchstone files in data/hfss, calibrates the leg probes with the
% reference-dipole solve, runs the calibrated and uncalibrated trials, then draws cutphi (Fig. 10)
% and cuttheta (supplement S7).
if ~ismember('fullwave', skip),  run_fullwave();  end

% ---- Conjugate-ambiguity checks (Section IV-B numbers; studies/) ----------------------------------
if ~ismember('checks', skip),    cpr_conj_check(); cpr_conj_check_k4(); end
fprintf('run_all done. Figures: %s\n', fullfile(cpr_paths().results, 'figs', 'paper'));
end
