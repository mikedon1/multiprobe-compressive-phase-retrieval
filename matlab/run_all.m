function run_all(varargin)
%RUN_ALL  Run every simulation behind the paper and its supplementary material, then draw the figures.
%
%   run_all                              everything with the 'paper' settings (takes days on a laptop)
%   run_all('Preset', 'quick')           everything at reduced size (results go to results_quick/, and the
%                                        numbers do NOT match the paper)
%   run_all('Skip', {'nearfield', 'fullwave', 'checks'})   leave parts out
%
% All settings are at the top of cpr_settings.m; the base values of the simulations are in core/cpr_cfg.m.
% Everything is checkpointed, so an interrupted run resumes where it stopped. A single study can be run with
% run_study('pr_ratio') and read back with cpr_load('pr_ratio'). MANIFEST.md (repository root) maps every
% figure and number of the paper to its study, script and results file.
p = inputParser; p.addParameter('Skip', {}); p.addParameter('Preset', '');
p.parse(varargin{:});
skip = p.Results.Skip;
cpr_setup();
if ~isempty(p.Results.Preset), cpr_settings(p.Results.Preset); end

% ---- Section VI and the supplement: phaseless studies (definitions in studies/cpr_jobs.m) ----
studies = { ...
    'pr_nd', 'pr_K', 'pr_P', ...                          % Fig. 4: wavelengths, tones per reading, active ports
    'pr_examples', ...                                    % Fig. 5 and the azimuth cut: IS examples with the pattern kept
    'pr_ratio', 'pr_snr', 'pr_rn', ...                    % Fig. 6(a), 6(b), supplement r_n figure
    'pr_calerr', 'pr_calerr30', 'pr_calerrsc', 'pr_calerrics', 'pr_calerrrc', ...   % Fig. 6(c)
    'pr_basiscmp32', ...                                  % Fig. 7: spherical-wave vs. Fourier basis
    'pr_icbad', ...                                       % Fig. 3: IC with the Fourier basis, phase per cut
    'pr_gauss', 'pr_nosparse', 'pr_combiner', ...         % Section VI-A/B text: Gaussian weights, no sparsity, combiner
    'pr_jitter', 'pr_tonespan', 'pr_trunc', 'pr_aut2', ...% Section VI-C and supplement tables
    'pr_cables', 'pr_scprior', 'pr_repeat', 'pr_caldrift', 'pr_calgain', 'pr_iters', 'pr_pairs', ...   % supplement
    'pr_detbw', 'pr_smooth'};                              % supplement: detector response, structured (smooth, dispersive) calibration errors
for s = studies, run_study(s{1}); end

% ---- Figures (figures/); written to <results>/figs/paper when S.figures.print is true ----
cpr_paper_figs('IV');                 % pr2_nd, pr2_K, pr2_P, pr2_ratio, pr2_snr, pr2_rn (paper files pr3_*.pdf)
cpr_paper_figs('v9');                 % pr3_basiscmp (Fig. 7), pr3_calerr (Fig. 6(c))
cpr_cuts_vsh();                       % pr3_cut (Fig. 5)
cpr_cuts_vsh('phi');                  % pr3_cut_phi (supplement)
cpr_ic_bad_vsh();                     % IC_bad (Fig. 3)

% ---- Section VIII: comparison with the two-scan technique (nearfield/) ----
% two-scan technique on two spheres with frequency continuation, the four architectures on one sphere;
% writes ift_freq (Fig. 10) and ift_accel (supplement)
if ~ismember('nearfield', skip), run_nearfield(); end

% ---- Section VII: full-wave example with embedded-probe calibration (fullwave/) ----
% true ring (65-port HFSS model): ten trials, m/n 0.6, perturbed calibration; moment-method ring; cut figures (cutphi, cuttheta); summary
if ~ismember('fullwave', skip),  run_fullwave();  end

% ---- Conjugate-ambiguity checks (Section IV-C numbers; studies/) ----
if ~ismember('checks', skip),    cpr_conj_check(); cpr_conj_check_k4(); end
fprintf('run_all done. Figures: %s\n', fullfile(cpr_paths().results, 'figs', 'paper'));
end
