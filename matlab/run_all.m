function run_all(varargin)
%RUN_ALL  Run every simulation behind the paper and its supplementary material, then draw the figures.
%
%   run_all                       everything, with the current settings preset
%   run_all('Skip', {'nearfield', 'fullwave', 'checks'})   leave parts out
%
% Parts: the Section IV sweeps (pr_*) and the supplementary complex-measurement sweeps (cs_*),
% 'nearfield' (Section V), 'fullwave' (Section VI), 'checks' (conjugate-ambiguity checks), then
% the figures. With the 'paper' preset this takes days on a laptop (the circular architectures
% and the two-scan technique dominate); run_quick runs the same things at reduced size.
% Everything is checkpointed, so an interrupted run resumes where it stopped.
p = inputParser; p.addParameter('Skip', {}); p.parse(varargin{:});
skip = p.Results.Skip;
cpr_setup();
for meas = {'pr', 'cs'}
    for kind = {'ratio', 'P', 'K', 'rn', 'nd', 'snr', 'legs', 'size', 'examples'}
        run_study([meas{1} '_' kind{1}]);
    end
end
run_study('pr_gauss');                % on/off switches vs Gaussian weights (Section IV-A)
run_study('pr_calerr');                % sensitivity to cable phase errors (Section IV-C)
cpr_robustness_report();
cpr_paper_figs('IV');                        % Section IV figures (phaseless)
cpr_paper_figs('II');                        % supplementary complex-measurement figures
cpr_cuts(1, 'pr'); cpr_cuts(1, 'pr', 'phi'); cpr_cuts(1, 'cs');   % pattern cuts at the example NMSE levels
if ~ismember('nearfield', skip), run_nearfield(); end
if ~ismember('fullwave', skip),  run_fullwave();  end
if ~ismember('checks', skip),    cpr_conj_check(); cpr_conj_check_k4(); end
fprintf('run_all done. Figures: %s\n', fullfile(cpr_paths().results, 'figs', 'paper'));
end