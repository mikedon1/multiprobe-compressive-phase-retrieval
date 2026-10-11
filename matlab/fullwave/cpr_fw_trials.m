function cpr_fw_trials(seeds, tag, varargin)
%CPR_FW_TRIALS  Section VI trials on the 32-rotation multi-frequency HFSS data.
%
%   cpr_fw_trials(1:10, 'calop', 'CalFile', fullfile(cpr_paths().results, 'fw_cal_reference.mat'))
%   cpr_fw_trials(1:10, 'uncal')                         % without the leg calibration
%
% Phaseless IC base (cpr_settings.m), K = 4, P = 4, m/n = 0.6, 3 leg probes, 50 dB SNR, using
% fw_data_hfss.mat from cpr_fw_diag. Saves <results>/fw_trials_<tag>.mat after every trial
% (resumable: finished seeds are skipped). Options:
%   'LegMode'  'real' (default: every leg is its own solved HFSS port) | 'ideal' (legs reuse
%              the anchor ring probe at the rotation offset, the earlier idealization)
%   'CalFile'  mat file with G (nlegs x nf): each leg's complex response relative to the anchor
%              ring probe at each frequency (cpr_fw_cal_reference). Empty: no calibration.
%   'CalMode'  'operator' (default): G enters the recovery operator's leg entries while the
%              detected data keep the actual leg responses, as a switch-only combiner gives
%              them. 'divide': leg data are divided by G before combination, which needs a
%              branch-compensation network that the system does not have.
p = inputParser;
p.addParameter('LegMode', 'real'); p.addParameter('CalFile', ''); p.addParameter('CalMode', 'operator');
p.parse(varargin{:}); o = p.Results;
S = cpr_settings();
maxNumCompThreads(S.threads);
outDir = cpr_paths().results;
L = load(fullfile(outDir, 'fw_data_hfss.mat')); D = L.D;
nf = S.fullwave.nd;                                    % frequencies of the sweep to use
if ~isempty(o.CalFile)
    C = load(o.CalFile);
    if strcmp(o.CalMode, 'divide')
        for k = 1:size(C.G, 1)
            D.Xleg(k, :, :) = D.Xleg(k, :, :) ./ reshape(C.G(k, :), 1, 1, []);
        end
    else
        D.legCal = C.G;
    end
    fprintf('legs calibrated with %s (%s)\n', o.CalFile, o.CalMode);
end
if nf < size(D.Xring, 3)                               % quick preset: first nf frequencies
    D.Xring = D.Xring(:, :, 1:nf); D.Xleg = D.Xleg(:, :, 1:nf); D.f = D.f(1:nf);
    if isfield(D, 'legCal'), D.legCal = D.legCal(:, 1:nf); end
end
c = cpr_cfg('meas', 'phaseless', 'arch', 'IC', 'K', 4, 'P', 4, 'ratio', 0.6, 'nlegs', 3, 'snrdB', 50);
f = fullfile(outDir, sprintf('fw_trials_%s.mat', tag));
if exist(f, 'file'), R0 = load(f); R = R0.R; else, R = struct([]); end
for s = seeds
    if ~isempty(R) && any([R.seed] == s), continue; end
    r = cpr_run_fw(c, s, D, o.LegMode);
    e = struct('seed', s, 'legMode', o.LegMode, 'score', r.score, 'time', r.time, 'm', r.m, 'n', r.n, ...
        'R', r.R, 'legMismatch_dB', r.legMismatch_dB, 'numTry', r.info.numTry, 'nres_dB', r.info.nres_dB, 'Xh', r.Xh, 'cfg', c);
    if isempty(R), R = e; else, R(end+1) = e; end %#ok<AGROW>
    save(f, 'R');
    fprintf('fw %s seed %d: global %.2f | per-lambda %.2f dB, legMismatch %.1f dB, %d restarts, %.0f s, m %d\n', ...
        o.LegMode, s, r.score.global, r.score.perlambda, r.legMismatch_dB, r.info.numTry, r.time, r.m);
end
fprintf('fw trials %s done\n', tag);
end
