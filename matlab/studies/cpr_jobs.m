function J = cpr_jobs(study)
%CPR_JOBS  Job list of a study. Each job = cpr_cfg overrides + a random seed.
%   J = cpr_jobs('pr_ratio')
% Every job uses cpr_cfg() for everything it does not override, and the sweep grids and trial
% counts come from cpr_settings.m. Studies are named <meas>_<kind> with meas = cs (complex
% measurements, supplementary material) or pr (phaseless measurements, main paper):
%   <meas>_ratio     sampling ratio m/n, all four architectures
%   <meas>_P         active probes P (IS)
%   <meas>_K         tones per reading K (IS)
%   <meas>_rn        installed-probe fraction r_n (RS and RC)
%   <meas>_nd        number of wavelengths, recovered jointly or each alone (IS)
%   <meas>_snr       SNR, all four architectures (IC and RC use settings.trials.snrCircular trials)
%   <meas>_legs      leg probes (IC and RC)
%   <meas>_size      angular resolution nw (IS)
%   <meas>_examples  IS runs that keep the recovered pattern (pattern-cut figures)
%   pr_gauss         on/off switches vs complex Gaussian weights on the same probes (IS at several
%                    sampling ratios, IC at its base ratio): Section IV-A
%   pr_calerr        sensitivity to rms cable phase errors in the data (IS and IC): Section IV-C
%   smoke            two tiny jobs, to check that the code runs
S = cpr_settings();
J = struct('study', {}, 'label', {}, 'over', {}, 'seed', {}, 'cost', {}, 'saveXh', {});
sw = 1:S.trials.sweep;
archs = {'IS', 'RS', 'IC', 'RC'};
if strcmp(study, 'smoke')
    J = add(J, study, 'cs', {'meas','complex','nw',16,'nd',2,'K',2,'ratio',0.3}, 1);
    J = add(J, study, 'pr', {'meas','phaseless','nw',16,'nd',2,'K',2,'ratio',0.6,'minTry',1,'maxTry',1}, 1);
    return
end
tok = regexp(study, '^(cs|pr)_(\w+)$', 'tokens', 'once');
if isempty(tok), error('cpr_jobs: unknown study %s', study); end
isCs = strcmp(tok{1}, 'cs');
meas = {'meas', ifelse(isCs, 'complex', 'phaseless')};
switch tok{2}
    case 'ratio'
        for a = archs, for rt = S.sweeps.ratio
            J = add(J, study, sprintf('%s_r%.2f', a{1}, rt), [meas {'arch',a{1},'ratio',rt}], sw);
        end, end
    case 'P'
        Pl = S.sweeps.Ppr; if isCs, Pl = S.sweeps.Pcs; end
        for P = Pl, J = add(J, study, sprintf('P%d', P), [meas {'P',P}], sw); end
    case 'K'
        for K = S.sweeps.K, J = add(J, study, sprintf('K%d', K), [meas {'K',K}], sw); end
    case 'rn'
        for rn = S.sweeps.rnRS, J = add(J, study, sprintf('RS_rn%.4f', rn), [meas {'arch','RS','rnRS',rn}], sw); end
        for rn = S.sweeps.rnRC, J = add(J, study, sprintf('RC_rn%.4f', rn), [meas {'arch','RC','rnRC',rn}], sw); end
    case 'nd'
        % joint: the first nd wavelengths together; alone: each of them by itself (nd = 1, f0idx = j)
        ndl = S.sweeps.ndPr; if isCs, ndl = S.sweeps.ndCs; end
        for nd = ndl
            J = add(J, study, sprintf('joint_nd%d', nd), [meas {'nd',nd,'K',min(4,max(1,nd/2))}], sw);   % K < nd: K = nd is the worst case
        end
        for j = 1:max(ndl)
            J = add(J, study, sprintf('single_f%d', j), [meas {'nd',1,'K',1,'f0idx',j}], sw);
        end
    case 'snr'
        for a = archs, for sn = S.sweeps.snr
            seeds = sw; if ~isCs && any(strcmp(a{1}, {'IC', 'RC'})), seeds = 1:S.trials.snrCircular; end
            J = add(J, study, sprintf('%s_snr%d', a{1}, sn), [meas {'arch',a{1},'snrdB',sn}], seeds);
        end, end
    case 'legs'
        for a = {'IC', 'RC'}, for nl = S.sweeps.legs
            J = add(J, study, sprintf('%s_legs%d', a{1}, nl), [meas {'arch',a{1},'nlegs',nl}], sw);
        end, end
    case 'size'
        for nw = S.sweeps.size, J = add(J, study, sprintf('nw%d', nw), [meas {'nw',nw}], sw); end
    case 'gauss'                              % same seeds: probes, tones and times match between the two
        for w = {'switch', 'gauss'}
            for rt = S.sweeps.gaussRatio
                J = add(J, study, sprintf('IS_%s_r%.2f', w{1}, rt), [meas {'weights',w{1},'ratio',rt}], sw);
            end
            J = add(J, study, sprintf('IC_%s', w{1}), [meas {'arch','IC','weights',w{1}}], sw);
        end
    case 'calerr'
        J = add(J, study, 'IS_cal0', meas, sw);
        for e = S.sweeps.calErr.IS, J = add(J, study, sprintf('IS_cal%g', e), [meas {'calErrDeg',e}], sw); end
        J = add(J, study, 'IC_cal0', [meas {'arch','IC'}], sw);
        for e = S.sweeps.calErr.IC, J = add(J, study, sprintf('IC_cal%g', e), [meas {'arch','IC','calErrDeg',e}], sw); end
    case 'examples'
        rl = S.sweeps.examples.pr; if isCs, rl = S.sweeps.examples.cs; end
        for rt = rl
            J = add(J, study, sprintf('IS_ex_r%.2f', rt), [meas {'ratio',rt}], sw);
            J(end).saveXh = true;
        end
    otherwise
        error('cpr_jobs: unknown study %s', study);
end
end

function J = add(J, study, label, over, seeds)
c = cpr_cfg(over{:});
cost = c.nw^2/2 * c.nd * max(c.ratio, 0.05) * c.P * c.K;            % rough relative cost
if strcmp(c.meas, 'phaseless'), cost = cost * c.maxTry * 20; end
for s = seeds
    J(end+1) = struct('study', study, 'label', label, 'over', {over}, 'seed', s, ...
        'cost', cost, 'saveXh', false); %#ok<AGROW>
end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
