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
%   <meas>_L         truncation degree of the VSH basis (IS)
%   pr_legsdft       IC and RC with the earlier 2D DFT basis, without and with 8 leg probes
%   <meas>_examples  IS runs that keep the recovered pattern (pattern-cut figures)
%   pr_gauss         on/off switches vs complex Gaussian weights on the same probes (IS at several
%                    sampling ratios, IC at its base ratio): Section V-A
%   pr_calerr        sensitivity to rms cable phase errors in the data (IS and IC): Section V-C
%   pr_calerrsc      pr_calerr with self-calibration (30 EM updates); pr_calerr30 = its control without
%   pr_jitter        detector timing jitter (IS and IC); pr_nosparse  band-limit-only ablation (no sparsity prior)
%   pr_detbw         (NOT queued, 2026-10-25) finite detector bandwidth / integration time in the data (IS, IC; tones from the
%                    whole band or an 8-bin window)
%   pr_smooth        (NOT queued, 2026-10-25) structured calibration errors: lengths and gains smooth along the arch, cable dispersion
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
% restart cap for studies whose runs cannot reach the residual target (cable errors, timing jitter, no sparsity
% prior): every restart would run to the end, so they use 3 instead of 10 (user, 2026-10-19)
cap = {'maxTry', 3};
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
    case 'basiscmp32'
        % basis comparison FIGURE at the paper's settings (user, 2026-10-19): 32 x 16, both polarizations, exact
        % truth for both bases (cpr_run routes DFT with npol = 2 through cpr_run_pattern), the paper's solver
        % setting, plain prior for both. IS Fourier / IS spherical-wave and IC Fourier with 8 legs / IC
        % spherical-wave without legs over the ratio grid; IC Fourier without legs at 0.2, 0.6, 1.0.
        % 2026-10-19 (user): legs dropped from the paper -> no IC Fourier + 8 legs curve; IC Fourier (no legs) fails by
        % design, so it gets 4 points and the restart cap of 3
        com = {'npol',2,'nw',32,'groupPrior','none'};
        for rt = S.sweeps.ratio
            J = add(J, study, sprintf('IS_dft_r%.2f', rt), [meas com {'arch','IS','angBasis','dft','ratio',rt}], sw);
            J = add(J, study, sprintf('IS_vsh_r%.2f', rt), [meas com {'arch','IS','angBasis','vsh','ratio',rt}], sw);
            J = add(J, study, sprintf('IC_vsh_legs0_r%.2f', rt), [meas com {'arch','IC','angBasis','vsh','nlegs',0,'ratio',rt}], sw);
        end
        for rt = [0.2 0.3 0.6 1.0]
            J = add(J, study, sprintf('IC_dft_legs0_r%.2f', rt), [meas com {'arch','IC','angBasis','dft','nlegs',0,'ratio',rt}], sw);   % round 7: equal effort (10 restarts) with IC VSH
        end
    case 'basiscmp'
        % basis comparison figure (user, 2026-10-19): the setup of the earlier work (one component, 64 x 32),
        % identical readings and solver effort for every curve (10 restarts, 20 updates of 1000 iterations,
        % target -50 dB), plain prior for both bases. Curves: IS Fourier, IS spherical-wave, IC Fourier with 8
        % legs, IC spherical-wave without legs, all over the ratio grid; IC Fourier without legs at 3 ratios
        % (it fails at all of them, and each failure runs every restart).
        com = {'npol',1,'nw',64,'groupPrior','none','maxTry',10,'nitEM',20,'nitGAMP',1000,'solverSNRdB',48};
        for rt = S.sweeps.ratio
            J = add(J, study, sprintf('IS_dft_r%.2f', rt), [meas com {'arch','IS','angBasis','dft','ratio',rt}], sw);
            J = add(J, study, sprintf('IS_vsh_r%.2f', rt), [meas com {'arch','IS','angBasis','vsh','ratio',rt}], sw);
            J = add(J, study, sprintf('IC_dft_legs8_r%.2f', rt), [meas com {'arch','IC','angBasis','dft','nlegs',8,'ratio',rt}], sw);
            J = add(J, study, sprintf('IC_vsh_legs0_r%.2f', rt), [meas com {'arch','IC','angBasis','vsh','nlegs',0,'ratio',rt}], sw);
        end
        for rt = [0.2 0.6 1.0]
            J = add(J, study, sprintf('IC_dft_legs0_r%.2f', rt), [meas com {'arch','IC','angBasis','dft','nlegs',0,'ratio',rt}], sw);
        end
    case 'legsmatch32'
        % basis comparison at the paper's resolution (32 x 16), one component (the DFT path is single-
        % polarization), the paper's solver setting and the plain prior for BOTH bases; IC and RC with 0 and
        % 8 legs at m/n 0.6, IS at 0.3 as a reference
        for b = {'dft', 'vsh'}
            com = {'angBasis', b{1}, 'npol', 1, 'nw', 32, 'groupPrior', 'none'};
            for a = {'IC', 'RC'}, for nl = [0 8]
                J = add(J, study, sprintf('%s_%s_legs%d', a{1}, b{1}, nl), [meas com {'arch',a{1},'nlegs',nl,'ratio',0.6}], sw);
            end, end
            J = add(J, study, sprintf('IS_%s', b{1}), [meas com {'arch','IS','ratio',0.3}], sw);
        end
    case 'legsdft'
        % the 2D DFT basis of earlier versions (one component, 64 x 32 grid, the heavier circular
        % solver settings): IC and RC without and with leg probes, at the earlier base ratio 0.6
        dft = {'angBasis','dft','npol',1,'nw',64,'maxTry',30,'nitEM',20,'nitGAMP',1000,'solverSNRdB',48,'ratio',0.6};
        for a = {'IC', 'RC'}, for nl = [0 8]
            J = add(J, study, sprintf('%s_dft_legs%d', a{1}, nl), [meas dft {'arch',a{1},'nlegs',nl}], sw);
        end, end
    case 'iters'
        % iteration budget (user, 2026-10-19): runs end on the iteration limit, not the residual target. IS and IC
        % at the base values with 4x the iterations: 4x the EM updates, or 4x the GAMP iterations per update
        for a = {'IS', 'IC'}
            J = add(J, study, sprintf('%s_base', a{1}), [meas {'arch',a{1}}], sw);
            J = add(J, study, sprintf('%s_em24', a{1}), [meas {'arch',a{1},'nitEM',24}], sw);
            J = add(J, study, sprintf('%s_gamp800', a{1}), [meas {'arch',a{1},'nitGAMP',800}], sw);
        end
    case 'calerrsc'
        % self-calibration (user, 2026-10-19): the pr_calerr points recovered with cpr_selfcal (one cable length
        % correction per port, 30 EM updates, MAP prior = the expected rms error), IS and IC
        for a = {'IS', 'IC'}
            for e = S.sweeps.calErr.(a{1})
                J = add(J, study, sprintf('%s_cal%g_sc', a{1}, e), [meas cap {'arch',a{1},'calErrDeg',e,'selfCal',true, ...
                    'nitEM',30,'selfCalPriorDeg',e}], sw);
            end
        end
    case 'trunc'
        % robustness (user, 2026-10-19): truth NOT in the span of the recovery basis (the inverse-crime check).
        % CRLH truth (degree-16 fit) recovered with L = 8 and 10 at m/n 0.3 and 0.6; the 5 x 5 array truth (degree-30
        % fit, energy above degree 12 at -37 dB) recovered on 64 x 32 with L = 12 (truth -45 dB outside the span) and
        % L = 16 (-83 dB). IS and IC.
        for a = {'IS', 'IC'}
            for L = [8 10], for rt = [0.3 0.6]
                J = add(J, study, sprintf('%s_crlh_L%d_r%.1f', a{1}, L, rt), [meas {'arch',a{1},'vshL',L,'ratio',rt}], sw);
            end, end
            for L = [12 16]
                J = add(J, study, sprintf('%s_array5_L%d', a{1}, L), [meas cap {'arch',a{1},'pattern','array5', ...
                    'patternLfit',30,'vshL',L,'nw',64}], sw);
            end
        end
    case 'pairs'
        % round 7 (user, 2026-10-20): random port sets vs deterministic pairs (Paulus-style combinations of probe signals,
        % cycling neighbour pairs of each polarization and the two ports of one probe), P = 2, same everything else
        for a = {'IS', 'IC'}, for rt = [0.12 0.3], for pr = {'random', 'pairs'}
            J = add(J, study, sprintf('%s_%s_r%.2f', a{1}, pr{1}, rt), [meas {'arch',a{1},'P',2,'ratio',rt,'portRule',pr{1}}], sw);
        end, end, end
    case 'combiner'
        % round 7 (H1): ideal passive P-way combiner (1/sqrt(P) branch transfers) vs the unit-transfer idealization, IS
        % at 50 dB SNR (fixed detector noise) over the active-port count
        for cn = {'unit', 'passive'}, for P = [1 2 4 8 16]
            J = add(J, study, sprintf('IS_%s_P%d', cn{1}, P), [meas {'combNorm',cn{1},'P',P,'snrdB',50}], sw);
        end, end
    case 'scprior'
        % round 7 (C1/C2): self-calibration with a mismatched prior width, IC and two-pass IC at 5 deg cable error
        for a = {{'arch','IC'}, {'arch','IC','icShiftFrac',1/2}}, for pr = [2.5 10]
            lab = 'IC'; if numel(a{1}) > 2, lab = 'ICS'; end
            J = add(J, study, sprintf('%s_cal5_sc_prior%g', lab, pr), [meas cap a{1} {'calErrDeg',5,'selfCal',true, ...
                'nitEM',30,'selfCalPriorDeg',pr}], sw);
        end, end
    case 'repeat'
        % round 7 (critique section 8): repeated trials at a hard point. IS and IC at m/n 0.12 (the transition) and IC
        % at 5 deg cable error with self-calibration and its 30-update control, seeds 1-5
        rs = 1:5;
        for a = {'IS', 'IC'}
            J = add(J, study, sprintf('%s_r0.12', a{1}), [meas {'arch',a{1},'ratio',0.12}], rs);
        end
        J = add(J, study, 'IC_cal5_em30', [meas cap {'arch','IC','calErrDeg',5,'nitEM',30}], rs);
        J = add(J, study, 'IC_cal5_sc', [meas cap {'arch','IC','calErrDeg',5,'selfCal',true,'nitEM',30,'selfCalPriorDeg',5}], rs);
    case 'caldrift'
        % NOT queued (ask): thermal drift (every cable longer by the same fraction) in the data, IC, 30 EM updates:
        % no self-calibration, per-port self-calibration, and the one-parameter drift self-calibration
        for e = S.sweeps.calErr.IC
            dr = {'arch','IC','calErrDeg',e,'calErrModel','drift','nitEM',30};
            J = add(J, study, sprintf('IC_drift%g_em30', e), [meas cap dr], sw);
            J = add(J, study, sprintf('IC_drift%g_sclen', e), [meas cap dr {'selfCal',true,'selfCalPriorDeg',e}], sw);
            J = add(J, study, sprintf('IC_drift%g_scdrift', e), [meas cap dr {'selfCal',true,'selfCalModel','drift'}], sw);
        end
    case 'calgain'
        % NOT queued (ask): gain errors per port (0.1-1 dB rms) plus 5 deg length errors, two-pass IC, 30 EM updates:
        % self-calibration of lengths only vs lengths and gains
        for gdb = [0.1 0.2 0.5 1]
            cg = {'arch','IC','icShiftFrac',1/2,'calErrDeg',5,'calGainDb',gdb,'nitEM',30,'selfCal',true,'selfCalPriorDeg',5};
            J = add(J, study, sprintf('ICS_gain%g_sclen', gdb), [meas cap cg], sw);
            J = add(J, study, sprintf('ICS_gain%g_scgain', gdb), [meas cap cg {'selfCalGain',true,'selfCalPriorGainDb',gdb}], sw);
        end
    case 'aut2'
        % NOT queued (ask the user): the electrically larger second example antenna (array_data). 4 x 4 array
        % (k r0 = 7.3) on the paper grid with L = 16, IS and IC over the ratio grid; 5 x 5 array (k r0 = 9.7) with
        % L = 20 on the 64 x 32 grid (the margin holds) and on the 32 x 16 grid (2L+1 > 32, the margin vanishes)
        a4 = {'pattern','array4','patternLfit',30,'vshL',16};
        for a = {'IS', 'IC'}, for rt = S.sweeps.ratio
            J = add(J, study, sprintf('%s_array4_r%.2f', a{1}, rt), [meas a4 {'arch',a{1},'ratio',rt}], sw);
        end, end
        a5 = {'pattern','array5','patternLfit',30,'vshL',20};
        for a = {'IS', 'IC'}, for g = [32 64]
            J = add(J, study, sprintf('%s_array5_nw%d', a{1}, g), [meas cap a5 {'arch',a{1},'nw',g}], sw);
        end, end
    case 'tonespan'
        % tone window (user, 2026-10-19, item 4c): K = 4 tones from a window of toneSpan adjacent bins (detector
        % bandwidth (span - 1) x 50 MHz: 150/250/350/550 MHz; 16 = the whole band, 750 MHz, as in every other study), IS and IC at the base values,
        % and the same with 30 ps timing jitter (cap 3: jitter runs cannot reach the target)
        for a = {'IS', 'IC'}
            for sp = [4 6 8 12 16]
                J = add(J, study, sprintf('%s_span%d', a{1}, sp), [meas {'arch',a{1},'toneSpan',sp}], sw);
            end
            for sp = [4 16]
                J = add(J, study, sprintf('%s_span%d_jit30', a{1}, sp), [meas cap {'arch',a{1},'toneSpan',sp,'jitterPs',30}], sw);
            end
        end
    case 'calerrics'
        % two-pass IC (user, 2026-10-19): one full phi rotation with half the readings, turn the arch one theta step,
        % a second full phi rotation with the other half (icShiftFrac 1/2: half of every rotation's readings use the
        % turned ring). Every cable-length error is visible (ic_shift_observability.m). With and without
        % self-calibration over the cable-error grid, and one point without errors (same accuracy as plain IC?)
        ics = {'arch','IC','icShiftFrac',1/2};
        J = add(J, study, 'ICS_cal0', [meas ics], sw);
        for e = S.sweeps.calErr.IC
            J = add(J, study, sprintf('ICS_cal%g_em30', e), [meas cap ics {'calErrDeg',e,'nitEM',30}], sw);
            J = add(J, study, sprintf('ICS_cal%g_sc', e), [meas cap ics {'calErrDeg',e, ...
                'selfCal',true,'nitEM',30,'selfCalPriorDeg',e}], sw);
        end
    case 'icbad'
        % Fig. 3(c) (user, 2026-10-19): IC with the 2D Fourier basis and no legs at the paper's settings (32 x 16,
        % both components, m/n 0.3, plain prior), keeping the recovered pattern. It fails by design (the block-
        % phase ambiguity), so every restart would run to the end: cap 3.
        J = add(J, study, 'IC_dft_legs0', [meas {'arch','IC','angBasis','dft','npol',2,'nw',32, ...   % round 7: 10 restarts
            'groupPrior','none','nlegs',0}], sw);
        J(end).saveXh = true;
    case 'calerrrc'
        % TODO (user, 2026-10-19), NOT queued: RC self-calibration curve, the architecture whose ports sweep every
        % theta row (no invisible smooth cable errors). RC with and without self-calibration, 30 EM updates each.
        for e = S.sweeps.calErr.IS
            J = add(J, study, sprintf('RC_cal%g_em30', e), [meas cap {'arch','RC','calErrDeg',e,'nitEM',30}], sw);
            J = add(J, study, sprintf('RC_cal%g_sc', e), [meas cap {'arch','RC','calErrDeg',e,'selfCal',true, ...
                'nitEM',30,'selfCalPriorDeg',e}], sw);
        end
    case 'calerr30'
        % control for pr_calerrsc (user, 2026-10-19): the same cable errors and the same 30 EM updates WITHOUT
        % self-calibration, so the self-calibration gain is not confounded with the larger iteration budget
        for a = {'IS', 'IC'}
            for e = S.sweeps.calErr.(a{1})
                J = add(J, study, sprintf('%s_cal%g_em30', a{1}, e), [meas cap {'arch',a{1},'calErrDeg',e,'nitEM',30}], sw);
            end
        end
    case 'jitter'
        % detector timing jitter (user, 2026-10-19): rms error of the sample times in the data only, IS and IC
        % at the base values (tones drawn from the whole band, beat frequencies up to 750 MHz)
        for a = {'IS', 'IC'}, for jp = [1 3 10 30 100]
            J = add(J, study, sprintf('%s_jit%gps', a{1}, jp), [meas cap {'arch',a{1},'jitterPs',jp}], sw);
        end, end
    case 'detbw'
        % (NOT queued, ask) finite detector response in the data only (cpr_detector): a first-order low-pass with 3 dB bandwidth B
        % (Inf = ideal) and a centered integration window T, IS and IC at the base values; tones from the whole band (beats up to
        % 750 MHz) or from a window of 8 bins (up to 350 MHz). The solver keeps the ideal |Ax| model. Restart cap as for jitter.
        for a = {'IS', 'IC'}
            for sp = [16 8]
                for bw = [Inf 4 2 1 0.5 0.25]
                    ov = {'arch',a{1},'toneSpan',sp}; if isfinite(bw), ov = [ov {'detBwGHz',bw}]; end %#ok<AGROW>
                    J = add(J, study, sprintf('%s_span%d_bw%g', a{1}, sp, bw), [meas cap ov], sw);
                end
            end
            for ti = [100 250 500 1000]
                J = add(J, study, sprintf('%s_span16_box%gps', a{1}, ti), [meas cap {'arch',a{1},'detModel','boxcar','detIntPs',ti}], sw);
            end
        end
    case 'smooth'
        % (NOT queued, ask) calibration errors with structure, in the data only. (1) length errors smooth along the arch
        % (calErrModel 'smooth') for IC, two-pass IC (ICS), RC and RS at 1 and 3 deg rms at the band centre, without correction and
        % with per-port self-calibration (30 EM updates); (2) a smooth gain error (0.5 dB) with 3 deg smooth lengths, two-pass IC:
        % no correction, lengths only, lengths and gains; (3) cable dispersion (quadratic in frequency) for IS and IC.
        shapes = {'IC', {'arch','IC'}; 'ICS', {'arch','IC','icShiftFrac',1/2}; 'RC', {'arch','RC'}; 'RS', {'arch','RS'}};
        for i = 1:size(shapes, 1)
            for e = [1 3]
                dr = [shapes{i, 2} {'calErrDeg',e,'calErrModel','smooth','nitEM',30}];
                J = add(J, study, sprintf('%s_smooth%g_em30', shapes{i, 1}, e), [meas cap dr], sw);
                J = add(J, study, sprintf('%s_smooth%g_sclen', shapes{i, 1}, e), [meas cap dr {'selfCal',true,'selfCalPriorDeg',e}], sw);
            end
        end
        cg = {'arch','IC','icShiftFrac',1/2,'calErrDeg',3,'calErrModel','smooth','calGainDb',0.5,'calGainModel','smooth','nitEM',30};
        J = add(J, study, 'ICS_sgain_em30', [meas cap cg], sw);
        J = add(J, study, 'ICS_sgain_sclen', [meas cap cg {'selfCal',true,'selfCalPriorDeg',3}], sw);
        J = add(J, study, 'ICS_sgain_scgain', [meas cap cg {'selfCal',true,'selfCalPriorDeg',3,'selfCalGain',true,'selfCalPriorGainDb',0.5}], sw);
        for a = {'IS', 'IC'}, for dd = [2 5]
            ov = {'arch',a{1},'calDispDeg',dd,'nitEM',30};
            J = add(J, study, sprintf('%s_disp%g_em30', a{1}, dd), [meas cap ov], sw);
            J = add(J, study, sprintf('%s_disp%g_sclen', a{1}, dd), [meas cap ov {'selfCal',true,'selfCalPriorDeg',dd}], sw);
        end, end
    case 'nosparse'
        % band-limit-only ablation (user, 2026-10-19; answers "how would an approach without compressive sensing
        % perform?"): the same spherical-wave basis with a Gaussian prior on every coefficient, IS and IC over
        % the ratio grid. Compare with pr_ratio (sparsity prior).
        for a = {'IS', 'IC'}, for rt = S.sweeps.ratio
            J = add(J, study, sprintf('%s_nosparse_r%.2f', a{1}, rt), [meas {'arch',a{1},'ratio',rt,'sparsePrior',false}], sw);   % round 7: equal effort with pr_ratio
        end, end
    case 'cables'
        % do random cable lengths help (user, 2026-10-19)? phase-matched cables (all lengths 0) vs U(0, 50) cm,
        % IS and IC at the base values, K = 1 (conjugate ambiguity exact for matched cables) and K = 4; one seed
        for a = {'IS', 'IC'}, for cm = [0 50], for K = [1 4]
            J = add(J, study, sprintf('%s_cable%d_K%d', a{1}, cm, K), [meas {'arch',a{1},'cableMaxCm',cm,'K',K}], sw);
        end, end, end
    case 'L'
        for L = S.sweeps.L, J = add(J, study, sprintf('L%d', L), [meas {'vshL',L}], sw); end
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
        J = add(J, study, 'IS_cal0', [meas cap], sw);
        for e = S.sweeps.calErr.IS, J = add(J, study, sprintf('IS_cal%g', e), [meas cap {'calErrDeg',e}], sw); end
        J = add(J, study, 'IC_cal0', [meas cap {'arch','IC'}], sw);
        for e = S.sweeps.calErr.IC, J = add(J, study, sprintf('IC_cal%g', e), [meas cap {'arch','IC','calErrDeg',e}], sw); end
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
