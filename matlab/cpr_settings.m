function S = cpr_settings(preset)
%CPR_SETTINGS  All user-adjustable settings in one place.
%
%   S = cpr_settings()          current preset (default 'paper', or the CPR_PRESET environment variable)
%   cpr_settings('quick')       switch to the quick preset for the rest of this MATLAB session
%   cpr_settings('paper')       switch back to the settings used for the paper's results
%
% 'paper' reproduces the paper (a full run takes days on a laptop: the circular architectures and
% the near-field comparison dominate). 'quick' runs EVERY study on a smaller pattern with few
% restarts and a few sweep points, in about 1.5 hours. Its numbers will NOT match the
% paper; it is meant to check that the code runs and to give the trends quickly.
% Other base parameters (P = 4 probes, K = 4 tones, sampling ratio, cable model, ...) are in
% core/cpr_cfg.m.
persistent active
if nargin > 0, active = preset; end
if isempty(active), active = getenv('CPR_PRESET'); end
if isempty(active), active = 'paper'; end

% =============================================================================================
% 1. PATTERN RESOLUTION AND FREQUENCIES
S.nw = 32;            % theta samples over 0-360 deg; phi uses nw/2 samples over 0-180 deg
                      %   (32 x 16 is near the mode-count Nyquist rate for vshL = 12, so m/n is a true
                      %   compression ratio; the DFT basis of earlier versions used 64)
S.nd = 16;            % number of wavelengths (50 MHz steps from 1.90 GHz; at most 32)
S.angBasis = 'vsh';   % angular sparse basis: 'vsh' (far-field vector spherical harmonics, the exact pattern
                      %   evaluated at the scan points) | 'dft' (2D DFT, earlier versions; needs leg probes)
S.vshL = 12;          % truncation degree of the 'vsh' basis
S.npol = 2;           % polarizations measured and recovered: 2 = dual-polarized probes (both tangential
                      %   components, shared mode coefficients) | 1 = one component (S.component)
S.nlegsCircular = 0;  % IC/RC leg probes: none are needed with the 'vsh' basis (the 'dft' basis needs about 8)
S.groupPrior = 'modesep'; % phaseless sparsity prior: 'modesep' (structured: each coefficient's rate = a factor
                      %   of its VSH mode, shared over frequency, times a factor of its frequency index; about
                      %   +1 dB over one rate, pilot P1) | 'none' (one rate for all coefficients)
S.component = 'across'; % recovered far-field component: 'across' (linear, across the scan circle, what a
                      %   linearly polarized arch probe measures), 'along', or 'cprh' (circular, earlier versions)
S.snrModel = 'detector'; % 'detector': fixed detector noise, SNR referred to one probe at one tone with the full
                      %   transmit power (tones share the power); 'normalized': same per-row SNR for every configuration

% 2. PHASE-RETRIEVAL SOLVER EFFORT   [spherical IS/RS , circular IC/RC]
%    One setting for all architectures with the 'vsh' basis (test T8); the 'dft' basis of earlier
%    versions used [10 30], [6 20], [200 1000], [40 48].
S.restarts  = [10 10];     % maximum random restarts (the solver stops at the first one that
                           %   reaches its residual target, so easy runs use far fewer)
S.nitEM     = [6 6];       % expectation-maximization updates per restart
S.nitGAMP   = [200 200];   % GAMP iterations per update
S.solverSNR = [40 40];     % stopping setting SNR_s in dB (residual target -(SNR_s + 2) dB)
S.ratioCircular = 0.3;     % sampling ratio m/n of the circular architectures (phaseless). Same as IS/RS: with the 'vsh'
                           %   basis all four perform alike (user, 2026-10-19); the 'dft' basis used 0.6

% 3. NUMBER OF TRIALS (random seeds) PER POINT
S.trials.sweep      = 1;   % points of the Section V sweeps and the supplementary sweeps
S.trials.snrCircular = 1;  % IC and RC points of the SNR sweep (1 like every sweep, user 2026-10-19; was 3)
S.trials.nearfield  = 1;   % Section VII compressive phase retrieval runs
S.trials.fullwave   = 10;  % Section VII full-wave trials (one trial: 0.5 to 3 h)
S.fullwaveExtras    = true; % Section VII: also m/n 0.6, half-ring calibration, perturbed calibration, MoM arch calibration (false: main trials only)

% 4. SWEEP GRIDS (the x-axis values of each figure)
S.sweeps.ratio  = [0.05 0.08 0.12 0.2 0.3 0.4 0.6 0.8 1.0];  % sampling ratio m/n (n = 32 x 16 x 16 x 2)
S.sweeps.Ppr    = [1 2 4 8 16];         % active probes, phaseless measurements
S.sweeps.K      = [1 2 4 8 16];         % tones per reading
S.sweeps.rnRS   = [1/16 1/8 1/4 1/2 1];  % installed-probe fraction, RS (n_phi = 16: 1, 2, 4, 8, 16 probes per row;
                                         %   1/32 would round to the same single probe as 1/16)
S.sweeps.rnRC   = [1/16 1/8 1/4 1/2 1];     % installed-probe fraction, RC
S.sweeps.ndPr   = [1 2 4 8 16];         % wavelengths, phaseless
S.sweeps.snr    = [10 20 30 40 50 60 80];    % SNR in dB
S.sweeps.examples = struct('pr', [0.1 0.15 0.2 0.3 0.4]);   % pattern-cut examples
S.sweeps.gaussRatio = [0.08 0.1 0.15 0.2 0.3];   % sampling ratios of the switches-vs-Gaussian-weights comparison (IS)
S.sweeps.calErr = struct('IS', [1 2 5 10 20], 'IC', [1 2 5 10 20]);   % same grid for both (user 2026-10-19)   % rms cable phase errors in degrees

% 5. FULL-WAVE EXAMPLE (Section VII, trials above) AND NEAR-FIELD COMPARISON (Section VIII)
S.nearfield.snr         = 50;     % dB
S.nearfield.twoScanStarts = 400;  % maximum random starts of the two-scan technique per frequency

% 6. PATHS AND COMPUTE
S.gampPath   = '';    % folder with GAMPmatlab (see cpr_setup.m); empty: $GAMPMATLAB or external/GAMPmatlab
S.resultsDir = '';    % empty: <repository>/results
S.threads    = 2;     % computational threads per MATLAB process

% 7. FIGURES
S.figures.print = true;   % false: only draw the figures in MATLAB (they stay open on screen);
                          % true: also write each one as PDF and PNG to <results>/figs/paper
S.figures.axes  = 'paper'; % 'paper': the axis limits and ticks used in the paper; 'auto': automatic limits and
                          %   ticks taken from the data (the quick preset uses this, since its sweeps are shorter)
% =============================================================================================

if strcmp(active, 'quick')                          % ---- QUICK PRESET (overrides) ----
    S.nw = 32;  S.nd = 8;
    S.restarts = [3 5];  S.nitEM = [6 6];  S.nitGAMP = [200 300];  S.ratioCircular = 0.8;
    S.trials.snrCircular = 1;  S.trials.fullwave = 1;  S.fullwaveExtras = false;
    S.sweeps.ratio = [0.1 0.3 0.6];
    S.sweeps.Ppr = [2 4];
    S.sweeps.K = [1 4];
    S.sweeps.rnRS = [1/8 1/2];   S.sweeps.rnRC = [1/4 1];
    S.sweeps.ndPr = [1 2 4 8];
    S.sweeps.snr = [20 50 80];
    S.sweeps.examples = struct('pr', 0.3);   % one example per figure
    S.sweeps.gaussRatio = [0.2 0.3];
    S.sweeps.calErr = struct('IS', [2 10], 'IC', 5);
    S.figures.axes = 'auto';
    S.nearfield.twoScanStarts = 10;
elseif ~strcmp(active, 'paper')
    error('cpr_settings: unknown preset "%s" (use ''paper'' or ''quick'')', active);
end
S.preset = active;
if ~isempty(getenv('CPR_FIGAXES')), S.figures.axes = getenv('CPR_FIGAXES'); end   % 'auto' | 'paper' (session override)
end
