function S = cpr_settings(preset)
%CPR_SETTINGS  All user-adjustable settings in one place.
%
%   S = cpr_settings()          current preset (default 'paper', or the CPR_PRESET environment variable)
%   cpr_settings('quick')       switch to the quick preset for the rest of this MATLAB session
%   cpr_settings('paper')       switch back to the settings used for the paper's results
%
% 'paper' reproduces the paper (a full run takes days on a laptop: the circular architectures and
% the near-field comparison dominate). 'quick' runs EVERY study on a smaller pattern with few
% restarts and a few sweep points, in minutes to about an hour. Its numbers will NOT match the
% paper; it is meant to check that the code runs and to give the trends quickly.
% Other base parameters (P = 4 probes, K = 4 tones, sampling ratio, cable model, ...) are in
% core/cpr_cfg.m.
persistent active
if nargin > 0, active = preset; end
if isempty(active), active = getenv('CPR_PRESET'); end
if isempty(active), active = 'paper'; end

% =============================================================================================
% 1. PATTERN RESOLUTION AND FREQUENCIES
S.nw = 64;            % theta samples over 0-360 deg; phi uses nw/2 samples over 0-180 deg
S.nd = 16;            % number of wavelengths (50 MHz steps from 1.90 GHz; at most 32)

% 2. PHASE-RETRIEVAL SOLVER EFFORT   [spherical IS/RS , circular IC/RC]
S.restarts  = [10 30];     % maximum random restarts (the solver stops at the first one that
                           %   reaches its residual target, so easy runs use far fewer)
S.nitEM     = [6 20];      % expectation-maximization updates per restart
S.nitGAMP   = [200 1000];  % GAMP iterations per update
S.solverSNR = [40 48];     % stopping setting SNR_s in dB (residual target -(SNR_s + 2) dB)
S.ratioCircular = 0.6;     % sampling ratio m/n of the circular architectures (phaseless): they need more measurements

% 3. NUMBER OF TRIALS (random seeds) PER POINT
S.trials.sweep      = 1;   % points of the Section IV sweeps and the supplementary sweeps
S.trials.snrCircular = 3;  % IC and RC points of the SNR sweep (mean in dB over the trials)
S.trials.nearfield  = 1;   % Section V compressive phase retrieval runs
S.trials.fullwave   = 10;  % Section VI full-wave trials

% 4. SWEEP GRIDS (the x-axis values of each figure)
S.sweeps.ratio  = [0.03 0.05 0.08 0.12 0.2 0.3 0.4 0.6 0.8];  % sampling ratio m/n
S.sweeps.Pcs    = [1 2 4 8 16 32];      % active probes, complex measurements
S.sweeps.Ppr    = [1 2 4 8 16];         % active probes, phaseless measurements
S.sweeps.K      = [1 2 4 8 16];         % tones per reading
S.sweeps.rnRS   = [1/32 1/16 1/8 1/4 1/2];   % installed-probe fraction, RS
S.sweeps.rnRC   = [1/16 1/8 1/4 1/2 1];     % installed-probe fraction, RC
S.sweeps.ndCs   = [1 2 4 8 16 32];      % wavelengths, complex (joint vs one at a time)
S.sweeps.ndPr   = [1 2 4 8 16];         % wavelengths, phaseless
S.sweeps.snr    = [10 20 30 40 50 60 80];    % SNR in dB
S.sweeps.legs   = [0 1 2 4 8 16];       % leg probes (IC, RC)
S.sweeps.size   = [16 32 64 128];       % angular resolution nw
S.sweeps.examples = struct('cs', [0.03 0.05 0.08 0.2], 'pr', [0.1 0.15 0.2 0.3 0.4]);   % pattern-cut examples
S.sweeps.gaussRatio = [0.08 0.1 0.15 0.2 0.3];   % sampling ratios of the switches-vs-Gaussian-weights comparison (IS)
S.sweeps.calErr = struct('IS', [1 2 5 10 20], 'IC', [2 5 10]);   % rms cable phase errors in degrees

% 5. NEAR-FIELD COMPARISON (Section V) AND FULL-WAVE EXAMPLE (Section VI)
S.nearfield.snr         = 50;     % dB
S.nearfield.twoScanStarts = 400;  % maximum random starts of the two-scan technique per frequency
S.fullwave.nd           = 16;     % frequencies of the HFSS sweep to use (2.375-3.125 GHz in 50 MHz steps)
S.fullwave.nitGAMP      = 200;    % GAMP iterations per EM update in Section VI (the Section VI runs used 200;
                                  %   the circular value S.nitGAMP(2) applies to Sections IV and V)

% 6. PATHS AND COMPUTE
S.gampPath   = '';    % folder with GAMPmatlab (see cpr_setup.m); empty: $GAMPMATLAB or external/GAMPmatlab
S.resultsDir = '';    % empty: <repository>/results
S.threads    = 2;     % computational threads per MATLAB process
% =============================================================================================

if strcmp(active, 'quick')                          % ---- QUICK PRESET (overrides) ----
    S.nw = 32;  S.nd = 8;
    S.restarts = [3 5];  S.nitEM = [6 6];  S.nitGAMP = [200 300];  S.ratioCircular = 0.8;
    S.trials.snrCircular = 1;  S.trials.fullwave = 2;
    S.sweeps.ratio = [0.1 0.3 0.6];
    S.sweeps.Pcs = [2 4];    S.sweeps.Ppr = [2 4];
    S.sweeps.K = [1 4];
    S.sweeps.rnRS = [1/8 1/2];   S.sweeps.rnRC = [1/4 1];
    S.sweeps.ndCs = [1 2 4 8];   S.sweeps.ndPr = [1 2 4 8];
    S.sweeps.snr = [20 50 80];
    S.sweeps.legs = [0 4 8];
    S.sweeps.size = [16 32];
    S.sweeps.examples = struct('cs', [0.05 0.2], 'pr', [0.15 0.3]);
    S.sweeps.gaussRatio = [0.2 0.3];
    S.sweeps.calErr = struct('IS', [2 10], 'IC', 5);
    S.nearfield.twoScanStarts = 10;
    S.fullwave.nd = 8;  S.fullwave.nitGAMP = 200;
elseif ~strcmp(active, 'paper')
    error('cpr_settings: unknown preset "%s" (use ''paper'' or ''quick'')', active);
end
S.preset = active;
end
