function c = cpr_cfg(varargin)
%CPR_CFG  Parameters of one simulation (single source of truth for the base values).
%   c = cpr_cfg()                      base configuration
%   c = cpr_cfg('K', 8, 'arch', 'IC')  base values with the listed fields changed
% Pattern size, restarts, solver effort and SNR_s come from cpr_settings.m; the remaining base
% parameters are below. Every sweep changes exactly one field.
S = cpr_settings();

% ---- measurement -------------------------------------------------------------------------------
c.meas       = 'phaseless'; % 'complex' (coherent receiver) | 'phaseless' (fast square-law detector)
c.arch       = 'IS';        % 'IS' | 'RS' | 'IC' | 'RC'
c.P          = 4;           % active probes per reading (on/off switches)
c.K          = 4;           % tones per reading (random subset of the nd wavelengths)
c.ratio      = 0.3;         % sampling ratio m/n (phaseless IS/RS; circular and complex bases below)
%   complex  : one reading gives K per-tone complex values (same probe set)
%   phaseless: one reading gives K time samples of |sum of tones|^2 (amplitude used)
c.rnRS       = 1/16;        % RS installed-probe fraction r_n
c.rnRC       = 1/2;         % RC installed-probe fraction r_n
c.nlegs      = 0;           % leg probes (IC/RC: S.nlegsCircular, set below)
c.legDir     = 'phi';       % 'phi': leg q sits at the anchor row (theta = 90 deg) of the rotation q steps ahead
                            %   (as in the full-wave model); 'theta': legs at theta rows nw-nlegs+1..nw of the next rotation
c.cableMaxCm = 50;          % cable lengths ~ U(0, cableMaxCm) cm, drawn at random on purpose and then fixed;
                            %   phase 2*pi*L/lambda per tone; calibrated (0 = phase-matched cables)
c.snrdB      = 100;         % near-noiseless base; the SNR sweep varies it
c.calErrDeg  = 0;           % cable calibration error, rms phase at the band centre (deg; data use the true phases)
c.calErrModel = 'length';   % 'length': a fixed length error per port (phase error grows with frequency);
                            %   'tone': independent per port and tone (earlier versions)
c.selfCal    = false;       % phaseless VSH runs: estimate a cable length correction per port during recovery
c.selfCalStart = 2;         %   first EM update after which the calibration is re-estimated
c.selfCalGN  = 5;           %   Gauss-Newton steps per re-estimation
c.selfCalPriorDeg = NaN;    %   expected rms calibration error (deg at the band centre) for a MAP prior; NaN: none
c.selfCalModel = 'length';  %   'length': one length per port | 'drift': one fractional length change for all cables
c.selfCalGain = false;      %   also estimate one gain per port
c.selfCalPriorGainDb = NaN; %   expected rms gain error (dB) for its MAP prior; NaN: none
c.calGainDb  = 0;           % fixed gain error per port in the data, rms dB (0 = none)
c.jitterPs   = 0;           % phaseless: rms error of the detector sample times (ps), random per sample, in the data
                            %   only (the solver assumes the nominal times); a fixed offset would be a common delay
c.icShift    = 1;           % IC calibration step: theta samples the ring is turned for part of the readings
c.icShiftFrac = 0;          %   fraction of the readings taken with the turned ring (0 = plain IC)
c.icShiftMode = 'perRot';   %   'perRot': last fraction of every rotation's readings; 'end': last fraction of all
                            %   readings (one turn of the ring at the end of the scan)
c.combNorm   = 'unit';      % combiner: 'unit' (idealized buffered, unit branch transfers) | 'passive' (1/sqrt(P))
c.toneSpan   = Inf;         % phaseless: each reading draws its K tones from a window of toneSpan adjacent frequency
                            %   bins at a random position (detector bandwidth (toneSpan - 1)*df); Inf = whole band
c.sparsePrior = true;       % false: band-limit-only ablation, a Gaussian prior on every coefficient (rate 1)
c.weights    = 'switch';    % 'switch' (on/off, the paper's system) | 'gauss' (benchmark: each active probe
%                             also gets a CN(0,1) weight, as an attenuator and phase shifter would give;
%                             drawn from a separate random stream, so probes, tones and times match the
%                             'switch' run of the same seed)

% ---- pattern (CRLH leaky-wave antenna, openEMS) ------------------------------------------------
c.nw         = S.nw;        % theta samples, 0-360 deg (nh = nw/2 phi samples, 0-180 deg)
c.nd         = S.nd;        % wavelengths (50 MHz steps)
c.f0idx      = 1;           % first frequency index (1 -> 1.90 GHz)
c.pattern    = 'crlh';      % example antenna: 'crlh' | 'array4' / 'array5' (array_data: n x n dipole array, 6 cm
                            %   spacing, the electrically larger second example, 2026-10-19)
c.patternLfit = 16;         % degree of the spherical-wave fit that defines the exact truth (30 for the arrays)

% ---- sparse basis: 2D DFT over angle x (DCT or DFT) along wavelength, with delay ---------------
c.angBasis   = S.angBasis;  % angular basis: 'vsh' (far-field vector spherical harmonics up to degree vshL;
                            %   the truth is the exact pattern at the scan points) | 'dft' (2D DFT over the grid)
c.vshL       = S.vshL;      % maximum degree of the 'vsh' basis
c.basis      = 'dct';       % 'dct' | 'dft'
c.gauge      = true;        % re-choose the representation between EM rounds:
%   phaseless K=1 : per-wavelength phase (unobservable) -> neighbor alignment
%   otherwise     : one global delay tau (pattern unchanged, basis re-expressed)
c.sqMode     = 'exact';     % |H|^2 products: 'exact' | 'uniform' (per-row average, faster)

% ---- solver (spherical values; the circular architectures are set below) -----------------------
c.solver     = 'cpr_prgamp'; % phaseless solver
c.minTry     = 1;           % restarts always run before early stopping is allowed
c.maxTry     = S.restarts(1);        % maximum restarts
c.nitEM      = S.nitEM(1);           % EM updates per restart
c.nitGAMP    = S.nitGAMP(1);         % GAMP iterations per update
c.solverSNRdB = S.solverSNR(1);      % accept a restart if its residual is below -(solverSNRdB + 2) dB
c.stopDB     = NaN;         % explicit acceptance target in dB; NaN: -(solverSNRdB + 2)
c.logRestarts = false;      % phaseless: score every restart (restart studies only)
c.npol       = S.npol;      % polarizations recovered jointly (2: dual-polarized probes, both components)
c.groupPrior = S.groupPrior; % phaseless sparsity prior: 'modesep' (structured, VSH basis) | 'none' (one rate)

% ---- per-measurement-type and per-architecture bases -------------------------------------------
meas = c.meas; arch = c.arch;
for i = 1:2:numel(varargin)
    if strcmp(varargin{i}, 'meas'), meas = varargin{i+1}; end
    if strcmp(varargin{i}, 'arch'), arch = varargin{i+1}; end
end
if strcmp(meas, 'complex')
    c.ratio = 0.08;         % sampling ratio
    c.nlegs = 0;            % complex measurements need no legs (no block-phase ambiguity)
elseif any(strcmp(arch, {'IC', 'RC'}))
    c.ratio = S.ratioCircular;  c.nlegs = S.nlegsCircular;  c.K = 4;
    c.maxTry = S.restarts(2);  c.nitEM = S.nitEM(2);  c.nitGAMP = S.nitGAMP(2);
    c.solverSNRdB = S.solverSNR(2);
end
for i = 1:2:numel(varargin)
    if ~isfield(c, varargin{i}), error('cpr_cfg: unknown field %s', varargin{i}); end
    c.(varargin{i}) = varargin{i+1};
end
c.nh = c.nw/2;
if ~strcmp(c.angBasis, 'vsh'), c.groupPrior = 'none'; end   % the structured prior is defined on VSH modes
c.K = min(c.K, c.nd);
end
