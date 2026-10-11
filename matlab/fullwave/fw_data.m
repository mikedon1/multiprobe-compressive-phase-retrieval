function D = fw_data(src, varargin)
%FW_DATA  "True pattern" for the Section VI full-wave IC geometry .
%
%   D = fw_data('analytic', 'NF', 16)          analytic dipole stand-in (tests)
%   D = fw_data('hfss', 'Glob', '...\sweep_redesign_mf_rot*.s20p')
%   D = fw_data('hfss17', 'Glob', '...\sweep_redesign_v1_rot*.s17p')   single-freq, NO legs
%
% Geometry: 16 ring probes on theta = 20..90 deg,
% 32 AUT rotations psi = 0:11.25:348.75 deg, AUT = 13.625 cm dipole tilted 10 deg,
% R = 100 cm, probes measure E_theta. Legs: copies of the anchor probe (ring index
% 16, theta = 90 deg) at psi offsets +k*11.25 deg, k = 1..3.
%
% D.Xring  16 x 32 x nf   ring responses (recovery target)
% D.Xleg   nl x 32 x nf   each leg's OWN response at each rotation (real ports or
%                         analytic at the leg's own position); [] if unavailable
% D.f      1 x nf (Hz);  D.anchor = 16;  D.legShift = s such that, ideally,
%          leg k at rotation L sees what the anchor sees at rotation L + s*k
%          (measured from the data here, not assumed)
% D.src    description
%
% HFSS port order: AUT = port 1; ring = ports 2..17 (theta ascending, VERIFIED
% against the 17-port sweep by cpr_fw_check_ports.m: best match is the identity);
% legs = ports 18..20 = k = 1..3 -- confirmed against the HFSS geometry: each leg port sheet sits at
% psi = 11.25*k exactly. Coupling note: port 20 (leg 3, the
% end leg) matches the anchor to ~1.4%, ports 18/19 (legs 1, 2) to ~12%.
p = inputParser;
p.addParameter('NF', 16);
p.addParameter('F0', 2.375e9);
p.addParameter('DF', 50e6);
p.addParameter('Glob', '');
p.parse(varargin{:});
o = p.Results;
nTh = 16; nPh = 32; dpsi = 11.25; anchor = 16; nl = 3;
switch src
    case 'analytic'
        f = o.F0 + (0:o.NF-1)*o.DF;
        lam = 29979245800 ./ f;                       % cm
        Xr = zeros(nTh, nPh, o.NF); Xl = zeros(nl, nPh, o.NF);
        for t = 1:o.NF
            args = {'NTheta', nTh, 'NPhi', nPh, 'ThetaDeg', [20 90], 'Rcm', 100, 'LambdaCm', lam(t)};
            Xr(:,:,t) = dipole_ic_pattern(13.625/lam(t), 10, args{:});
            for k = 1:nl
                Xk = dipole_ic_pattern(13.625/lam(t), 10, args{:}, 'ProbePsiDeg', k*dpsi);
                Xl(k,:,t) = Xk(anchor, :);
            end
        end
        D.src = sprintf('analytic dipole, %d freq %.3f-%.3f GHz', o.NF, f(1)/1e9, f(end)/1e9);
    case 'hfss'
        [Xp, info] = load_hfss_mf(o.Glob);
        if size(Xp, 1) ~= nTh + nl, error('fw_data:ports', 'expected %d probe ports, got %d', nTh+nl, size(Xp,1)); end
        if size(Xp, 2) ~= nPh, error('fw_data:rot', 'expected %d rotations, got %d', nPh, size(Xp,2)); end
        Xr = Xp(1:nTh, :, :); Xl = Xp(nTh+1:end, :, :); f = info.freq;
        D.src = sprintf('HFSS real legs, %d freq', numel(f));
    case 'hfss17'
        [Xp, info] = load_hfss_mf(o.Glob);
        Xr = Xp; Xl = []; f = info.freq;
        D.src = sprintf('HFSS 17-port (no leg ports), %d freq', numel(f));
    otherwise
        error('fw_data:src', 'unknown source %s', src);
end
s = max(abs(Xr(:))); Xr = Xr/s; if ~isempty(Xl), Xl = Xl/s; end
D.Xring = Xr; D.Xleg = Xl; D.f = f(:).'; D.anchor = anchor; D.nlegsAvail = size(Xl, 1);
% leg direction from data: compare leg 1 with the anchor at L-1 and L+1
D.legShift = -1;
if ~isempty(Xl)
    a = reshape(Xr(anchor, :, :), nPh, []); l1 = reshape(Xl(1, :, :), nPh, []);   % rotations x freq
    em = norm(l1 - circshift(a, 1, 1), 'fro'); ep = norm(l1 - circshift(a, -1, 1), 'fro');
    if ep < em, D.legShift = +1; end
    D.legResidual = min(em, ep)/norm(l1, 'fro');
    D.legResidualOther = max(em, ep)/norm(l1, 'fro');
end
end
