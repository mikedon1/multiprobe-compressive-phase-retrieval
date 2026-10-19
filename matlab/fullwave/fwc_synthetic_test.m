function fwc_synthetic_test(outDir)
%FWC_SYNTHETIC_TEST  End-to-end check of the embedded-probe calibration chain with exact dipole fields.
% Builds an HFSS-like case in HFSS axes: AUT = three Hertzian dipoles along the AUT's tilted axis
% (10 deg from +y toward +z) plus one offset dipole (breaks symmetry), rotated about HFSS +y by psi
% (right-handed, as AEDT rotates); 32 probe dipoles on the dual ring (theta-hat at R = 1 m, x-directed at
% R = 1.01 m, theta = 20..90 deg from +y in the x = 0 plane). Writes the AUT far field and each probe's
% near field on two spheres as AEDT-style data-table CSVs, reads them back with read_aedt_table, runs
% fwc_aut_modes, fwc_probe_resp and fwc_predict, and compares with the exact couplings p_b . E_a(r_b)
% at four rotations for both rotation senses. Expected: about -80 dB or better for sense +1, poor
% for -1. (No probe-probe coupling here: the test checks frames, readers, fits and the rotation
% sense, not the physics of the embedded responses.)
if nargin < 1, outDir = fullfile(tempdir, 'fwc_synth'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
maxNumCompThreads(1);
c0 = 299792458; freq = [2.5e9 3.125e9];
% AUT dipoles (m), HFSS axes
u = [0, cosd(10), sind(10)];
ra = [-0.03*u; 0*u; 0.03*u; 0.01 0.005 -0.02];
pa = [u; 1.6*u; u; 0.3*[1 0.5i 0]];
% ring
thr = 20 + (0:15)*70/15;
rb = [zeros(16,1), cosd(thr(:)), sind(thr(:))];
pb = [rb(:,1)*0, -sind(thr(:)), cosd(thr(:)); repmat([1 0 0], 16, 1)];
rb = [rb * 1.00; rb * 1.01];
psi = (0:3) * 11.25;
% AUT far field (rE up to a constant), HFSS grid of 3D_ATK
thF = -180:2:180; phF = 0:5:180;
[TH, PH] = ndgrid(thF*pi/180, phF*pi/180);
n = [sin(TH(:)).*cos(PH(:)), sin(TH(:)).*sin(PH(:)), cos(TH(:))];
tt = [cos(TH(:)).*cos(PH(:)), cos(TH(:)).*sin(PH(:)), -sin(TH(:))];
pv = [-sin(PH(:)), cos(PH(:)), zeros(numel(TH),1)];
D = struct();
for f = 1:numel(freq)
    k = 2*pi*freq(f)/c0; F = zeros(size(n));
    for i = 1:size(ra,1)
        F = F + cross(cross(n, repmat(pa(i,:), size(n,1), 1), 2), n, 2) .* exp(1i*k*(n*ra(i,:).'));
    end
    D.rETheta(:, :, f) = reshape(sum(F.*tt, 2), size(TH));
    D.rEPhi(:, :, f) = reshape(sum(F.*pv, 2), size(TH));
end
write_tab(fullfile(outDir, 'aut_ff.csv'), thF, phF, freq, D, 'V');
% probe near fields on two spheres
radii = [0.12 0.17]; thN = 0:5:180; phN = 0:5:355;
[TH, PH] = ndgrid(thN*pi/180, phN*pi/180);
n = [sin(TH(:)).*cos(PH(:)), sin(TH(:)).*sin(PH(:)), cos(TH(:))];
for p = 1:32
    for s = 1:2
        D = struct();
        for f = 1:numel(freq)
            k = 2*pi*freq(f)/c0;
            E = dipole_field(pb(p,:).', rb(p,:).', radii(s)*n, k);
            D.NearEX(:, :, f) = reshape(E(:,1), size(TH));
            D.NearEY(:, :, f) = reshape(E(:,2), size(TH));
            D.NearEZ(:, :, f) = reshape(E(:,3), size(TH));
        end
        write_tab(fullfile(outDir, sprintf('nf%d_p%02d.csv', s, p)), thN, phN, freq, D, 'V/m');
    end
end
% exact couplings
Sx = zeros(32, numel(psi), numel(freq));
for f = 1:numel(freq)
    k = 2*pi*freq(f)/c0;
    for j = 1:numel(psi)
        Ry = [cosd(psi(j)) 0 sind(psi(j)); 0 1 0; -sind(psi(j)) 0 cosd(psi(j))];
        for p = 1:32
            E = zeros(1,3);
            for i = 1:size(ra,1)
                E = E + dipole_field((Ry*pa(i,:).'), Ry*ra(i,:).', rb(p,:), k);
            end
            Sx(p, j, f) = E * pb(p,:).';
        end
    end
end
% pipeline
A = fwc_aut_modes(read_aedt_table(fullfile(outDir, 'aut_ff.csv')), 12);
fprintf('AUT fit residual: %s dB\n', mat2str(round(20*log10(A.relres), 1)));
Lreg = 20; R = zeros(32, size(A.a, 1), numel(freq));
for p = 1:32
    Ts = {read_aedt_table(fullfile(outDir, sprintf('nf1_p%02d.csv', p))), ...
          read_aedt_table(fullfile(outDir, sprintf('nf2_p%02d.csv', p)))};
    [Rp, info] = fwc_probe_resp(Ts, radii, Lreg, 12, 0.145);
    R(p, :, :) = permute(Rp, [3 1 2]);
    if p == 1 || p == 17
        fprintf('port %d: regular fit %s dB, radius consistency %s dB\n', p, ...
            mat2str(round(20*log10(info.fitres.'), 1)), mat2str(round(20*log10(info.radres.'), 1)));
    end
end
for sense = [1 -1]
    Sp = fwc_predict(A.a, A.modes, R, psi*pi/180, sense);
    [e, kap] = fwc_fit_scale(Sp, Sx);
    fprintf('sense %+d: NMSE per frequency %s dB, |kappa| %s\n', sense, mat2str(round(e, 1)), mat2str(abs(kap), 4));
end
end

function write_tab(file, x, phi, freq, D, unit)
names = fieldnames(D); lab = {'"Theta [deg]"'}; M = x(:);
for q = 1:numel(names)
    for part = {'re', 'im'}
        for f = 1:numel(freq)
            for j = 1:numel(phi)
                lab{end+1} = sprintf('"%s(%s) [%s] - Freq=''%gGHz'' Phi=''%gdeg''"', part{1}, names{q}, unit, freq(f)/1e9, phi(j)); %#ok<AGROW>
                v = D.(names{q})(:, j, f);
                if strcmp(part{1}, 're'), M(:, end+1) = real(v); else, M(:, end+1) = imag(v); end %#ok<AGROW>
            end
        end
    end
end
fid = fopen(file, 'w'); fprintf(fid, '%s\n', strjoin(lab, ',')); fclose(fid);
writematrix(M, file, 'WriteMode', 'append');
end

function E = dipole_field(p, r0, P, k)
R = P - r0(:).'; d = sqrt(sum(R.^2, 2)); n = R ./ d;
g = exp(-1i*k*d); pn = n * p(:);
nxp = cross(n, repmat(p(:).', size(n,1), 1), 2);
E = (1/(4*pi)) * (k^2 * cross(nxp, n, 2) .* (g./d) + (3*n.*pn - p(:).') .* ((1./d.^3 + 1i*k./d.^2) .* g));
end
