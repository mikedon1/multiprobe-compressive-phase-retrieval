function [Rring, info] = fwc_ring_calibration(V, dirs, varargin)
%FWC_RING_CALIBRATION  Embedded probe responses of the FULL ring alone (64 ports) from the near-field exports of the ring-alone HFSS solve.
%   [Rring, info] = fwc_ring_calibration(V, dirs)
%   [Rring, info] = fwc_ring_calibration(V, dirs, 'AllowPartial', true)   missing ports fall back on the arch response (tests only)
% V     fwc_validate_hfss result (V.map: the 32 arch ports as (pol, ring_index) in the order of V.R; V.A.modes)
% dirs  folder or cell array of folders, each with ringalone_fullring_map.csv and ringalone_fullring_nf/ (fwc_check_ring_resp)
% Rring 64 x nmode x nf in the row order of fwc_ring_true: rows 1..32 near half in the order of V.map, rows 33..64 far half in the same order.
%       A near port uses its measured response; a far port uses its own measured response if exported, else (-1)^m times the near one (the ring
%       alone has an exact 180-deg symmetry, confirmed to 0.1 dB by the first exports). Use with fwc_run_trials(..., 'Ring', 'true',
%       'RingData', T, 'RingCal', 'ring') where T = fwc_ring_true(...) with T.Rring = Rring.
% info  .nNear, .nFar (measured ports used), .missing (arch ports without a measured near response), .fromDirs
p = inputParser; p.addParameter('AllowPartial', false); p.parse(varargin{:}); allowPartial = p.Results.AllowPartial;
if ischar(dirs) || isstring(dirs), dirs = {char(dirs)}; end
np = height(V.map); modes = V.A.modes; sgn = reshape((-1) .^ modes(:, 2), 1, []);
nmode = size(V.R, 2); nf = size(V.R, 3);
near = nan(np, 1); far = nan(np, 1); Rn = zeros(np, nmode, nf); Rf = zeros(np, nmode, nf); src = strings(np, 2);
for d = 1:numel(dirs)
    o = fwc_check_ring_resp(V, dirs{d});                          % prints the comparison with the arch calibration
    for q = 1:height(o.map)
        k = find(V.map.pol == o.map.pol(q) & V.map.ring_index == o.map.ring_index(q));
        if numel(k) ~= 1, error('fwc_ring_calibration: port %s %d not unique', o.map.pol(q), o.map.ring_index(q)); end
        if o.map.half(q) == "near", Rn(k, :, :) = o.R(q, :, :); near(k) = 1; src(k, 1) = string(dirs{d});
        else, Rf(k, :, :) = o.R(q, :, :); far(k) = 1; src(k, 2) = string(dirs{d}); end
    end
end
missing = find(isnan(near));
if ~isempty(missing) && ~allowPartial
    error('fwc_ring_calibration: %d near ports have no export (first: %s %d). Pass ''AllowPartial'', true to use the arch response there.', ...
        numel(missing), V.map.pol(missing(1)), V.map.ring_index(missing(1)));
end
Rarch = V.R;
for k = missing(:).', Rn(k, :, :) = Rarch(k, :, :); end          % partial test only
for k = find(isnan(far)).', Rf(k, :, :) = Rn(k, :, :) .* sgn; end % far half by symmetry
Rring = [Rn; Rf];
info = struct('nNear', nnz(~isnan(near)), 'nFar', nnz(~isnan(far)), 'missing', missing, 'fromDirs', {dirs});
fprintf('fwc_ring_calibration: %d near ports measured, %d far ports measured, %d near ports from the arch response\n', info.nNear, info.nFar, numel(missing));
end
