function err = cpr_check_component()
%CPR_CHECK_COMPONENT  Check the 'across' component of pattern_crlh against a direct computation.
% On the 2-degree source grid (nw = 180), the scan-grid value at scan angles (theta_s, phi_s) must
% equal E . n, the far field of the (90-degree rotated) antenna in direction
% (theta_s, phi_s) for theta_s <= 180 or (360 - theta_s, phi_s + 180) otherwise, projected on the
% fixed normal n = (-sin phi_s, cos phi_s, 0) of the scan circle. Returns the relative error.
S0 = cpr_settings(); if isfield(S0, 'component') && ~strcmp(S0.component, 'across')
    error('cpr_check_component: set cpr_settings component to ''across''');
end
D = crlh_data(); k = 1;
[X, ~] = pattern_crlh(180, 1, k);                    % 180 x 90, scan theta 0:2:358, phi 0:2:178
Et = D.E_theta{k}; Ep = D.E_phi{k};                  % 91 x 181, theta 0:2:180, phi 0:2:360
th = D.theta(:); ph = D.phi(:).';
Y = zeros(180, 90);
for j = 1:90
    phs = (j-1)*2;
    n = [-sind(phs), cosd(phs), 0];
    for i = 1:180
        ths = (i-1)*2;
        if ths <= 180, t = ths; p = phs; else, t = 360 - ths; p = mod(phs + 180, 360); end
        ps = mod(p - 90, 360);                       % source phi after the 90-degree rotation
        it = round(t/2) + 1; ip = round(ps/2) + 1;
        thv = [cosd(t)*cosd(p), cosd(t)*sind(p), -sind(t)];   % theta-hat at (t, p)
        phv = [-sind(p), cosd(p), 0];                         % phi-hat at (t, p)
        Evec = Et(it, ip)*thv + Ep(it, ip)*phv;
        Y(i, j) = Evec * n.';
    end
end
a = (Y(:)'*X(:)) / (Y(:)'*Y(:));                     % X is normalized: compare up to a complex scale
err = norm(X(:) - a*Y(:)) / norm(X(:));
fprintf('across component vs direct E . n: relative error %.2e\n', err);
end
