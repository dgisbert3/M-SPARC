function [excScan, VxcScan1, VxcScan2, VxcScan3] = xcSScan(rho, normDrho, tau)
 % all three input variables have 3 cols
 % col1: n; col2: n up; col3: n dn
 % excScan has one column; all VxcScan1 and VxcScan3 have two columns;
 % VxcScan2 has three columns: Vc2, Vx2_up, Vx2_dn, corresonding to
 % normDrho
[s_up, alpha_up, DsDn_up, DsDDn_up, DalphaDn_up, DalphaDDn_up, DalphaDtau_up]...
    = basicMGGAvariables(2*rho(:, 2), 2*normDrho(:, 2), 2*tau(:, 2));
[ex_up, Vx1_up, Vx2_up, Vx3_up]...
    = exchangeSCAN(2*rho(:, 2), s_up, alpha_up, DsDn_up, DsDDn_up, DalphaDn_up, DalphaDDn_up, DalphaDtau_up);
[s_dn, alpha_dn, DsDn_dn, DsDDn_dn, DalphaDn_dn, DalphaDDn_dn, DalphaDtau_dn]...
    = basicMGGAvariables(2*rho(:, 3), 2*normDrho(:, 3), 2*tau(:, 3));
[ex_dn, Vx1_dn, Vx2_dn, Vx3_dn]...
    = exchangeSCAN(2*rho(:, 3), s_dn, alpha_dn, DsDn_dn, DsDDn_dn, DalphaDn_dn, DalphaDDn_dn, DalphaDtau_dn);
[s, alpha, DsDn, DsDDn, DalphaDn, DalphaDDn, DalphaDtau]...
    = basicMGGAvariables(rho(:, 1), normDrho(:, 1), tau(:, 1));
[zeta, DzetaDnup, DzetaDndn] = polarizeVarZeta(rho);
[ec, Vc1, Vc2, Vc3] = correlationSCAN(rho(:, 1), rho(:, 2), rho(:, 3), s, alpha, zeta, DsDn, DsDDn, DalphaDn, DalphaDDn, DalphaDtau, DzetaDnup, DzetaDndn);
excScan = (ex_up.*rho(:, 2) + ex_dn.*rho(:, 3)) ./ rho(:, 1) + ec;
VxcScan1 = [Vx1_up, Vx1_dn] + Vc1;
VxcScan2 = [Vc2, Vx2_up, Vx2_dn];
VxcScan3 = [Vx3_up + Vc3, Vx3_dn + Vc3];
end

function [s, alpha, DsDn, DsDDn, DalphaDn, DalphaDDn, DalphaDtau] ...
    = basicMGGAvariables(rho, normDrho, tau)
%     normDrho = (Drho(:,1).^2 + Drho(:,2).^2 + Drho(:,3).^2).^0.5;
    %% value
    s = normDrho ./ (2*(3*pi^2)^(1/3)*rho.^(4/3));
    tauw = normDrho.^2 ./ (8*rho);
    tauUnif = 3/10*(3*pi^2)^(2/3) * rho.^(5/3);
    alpha = (tau - tauw) ./ tauUnif; % when there is no spin, ds = 1
    %% variation
    DsDn = -2*normDrho ./ (3*(3*pi^2)^(1/3).*rho.^(7/3));
    DsDDn = 1 ./ (2*(3*pi^2)^(1/3)*rho.^(4/3));
    DtauwDn = -normDrho.^2 ./ (8*rho.^2);
    DtauwDDn = normDrho ./ (4*rho);
    DtauUnifDn = (3*pi^2)^(2/3)/2 * rho.^(2/3);
    DalphaDn = (-DtauwDn.*tauUnif - (tau - tauw).*DtauUnifDn) ./ (tauUnif.^2);
    DalphaDDn = (-DtauwDDn) ./ tauUnif;
    DalphaDtau = 1 ./ tauUnif;
end

function [ex, Vx1, Vx2, Vx3] = exchangeSCAN(rho, s, alpha, DsDn, DsDDn, DalphaDn, DalphaDDn, DalphaDtau)
    %% solve epsilon_x
    epsixUnif = -3/(4*pi)*(3*pi^2*rho).^(1/3);
    % compose h_x^1
    k1 = 0.065;
    muak = 10/81;
    b2 = sqrt(5913/405000);
    b1 = 511/13500/(2*b2);
    b3 = 0.5;
    b4 = muak*muak/k1 - 1606/18225 - b1*b1;
    % compose x
    x = muak*s.^2 .* (1 + b4*s.^2/muak .* exp(-abs(b4)*s.^2/muak))...
        + (b1*s.^2 + b2*(1 - alpha).*exp(-b3*(1 - alpha).^2)).^2;
    hx1 = 1 + k1 - k1./(1 + x/k1);
    % interpolate and extrapolate h_x to get F_x
    hx0 = 1.174;
    % switching function f_x
    c1x = 0.667;
    c2x = 0.8;
    dx = 1.24;
    alphaG1 = (alpha > 1);
    alphaE1 = (alpha == 1);
    alphaL1 = (alpha < 1);
    fx = zeros(size(rho, 1), size(rho, 2));
    fx(alphaG1) = -dx*exp(c2x ./ (1 - alpha(alphaG1)));
    fx(alphaE1) = 0.0;
    fx(alphaL1) = exp(-c1x*alpha(alphaL1) ./ (1 - alpha(alphaL1)));
    a1 = 4.9479;
    gx = 1 - exp(-a1*s.^(-0.5));
    Fx = (hx1 + fx.*(hx0 - hx1)).*gx;
    % get epsilon_x
    ex = epsixUnif.*Fx;
    %% solve variation of F_x
    s2 = s.*s;
    term1 = 1 + (b4*s2)/muak.*exp(-abs(b4)*s2/muak);
    term2 = s2.*(b4/muak*exp(-abs(b4)*s2/muak) + b4*s2/muak.*exp(-abs(b4)*s2/muak)*(-abs(b4)/muak));
    term3 = 2*(b1*s2 + b2*(1 - alpha).*exp(-b3*(1 - alpha).^2));
    term4 = b2*(-exp(-b3*(1-alpha).^2) + (1-alpha).*exp(-b3*(1-alpha).^2).*(2*b3*(1-alpha)));
    DxDs = 2*s.*(muak*(term1 + term2) + b1*term3);
    DxDalpha = term3.*term4;
    DxDn = DsDn.*DxDs + DalphaDn.*DxDalpha;
    DxDDn = DsDDn.*DxDs + DalphaDDn.*DxDalpha;
    DxDtau = DalphaDtau.*DxDalpha;
    
    DgxDn = -exp(-a1*s.^(-0.5)).*(a1/2*s.^(-1.5)).*DsDn;
    DgxDDn = -exp(-a1*s.^(-0.5)).*(a1/2*s.^(-1.5)).*DsDDn;
    Dhx1Dx = 1 ./ (1 + x/k1).^2;
    Dhx1Dn = DxDn .* Dhx1Dx;
    Dhx1DDn = DxDDn .* Dhx1Dx;
    Dhx1Dtau = DxDtau .* Dhx1Dx;
    DfxDalpha = zeros(size(rho, 1), size(rho, 2));
    DfxDalpha(alphaG1) = -dx*exp(c2x./(1-alpha(alphaG1))).*(c2x./(1-alpha(alphaG1)).^2);
    DfxDalpha(alphaE1) = 0.0;
    DfxDalpha(alphaL1) = exp(-c1x*alpha(alphaL1)./(1-alpha(alphaL1))).*(-c1x./(1-alpha(alphaL1)).^2);
    DfxDn = DfxDalpha.*DalphaDn;
    DfxDDn = DfxDalpha.*DalphaDDn;
    DfxDtau = DfxDalpha.*DalphaDtau;
    
    DFxDn = (hx1+fx.*(hx0-hx1)).*DgxDn + gx.*(1-fx).*Dhx1Dn + gx.*(hx0-hx1).*DfxDn;
    DFxDDn = (hx1+fx.*(hx0-hx1)).*DgxDDn + gx.*(1-fx).*Dhx1DDn + gx.*(hx0-hx1).*DfxDDn;
    DFxDtau = gx.*(1-fx).*Dhx1Dtau + gx.*(hx0-hx1).*DfxDtau;
    %% solve variant of n*epsilon_x^{unif}*F_x
    DepsixUnifDn = -(3*pi^2)^(1/3)/(4*pi) * rho.^(-2/3);
    Vx1 = (epsixUnif+rho.*DepsixUnifDn).*Fx + rho.*epsixUnif.*DFxDn;
    Vx2 = rho.*epsixUnif.*DFxDDn;
    Vx3 = rho.*epsixUnif.*DFxDtau;
end

function [zeta, DzetaDnup, DzetaDndn] = polarizeVarZeta(rho)
    zeta = (rho(:, 2) - rho(:, 3)) ./ rho(:, 1);
    DzetaDnup = 2*rho(:, 3) ./ (rho(:, 1).^2);
    DzetaDndn = -2*rho(:, 2) ./ (rho(:, 1).^2);
end

function [ec, Vc1, Vc2, Vc3] = correlationSCAN(rho, rhoup, rhodn, s, alpha, zeta, DsDn, DsDDn, DalphaDn, DalphaDDn, DalphaDtau, DzetaDnup, DzetaDndn)
% all input variables of this function only have one column. 
% ec has only one column; Vc1 has two columns; Vc2 and Vc3 have only one column.
    phi = ((1 + zeta).^(2/3) + (1 - zeta).^(2/3)) / 2;
    rs = (0.75./(pi*rho)).^(1/3) ;
    %% compute epsilon_c^0 when alpha \approx 0
    b1c = 0.0285764;
    b2c = 0.0889;
    b3c = 0.125541;
    ecLDA0 = -b1c ./ (1 + b2c*(rs.^0.5) + b3c*rs);
    dx = 0.5*((1 + zeta).^(4/3) + (1 - zeta).^(4/3)); % no spin, it should be 1
%     cx = -3/(4*pi)*(9*pi/4)^(1/3)*dx;
    cx0 = -3/(4*pi)*(9*pi/4)^(1/3);
    Gc = (1 - 2.3631*(dx - 1)) .* (1 - zeta^12);
    w0 = exp(-ecLDA0/b1c) - 1;
    betaConst = 0.06672455060314922;
%     betaRsInf = 0.066725*0.1/0.1778;
    betaRsInf = betaConst*0.1/0.1778;
    f0 = -0.9;
%     xiInf = (3*pi^2/16)^(2/3) * (betaRsInf*phi/(cx - f0)); % should be 0.128026 if zeta = 0
    xiInf0 = (3*pi^2/16)^(2/3) * (betaRsInf*1/(cx0 - f0));
%     gInf = (1 + 4*xiInf*s.^2).^(-0.25);
    gInf0s = (1 + 4*xiInf0*s.^2).^(-0.25);
    H0 = b1c*log(1 + w0.*(1 - gInf0s));
    ec0 = (ecLDA0 + H0).*Gc;
    %% compute epsilon_c^1 when alpha \approx 1
    sqr_rs = rs .^ 0.5;
    rsm1_2 = 1 ./ sqr_rs;
%     beta = 0.066725* (1 + 0.1*rs) ./ (1 + 0.1778*rs);
    beta = betaConst* (1 + 0.1*rs) ./ (1 + 0.1778*rs);
    % epsilon_c LSDA1
    % correlation parameters
	p = 1 ;
% 	AA = 0.031091 ;
	AA = 0.0310907 ;
	alpha1 = 0.21370 ;
	beta1 = 7.5957 ;
	beta2 = 3.5876 ;
	beta3 = 1.6382 ;
	beta4 = 0.49294 ;
    % conpute
    ec0_q0 = -2.0 * AA * (1.0 + alpha1*rs);
	ec0_q1 = 2.0 * AA *(beta1*sqr_rs + beta2*rs + beta3*rs.*sqr_rs + beta4*rs.*rs);
	ec0_q1p = AA * (beta1*rsm1_2 + 2.0*beta2 + 3.0*beta3*sqr_rs + 4.0*beta4*rs);
	ec0_den = 1.0 ./ (ec0_q1.*ec0_q1 + ec0_q1);
	ec0_log = -log(ec0_q1.*ec0_q1 .* ec0_den);
	ecrs0 = ec0_q0 .* ec0_log;
    ec_lsda1 = ecrs0;
% 	Dec_lsda1Drs = -2.0 * AA * alpha1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;
    
    r = 0.031091;
	w1 = exp(-ec_lsda1 ./ (r*phi.^3)) - 1;
    A = beta ./ (r*w1);
    t = (3*pi^2/16)^(1/3) * s./(phi*sqr_rs);
    g = (1 + 4*A.*t.*t).^(-0.25);
    H1 = r*phi.^3*log(1 + w1.*(1 - g));
    ec1 = ec_lsda1 + H1;
    %% interpolate and extrapolate epsilon_c
    c1c = 0.64;
    c2c = 1.5;
    dc = 0.7;
    alphaG1 = (alpha > 1);
    alphaE1 = (alpha == 1);
    alphaL1 = (alpha < 1);
    fc = zeros(size(rho, 1), size(rho, 2));
    fc(alphaG1) = -dc*exp(c2c ./ (1 - alpha(alphaG1)));
    fc(alphaE1) = 0.0;
    fc(alphaL1) = exp(-c1c*alpha(alphaL1) ./ (1 - alpha(alphaL1)));
    ec = ec1 + fc.*(ec0 - ec1);
    %% compute variation of epsilon_c^0
%     DzetaDn = 0; % no spin
    DrsDn = -4*pi/9*(4*pi/3*rho).^(-4/3);
    DdxDnup = (4/3*(1 + zeta).^(1/3) - 4/3*(1 - zeta).^(1/3)).*DzetaDnup; 
    DdxDndn = (4/3*(1 + zeta).^(1/3) - 4/3*(1 - zeta).^(1/3)).*DzetaDndn; 
    DGcDnup = -2.3631*DdxDnup.*(1 - zeta.^12) + (1 - 2.3631*(dx - 1)).*(12*zeta.^11.*DzetaDnup);
    DGcDndn = -2.3631*DdxDndn.*(1 - zeta.^12) + (1 - 2.3631*(dx - 1)).*(12*zeta.^11.*DzetaDndn);
    DgInf0sDs = -0.25*(1 + 4*xiInf0*s.*s).^(-1.25) .* (4*xiInf0*2*s);
    DgInf0sDn = DgInf0sDs .* DsDn;
    DgInf0sDDn = DgInf0sDs .* DsDDn;
    DecLDA0Dn = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2).*DrsDn;
    Dw0Dn = (w0 + 1) .* (-DecLDA0Dn/b1c);
    DH0Dn = b1c*(Dw0Dn.*(1 - gInf0s) - w0.*DgInf0sDn) ./ (1 + w0.*(1 - gInf0s));
    DH0DDn = b1c* (-w0.*DgInf0sDDn) ./ (1 + w0.*(1 - gInf0s));
%     Dec0Dn = (DecLDA0Dn + DH0Dn).*Gc + (ec0 + H0).*DGcDn;
    Dec0Dnup = (DecLDA0Dn + DH0Dn).*Gc + (ec0 + H0).*DGcDnup;
    Dec0Dndn = (DecLDA0Dn + DH0Dn).*Gc + (ec0 + H0).*DGcDndn;
    Dec0DDn = DH0DDn.*Gc;
    %% compute variation of epsilon_c^1
    Dec_lsda1Dn = - (rs./rho/3).*(-2*AA*alpha1*log(1+1./(2*AA*( beta1*sqr_rs + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1))))) ...
		- ((-2*AA*(1+alpha1*rs)).*(AA*( beta1*(rs.^-0.5)+ 2*beta2 + 3*beta3*(rs.^0.5) + 2*(p+1)*beta4*(rs.^p) ))) ...
		./((2*AA*( beta1*sqr_rs + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) ) ...
		.*(2*AA*( beta1*sqr_rs + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) )+(2*AA*( beta1*(rs.^0.5) + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) )) ) ;% from LDA_PW(S)
    DbetaDn = 0.066725*(0.1*(1+0.1778*rs) - 0.1778*(1+0.1*rs)) ./ ((1+0.1778*rs).^2) .* DrsDn;
%     DphiDn = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDn; % no spin, should be 0
    DphiDnup = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDnup;
    DphiDndn = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDndn;
%     DtDn = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDn.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDnup = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDnup.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDndn = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDndn.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDDn = t.*DsDDn./s;
%     Dw1Dn = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dn - r*ec_lsda1.*(3.*phi.^2.*DphiDn)) ./ ((r*phi.^3).^2));
    Dw1Dnup = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dn - r*ec_lsda1.*(3.*phi.^2.*DphiDnup)) ./ ((r*phi.^3).^2));
    Dw1Dndn = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dn - r*ec_lsda1.*(3.*phi.^2.*DphiDndn)) ./ ((r*phi.^3).^2));
%     DADn = (w1.*DbetaDn - beta.*Dw1Dn) ./ (r*w1.^2);
    DADnup = (w1.*DbetaDn - beta.*Dw1Dnup) ./ (r*w1.^2);
    DADndn = (w1.*DbetaDn - beta.*Dw1Dndn) ./ (r*w1.^2);
%     DgDn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADn.*t.*t + 2*A.*t.*DtDn));
    DgDnup = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADnup.*t.*t + 2*A.*t.*DtDnup));
    DgDndn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADndn.*t.*t + 2*A.*t.*DtDndn));
    DgDDn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*2*A.*t.*DtDDn);
%     DH1Dn = r*(3*phi.^2.*DphiDn.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dn.*(1 - g) - w1.*DgDn) ./ (1 + w1.*(1 - g)));
    DH1Dnup = r*(3*phi.^2.*DphiDnup.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dnup.*(1 - g) - w1.*DgDnup) ./ (1 + w1.*(1 - g)));
    DH1Dndn = r*(3*phi.^2.*DphiDndn.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dndn.*(1 - g) - w1.*DgDndn) ./ (1 + w1.*(1 - g)));
    DH1DDn = r* (phi.^3.*(-w1.*DgDDn) ./ (1 + w1.*(1 - g)));
%     Dec1Dn = Dec_lsda1Dn + DH1Dn;
    Dec1Dnup = Dec_lsda1Dn + DH1Dnup;
    Dec1Dndn = Dec_lsda1Dn + DH1Dndn;
    Dec1DDn = DH1DDn;
    %% variant of fc and ec
    DfcDalpha = zeros(size(rho, 1), size(rho, 2));
    DfcDalpha(alphaG1) = fc(alphaG1).*(c2c./(1 - alpha(alphaG1)).^2);
    DfcDalpha(alphaE1) = 0.0;
    DfcDalpha(alphaL1) = fc(alphaL1).*(-c1c./(1 - alpha(alphaL1)).^2);
    DfcDn = DfcDalpha.*DalphaDn;
    DfcDDn = DfcDalpha.*DalphaDDn;
    DfcDtau = DfcDalpha.*DalphaDtau;
%     DepsiloncDn = Dec1Dn + fc.*(Dec0Dn -Dec1Dn) + DfcDn.*(ec0 - ec1);
    DepsiloncDnup = Dec1Dnup + fc.*(Dec0Dnup -Dec1Dnup) + DfcDn.*(ec0 - ec1);
    DepsiloncDndn = Dec1Dndn + fc.*(Dec0Dndn -Dec1Dndn) + DfcDn.*(ec0 - ec1);
    DepsiloncDDn = Dec1DDn + fc.*(Dec0DDn -Dec1DDn) + DfcDDn.*(ec0 - ec1);
    DepsiloncDtau = DfcDtau.*(ec0 - ec1);
    Vc1 = [ec + rhoup.*DepsiloncDnup, ec + rhodn.*DepsiloncDndn];
    Vc2 = rho.*DepsiloncDDn;
    Vc3 = rho.*DepsiloncDtau;
end