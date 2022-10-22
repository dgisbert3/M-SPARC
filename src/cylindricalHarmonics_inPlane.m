function [Cm,Sm] = cylindricalHarmonics_inPlane(X,Y,m)
% @brief    cylindricalHarmonics_inPlane(X,Y,m) calculates the cylindrical harmonics
%           Cm,Sm (real) at the given positions (X,Y) in plane
%
% @param X      The x coordinates of the positions in the x-y plane
% @param Y      The y coordinates of the positions in the x-y plane
% @param m      The index for the sum in angular direction (cos and sin)
%
% @authors  David Codony <dgisbert3@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2022 Material Physics & Mechanics Group, Georgia Tech
%============================================================================================

isInt = m - round(m);
if (isInt ~= 0)||(m < 0)
	error('cylindricalHarmonics_inPlane(X,Y,Lz,m): <m> must be a positive integer');
end

r = sqrt(X.^2 + Y.^2);

cosXY = X./r;
sinXY = Y./r;

c  = 0.398942280401433; % 1/sqrt(2*pi)

switch m
    case 1
        Cm = c*cosXY;                                                          % c*cos(theta);
        Sm = c*sinXY;                                                          % c*sin(theta);
    case 2
        Cm = c*(2*cosXY.^2-1);                                                 % c*cos(2*theta);
        Sm = c*(2.*cosXY.*sinXY);                                              % c*sin(2*theta);
    case 3
        Cm = c*(4*cosXY.^2-3).*cosXY;                                          % c*cos(3*theta);
        Sm = c*(3-4*sinXY.^2).*sinXY;                                          % c*sin(3*theta);
    case 4
        Cm = c*(1-8*cosXY.^2+8*cosXY.^4);                                      % c*cos(4*theta);
        Sm = c*cosXY.*(4*sinXY-8.*sinXY.^3);                                   % c*sin(4*theta);
    case 5
        cs = cosXY.^2.*sinXY.^2;
        Cm = c*cosXY.*(cosXY.^4-10*cs+5*sinXY.^4);                             % c*cos(5*theta);
        Sm = c*sinXY.*(sinXY.^4-10*cs+5*cosXY.^4);                             % c*sin(5*theta);
    case 6
        Cm = c*(-1+18.*cosXY.^2-48.*cosXY.^4+32*cosXY.^6);                     % c*cos(6*theta);
        Sm = c*cosXY.*sinXY.*(6*cosXY.^4-20.*cosXY.^2.*sinXY.^2+6.*sinXY.^4);  % c*sin(6*theta);
    otherwise
        theta = m*atan2(Y,X);
        Cm = c*cos(theta);
        Sm = c*sin(theta);
end
end