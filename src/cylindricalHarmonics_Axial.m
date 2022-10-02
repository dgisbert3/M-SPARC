function [BCn,BSn] = cylindricalHarmonics_Axial(Z,Lz,n)
% @brief    cylindricalHarmonics(X,Y,Z,m,n,kind) calculates the cylindrical harmonics
%           BCn,BSn (real) at the given positions (X,Y,Z)
%
% @param Z      The z (axial) coordinates of the grid 
% @param Lz     The z (axial) length of the periodic domain 
% @param n      The index for the sum in axial direction (modified Bessel function of second kind B)
%
% @authors  David Codony <dgisbert3@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2022 Material Physics & Mechanics Group, Georgia Tech
%============================================================================================

if (Lz <= 0)
	error('cylindricalHarmonics_Axial(Z,Lz,n): <Lz> must be positive');
end

isInt = n - round(l);
if (isInt ~= 0)||(n < 0)
	error('cylindricalHarmonics_Axial(Z,Lz,n): <n> must be a positive integer');
end


theta = 2*pi*n*Z/Lz ;

c = 0.564189583547756/Lz^0.5; % 1/sqrt(pi*Lz)

BCn = c * cos(theta);
BSn = c * sin(theta);
end