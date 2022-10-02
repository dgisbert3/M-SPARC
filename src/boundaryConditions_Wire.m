function phi = boundaryConditions_Wire(rho,S)
    % find boundary conditions for periodic in 1D, dirichlet in 2D

    % Implementation by David Codony, after expansion of the log term, and
    % introducing higher order (modified Bessel of second kind) terms
    
	% Calculate phi
	phi = zeros(S.Nx+2*S.FDn, S.Ny+2*S.FDn, S.Nz+2*S.FDn);

	cellsize  = [S.L1, S.L2, S.L3];
	gridsizes = [S.Nx, S.Ny, S.Nz];
	meshsizes = [S.dx, S.dy, S.dz];
	bcs       = [S.BCx,S.BCy,S.BCz];

	% find which direction has Periodic BC
	dir_Z = find(bcs == 0); 
	dir_X = mod(dir_Z, 3) + 1;
	dir_Y = mod(dir_X, 3) + 1;

	% once we find the direction, we assume that direction is the Z
	% direction, the other two directions are then called X, Y
	NX = gridsizes(dir_X); NY = gridsizes(dir_Y); NZ = gridsizes(dir_Z);
	LX = cellsize(dir_X);  LY = cellsize(dir_Y);  LZ = cellsize(dir_Z);
	dX = meshsizes(dir_X); dY = meshsizes(dir_Y); dZ = meshsizes(dir_Z);
	% A_XY = LX * LY; % area of (x',y') surface

	% reshape rho to 3D 
	rho = reshape(rho+S.b, S.Nx, S.Ny, S.Nz); % note here after rho = rho + b

	% permute rho so that the new Z direction has Periodic BC
	new_order = [dir_X, dir_Y, dir_Z]; % a permutation of [1,2,3]
	rho = permute(rho,new_order);
	phi = permute(phi,new_order);

	% sum over Z direction, \int (\rho) dz / LZ
	rho_XY_av = sum(rho,3) * (dZ/LZ); % NX * NY
	rho_XY_av = rho_XY_av(:);
	X = (0 : NX-1) * dX - L1/2; % Centered
	Y = (0 : NY-1) * dY - L2/2; % Centered
	Z = (0 : NZ-1) * dZ - L3/2; % Centered

	[XX,YY] = ndgrid(X,Y);
    RR = sqrt(XX.^2+YY.^2);

	I_ex = 1:NX+2*S.FDn;
	J_ex = 1:NY+2*S.FDn;
	[II_ex,JJ_ex] = ndgrid(I_ex,J_ex);

	% flag for points inside the domain
	isIn = zeros(NX+2*S.FDn,NY+2*S.FDn);
	isIn(S.FDn+1:NX+S.FDn, S.FDn+1:NY+S.FDn) = 1;

	% positions of nodes outside the domain
	isbc = find(~isIn);

	XX_bc = (II_ex(isbc)-1-S.FDn) * dX - L1/2; % Centered
	YY_bc = (JJ_ex(isbc)-1-S.FDn) * dY - L2/2; % Centered

    % First contributions: In-plane
    %
    % 4*pi*sum_{m = 1:m_cut}[
    %    1/(m*RR^m) (Cm)* int [ (RR)'^m * (Cm)' * (rho)' dX'dY' ]
    %  + 1/(m*RR^m) (Sm)* int [ (RR)'^m * (Sm)' * (rho)' dX'dY' ]
    %                       ]
    m_cut = 6; % HARDCODED

    V_XY = 0;
    for m = 1:m_cut
        RRm = RR.^m;
        [Cm,Sm] = cylindricalHarmonics_inPlane(XX,YY,LZ,m); % Should be computed in initialization
        int_Cm = sum(sum( repmat( RRm .*Cm ,1,1,NZ ) .* rho )) * dX * dY;
        int_Sm = sum(sum( repmat( RRm .*Sm ,1,1,NZ ) .* rho )) * dX * dY;
        V_XY = V_XY + (Cm*int_Cm + Sm*int_Sm) ./ (m*RRm);
    end
    
    % Second contributions: Axial
    %
    % 4*pi*sum_{n = 1:n_cut}[
    %    (BCn)* int [ bessel(x-x') * (BCn)' * (rho)' dX'dY' ]
    %  + (BSn)* int [ bessel(x-x') * (BSn)' * (rho)' dX'dY' ]
    %                       ]
    n_cut = 6; % HARDCODED

    diff_RR = sqrt(colminusrow(XX_bc,XX(:)').^2 + colminusrow(YY_bc,YY(:)').^2);
    for n = 1:n_cut
        [BCn,BSn] = cylindricalHarmonics_Axial(Z,LZ,n); % Should be computed in initialization
        Besseln = besselk(0,2*pi*n*diff_RR/LZ);
        % HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE
        int_Cm = sum( repmat( RRm .*Cm ,1,1,NZ ) .* rho ,2) * dX * dY;
        int_Sm = sum( repmat( RRm .*Sm ,1,1,NZ ) .* rho ,2) * dX * dY;
        V_XY = V_XY + (Cm*int_Cm + Sm*int_Sm) ./ (m*RRm);
    end

	V_XY_full = zeros(NX+2*S.FDn,NY+2*S.FDn);
	V_XY_full(isbc) = V_XY_av;

	phi = bsxfun(@plus, phi, V_XY_full);

	% permute phi back to original directions
	phi(neworder) = real(phi);
end


% tool function colminusrow
function xmy = colminusrow(x,y)
% A column vector x minus a row vector.
% In Matlab versions after R2018b, it's just x - y
if (size(x,2) ~= 1)
	error('ERROR: the first vector must be a column vector');
end
if (size(y,1) ~= 1)
	error('ERROR: the second vector must be a row vector');
end

[xx,yy] = ndgrid(x,y);
xmy = xx - yy;

end