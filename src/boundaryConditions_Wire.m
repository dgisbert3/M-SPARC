function phi = boundaryConditions_Wire(rho,S,m_cut,n_cut)
    % find boundary conditions for periodic in 1D, dirichlet in 2D

    % Implementation by David Codony, after expansion of the log term, and
    % introducing higher order (modified Bessel of second kind) terms
    
    if nargin<3
        m_cut = 6; % Default value for in-plane terms
    end
    if nargin<4
        n_cut = 6; % Default value for axial terms
    end
    
	% Calculate phi
	phi = zeros(S.Nx+2*S.FDn, S.Ny+2*S.FDn, S.Nz+2*S.FDn);

	cellsize  = [S.L1, S.L2, S.L3];
	gridsizes = [S.Nx, S.Ny, S.Nz];
	meshsizes = [S.dx, S.dy, S.dz];
	bcs       = [S.BCx,S.BCy,S.BCz];
    ef_type   = [S.EF_TYPEx,S.EF_TYPEy,S.EF_TYPEz];
    ef_value  = [S.EFx,S.EFy,S.EFz];
    
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
	[~, reverse_order] = sort(new_order); % reverse order to get back
    rho = permute(rho,new_order);
	phi = permute(phi,new_order);
    % rho = permute(rho,reverse_order); % this will recover original

	% sum over Z direction, \int (\rho) dz / LZ
	rho_XY_av = sum(rho,3) * (dZ/LZ); % NX * NY
	rho_XY_av = rho_XY_av(:);
	X = (0 : NX-1) * dX - LX/2 - dX/2; % (Almost) centered
	Y = (0 : NY-1) * dY - LY/2 - dY/2; % (Almost) centered
	Z = (0 : NZ-1) * dZ - LZ/2 - dZ/2; % (Almost) centered
    
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

	XX_bc = (II_ex(isbc)-1-S.FDn) * dX - LX/2 - dX/2; % (Almost) centered
	YY_bc = (JJ_ex(isbc)-1-S.FDn) * dY - LY/2 - dY/2; % (Almost) centered
    RR_bc = sqrt(XX_bc.^2+YY_bc.^2);
    
    diff_RR = sqrt(colminusrow(XX_bc,XX(:)').^2 + colminusrow(YY_bc,YY(:)').^2);
    
    % First contributions: In-plane
    %
    % 4*pi*sum_{m = 1:m_cut}[
    %    1/(m*RR^m) (Cm)* int [ (RR)'^m * (Cm)' * (rho)' dX'dY' ]
    %  + 1/(m*RR^m) (Sm)* int [ (RR)'^m * (Sm)' * (rho)' dX'dY' ]
    %                       ]

    V_XY = 0;
    for m = 1:m_cut
        RRm = RR.^m;
        [Cm ,Sm ] = cylindricalHarmonics_inPlane(XX   ,YY   ,m); % Should be computed in initialization
        [Cmp,Smp] = cylindricalHarmonics_inPlane(XX_bc,YY_bc,m); % Should be computed in initialization
        int_Cm = sum( RRm(:) .*Cm(:) .* rho_XY_av ) * dX * dY;
        int_Sm = sum( RRm(:) .*Sm(:) .* rho_XY_av ) * dX * dY;
        V_XY = V_XY + 4*pi* (Cmp*int_Cm + Smp*int_Sm) ./ (m*RR_bc.^m);
    end
    
    % Second contributions: Axial
    %
    % 4*pi*sum_{n = 1:n_cut}[
    %    (BCn)* int [ bessel(x-x') * (BCn)' * (rho)' dX'dY' ]
    %  + (BSn)* int [ bessel(x-x') * (BSn)' * (rho)' dX'dY' ]
    %                    
    
    V_XYZ = 0;
    for n = 1:n_cut
        [BCn ,BSn ] = cylindricalHarmonics_Axial(Z   ,LZ,n); % Should be computed in initialization
        Besseln = besselk(0,2*pi*n*diff_RR/LZ); %(xy)(x'y')
        BesselxRho = bsxfun(@times,Besseln,permute(reshape(rho,[],NZ),[3,1,2])); %(xy)(x'y')(z')
        
        int_BCn = sum(sum(bsxfun(@times,BesselxRho,permute(BCn,[2,3,1])),3),2) * dX * dY * dZ; %(xy)
        int_BSn = sum(sum(bsxfun(@times,BesselxRho,permute(BSn,[2,3,1])),3),2) * dX * dY * dZ; %(xy)

        V_XYZ = V_XYZ + 4*pi*bsxfun(@times,int_BCn,BCn') + 4*pi*bsxfun(@times,int_BSn,BSn'); %(xy)(z)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%   ELECTRIC  FIELD   %%%%%%%%%%%%%%%%%%%%%%%%% 
    E_X = 0;
    E_Y = 0;
    if ef_type(dir_X)==1                          % Nearby capacitor
        % Compute voltage drop at both sides
        X_drop = [0 , NX-1] * dX - LX/2 - dX/2; % (Almost) centered
        Y_drop = (0 : NY-1) * dY - LY/2 - dY/2; % (Almost) centered

        [XX_drop,YY_drop] = ndgrid(X_drop,Y_drop);
        RR_bc = sqrt(XX_drop.^2+YY_drop.^2);
    
        [Cm ,Sm ] = cylindricalHarmonics_inPlane(XX     ,YY     ,1); % Should be computed in initialization
        [Cmp,Smp] = cylindricalHarmonics_inPlane(XX_drop,YY_drop,1);
        int_Cm = sum( RR(:) .*Cm(:) .* rho_XY_av ) * dX * dY;
        int_Sm = sum( RR(:) .*Sm(:) .* rho_XY_av ) * dX * dY;
        V_XX = (Cmp*int_Cm + Smp*int_Sm) ./ RR_bc;
        
        Vdrop_X = 4*pi*mean(V_XX,2);

        % Revert macroscopic electric field Z that naturally arises
        E_X = E_X + diff(Vdrop_X)/LX;             
    end
    
    if ef_type(dir_Y)==1                          % Nearby capacitor
        % Compute voltage drop at both sides
        X_drop = (0 : NX-1) * dX - LX/2 - dX/2; % (Almost) centered
        Y_drop = [0 , NY-1] * dY - LY/2 - dY/2; % (Almost) centered

        [XX_drop,YY_drop] = ndgrid(X_drop,Y_drop);
        RR_bc = sqrt(XX_drop.^2+YY_drop.^2);
    
        [Cm ,Sm ] = cylindricalHarmonics_inPlane(XX     ,YY     ,1); % Should be computed in initialization
        [Cmp,Smp] = cylindricalHarmonics_inPlane(XX_drop,YY_drop,1);
        int_Cm = sum( RR(:) .*Cm(:) .* rho_XY_av ) * dX * dY;
        int_Sm = sum( RR(:) .*Sm(:) .* rho_XY_av ) * dX * dY;
        V_YY = (Cmp*int_Cm + Smp*int_Sm) ./ RR_bc;
        
        Vdrop_Y = 4*pi*mean(V_YY,1);

        % Revert macroscopic electric field Z that naturally arises
        E_Y = E_Y + diff(Vdrop_Y)/LY;             
    end   
    
	E_X = E_X + ef_value(dir_X);                  % Add user-defined field
	E_Y = E_Y + ef_value(dir_Y);                  % Add user-defined field
    
    % Add linear fields
    V_XY = V_XY - (XX_bc(:)-mean(XX_bc)) * E_X;
    V_XY = V_XY - (YY_bc(:)-mean(YY_bc)) * E_Y;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    V_XY_full       = zeros(NX+2*S.FDn,NY+2*S.FDn);
	V_XY_full(isbc) = V_XY;

	phi = bsxfun(@plus, phi, real(V_XY_full));
    
    V_XYZ_full         = zeros((NX+2*S.FDn)*(NY+2*S.FDn),NZ);
	V_XYZ_full(isbc,:) = V_XYZ;
    
    % Repeat periodic points in Z
    V_XYZ_full = cat(2,V_XYZ_full(:,end-S.FDn+1:end),V_XYZ_full,V_XYZ_full(:,1:S.FDn));

	phi = phi + reshape(V_XYZ_full,[NX+2*S.FDn,NY+2*S.FDn,NZ+2*S.FDn]);
    
	% permute phi back to original directions
	phi = permute(real(phi),reverse_order);
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