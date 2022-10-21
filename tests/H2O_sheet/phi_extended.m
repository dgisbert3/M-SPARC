function phi = phi_extended(S)
II = (1+S.FDn):(S.Nx+S.FDn);
JJ = (1+S.FDn):(S.Ny+S.FDn);
KK = (1+S.FDn):(S.Nz+S.FDn);

phi           = reshape(S.phi_extended,S.Nx+S.FDn*2,S.Ny+S.FDn*2,S.Nz+S.FDn*2);
phi(II,JJ,KK) = reshape(S.phi,S.Nx,S.Ny,S.Nz);

dphidx = reshape(gradOperator_extended(S,1)*phi(:),S.Nx,S.Ny,S.Nz); % Test
dphidy = reshape(gradOperator_extended(S,1)*phi(:),S.Nx,S.Ny,S.Nz); % Test
dphidz = reshape(gradOperator_extended(S,1)*phi(:),S.Nx,S.Ny,S.Nz); % Test

% %% Copy-paste periodic directions
% %% This is not strictly necessary
if S.BCx == 0
    phi(1:S.FDn         ,:,:) = phi((1:S.FDn         )+S.Nx,:,:);
    phi(S.Nx+S.FDn+1:end,:,:) = phi((S.Nx+S.FDn+1:end)-S.Nx,:,:);
end
if S.BCy == 0
    phi(:,1:S.FDn         ,:) = phi(:,(1:S.FDn         )+S.Ny,:);
    phi(:,S.Ny+S.FDn+1:end,:) = phi(:,(S.Ny+S.FDn+1:end)-S.Ny,:);
end
if S.BCz == 0
    phi(:,:,1:S.FDn         ) = phi(:,:,(1:S.FDn         )+S.Nz);
    phi(:,:,S.Nz+S.FDn+1:end) = phi(:,:,(S.Nz+S.FDn+1:end)-S.Nz);
end

% dphidx2 = reshape(gradOperator_extended(S,1)*phi(:),S.Nx,S.Ny,S.Nz); % Test
% dphidy2 = reshape(gradOperator_extended(S,1)*phi(:),S.Nx,S.Ny,S.Nz); % Test
% dphidz2 = reshape(gradOperator_extended(S,1)*phi(:),S.Nx,S.Ny,S.Nz); % Test
% 
% norm(dphidx(:)-dphidx2(:)) % Test
% norm(dphidy(:)-dphidy2(:)) % Test
% norm(dphidz(:)-dphidz2(:)) % Test

phi           = phi(:);
end