function Grad = gradOperator_extended(S,dir)
% gradOperator_extended     Computes the gradient operator in direction dir. 
% The result is a rectangular sparse matrix, which has as many rows as grid
% points, and as many columns as grid points including ghost points along 
% every direction. Multiplying this
% operator by an extended version of phi, where ghost points are not
% assumed to be zero, gives the correct gradient.
%
% by David Codony

assert( S.cell_typ == 1 , 'gradOperator_Ghosts(S,dir): Implemented only for S.cell_typ == 1');
Nx = S.Nx;
Ny = S.Ny;
Nz = S.Nz;
N = S.N;
n0 = S.FDn;
w1 = S.w1;
dx = S.dx;
dy = S.dy;
dz = S.dz;

% Initial number of non-zeros: including ghost nodes
nnz_count = 2*n0*N;

% Row numbers and non-zero values
I = zeros(nnz_count,1) ;
V = zeros(nnz_count,1) ;

% Indices of the columns
II = zeros(nnz_count,1);
JJ = zeros(nnz_count,1) ;
KK = zeros(nnz_count,1) ;

switch dir
    case 1
        % Gradient along x_direction
        row_count = 1;
        count = 1 ;
        for kk=1:Nz
            for jj=1:Ny
                for ii=1:Nx
                    % off-diagonal elements
                    for p=1:n0
                        % ii+p
                        I(count) = row_count; II(count) = ii+p; JJ(count) = jj; KK(count) = kk;
                        V(count) = w1(p+1)/dx;
                        count = count + 1;
                        % ii-p
                        I(count) = row_count; II(count) = ii-p ; JJ(count) = jj; KK(count) = kk;
                        V(count) = -w1(p+1)/dx;
                        count = count + 1;
                    end
                    row_count = row_count+1;
                end
            end
        end
    case 2
        % Gradient along y_direction
        row_count = 1;
        count = 1 ;
        for kk=1:Nz
            for jj=1:Ny
                for ii=1:Nx
                    % off-diagonal elements
                    for p=1:n0
                        % ii+p
                        I(count) = row_count; II(count) = ii; JJ(count) = jj+p; KK(count) = kk;
                        V(count) = w1(p+1)/dy;
                        count = count + 1;
                        % ii-p
                        I(count) = row_count; II(count) = ii ; JJ(count) = jj-p; KK(count) = kk;
                        V(count) = -w1(p+1)/dy;
                        count = count + 1;
                    end
                    row_count = row_count+1;
                end
            end
        end
    case 3
        % Gradient along z_direction
        row_count = 1;
        count = 1 ;
        for kk=1:Nz
            for jj=1:Ny
                for ii=1:Nx
                    % off-diagonal elements
                    for p=1:n0
                        % ii+p
                        I(count) = row_count; II(count) = ii; JJ(count) = jj; KK(count) = kk+p;
                        V(count) = w1(p+1)/dz;
                        count = count + 1;
                        % ii-p
                        I(count) = row_count; II(count) = ii ; JJ(count) = jj; KK(count) = kk-p;
                        V(count) = -w1(p+1)/dz;
                        count = count + 1;
                    end
                    row_count = row_count+1;
                end
            end
        end
end

% Map periodic indices 
if S.BCx == 0
    II = mod(II+(Nx-1),Nx)+1;
end
if S.BCy == 0
    JJ = mod(JJ+(Ny-1),Ny)+1;
end
if S.BCz == 0
    KK = mod(KK+(Nz-1),Nz)+1;
end

% Extend indices in every dimension
II = II + n0;
JJ = JJ + n0;
KK = KK + n0;

Nx = Nx + 2*n0;
Ny = Ny + 2*n0;
Nz = Nz + 2*n0;
    
% Getting linear indices of the columns
J = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;
        
Grad = sparse(I,J,V,N,Nx*Ny*Nz);
end


