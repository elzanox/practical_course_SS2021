Lx = 0.5;     % [m] length of the plate in x-direction
Ly = 0.4;     % [m] length of the plate in y-direction
k  = 57.7;    % [W/m/K] thermal conductivity of steel
c  = 448;     % [J/kg/K] heat capacity of steel
rho = 7840;   % [kg/m3] mass density of steel
H  = 0.003;   % [m] thickness of the plate

a  = 0.05;     % [m] width parameter of the heat load
x0 = 0.7*Lx;   % [m] x-coordinate of the center of the heat load
y0 = 0.8*Ly;   % [m] y-coordinate of the center of the heat load
B0 =@(x,y) 1/pi/a^2*exp(-((x-x0).^2+(y-y0).^2)/a^2);  % [1/m2] applied heat load

Nx = 81;      % number of grid points in the x-direction
Ny = 65;      % number of grid points in the y-direction
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
dx = x(2) - x(1);   % note that the grid spacings can be slightly smaller than the spec
dy = y(2) - y(1);

nn = Nx*Ny;
node_nmbrs = 1:nn;
node_nmbrs = reshape(node_nmbrs, Nx, Ny);

%% construct A-matrix

A = sparse(nn,nn);

for ii = 2:Nx-1
    for jj = 2:Ny-1
        node11 = node_nmbrs(ii,jj);
        node21 = node_nmbrs(ii+1,jj);
        node12 = node_nmbrs(ii,jj+1);
        node01 = node_nmbrs(ii-1,jj);
        node10 = node_nmbrs(ii,jj-1);
        
        A(node11, node11) = -2/dx^2 -2/dy^2;
        A(node11, node21) =  1/dx^2;
        A(node11, node12) =  1/dy^2;
        A(node11, node01) =  1/dx^2;
        A(node11, node10) =  1/dy^2;
    end
end

for ii = 2:Nx-1
    node11 = node_nmbrs(ii,1);
    node21 = node_nmbrs(ii+1,1);
    node12 = node_nmbrs(ii,2);
    node01 = node_nmbrs(ii-1,1);
    
    A(node11, node11) = -2/dx^2 -2/dy^2;
    A(node11, node21) =  1/dx^2;
    A(node11, node12) =  2/dy^2;
    A(node11, node01) =  1/dx^2;
end

for ii = 2:Nx-1
    node11 = node_nmbrs(ii,end);
    node21 = node_nmbrs(ii+1,end);
    node12 = node_nmbrs(ii,end-1);
    node01 = node_nmbrs(ii-1,end);
    
    A(node11, node11) = -2/dx^2 -2/dy^2;
    A(node11, node21) =  1/dx^2;
    A(node11, node12) =  2/dy^2;
    A(node11, node01) =  1/dx^2;
end

for jj = 2:Ny-1
    node11 = node_nmbrs(1,jj);
    node21 = node_nmbrs(2,jj);
    node12 = node_nmbrs(1,jj+1);
    node10 = node_nmbrs(1,jj-1);
        
    A(node11, node11) = -2/dx^2 -2/dy^2;
    A(node11, node21) =  2/dx^2;
    A(node11, node12) =  1/dy^2;
    A(node11, node10) =  1/dy^2;
end

for jj = 2:Ny-1
    node11 = node_nmbrs(end,jj);
    node21 = node_nmbrs(end-1,jj);
    node12 = node_nmbrs(end,jj+1);
    node10 = node_nmbrs(end,jj-1);
        
    A(node11, node11) = -2/dx^2 -2/dy^2;
    A(node11, node21) =  2/dx^2;
    A(node11, node12) =  1/dy^2;
    A(node11, node10) =  1/dy^2;
end

node11 = node_nmbrs(1,1);
node21 = node_nmbrs(2,1);
node12 = node_nmbrs(1,2);
A(node11, node11) = -2/dx^2 - 2/dy^2;
A(node11, node21) =  2/dx^2;
A(node11, node12) =  2/dy^2;

node11 = node_nmbrs(end,1);
node21 = node_nmbrs(end-1,1);
node12 = node_nmbrs(end,2);
A(node11, node11) = -2/dx^2 - 2/dy^2;
A(node11, node21) =  2/dx^2;
A(node11, node12) =  2/dy^2;

node11 = node_nmbrs(1,end);
node21 = node_nmbrs(2,end);
node12 = node_nmbrs(1,end-1);
A(node11, node11) = -2/dx^2 - 2/dy^2;
A(node11, node21) =  2/dx^2;
A(node11, node12) =  2/dy^2;

node11 = node_nmbrs(end,end);
node21 = node_nmbrs(end-1,end);
node12 = node_nmbrs(end,end-1);
A(node11, node11) = -2/dx^2 - 2/dy^2;
A(node11, node21) =  2/dx^2;
A(node11, node12) =  2/dy^2;

A = k*H*A;

%% construct the forcing vector

B = zeros(nn, 1);
for ii = 1:Nx
    for jj = 1:Ny
        B(node_nmbrs(ii,jj)) = B0(x(ii), y(jj));
    end
end

save('Week3_spatial_discretization.mat', 'A', 'B', 'x', 'y', 'node_nmbrs')