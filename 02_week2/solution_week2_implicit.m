clear all
close all
clc

Lx = 0.5;     % [m] length of the plate in x-direction
Ly = 0.75;    % [m] length of the plate in y-direction
h  = 5;       % [W/m/K] thermal conductance of the boundary
k  = 400;     % [W/m/K] thermal conductivity of copper
H  = 0.01;    % [m] thickness of the plate
a  = 0.05;    % [m] with parameter of the heat load
Q0 = 100;     % [W/m2] applied heat load 
x0 = 0.3*Lx;  % [m] x-coordinate of the center of the heat load
y0 = 0.3*Ly;  % [m] y-coordinate of the center of the heat load
Q  = @(x,y) Q0*exp(-((x-x0).^2+(y-y0).^2)/a^2);  % [W/m2] applied heat load

Nx = 32;      % Number of grid points in the x- and y-direction
Ny = 42;
x  = linspace(0,Lx,Nx);   % grids in the x- and y-direction
y  = linspace(0,Ly,Ny);
dx = x(2) - x(1);         % grid spacings in the x- and y-direction
dy = y(2) - y(1);

% Step I: create the matrix with the node numbers
% Initialize the node numbers matrix
node_nmbrs = zeros(Nx+2, Ny+2);

% Step I.A internal nodes
for ii = 1:Nx
    for jj = 1:Ny
        node_nmbrs(ii+1,jj+1) = (ii-1)*(Ny) + jj;
    end
end

% Total number of nodes in the implicit formulation
nn = (Nx + 2) * (Ny + 2);

A = zeros(nn, nn); % Build the stiffness matrix

% Loop through internal nodes to fill A
for ii = 1:Nx
    for jj = 1:Ny
        node11 = node_nmbrs(ii+1,jj+1);
        node21 = node_nmbrs(ii+2,jj+1);
        node12 = node_nmbrs(ii+1,jj+2);
        node01 = node_nmbrs(ii  ,jj+1);
        node10 = node_nmbrs(ii+1,jj  );
        
        % Implement equations for the internal nodes
        A(node11, node11) = -2*(1/dx^2 + 1/dy^2);
        A(node11, node21) = 1/dx^2;
        A(node11, node12) = 1/dy^2;
        A(node11, node01) = 1/dx^2;
        A(node11, node10) = 1/dy^2;
    end
end

% Boundary conditions at x = 0
for jj = 1:Ny
    node11 = node_nmbrs(1,jj+1);
    
    % Implement boundary conditions
    A(node11, node11) = 1;
end

% Boundary conditions at x = Lx
for jj = 1:Ny
    node11 = node_nmbrs(Nx+2,jj+1); 
    
    % Implement boundary conditions
    A(node11, node11) = 1;
end

% Boundary conditions at y = 0
for ii = 1:Nx
    node11 = node_nmbrs(ii+1,1);
    
    % Implement boundary conditions
    A(node11, node11) = 1;
end

% Boundary conditions at y = Ly
for ii = 1:Nx
    node11 = node_nmbrs(ii+1,Ny+2);
    
    % Implement boundary conditions
    A(node11, node11) = 1;
end

f = zeros(nn,1); % Forcing vector

% Loop through internal nodes to define the forcing vector
for ii = 1:Nx
    for jj = 1:Ny
        node11 = node_nmbrs(ii+1, jj+1);
        
        % Define the forcing vector
        f(node11) = Q(x(ii), y(jj));
    end
end

% Solve the constructed system of equations
T = A\f;

nodes_int = node_nmbrs(2:end-1, 2:end-1);  % Only the values at the internal nodes will be plotted

fig1 = figure(1);
surf(x, y, reshape(f(nodes_int), Ny, Nx).')
daspect([1, 1, max(f)])
axis tight
xlabel 'x [m]'
ylabel 'y [m]'
title 'Applied heat load [W/m^2]'
colorbar
saveas(fig1, 'Exercise1_heatload.jpg')

fig2 = figure(2);
surf(x, y, reshape(T(nodes_int), Ny, Nx).')
daspect([1, 1, 0.2])
axis tight
xlabel 'x [m]'
ylabel 'y [m]'
title 'Temperature field [K]'
colorbar
saveas(fig2, 'Exercise1_temperature.jpg')
