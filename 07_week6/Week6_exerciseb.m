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
Q  =@(x,y) Q0*exp(-((x-x0).^2+(y-y0).^2)/a^2);  % [W/m2] applied heat load

%% Divide the domain into elements and define grids for node positions
Mx = 32;      % Number of elements in the x and y directions
My = 42;
ne = Mx*My;   % Total number of elements

%% Step 1: Divide the domain into elements
Lex = Lx / Mx;         % element sizes in the x- and y-direction
Ley = Ly / My;

%% Step 2: Choose the element type
% we choose linear elements

%% Step 3: Assign numbers to all nodes

Nx = Mx+1;    % number of nodes in the x and y directions
Ny = My+1;
nn = Nx*Ny;   % Total number of nodes

xgrid  = linspace(0,Lx,Nx);   % grids in the x- and y-direction
ygrid  = linspace(0,Ly,Ny);

% TODO: Set up the node_nmbrs matrix

%% Step 4: Build the element list
elem_list = zeros(Mx*My, 4);
elem_nmbrs = zeros(Mx, My);
e = 0;
for ii = 1:Mx
    for jj = 1:My
        e = e+1;
        elem_list(e, :) = [node_nmbrs(ii,jj), TODO];
        % TODO: add the correct element numbers (in the right order) to the element list
        elem_nmbrs(ii,jj) = e;
        % We also build the matrix elem_nmbrs which you will need for the
        % Robin BC later on. 
    end
end

%% Step 5: Compute the element matrices
% pkg load symbolic % uncomment this line when you use Octave. 
                    % This will lead to a few warnings that you can ignore. 
syms x y
xi  = x/Lex; 
eta = y/Ley;
Ne = [(1-xi)*(1-eta), xi*(1-eta), xi*eta, (1-xi)*eta];

% This is how you can build the element mass matrix 
% (computations by hand get quite cumbersome in 2D)
Ee = double(int(int(Ne.'*Ne, x, [0, Lex]), y, [0,Ley])); 

dxNe = diff(Ne, x); 
dyNe = diff(Ne, y);

% TODO: Define the element stiffness matrices (for the Laplacian)
Ae   = TODO;  

% TODO: Define the element load vector
fe =  TODO;

% TODO: Compute the element matrices you need to apply the Robin BC on each
% of the four edges of the element
Eebot   = TODO;
Eetop   = TODO;
Eeleft  = TODO;
Eeright = TODO;

%% Step 6: Assemble stiffness matrix and load vector
A = zeros(nn,nn);
for e = 1:ne
    nodes = elem_list(e, :);
    A(nodes,nodes) = A(nodes,nodes) + TODO;
end

% element load vector
f = zeros(nn,1);
for ii = 1:Mx
    for jj = 1:My
        e = elem_nmbrs(ii,jj);
        nodes = elem_list(e,:);
        f(nodes) = f(nodes) + TODO;
    end
end

%% Step 7: Include the Robin BCs

% lower edge
for ii = 1:Mx
    e = elem_nmbrs(ii,1);
    nodes = elem_list(e,:);
    A(nodes,nodes) = A(nodes,nodes) + TODO;
end

% TODO: top edge

% TODO: left edge

% TODO: right edge

%% Step 8: Include Dirichlet BCs (not needed in this exercise)

%% Compute and plot the steady state solution
T = -A \ f;

surf(xgrid,ygrid,T(node_nmbrs).')
daspect([1, 1, 0.1])
axis tight
xlabel 'x [m]'
ylabel 'y [m]'
title 'Temperature field [K]'