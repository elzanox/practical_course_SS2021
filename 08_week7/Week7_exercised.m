clear all
close all
clc

Lx = 0.4;     % [m] length of the plate in x-direction
Ly = 0.5;     % [m] length of the plate in y-direction
H  = 0.001;   % [m] thickness
E  = 69e9;    % [Pa] Young's modulus
nu = 0.3;     % [-] Poisson's ratio
a  = 0.1;     % [m] width parameter
fleft =@(y) 100*exp(-(y-Ly/2).^2/a^2); % [N/m] force per unit length applied at the left edge
fbot  =@(x) 100*exp(-(x-Lx/2).^2/a^2); % [N/m] force per unit length applied at the bottom edge

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

% TODO: define the node_nmbrs matrix
node_nmbrs = TODO;

%% Step 4: Build the element list
elem_list = zeros(Mx*My, 4);
elem_nmbrs = zeros(Mx, My);
e = 0;
for ii = 1:Mx
    for jj = 1:My
        e = e+1;
        % TODO: build the element list
        elem_list(e, :) = TODO;
        % We also build the matrix elem_nmbrs which you will need for the BCs later on
        elem_nmbrs(ii,jj) = e;
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
Ee = double(int(int(Ne.'*Ne, x, [0, Lex]), y, [0,Ley])); 

dxNe = diff(Ne, x); 
dyNe = diff(Ne, y);

% TODO: Define the element stiffness matrices Ae11, Ae12, Ae21, and Ae22
Ae11 = TODO;
Ae12 = TODO;
Ae21 = TODO; 
Ae22 = TODO;

% TODO: Compute the element load vectors for the left and bottom edge
feleft  = TODO;
febot   = TODO;

%% Step 6: Assemble stiffness matrix and load vector
A = zeros(2*nn,2*nn);
for e = 1:ne
    nodes = elem_list(e, :);
    % TODO: assemble the stiffness matrix using the element matrices Ae11, Ae12, Ae21, Ae22.
end

% build the global load vector
f = zeros(2*nn,1);
for jj = 1:My                     % left edge
    e = elem_nmbrs(1,jj);
    nodes = elem_list(e,:);
    % TODO: build load vector for the left edge
end

% TODO: bottom edge

%% Step 7: Include the Robin BCs (not needed in this exercise)

%% Step 8: Include Dirichlet BCs

% TODO: select the nodes on the right edge
nodes_right = TODO;
% TODO: select the numbers of the dofs (so not the nodes!)
cdofs = TODO;
% The remaining DOFs are the free dofs
fdofs = setdiff(1:2*nn, cdofs);

%% Compute and plot the steady state solution

u = zeros(2*nn,1);
u(fdofs) = A(fdofs,fdofs) \ f(fdofs);

% TODO: compute the x and y components of the resulting displacement field
ux = TODO;
uy = TODO;

normu = sqrt(ux.^2 + uy.^2);
[X,Y] = ndgrid(xgrid, ygrid);
scale = 1e5;
surf(X+scale*ux(node_nmbrs),Y+scale*uy(node_nmbrs),normu(node_nmbrs))
view(2);
xlabel 'x [m]'
ylabel 'y [m]'
title '(scaled) displacement [m]'
colorbar
axis equal