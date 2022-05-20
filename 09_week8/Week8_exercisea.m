clear all
close all
clc

L  = 0.4;     % [m] length of the beam
rho = 2700;   % [kg/m3] mass density
E  = 69e9;    % [Pa] Young's modulus
W  = 0.02;    % [m] width of the beam
H  = 0.01;    % [m] height of the beam

A = H*W;      % [m2] cross sectional area
I = W*H^3/12; % [m4] second moment of area


M = 100;      % [-] number of elements
N = TODO;      % [-] number of nodes
xgrid  = linspace(0,L,N);
Le = L/M;     % length of one element

%% Compute the element matrices
% pkg load symbolic % uncomment this line when you use Octave. 
                    % This will lead to a few warnings that you can ignore. 
syms x
xi  = x/Le; 
% TODO define the element shape functions in terms of xi
Ne = TODO

% TODO: Define the element mass matrix
Me = TODO; 

% TODO: Define the element stiffness matrices Ke
Ke = TODO;

% TODO: Assemble stiffness matrix and load vector
MM = zeros(2*N,2*N);
KK = zeros(2*N,2*N);
for e = 1:M
    TODO
end

% TODO: Determine constrained and free DOFs
cdofs = TODO;
fdofs = setdiff(1:2*N, cdofs);

% TODO: Compute the free-free partitions of M and K
Mff = MM(fdofs, fdofs);
Kff = KK(fdofs, fdofs);

% TODO: Solve the eigenmodes and compute the eigenfrequencies
angular_eigen_frequencies = TODO % [rad/s]

% TODO: Define a matrix with all DOFs of the eigenmodes

% TODO: Plot the eigenmodes
figure(1)
plot(xgrid, TODO)
xlabel 'x [m]'
ylabel 'eigenmode [-]'