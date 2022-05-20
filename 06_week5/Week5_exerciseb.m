clear all
close all
clc

L   = 0.3;   % [m] length of the rod
Acs = 0.01;  % [m2] cross sectional area of the rod
k   = 237;   % [W/m/K] thermal conductivity of the rod
T0  = 293;   % [K] reference temperature (not used)
h   = 3;     % [W/K] cooling coefficient at the right boundary
Q0  = 100;   % [W/m] intensity of the applied heat load 
a   = 0.1;   % [m] width parameter for the applied heat load
Q   =@(x) Q0*exp(-(x-L/2).^2/a^2); % [W/m] applied heat load

M = 10;    % number of elements
N = TODO;   % number of nodes

Le = TODO;  % length of one element
x = TODO;  % spatial grid with the node positions

Ae = TODO;      % Element stiffness matrix
fe = TODO;        % Element load vector

%% Assemble the matrices
A = sparse(N,N);        % Initialize zero matrices
f = zeros(N,1);

for e = 1:M % you need a loop over all elements
    TODO
end
% dont forget the Robin BC...

% Apply the dirichlet BCs

% Compute the steady state temperature field T

% plot the result
plot(x,T)
xlabel 'x [m]'
ylabel 'T [K]'