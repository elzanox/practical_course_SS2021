clear all   % clear all variables that existed before
close all   % close all figures
clc         % clear the command window

L   = 0.3;   % [m] length of the rod
Acs = 0.01;  % [m2] cross sectional area of the rod
k   = 237;   % [W/m/K] thermal conductivity of the rod
T0  = 293;   % [K] reference temperature (not used)
h   = 3;     % [W/K] cooling coefficient at the right boundary
Q0  = 100;   % [W/m] intensity of the applied heat load 
a   = 0.1;   % [m] width parameter for the applied heat load
Q   =@(x) Q0*exp(-(x-L/2).^2/a^2); % [W/m] applied heat load

N = 11;
x = linspace(0,L,N).';      % spatial grid
dx = x(2) - x(1);            % grid spacing (is constant)

% make a figure of the applied heat load (just for understanding)
fig = figure(1);
plot(x, Q(x))
xlabel 'x [m]'
ylabel 'Q(x) [W/m]'

%% TODO: Construct the stiffness matrix A here
% You can choose an explicit or an implicit implementation of the boundary conditions,
% but note that the stiffness matrix will become larger in the implicit implementation
A = sparse(N+2,N+2);    % creates a NxN sparse matrix (so all elements are zero)
A(1,2) = 1;

for ii = 2:N+1
  A(ii,ii-1) = 1/dx^2;
  A(ii,ii) = -2/dx^2;
  A(ii,ii+1) = 1/dx^2;
end

A(end, end) = 1;
A(end, end-2) = -1;

%% TODO: Construct the load vector f here
% note that the size of b should match the size of the stiffness matrix A
f = zeros(N,1);     % creates a zero column vector of length N

%% TODO: Compute the resulting temperature field Tr
% note that Tr should be a vector of length N

% save a figure of the resulting temperature field
fig = figure(2);
plot(x,T)
xlabel 'x [m]'
ylabel 'T(x) [K]'
print('Week1_T.jpg')