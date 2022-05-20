clear all
close all
clc

Lx = 0.5;     % [m] length of the plate in x-direction
Ly = 0.4;     % [m] length of the plate in y-direction
k  = 57.7;    % [W/m/K] thermal conductivity of steel
c  = 448;     % [J/kg/K] heat capacity of steel
rho = 7840;   % [kg/m3] mass density of steel
H  = 0.003;   % [m] thickness of the plate
Tsim = 60;     % [s] duration of the simulation     
P0 = 50;       % [W] heating power
t1 = 10;       % [s] the first heat load is active until t=t1
P =@(t) P0*(tanh(t1-t)+1)/2; % applied heating power [W]

load('Week3_spatial_discretization')

nn = length(B);

%% Compute reference solution

NT = 769;
Tref = crank_nicolson(NT, nn, A, B, P, Tsim);

% TODO: compute the norm of the reference solution (used to compute the reltaive error later)
norm_Tref = 1;

%% Do the convergence study

e_fe = zeros(1,5); % Errors in the Forward Euler scheme
e_cn = zeros(1,5); % Errors in the Crank-Nicolson scheme
e_be = zeros(1,5); % Errors in the Backward Euler scheme
dt   = zeros(1,5); % time increments for the different values of NT
for ii = 1:5
    NT = 12*2^ii + 1;
    
    T = forward_Euler(NT, nn, A, B, P, Tsim);
    % TODO compute the error in T using Tref
    
    T = crank_nicolson(NT, nn, A, B, P, Tsim);
    % TODO compute the error in T using Tref
    
    T = backward_Euler(NT, nn, A, B, P, Tsim);
    % TODO compute the error in T using Tref
    
    dt(ii) = 5*2^-ii;
end

%% plotting

loglog(dt, e_fe/norm_Tref, dt, e_cn/norm_Tref, dt, e_be/norm_Tref)
ylim([0, 1])
xlabel 'time increment [s]'
ylabel 'relative error [-]'
legend('Forward Euler','Crank-Nicolson','Backward Euler')

function Tfe = forward_Euler(NT, nn, A, B, P, Tsim)
time = linspace(0,Tsim,NT);
dt = time(2) - time(1);

Tfe = zeros(nn, NT); 
% TODO insert your Forward Euler code
end

function Tcn = crank_nicolson(NT, nn, A, B, P, Tsim)
time = linspace(0,Tsim,NT);
dt = time(2) - time(1);

I = speye(nn, nn);

Tcn = zeros(nn, NT); 
% TODO insert your Crank-Nicolson code
end

function Tbe = backward_Euler(NT, nn, A, B, P, Tsim)
time = linspace(0,Tsim,NT);
dt = time(2) - time(1);

I = speye(nn, nn);

Tbe = zeros(nn, NT); 
% TODO insert your Backward Euler code
end