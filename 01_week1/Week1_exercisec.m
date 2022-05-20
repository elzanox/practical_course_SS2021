clear all   % clear all variables that existed before
close all   % close all figures
clc         % clear the command window

%% main loop: compute finite difference approximations on various grids

N_list = [11, 31, 101, 301, 1001]; % list of considered number of grid points
for k = 1:length(N_list)  % loop over the number of considered grid points
  N = N_list(k);
  
  % finite difference approximation
  [Tr_fd, dx(k)] = finite_difference_approximation(N);
  
  % TODO: obtain a reference solution and ...
  % ... compute the error in the finite difference approximation
  eLinf_rel(k) = TODO;
end

% plot the relative L2 error vs the grid spacing
loglog(dx, eLinf_rel)
grid on
xlabel 'grid spacing [m]'
ylabel 'relative L^\infty error [-]'
print('Week1_Tr_convergence2.jpg')

function [T, dx] = finite_difference_approximation(N)
  L   = 0.3;   % [m] length of the rod
  Acs = 0.01;  % [m2] cross sectional area of the rod
  k   = 237;   % [W/m/K] thermal conductivity of the rod
  T0  = 293;   % [K] reference temperature (not used)
  h   = 3;     % [W/K] cooling coefficient at the right boundary
  Q0  = 100;   % [W/m] intensity of the applied heat load 
  a   = 0.1;   % [m] width parameter for the applied heat load
  Q   =@(x) Q0*exp(-(x-L/2).^2/a^2); % [W/m] applied heat load
  
  x = linspace(0,L,N).';     % create the (uniform) spatial grid
  dx = x(2) - x(1);          % grid spacing (is thus constant)
  
  %% TODO: Insert the code of your finite difference approximation from b. here
end

function T = exact_solution(N)
  %% TODO: You can define a second function to compute the exact solution of the BVP
  % (but it is also possible to solve the exercise without that)
  
end