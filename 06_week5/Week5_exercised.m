clear all
close all
clc

M_list = [10, 30, 100, 300, 1000];

eP1 = zeros(1,length(M_list));
eP2 = zeros(1,length(M_list));
Le  = zeros(1,length(M_list));
for kk = 1:length(M_list)
    M = M_list(kk);
    
    [eP1(kk), Le(kk)] = error_in_P1solution(M);
    [eP2(kk), ~     ] = error_in_P2solution(M);
end

figure(1)
loglog(Le, eP1, Le, eP2)
grid on
xlabel 'element size Le [m]'
ylabel 'relative Linfty error [-]'
legend('P1 elements', 'P2 elements')

figure(2)
loglog(M_list+1, eP1, 2*M_list+1, eP2)
grid on
xlabel 'number of nodes [-]'
ylabel 'relative Linfty error [-]'
legend('P1 elements', 'P2 elements')

function [out,Le] = error_in_P1solution(M)
L   = 0.3;   % [m] length of the rod
Acs = 0.01;  % [m2] cross sectional area of the rod
k   = 237;   % [W/m/K] thermal conductivity of the rod
T0  = 293;   % [K] reference temperature (not used)
h   = 3;     % [W/K] cooling coefficient at the right boundary
Q0  = 100;   % [W/m] intensity of the applied heat load 
a   = 0.1;   % [m] width parameter for the applied heat load
Q   =@(x) Q0*exp(-(x-L/2).^2/a^2); % [W/m] applied heat load

% TODO: Insert your code to compute the solution T based on P1-elements

% determine the error
Tref = exact_solution(M);
err = T - Tref;
out = TODO;  % compute the relative L2 error.

end


function [out,Le] = error_in_P2solution(M)
L   = 0.3;   % [m] length of the rod
Acs = 0.01;  % [m2] cross sectional area of the rod
k   = 237;   % [W/m/K] thermal conductivity of the rod
T0  = 293;   % [K] reference temperature (not used)
h   = 3;     % [W/K] cooling coefficient at the right boundary
Q0  = 100;   % [W/m] intensity of the applied heat load 
a   = 0.1;   % [m] width parameter for the applied heat load
Q   =@(x) Q0*exp(-(x-L/2).^2/a^2); % [W/m] applied heat load

% TODO insert your code that computes the solution T based on P2 elements

% determine the error
Tref = exact_solution(2*M);
err = T - Tref;
out = TODO; % compute the relative L2 error

end

function T = exact_solution(M)
  %% This function computes the analytic solution of the considered BVP
  
  L   = 0.3;   % [m] length of the rod
  Acs = 0.01;  % [m2] cross sectional area of the rod
  k   = 237;   % [W/m/K] thermal conductivity of the rod
  h   = 3;     % [W/K] cooling coefficient at the right boundary
  Q0  = 100;   % [W/m] intensity of the applied heat load 
  a   = 0.1;   % [m] width parameter for the applied heat load
  
  Tpart =@(x) -Q0*a^2/2/k/Acs*exp(-(x-L/2).^2/a^2) - sqrt(pi)*a*Q0*(x-L/2)/2/Acs/k.*erf((x-L/2)/a);
  dTpartdx =@(x) -sqrt(pi)*a*Q0/2/Acs/k.*erf((x-L/2)/a);
  Mat = [0 1; k*Acs+h*L, h];
  vec = [Tpart(0); k*Acs*dTpartdx(L)+h*Tpart(L)];
  coeff = -Mat\vec;
  
  xgrid = linspace(0,L,M+1).';    % grid on which the solution will be evaluated
  T = Tpart(xgrid) + coeff(1)*xgrid + coeff(2);
end