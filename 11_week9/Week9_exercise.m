clear all
close all
clc

load('Week9_rod_model')

w = 1e-6;               % weight
u0 = zeros(N,1);        % initial guess 
beta = 100;             % initial step size

% compute cost functional J0 at the initial guess
T0 = -A\(E*u0);
dT0meas = Emeas*T0 - Tmeas;
J0 = 0.5*(dT0meas.')*dT0meas + 0.5*w*u0.'*E*u0;

tolJ = 1e-3;
tolu = 1e-3;
max_iters = 10000;
for ii = 1:max_iters
    
    % TODO: Implement the basic gradient descent algorithm from the lecture
    
    if TODO % put your convergence conditions here. The break statement will make you exit the optimization loop. 
        break;
    end
    
    
end


disp(['number of iterations: ', num2str(ii)])

figure(1)
plot(x,u1)
xlabel 'x'
ylabel 'u_{opt} [W/m]'
title 'heat load'

figure(2)
plot(x,T1)
hold on
plot(xmeas,Tmeas, 'x')
xlabel 'x'
ylabel 'T_{opt}'
title 'Temperature field'