clear all
close all
clc

load('Week12_string_model')
R = eye(2);

%% Parameters for time discretization
NT = 1000;
T = 0.25;
tgrid = linspace(0,T,NT);
dt = tgrid(2) - tgrid(1);
tgrid2 = tgrid(1:end-1) + diff(tgrid)/2;

U0 = zeros(2,NT-1);  % initial guess for the control

%% Create the matrices for the first order system

E = TODO;
A = TODO;
B = TODO;
Q = TODO;
Xinit = TODO;

X0 = compute_X(E, A, Xinit, B, U0, NT, dt);
J0 = cost_function(Q,R,X0,U0,NT,dt);
Jinit = J0

%% optimization

max_iters = 500; tolJ = 1e-5; tolu = 1e-4;
for ii = 1:max_iters
    % compute the adjoint state and
    Phi = compute_phi(E, A, Q, X0, NT, dt);
    g = TODO;
    
    G = TODO;   % linear term G
    
    % quadratic term
    Xnabla = TODO;
    H = hessian(Q,R,Xnabla,g,NT,dt);
    
    % determine the step size
    step  = TODO;
    
%     % uncomment this part to check your quadratic approximation during the
%     % optimization (this is time consuming, so only for debugging)
%     steps = linspace(0,2*step,20);
%     for kk = 1:length(steps)
%         U1 = U0 - steps(kk)*g;
%         X1 = compute_X(E, A, Xinit, B, U1, NT, dt);
%         J(kk) = cost_function(Q,R,X1,U1,NT,dt);
%     end
%     figure
%     plot(steps, J, steps, J0-G*steps+H/2*steps.^2)
%     xlabel 'stepsize \beta'
%     ylabel 'cost function'
%     legend('cost function', 'quadratic approximation')
%     pause
    
    % find the new iterate for the control U1 and compute the corresponding
    % state X1 and cost functional J1
    U1 = TODO;
    X1 = TODO;
    J1 = TODO;
    
    if abs(J1-J0)/abs(J0) < tolJ && norm(U1-U0)/norm(U0) < tolu
        break;
    end
    
    % update U, X, and J for the new iteration
    U0 = U1;
    X0 = X1;
    J0 = J1;
end
Jfinal = J1

figure(3)
plot(tgrid2, U1)
title 'optimal control'
xlabel 't [s]'
ylabel 'u_{opt}'

figure(4)
for ii = 1:10:NT
    plot(x,[0; X1(1:Nf,ii); 0])
    title(['t = ', num2str(tgrid(ii))]) 
    xlabel 'x'
    ylabel 'u(t,x)'
    ylim([-1,1])
    pause(0.1)
end