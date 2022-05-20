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

U0 = zeros(2,NT-1);  % control variable

%% Create the matrices for the first order system

E = TODO;
A = TODO;
B = TODO;
Q = TODO;
Xinit = TODO;

%% forward dynamics
X0 = compute_X(E, A, Xinit, B, U0, NT, dt);
showX = 0;
if showX
    figure(1)
    for ii = 1:10:NT
        plot(x,[0; X0(1:Nf,ii); 0])
        title(['t = ', num2str(tgrid(ii))])
        xlabel 'x'
        ylabel 'u(t,x)'
        ylim([-1,1])
        pause(0.1)
    end
end

%% cost function
J0 = cost_function(Q,R,X0,U0,NT,dt);

%% adjoint equation
Phi = compute_phi(E, A, Q, X0, NT, dt);
showPhi = 0;
if showPhi
    figure(2)
    for ii = 1:10:NT
        plot(x,[0; Phi(1:Nf,ii); 0])
        title(['t = ', num2str(tgrid(ii))])
        xlabel 'x'
        ylabel '\phi(t,x)'
        ylim([-1.5,1.5])
        pause(0.1)
    end
end

%% gradient
g = TODO;
plot(tgrid2, g.')
xlabel 'time [s]'
ylabel 'gradient'

%% quadratic approximation
G = 0; % TODO: compute the linear coefficient

% TODO: compute Xnabla and the quadratic coefficient H
Xnabla = TODO;
H = hessian(Q,R,Xnabla,g,NT,dt); 

step  = TODO;

steps = linspace(0,2*step,20);
for kk = 1:length(steps)
    U1 = U0 - steps(kk)*g;
    X1 = compute_X(E, A, Xinit, B, U1, NT, dt);
    J(kk) = cost_function(Q,R,X1,U1,NT,dt);
end
figure
plot(steps, J, steps, J0-G*steps+H/2*steps.^2)
xlabel 'stepsize \beta'
ylabel 'cost function'
legend('cost function', 'quadratic approximation')