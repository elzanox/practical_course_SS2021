clear all
close all
clc

load('Week11_bar_model')

T  = 8*60;     % length of the considered time interval
NT = 201;      % number of temporal grid points
tgrid = linspace(0,T,NT);    % temporal grid
dt = tgrid(2) - tgrid(1);

tgrid2 = tgrid(1:end-1) + diff(tgrid)/2;  % intermediate points in time grid
u0 = sin(pi*tgrid2/T/2);

% compute the state, cost functional, and adjoint state for the selected u0
X = computeX(E,A,B,Xinit,u0,dt,N,NT);
J = evalJ(u0,X,Q,NT,dt);
Phi = computePhi(E,A,Q,X,dt,N,NT);

%% TODO: complete part c.
g = TODO;   % TODO: Compute the gradient
    
G = TODO;                % TODO: Compute the linear coefficient

% TODO: Compute Xnabla, the change in state resulting from the gradient
Xnabla = TODO;  

% TODO: Use Xnabla to compute H, the quadratic coefficient
H      = evalHessian(g,Xnabla,Q,NT,dt);

% TODO: Find the (optimal) stepsize
step = TODO;

% Compare the quadratic approximation to the true function
step_list = linspace(0,2*step,20);
for kk = 1:length(step_list)
    u1 = u0 - step_list(kk)*g;
    X1 = computeX(E,A,B,Xinit,u1,dt,N,NT);
    Jlist(kk) = evalJ(u1,X1,Q,NT,dt);
end
plot(step_list, Jlist, step_list, J-G*step_list+H/2*step_list.^2)

% TODO: Use the answer from a to complete the function that computes the
% state x
function X = computeX(E,A,B,Xinit,u,dt,N,NT)
X = zeros(N, NT);
X(:,1) = TODO;
for ii = 2:NT
    X(:,ii) = TODO;
end
end

% TODO: Use your answer from b to make the function that computes the
% adjoint state Phi
function Phi = computePhi(E,A,Q,X,dt,N,NT)
Phi = zeros(N, NT-1);
Phi(:,end) = TODO;
for ii = NT-2:-1:1
    Phi(:,ii) = TODO;
end
end

% TODO: Complete the function that computes the cost functional J
function J = evalJ(u,X,Q,NT,dt)
J = 0;
for ii = 1:NT-1
    J = J + TODO;
end
end

% TODO: Complete the function that computes the Hessian
function H = evalHessian(u,Xnabla,Q,NT,dt)
H = 0;
for ii = 1:NT-1
    H = H + TODO;
end
end