clear all
close all
clc

load('Week10_plate_model')

w = 1e-6;               % weight
u0 = zeros(nn,1);        % initial guess 
beta = 100;             % initial step size

% compute cost functional J0 at the initial guess
T0 = -A\(E*u0);
em = Em*T0 - Tm;
J0 = 0.5*(em.')*em + 0.5*w*u0.'*E*u0;

gA =  TODO;      % gradient w.r.t. weighted inner product
gB =  TODO;      % gradient w.r.t. Eucliden inner product

GA = TODO;       % linear terms in quadratic approximations
GB = TODO;

TnablaA = TODO;  % temperature fields resulting from gradients gA and gB
TnablaB = TODO;

HA = TODO;  % second order terms in quadratic approximations
HB = TODO;

betaoptA = TODO;  % optimal step sizes betaopt
betaoptB = TODO;

% for validation: compute the cost functional for a range of values for
% beta and compare with the quadratic approximation constructed above
% Note: because the considered functional is quadratic, the quadratic 
% approximation should be exactly equal to the values of the cost
% functional
betaA = linspace(0, 2*betaoptA, 20);
betaB = linspace(0, 2*betaoptB, 20);

JtestA = zeros(size(betaA));
for kk = 1:length(betaA)
    u1 = u0 - gA*betaA(kk);
    T0 = -A\(E*u1);
    em = Em*T0 - Tm;
    JtestA(kk) = 0.5*(em.')*em + 0.5*w*u0.'*E*u0;
end

JtestB = zeros(size(betaB));
for kk = 1:length(betaA)
    u1 = u0 - gB*betaB(kk);
    T0 = -A\(E*u1);
    em = Em*T0 - Tm;
    JtestB(kk) = 0.5*(em.')*em + 0.5*w*u0.'*E*u0;
end

% plots to compare the quadratic approximation and the cost functional
figure(1)
plot(betaA, JtestA, betaA, J0 -GA*betaA + 0.5*HA*betaA.^2)
xlabel 'step size \beta'
ylabel 'J(u_0 + \beta \nabla J(u_0))'
title 'gradient w.r.t. the weighted inner product'

figure(2)
plot(betaB, JtestB, betaB, J0 -GB*betaB + 0.5*HB*betaB.^2)
xlabel 'step size \beta'
ylabel 'J(u_0 + \beta \nabla J(u_0))'
title 'gradient w.r.t. the Euclidean inner product'