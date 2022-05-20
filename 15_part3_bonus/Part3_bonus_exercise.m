clear all
close all
clc

N = 2;    % state dimension
I = 64;  % number of data points
NT = 101; % number of time points
weights = [0.01,0.01];  % weights [w1, w2]

sigma =@(x) max(x, 0); % nonlinearity
dsigma =@(x) (x > 0);  % derivative of nonlinearity

V0 = zeros(2,2,NT-1);   % initial guess for Vk
b0 = zeros(2,NT-1);     % initial guess for bk

time = linspace(0,1,NT); % time grid
dt = time(2) - time(1);

theta = 2*pi*rand(1,I); % generate initial points x_in^i
r     = [0.4*rand(1,I/2), 0.8+0.4*rand(1,I/2)];
Xin    = [r.*cos(theta); r.*sin(theta)];
% Xin     = [-0.5+0.2*rand(1,I/4), 0.5+0.2*rand(1,I/4), 0.2*rand(1,I/2);
%          zeros(1,I)];

figure(1)  % visiualize the data
plot(Xin(1,1:I/2), Xin(2,1:I/2), 'b.')
hold on
plot(Xin(1,I/2+(1:I/2)), Xin(2,I/2+(1:I/2)), 'r.')
xlabel 'x_1'
ylabel 'x_2'
title 'X_0'
Xin = Xin(:);  % convert the initial data to a column vector

% create final poinst Yout^i
Yout = [0*ones(1,I/2), 2*ones(1,I/2);  % final condition in matrix form
    ones(1,I)];
Yout = Yout(:);  % convert the final data to a column vector


% compute the state and cost function at the initial point
XX0 = computeX(Xin, V0, b0, N, I, NT, sigma, dt); % time history

J0 = evalJ(Yout, XX0, V0, b0, NT, weights, dt);
Jstart = J0



%% use the commented part below to verify your implementation of the gradient
% Phi = computePhi(Yout, XX0, V0, b0, N, I, NT, dsigma, weights, dt);
% [gradV, gradb] = computeGradients(Phi, XX0, V0, b0, N, I, NT, sigma, weights, dt);

% %% test gradb
% db = -gradb;
% atest = linspace(0, 0.000001, 10);
% Jtest = zeros(size(atest));
% for ii = 1:length(atest)
%     b1 = b0 + atest(ii)*db;
%     XX1 = computeX(Xin, V0, b1, N, I, NT, sigma, dt);
%     Jtest(ii) = evalJ(Yout, XX1, V0, b1, NT, weights, dt);
% end
% Jappr = Jtest(1) + sum(sum(gradb.*db))*atest;
% figure(4)
% % plot(atest,Jtest);
% plot(atest, Jtest, atest, Jappr)
% %
% %% test gradV
% dV = -gradV;
% atest = linspace(0, 0.0001, 10);
% Jtest = zeros(size(atest));
% for ii = 1:length(atest)
%     V1 = V0 + atest(ii)*dV;
%     XX1 = computeX(Xin, V1, b0, N, I, NT, sigma, dt);
%     Jtest(ii) = evalJ(Yout, XX1, V1, b0, NT, weights, dt);
% end
% Jappr = Jtest(1) + sum(sum(sum(gradV.*dV)))*atest;
% figure(5)
% % plot(atest,Jtest);
% plot(atest, Jtest, atest, Jappr)

%% optimization algo
tic
Niters = 1000;  % 1000 iterations takes about 2 minutes on my laptop.
gamma = 0.001;
tries_hist = [];
Jhist = J0;
gamma_hist = [];
for ii = 1:Niters
    % TODO: Compute the gradients
    
    gradV = gradV / max(max(max(abs(gradV))));
    gradb = gradb / max(max(abs(gradb)));
    
    J1 = Inf;
    tries = 0;
    gamma = gamma*4;
    while J1 > J0
        V1 = V0 - gamma*gradV;
        b1 = b0 - gamma*gradb;
        % TODO: Compute J1 = J(V1, b1)
        gamma = gamma / 2;
        tries = tries + 1;
    end
    % you can keep a track of the history of the number of tries and the
    % learning rate if you want to
    %     tries_hist = [tries_hist, tries];
    %     gamma_hist = [gamma_hist, gamma];
    
    J0 = J1;
    V0 = V1;
    b0 = b1;
    XX0 = XX1;
    
    Jhist = [Jhist, J0]; % keep track of the cost function values in the different iterations
end
toc

Jend = J0

Xout = XX1(:, end);
Xout = reshape(Xout, 2, I);
figure(5)
plot(Xout(1,1:I/2), Xout(2,1:I/2), 'b.')
hold on
plot(Xout(1,I/2+(1:I/2)), Xout(2,I/2+(1:I/2)), 'r.')
xlabel 'x_1'
ylabel 'x_2'
title 'X_{out}'

figure(6)
plot(XX1(1:N:I, :).', XX1(2:N:I, :).', 'b-')
hold on
plot(XX1(I+(1:N:I), :).', XX1(I+(2:N:I), :).', 'r-')

figure(7)
plot(time,XX1(1:N:I, :).', 'b-')
hold on
plot(time,XX1(I+(2:N:I), :).', 'r-')
ylabel 'x_{i,1}'

figure(8)
plot(time, XX1(2:N:I, :).', 'b-')
hold on
plot(time, XX1(I+(2:N:I), :).', 'r-')
ylabel 'x_{i,2}'