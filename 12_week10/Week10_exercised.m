clear all
close all
clc

load('Week10_plate_model')

w = 1e-9;               % weight
u0 = zeros(nn,1);        % initial guess
beta = 100;             % initial step size

% compute cost functional J0 at the initial guess
T0 = -A\(E*u0);
em0 = Em*T0 - Tm;
J0 = 0.5*(em0.')*em0 + 0.5*w*u0.'*E*u0;

tolJ = 1e-5;
tolu = 1e-5;
max_iters = 10000;
umin = 0;            % minimal and maximal values for u
umax = 100;
for ii = 1:max_iters
    
    g =  TODO;      % gradient w.r.t. weighted inner product
    gPr = TODO;     % compute the projected gradient
    betaopt = TODO; % determine the optimal step size based on a quadratic approximation
    
    beta = 2*betaopt;
    J1 = Inf;
    while J1 > J0
        beta = beta/2;
        % TODO: do a line search to assure that the cost functional
        % decreases in every iteration
    end
    
    if %TODO: check convergence
        break;
    end
    
    % TODO: update u0, J0, etc.
end


disp(['number of iterations: ', num2str(ii)])

figure(1)
surf(xgrid, ygrid, u1(node_nmbrs).')
xlabel 'x [m]'
ylabel 'y [m]'
ylabel 'u_{opt}'

figure(2)
surf(xgrid, ygrid, T1(node_nmbrs).')
hold on
plot3(xm, ym, Tm, 'rx', 'MarkerSize',10,'linewidth',2)
xlabel 'x [m]'
ylabel 'y [m]'
ylabel 'T_{opt}'