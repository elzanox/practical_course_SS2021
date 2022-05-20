clear all
close all
clc

load('Week11_bar_model')

T  = 8*60;     % length of the considered time interval
NT = 201;      % number of temporal grid points
tgrid = linspace(0,T,NT);    % temporal grid
dt = tgrid(2) - tgrid(1);

tgrid2 = tgrid(1:end-1) + diff(tgrid)/2;  % intermediate points in time grid
u0 = sin(2*pi*tgrid2/T);

X0 = computeX(E,A,B,Xinit,u0,dt,N,NT);
J0 = evalJ(u0,X0,Q,NT,dt);
Jinit = J0

tolJ = 1e-5;
tolu = 1e-3;
max_iters = 2000;
for ii = 1:max_iters
    Phi = TODO;
    g = TODO;
    
    G = TODO;
    Xnabla = TODO;
    H      = TODO;
    step = TODO;

    % You can uncomment the part below to check your implementation
%     step_list = linspace(0,2*step,20);
%     for kk = 1:length(step_list)
%         u1 = u0 - step_list(kk)*g;
%         X1 = computeX(E,A,B,Xinit,u1,dt,N,NT);
%         Jlist(kk) = evalJ(u1,X1,Q,NT,dt);
%     end
%     plot(step_list, Jlist, step_list, J0-G*step_list+H/2*step_list.^2)
    
    u1 = TODO;
    X1 = TODO;
    J1 = TODO;
    
    
    if J0 - J1 < tolJ*J0 && norm(u1-u0) < tolu*norm(u0)
        break
    end
    
    u0 = TODO;
    X0 = TODO;
    J0 = TODO;
end

% plotting
Jfinal = J1

figure()
plot(tgrid2, u1)
xlabel 't [s]'
ylabel 'u(t) [W]'

figure()
plot(tgrid, X1(end,:))
xlabel 't [s]'
ylabel 'T(L,t) [K]'

figure()
Xmax = max(max(X1));
Xmin = min(min(X1));
for ii = 1:NT
    plot(x,X1(:,ii),x,ones(size(x)))
    xlabel 'x [m]'
    ylabel 'T(x,t) [K]'
    ylim([Xmin, Xmax])
    title(['t = ', num2str(tgrid(ii)), ' [s]'])
    pause(0.1)
end

% Use your answers from part c to complete the functions below
function X = computeX(E,A,B,Xinit,u,dt,N,NT)
X = zeros(N, NT);
X(:,1) = TODO;
for ii = 2:NT
    X(:,ii) = TODO;
end
end

function Phi = computePhi(E,A,Q,X,dt,N,NT)
Phi = zeros(N, NT-1);
Phi(:,end) = TODO;
for ii = NT-2:-1:1
    Phi(:,ii) = TODO;
end
end

function J = evalJ(u,X,Q,NT,dt)
J = 0;
for ii = 1:NT-1
    J = J + TODO;
end
end

function H = evalHessian(u,Xnabla,Q,NT,dt)
H = 0;
for ii = 1:NT-1
    H = H + TODO;
end
end