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

%% TODO: complete question a

% TODO: Compute the state x using the Crank-Nocolson scheme
X = zeros(N, NT);
X(:,1) = TODO;
for ii = 2:NT
    X(:,ii) = TODO;
end

% TODO: Compute the cost functional J in the considered point
J = 0;
for ii = 1:NT-1
    J = J + TODO;
end

% plotting for part a
show_a = 0;  % set to 1 to visualize the solution for part a
if show_a
    Xmin = min(min(X));
    Xmax = max(max(X));
    for ii = 1:NT
        plot(x,X(:,ii),x,ones(size(x)))
        xlabel 'x [m]'
        ylabel 'T(x,t) [K]'
        ylim([Xmin, Xmax])
        title(['t = ', num2str(tgrid(ii)), ' [s]'])
        pause(0.1)
    end
end

%% TODO: compute the adjoint state
Phi = zeros(N, NT-1);
Phi(:,end) = TODO;
for ii = NT-2:-1:1
    Phi(:,ii) = TODO;
end

% TODO: use Phi to compute the gradient
g = TODO; 

show_b = 0;   % set to 1 so visualize the solution for b
if show_b
    figure()
    Phimin = min(min(Phi));
    Phimax = max(max(Phi));
    for ii = 1:NT-1
        plot(x,Phi(:,ii))
        xlabel 'x [m]'
        ylabel '\varphi(x,t) [K]'
        ylim([Phimin, Phimax])
        title(['t = ', num2str(tgrid2(ii)), ' [s]'])
        pause(0.1)
    end
end

figure
plot(tgrid2, g)
xlabel 't [s]'
ylabel 'gradient of J'