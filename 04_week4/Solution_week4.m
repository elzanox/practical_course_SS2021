clear all
close all
clc

load('Week4_velocity_field')   % load the given velocity field. 
% The provided file contains the grid vectors x and y and 
% the x and y components of the velocity field vx and vy. 
% The velocity fields are given as matrices which contain the values in the
% grid points. 
% Mengalikan komponen x dan y medan kecepatan dengan faktor 4
%vx = 40 * vx;
%vy = 40 * vy;
%% plot the given velocity field
figure(1)
normv = sqrt(vx.^2 + vy.^2);   % the colors in the figure show the magnitude normv
surf(x,y,0*normv.', normv.')
colormap jet
shading interp
hold on
[X,Y] = ndgrid(x,y);
quiver(X,Y,vx,vy)              % the arrows in the figure show the direction
xlabel 'x [m]'
ylabel 'y [m]'
title 'velocity field'
cbar = colorbar;
ylabel(cbar, 'magnitude of the velocity [m/s]')
axis equal
axis tight
view(2)

%% define other parameters
Lx = max(x);             % determine the parameters of the spatial 
Ly = max(y);                   % discretization from the given grid vectors
Nx = length(x);
Ny = length(y);
dx = x(2) - x(1);
dy = y(2) - y(1);

kappa  = 1.25e-5;    % [m^2/s] diffusivity of chlorine in water
c0 = 10;             % [mg/l] magnitude of the initial concentration of chlorine
x0 = 0.1*Lx;         % [m] x-coordinate of the center of the initial concentration
y0 = 0.7*Ly;         % [m] y-coordinate of the center of the initial concentration
a = 0.03*Ly;         % [m] width parameter for the initial concentration
Tsim = 8*60;         % [s] length of the simulation (=8 min)

NT = 100;            % number of time steps
time = linspace(0,Tsim,NT);
dt = time(2) - time(1);

%% construct the node_nmbrs matrix
nn = Nx*Ny;
node_nmbrs = 1:nn;
node_nmbrs = reshape(node_nmbrs, Nx, Ny);

%% construct the A-matrix
A = sparse(nn,nn);
% Konstruksi A-matrix untuk node-node internal
for ii = 2:Nx-1
    for jj = 2:Ny-1
        node11 = node_nmbrs(ii,jj);
        node21 = node_nmbrs(ii+1,jj);
        node12 = node_nmbrs(ii,jj+1);
        node01 = node_nmbrs(ii-1,jj);
        node10 = node_nmbrs(ii,jj-1);

        A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
        %A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2 + (vx(ii)*4*dt/dx^2 + vy(jj)*4*dt/dy^2);
        A(node11, node21) = kappa/dx^2;
        A(node11, node01) = kappa/dx^2;
        A(node11, node12) = kappa/dy^2;
        A(node11, node10) = kappa/dy^2;
    end
end


% We simply neglect the velocity field in the equations for nodes on the
% boundary. Because the velocity on the boundary is zero, this will
% (probably) not lead to a significant error. 
for ii = 2:Nx-1                    % lower edge
    node11 = node_nmbrs(ii,1);
    node21 = node_nmbrs(ii+1,1);
    node12 = node_nmbrs(ii,2);
    node01 = node_nmbrs(ii-1,1);
    
    A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
    A(node11, node21) =    kappa/dx^2;
    A(node11, node12) =  2*kappa/dy^2;
    A(node11, node01) =    kappa/dx^2;
end

for ii = 2:Nx-1                   % top edge
    node11 = node_nmbrs(ii,end);
    node21 = node_nmbrs(ii+1,end);
    node12 = node_nmbrs(ii,end-1);
    node01 = node_nmbrs(ii-1,end);
    
    A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
    A(node11, node21) =    kappa/dx^2;
    A(node11, node12) =  2*kappa/dy^2;
    A(node11, node01) =    kappa/dx^2;
end

for jj = 2:Ny-1                 % left edge
    node11 = node_nmbrs(1,jj);
    node21 = node_nmbrs(2,jj);
    node12 = node_nmbrs(1,jj+1);
    node10 = node_nmbrs(1,jj-1);
        
    A(node11, node11) = -2*kappa/dx^2 -2*kappa/dy^2;
    A(node11, node21) =  2*kappa/dx^2;
    A(node11, node12) =    kappa/dy^2;
    A(node11, node10) =    kappa/dy^2;
end

for jj = 2:Ny-1              % right edge
    node11 = node_nmbrs(end,jj); 
    node21 = node_nmbrs(end-1,jj);
    node12 = node_nmbrs(end,jj+1);
    node10 = node_nmbrs(end,jj-1);
        
    A(node11, node11) = -2*kappa/dx^2 -2*kappa/dy^2;
    A(node11, node21) =  2*kappa/dx^2;
    A(node11, node12) =    kappa/dy^2;
    A(node11, node10) =    kappa/dy^2;
end

node11 = node_nmbrs(1,1);      % lower left corner
node21 = node_nmbrs(2,1);
node12 = node_nmbrs(1,2);
A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
A(node11, node21) =  2*kappa/dx^2;
A(node11, node12) =  2*kappa/dy^2;

node11 = node_nmbrs(end,1);    % lower right corner
node21 = node_nmbrs(end-1,1);
node12 = node_nmbrs(end,2);
A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
A(node11, node21) =  2*kappa/dx^2;
A(node11, node12) =  2*kappa/dy^2;

node11 = node_nmbrs(1,end);    % top left corner
node21 = node_nmbrs(2,end); 
node12 = node_nmbrs(1,end-1);
A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
A(node11, node21) =  2*kappa/dx^2;
A(node11, node12) =  2*kappa/dy^2;

node11 = node_nmbrs(end,end);  % top right corner
node21 = node_nmbrs(end-1,end);
node12 = node_nmbrs(end,end-1);
A(node11, node11) = -2*kappa/dx^2 - 2*kappa/dy^2;
A(node11, node21) =  2*kappa/dx^2;
A(node11, node12) =  2*kappa/dy^2;

%% Dirchlet boundary conditions
% We impose zero dirichlet BCs at the inlet and at the outlet
% To determine the free dofs in our problem, we thus need to find the node
% numbers of the nodes at the inlet and the outlet

% Identifikasi node-node di inlet dan outlet
nodes_inlet = node_nmbrs(1,:);
nodes_outlet = node_nmbrs(end,:);
% nodes_inlet adalah semua node pada tepi kiri, nodes_outlet adalah semua node pada tepi kanan
fdofs = setdiff(1:nn, [nodes_inlet, nodes_outlet]);

%% Initial condition
rho0 = zeros(nn,1);   % this vector will represent the intial condition.
                      % Use the node_nmbrs matrix to build this vector.
                      
% TODO: Define rho0
% Hitung kondisi awal konsentrasi
for ii = 1:Nx
    for jj = 1:Ny
        x_diff = x(ii) - x0;
        y_diff = y(jj) - y0;
        rho0(node_nmbrs(ii,jj)) = c0 * exp(-((x_diff^2 + y_diff^2) / (2*a^2)));
    end
end
%% Time integration
rho = zeros(nn, NT); % Solution matrix for the density of chlorine. 
                     % The columns of this matrix contain the solutions at different time instances
rho(:,1) = rho0;     % The first column is thus equal to the initial condition

I = speye(nn,nn);
% TODO: Implement a time discretization of scheme of your choice.
for kk = 2:NT
    rhokk = rho(:,kk-1);
    rhokk(nodes_inlet) = 0; % Kondisi batas Dirichlet di inlet (konsentrasi nol)
    rhokk(nodes_outlet) = 0; % Kondisi batas Dirichlet di outlet (konsentrasi nol)

    % Perbarui vektor rhokk dengan hasil perkalian matriks A dan rhokk
    rhokk = rhokk + dt * (A * rhokk);

    rho(:,kk) = rhokk;
end


%% plot the obtained solution

create_movie = 0;  % do you want to create a movie?  Set to zero if not. 
% Unfortunately, it is not (yet) possible to create a movie in Octave (this is not important for the exercise)
if create_movie
    movie_name = 'Week4_exercise_movie';  % file name
    vidfile = VideoWriter(movie_name,'MPEG-4');
    vidfile.FrameRate = 10;      % change this number to slow down or speed up the movie
    open(vidfile);
    fig = figure(2);
    set(fig,'color','w');
end

rho_max = max(max(rho)); % used to fix the color scale for all frames
for kk = 1:NT
    rhokk = rho(:,kk);
    surf(x,y,0*rhokk(node_nmbrs).',rhokk(node_nmbrs).')
    hold on
    quiver(X,Y,vx,vy)
    view(2)
    title(['t = ', num2str(time(kk)), ' [s]'])
    xlabel 'x [m]'
    ylabel 'y [m]'
    axis tight
    hcb = colorbar();
    ylabel(hcb,'concentration of Cl [mg/l]')
%     caxis([0,rho_max])  % I did not fix to color axis to see the solution clearer, but this is somewhat deceiving. 
    colormap default
    shading interp
    axis equal
    hold off
    
    if create_movie
        frame = getframe(fig);
        writeVideo(vidfile, frame);
    end
    
    pause(0.1);
end

if create_movie
    close(vidfile);
end