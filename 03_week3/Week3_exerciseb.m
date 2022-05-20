clear all
close all
clc

Lx = 0.5;     % [m] length of the plate in x-direction
Ly = 0.4;     % [m] length of the plate in y-direction
k  = 57.7;    % [W/m/K] thermal conductivity of steel
c  = 448;     % [J/kg/K] heat capacity of steel
rho = 7840;   % [kg/m3] mass density of steel
H  = 0.003;   % [m] thickness of the plate
Tsim = 60;     % [s] duration of the simulation     
P0 = 50;       % [W] heating power
t1 = 10;       % [s] the first heat load is active until t=t1
P =@(t) P0*(tanh(t1-t)+1)/2;  % applied heating power [W]

load('Week3_spatial_discretization')

NT = 49;       % number of grid points in time
time = linspace(0,Tsim,NT);
dt = time(2) - time(1);

nn = length(B);
I = speye(nn, nn); % define identity matrix of appropriate size

%% Forward Euler

Tfe = zeros(nn, NT); % every column of Tfe the solution vector at a different time instant

% TODO implement Forward Euler

%% Crank-Nicolson

Tcn = zeros(nn, NT); 

% TODO implement Crank-Nicolson

%% Backward Euler

Tbe = zeros(nn, NT); 

% TODO implement Backward Euler

%% Display the solution and make a movie if desired

% Choose which of the three solutions you want to display
T = Tfe;
T = Tcn;
T = Tbe;

create_movie = 0;  % do you want to create a movie?  Set to zero if not. 
% Unfortunately, it is not (yet) possible to create a movie in Octave (this is not important for the exercise)
if create_movie
    movie_name = 'Exercise4b_movie';  % file name
    vidfile = VideoWriter(movie_name,'MPEG-4');
    vidfile.FrameRate = 10;      % change this number to slow down or speed up the movie
    open(vidfile);
    fig = figure;
    set(fig,'color','w');
end

Tmax = max(max(T));
for kk = 1:NT
    Tkk = T(:,kk);
    surf(x,y,Tkk(node_nmbrs).')
    view(2)
    title(['t = ', num2str(time(kk)), ' [s]'])
    xlabel 'x [m]'
    ylabel 'y [m]'
    axis tight
    hcb = colorbar();
    ylabel(hcb,'temperature increase [K]')
    caxis([0,Tmax])
    shading interp
    axis equal
    
    if create_movie
        frame = getframe(fig);
        writeVideo(vidfile, frame);
    end
    
    pause(0.1);
end

if create_movie
    close(vidfile);
end
