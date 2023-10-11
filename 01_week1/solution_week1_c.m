clear all   % clear all variables that existed before
close all   % close all figures
clc         % clear the command window

L   = 0.3;   % [m] length of the rod
Acs = 0.01;  % [m2] cross sectional area of the rod
k   = 237;   % [W/m/K] thermal conductivity of the rod
T0  = 293;   % [K] reference temperature (not used)
h   = 3;     % [W/K] cooling coefficient at the right boundary
Q0  = 100;   % [W/m] intensity of the applied heat load
a   = 0.1;   % [m] width parameter for the applied heat load
Q   =@(x) Q0*exp(-(x-L/2).^2/a^2); % [W/m] applied heat load


N_list = [11, 31, 101, 301, 1001]; % list of considered number of grid points



for k = 1:length(N_list)  % loop over the number of considered grid points
    % make a figure for each N
    x = linspace(0, L, N_list(k) + 2).';      % spatial grid
    dx = x(2) - x(1);            % grid spacing (is constant)

    % Construct the stiffness matrix A here
    A = sparse(N_list(k) + 2, N_list(k) + 2);    % creates an NxN sparse matrix (so all elements are zero)
    A(1, 2) = 1;

    for ii = 2:N_list(k) + 1
        A(ii, ii - 1) = 1 / dx^2;
        A(ii, ii) = -2 / dx^2;
        A(ii, ii + 1) = 1 / dx^2;
    end

    A(end, end) = 1;
    A(end, end - 2) = -1;

    % Construct the load vector f here
    f = zeros(N_list(k) + 2, 1);     % creates a zero column vector of length N
    % f = Q(x)*dx^2;      % update f with the heat source term
    % f = Q(x)*dx;      % update f with the heat source term
    
    for ii = 1:N_list(k)+2
    f(ii) = Q(x(ii)) * dx;
    end
    % Compute the load vector by integrating the heat load over each cell
    % for ii = 2:N_list(k) + 1
    %     f(ii) = (Q(x(ii - 1)) + Q(x(ii))) / 2 * dx;
    % end

    % Compute the resulting temperature field Tr
    Tr = A \ f;

    % exact
    % T_exact = Ta + (Tb - Ta) * x / L;
    % T_exact = T0 + (Q0 / (2 * 237)) * (1 - exp(-(x - L/2).^2 / a^2)) - (h * L / (2 * 237)) * (x - L);
    T_exact = T0 + (Q0 / (2 * 237)) * (1 - exp(-(x - L/2).^2 / a^2)) - (h * L^2 / (2 * 237)) * (x - L);

    % Menghitung error relatif L^∄1�7
    % eLinf_rel(k) = max(abs(T - T_exact)); 
    % eLinf_rel(k) = max(abs(T_exact - T)) / max(abs(T_exact));
    % eLinf_rel(k) = max(abs(Tr_fd - Tr_exact)) / max(abs(Tr_exact));
    eLinf_rel(k) = norm(Tr - T_exact, inf) / norm(T_exact, inf); % Hitung kesalahan relatif L^∞

   
    figure(1); % Create a single figure
    subplot(2, 3, k); % Gambar pertama di jendela gambar utama
    plot(x, Q(x));
    xlabel 'x [m]'
    ylabel 'Q(x)'
    title(['Q(x) Profile for N = ', num2str(N_list(k))]);

    figure(2)
    subplot(2, 3, k); % Gambar pertama di jendela gambar utama
    plot(x, Tr);
    xlabel 'x [m]'
    ylabel 'T(x)'
    title(['T(x) [K] Profile for N = ', num2str(N_list(k))]);

    % plot the relative L2 error vs the grid spacing
    figure(3)
    subplot(2, 3, k); % Gambar pertama di jendela gambar utama
    loglog(dx, eLinf_rel)
    grid on
    xlabel 'grid spacing [m]'
    ylabel 'relative L^\infty error [-]'

end

pause;


%N = 11;
% N_list = [11, 31, 101, 301, 1001]; % list of considered number of grid points

% for k = 1:length(N_list)  % loop over the number of considered grid points
% % make a figure for each N
%   x = linspace(0,L,N_list(k)+2).';      % spatial grid
%   dx = x(2) - x(1);            % grid spacing (is constant)

% % make a figure of the applied heat load (just for understanding)

%   %% Construct the stiffness matrix A here
%   A = sparse(N_list(k)+2,N_list(k)+2);    % creates a NxN sparse matrix (so all elements are zero)
%   A(1,2) = 1;

%   for ii = 2:N_list(k)+1
%     A(ii,ii-1) = 1/dx^2;
%     A(ii,ii) = -2/dx^2;
%     A(ii,ii+1) = 1/dx^2;
%   end

%   A(end, end) = 1;
%   A(end, end-2) = -1;

%   %% Construct the load vector f here
%   f = zeros(N_list(k)+2,1);     % creates a zero column vector of length N
%   % f = Q(x)*dx^2;      % update f with the heat source term
%   % Compute the load vector by integrating the heat load over each cell
%   for ii = 2:N_list(k)+1
%       f(ii) = (Q(x(ii-1)) + Q(x(ii))) / 2 * dx;
%   end

%   %% Compute the resulting temperature field Tr
%   T = A \ f;
%   %% exact
%   % T_exact = Ta + (Tb - Ta) * x / L;
%   % T_exact = T0 + (Q0 / (2 * 237)) * (1 - exp(-(x - L/2).^2 / a^2)) - (h * L / (2 * 237)) * (x - L);
%   T_exact = T0 + (Q0 / (2 * 237)) * (1 - exp(-(x - L/2).^2 / a^2)) - (h * L^2 / (2 * 237)) * (x - L);

%   % Menghitung error relatif L^∄1�7
%   % eLinf_rel(k) = max(abs(T - T_exact)); 
%   % eLinf_rel(k) = max(abs(T_exact - T)) / max(abs(T_exact));
%   % eLinf_rel(k) = max(abs(Tr_fd - Tr_exact)) / max(abs(Tr_exact));
%   eLinf_rel(k) = norm(T - T_exact, inf) / norm(T_exact, inf); % Hitung kesalahan relatif L^∞
  
%   disp(eLinf_rel)
  
%   fig = figure();
%   %set(gcf, 'Position', [0, 0, 640, 480]); % Mengatur ukuran jendela gambar
%   subplot(2, 2, 1); % Gambar pertama di jendela gambar utama
%   plot(x, Q(x));
%   xlabel 'x [m]'
%   ylabel 'Q(x)'
%   title(['Q(x) Profile for N = ', num2str(N_list(k))]);
%   % %fig = figure();
%   % subplot(2, 2, 2); % Gambar pertama di jendela gambar utama
%   % plot(x, T);
%   % xlabel 'x [m]'
%   % ylabel 'T(x) [K]'
%   % title(['T(x) [K] Profile for N = ', num2str(N_list(k))]);
%   % % plot the relative L2 error vs the grid spacing
%   % subplot(2, 2, 3); % Gambar pertama di jendela gambar utama
%   % loglog(dx, eLinf_rel)
%   % grid on
%   % xlabel 'grid spacing [m]'
%   % ylabel 'relative L^\infty error [-]'
%   % print('Week1_Tr_convergence2.jpg')
% end
% pause;
%fig = figure(1);
%plot(x, Q(x))
%xlabel 'x [m]'
%ylabel 'Q(x) [W/m]'



% save a figure of the resulting temperature field
%fig = figure(2);
%plot(x,T)
%xlabel 'x [m]'
%ylabel 'T(x) [K]'
%print('Week1_T.jpg')

