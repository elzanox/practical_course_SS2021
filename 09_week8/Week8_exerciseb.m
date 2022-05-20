clear all
close all
clc

L  = 0.4;     % [m] length of the beam
rho = 2700;   % [kg/m3] mass density
E  = 69e9;    % [Pa] Young's modulus
W  = 0.02;    % [m] width of the beam
H  = 0.01;    % [m] height of the beam

A = H*W;      % [m2] cross sectional area
I = W*H^3/12; % [m4] second moment of area

%TODO: Compute the analytic expression for the eigenvalues
om_ana = TODO;

M_list = [10 30 100 300 1000 3000, 10000];   % considered values for M
omega_FEM = zeros(6, length(M_list)); % matrix to store the eigenfrequencies for different values of M
err_FEM = zeros(6, length(M_list));   % matrix to store the errors for different values of M

for kk = 1:length(M_list)
    M = M_list(kk);
    omega_FEM(:,kk) = get_eigenfrequencies_FE(M);
    % note: we assume that get_eigenfrequencies_FE returns a column vector
    % when it returns a row vector you can add a transpose .' to the RHS
    err_FEM(:,kk) = abs(omega_FEM(:,kk) - om_ana.');
    % note: the last line assumes that om_ana is a row vector
    % when you defined om_ana as a column vector, remove the transpose .'
end

loglog(M_list,err_FEM.', '--x')
xlabel 'number of elements M [-]'
ylabel 'error in eigenfrequency [rad/s]'
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5', 'mode 6')

function angular_eigen_frequencies = get_eigenfrequencies_FE(M)
L  = 0.4;     % [m] length of the beam
rho = 2700;   % [kg/m3] mass density
E  = 69e9;    % [Pa] Young's modulus
W  = 0.02;    % [m] width of the beam
H  = 0.01;    % [m] height of the beam

A = H*W;      % [m] cross sectional area
I = W*H^3/12; % [m4] second moment of area

% TODO: insert your code from b.
angular_eigen_frequencies = TODO; % [rad/s]
end