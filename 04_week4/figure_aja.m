clear all
close all
clc

load('Week4_velocity_field')   % load the given velocity field. 
% The provided file contains the grid vectors x and y and 
% the x and y components of the velocity field vx and vy. 
% The velocity fields are given as matrices which contain the values in the
% grid points. 

%% plot the given velocity field
figure(1)
normv = sqrt(vx.^2 + vy.^2);   % the colors in the figure show the magnitude normv
disp(normv)
surf(x,y,0*normv.', normv.')
% view(2)
% print('normal y', '-dpng', '-S800,600');
colormap jet
% shading interp
% hold on
% [X,Y] = ndgrid(x,y);
% quiver(X,Y,vx,vy)              % the arrows in the figure show the direction
xlabel 'x [m]'
ylabel 'y [m]'
title 'velocity field'
cbar = colorbar;
ylabel(cbar, 'magnitude of the velocity [m/s]')
axis equal
axis tight
view(2)

pause;