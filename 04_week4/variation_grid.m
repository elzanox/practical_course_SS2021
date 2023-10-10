clear all
close all
clc

% load the given velocity field.
load('Week4_velocity_field')

% plot the given velocity field
figure(1)
qx = 0.2;
qy = 0.0003;

% vary vy
vy = vy*qx;

% vary vx
vx = vx*qy;

% calculate the magnitude of the velocity
normv = sqrt(vx.^2 + vy.^2);

% plot the velocity field
surf(x,y,0*normv.', normv.');
colormap jet
shading interp
hold on
[X,Y] = ndgrid(x,y);
quiver(X,Y,vx,vy)
xlabel 'x [m]'
ylabel 'y [m]'
title 'velocity field'
cbar = colorbar;
ylabel(cbar, 'magnitude of the velocity [m/s]')
axis equal
axis tight
view(2)
pause;
