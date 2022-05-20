clear all
close all
clc

load('Week10_plate_model')

w = 1e-6;               % weight
u0 = zeros(nn,1);        % initial guess 
beta = 100;             % initial step size

% compute cost functional J0 at the initial guess
T0 = -A\(E*u0);
em = Em*T0 - Tm;
J0 = 0.5*(em.')*em + 0.5*w*u0.'*E*u0;

gA =  TODO;      % gradient w.r.t. weighted inner product
gB =  TODO;      % gradient w.r.t. Eucliden inner product

figure(1)
surf(xgrid, ygrid, gA(node_nmbrs).')
xlabel 'x [m]'
ylabel 'y [m]'
title 'gradient w.r.t. the weighted inner product'


figure(2)
surf(xgrid, ygrid, gB(node_nmbrs).')
xlabel 'x [m]'
ylabel 'y [m]'
title 'gradient w.r.t. the Euclidean inner product'