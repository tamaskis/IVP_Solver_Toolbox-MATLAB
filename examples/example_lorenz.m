%% example_lorenz.m
% IVP Solver Toolbox
%
% Example simulating the Lorenz system. This script is also used to create 
% the thumbnail image for the IVP solver toolbox.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-05
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% SOLVES LORENZ SYSTEM

% Lorenz parameters
rho = 28;
sigma = 10;
beta = 8/3;

% Lorenz equations in vector form
f = @(t,x) [sigma*(x(2)-x(1));
            x(1)*(rho-x(3))-x(2);
            x(1)*x(2)-beta*x(3)];

% solution
[t,y] = odeabm(f,[0,100],[0;1;1.05],0.001);

% plots figure
figure;
plot3(y(:,1),y(:,2),y(:,3));
view(45,20);
axis off;