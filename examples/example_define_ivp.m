%% example_define_ivp.m
% IVP Solver Toolbox
%
% Example of defining an IVP.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-07
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% SOLUTION

% parameters
b = 5;          % damping constant [N.s/m]
k = 1;          % spring constant [N/m]
m = 2;          % mass [kg]
x0 = 1;         % initial position [m]
xdot0 = 0;        % initial velocity [m/s]

% forcing function
F = @(t) cos(pi*t);

% initial condition
y0 = [x0;
      xdot0];

% ODE
f = @(t,y) [y(2);
            -(b/m)*y(2)-(k/m)*y(1)+(1/m)*F(t)];



%% ALTERNATE SOLUTION

% parameters
b = 5;          % damping constant [N.s/m]
k = 1;          % spring constant [N/m]
m = 2;          % mass [kg]
x0 = 1;         % initial position [m]
xdot0 = 0;      % initial velocity [m/s]

% forcing function
F = @(t) cos(pi*t);

% initial condition
y0 = [x0;
      xdot0];

% assigns function handle to ODE
f = @(t,y) f_extra(t,y,b,k,m,F);

% defines ODE
function dy = f_extra(t,y,b,k,m,F)
    
    % unpacks state vector
    x = y(1);
    xdot = y(2);
    
    % preallocates state vector derivative
    dy = zeros(size(y));
    
    % assembles state vector derivative
    dy(1) = xdot;
    dy(2) = -(b/m)*xdot-(k/m)*x+(1/m)*F(t);
    
end