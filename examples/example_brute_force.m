%% example_brute_force.m
% ODE Solver Toolbox
%
% Example of a "brute force" method for solving ODEs.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-04-16
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% SOLUTION

% car parameters
m = 1500;       % mass [kg]
A = 2.5;        % cross-sectional area [m^2]
CD = 0.25;      % drag coefficient [-]
F_eng = 4500;   % applied force [N]

% IVP parameters
x0 = 0;         % initial position [m]
v0 = 0;         % initial velocity [m/s]
xf = 300;       % final position [m]

% physical parameters
rho = 1.2;      % air density [kg/m^3]

% computational parameters
dt = 0.1;       % time step [s]

% preallocate arrays
t = zeros(10000,1);
x = zeros(size(t));
v = zeros(size(t));
a = zeros(size(t));

% solves initial value problem
n = 1;
while x(n) < xf
    
    % calculates drag force [N]
    FD = 0.5*CD*A*rho*v(n)^2;
    
    % calculates acceleration [m/s^2]
    a(n) = (F_eng-FD)/m;
    
    % integrates acceleration to find velocity [m/s]
    v(n+1) = v(n)+a(n)*dt;
    
    % integrates velocity to find position [m]
    x(n+1) = x(n)+v(n)*dt;
    
    % increments loop index
    n = n+1;
    
end
    
% trims arrays
t = t(1:(n-1));
x = x(1:(n-1));
v = v(1:(n-1));
a = a(1:(n-1));

% velocity and acceleration of car when position is x = 300 m
v300 = v(end);
a300 = a(end);

% prints results
fprintf("    velocity = %.3f m/s\n",v300)
fprintf("acceleration = %.3f m/s^2\n",a300)
fprintf("  drag force = %.3f N\n\n",FD)