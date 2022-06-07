%% example_waitbar.m
% IVP Solver Toolbox
%
% Examples demonstrating how to use waitbar.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-06
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% DEFINES LORENZ SYSTEM

% Lorenz parameters
rho = 28;
sigma = 10;
beta = 8/3;

% Lorenz equations in vector form
f = @(t,x) [sigma*(x(2)-x(1));
            x(1)*(rho-x(3))-x(2);
            x(1)*x(2)-beta*x(3)];



%% EXAMPLE #1 - WAITBAR WITH DEFAULT MESSAGE

% solves IVP using RK4, displaying waitbar with default message
[t,y] = solve_ivp(f,[0,200],[0;1;1.05],0.001,[],true);



%% EXAMPLE #2 - WAITBAR WITH CUSTOM MESSAGE

% solves IVP using ABM8, displaying waitbar with custom message
[t,y] = solve_ivp(f,[0,200],[0;1;1.05],0.001,'ABM8',...
    'Solving Lorenz system...');