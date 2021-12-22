%% example_waitbar.m
% ODE Solver Toolbox
%
% Examples demonstrating how to use waitbar.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-22
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to all "ODE Solver Toolbox" functions
addpath(genpath('..'));



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

% solves ODE, displaying waitbar with default message
[t,y] = ABM8(f,[0,200],[0;1;1.05],0.001,true);



%% EXAMPLE #2 - WAITBAR WITH CUSTOM MESSAGE

% solves ODE, displaying waitbar with custom message
[t,y] = ABM8(f,[0,200],[0;1;1.05],0.001,'Solving Lorenz system...');