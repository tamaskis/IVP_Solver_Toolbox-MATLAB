%% example_matrix_ode.m
% ODE Solver Toolbox
%
% Example for solving a matrix-valued ODE (the Riccati differential eqn.).
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



%% DEFINING SYSTEM

% state matrix
A = [1   1;
     2   1];

% input matrix
B = [1;
     1];

% cross-coupling weighting matrix
N = [0  0;
     0  0];

% state weighting matrix
Q = [2   1;
     1   1];

% input weighting matrix
R = [1   0;
     0   1];

% final condition
PT = [1   1;
      1   1];

% final time
T = 5;



%% SOLVING RICCATI DIFFERENTIAL EQUATION

% defines the Riccati differential equation (a matrix-valued ODE)
F = @(t,P) -(A.'*P+P*A-(P*B+N)/R*(B.'*P+N.')+Q);

% converts the matrix-valued ODE to a vector-valued ODE
f = @(t,y) odefun_mat2vec(F,t,y);

% final condition
yT = odeIC_mat2vec(PT);

% solves vector-valued ODE
[t,y] = ode45(f,[T,0],yT);

% transforms solution matrix for vector-valued ODE into solution array for
% matrix-valued ODE
P = odesol_vec2mat(y);

% solution for P0 (will be at end of array since P solved for backwards in
% time)
tic
P0 = P(:,:,end);
toc



%% SOLVING RICCATI DIFFERENTIAL EQUATION

% defines the Riccati differential equation (a matrix-valued ODE)
F = @(t,P) -(A.'*P+P*A-(P*B+N)/R*(B.'*P+N.')+Q);

% converts the matrix-valued ODE to a vector-valued ODE
f = odefun_mat2vec_new(F);

% final condition
yT = odeIC_mat2vec(PT);

% solves vector-valued ODE
[t,y] = ode45(f,[T,0],yT);

% transforms solution matrix for vector-valued ODE into solution array for
% matrix-valued ODE
P = odesol_vec2mat(y);

% solution for P0 (will be at end of array since P solved for backwards in
% time)
tic
P0 = P(:,:,end);
toc