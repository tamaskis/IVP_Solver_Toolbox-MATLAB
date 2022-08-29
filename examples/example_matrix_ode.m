%% example_matrix_ode.m
% IVP Solver Toolbox
%
% Example for solving a matrix-valued IVP (the Riccati differential
% equation).
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% DEFINING SYSTEM

% state matrix
A = [1   1;
     2   1];

% input matrix
B = [1;
     1];

% state weighting matrix
Q = [2   1;
     1   1];

% input weighting matrix
R = 1;

% cross-coupling weighting matrix
S = [0;
     0];

% final condition
PT = [1   1;
      1   1];

% final time
T = 5;



%% SOLVING RICCATI DIFFERENTIAL EQUATION USING IVP SOLVER

% defines the Riccati differential equation (a matrix-valued ODE)
F = @(t,P) -(A.'*P+P*A-(P*B+S)/R*(B.'*P+S.')+Q);

% solves matrix-valued IVP using a step size of h = 0.001
[~,P] = solve_ivp_matrix(F,[T,0],PT,0.001,[],'RK4');

% solution for P0 (will be at end of array since P solved for backwards in
% time)
P0 = P(:,:,end)



%% SOLVING RICCATI DIFFERENTIAL EQUATION USING ONE-STEP PROPAGATION

% time vector between t = 5 and t = 0 with a spacing of h = 0.001.
h = -0.001;
t = (5:h:0)';

% preallocate vector to store solution
P = zeros(2,2,length(t));

% store initial condition
P(:,:,1) = PT;

% solving using "RK4"
for i = 1:(length(t)-1)
    P(:,:,i+1) = RK4(F,t(i),P(:,:,i),h);
end

% solution for P0 using one-step propagation
P0_step = P(:,:,end);

% maximum absolute error between the two results (should be 0)
max(abs(P0-P0_step),[],'all')