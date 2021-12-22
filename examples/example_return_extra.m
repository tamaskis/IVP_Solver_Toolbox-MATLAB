%% example_return_extra.m
% ODE Solver Toolbox
%
% Example of returning extra parameters from an ODE solution.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-22
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% SOLUTION

% assign function handle and pass extra parameter
f = @(t,y) f_extra(t,y,5);

% solve ODE
y0 = [1;
      1];
t0 = 0;
tf = 10;
[t,y] = ode45(f,[t0,tf],y0);

% recover "x" at every time step
x = zeros(size(t));
for i = 1:length(t)
    [~,x(i)] = f(t(i),y(i,:)');
end

% plot x vs. t
figure;
plot(t,x,'linewidth',1.5);
grid on;
xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$x$','Interpreter','latex','FontSize',18);

% MATLAB function must be declared at end
function [f,x] = f_extra(t,y,a)
    x = 2*a-5*t;
    f = x*y;
end