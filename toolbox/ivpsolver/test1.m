%%
clear;clc;close all;
addpath(genpath('..'));
%%
f = @(t,y) y;
y2 = 3;
h = 0.01;
t = (2:h:10)';
y = zeros(size(t));
y(1) = y2;
for i = 1:(length(t)-1)
    y(i+1) = RK1_euler(f,t(i),y(i),h);
end
[t_ivp,y_ivp] = solve_ivp(f,[2,10],y2,h,'Euler');
[t_ivp,y_ivp] = solve_ivp(f,[10,2],y2,h,'Euler');
%max(abs(y_ivp-y))
t_ivp