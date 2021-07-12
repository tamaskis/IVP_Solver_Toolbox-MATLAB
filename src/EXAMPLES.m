%% Copyright (c) 2021 Tamas Kis

% Examples for using the AM_coefficients function.



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

%addpath('src/');



%% EXAMPLE

% finds the coefficients for the 3rd-order Adams-Moulton corrector
%f = @(t,x) 2*x;
f = @(t,x) func(t,x);
x0 = 1;

C = @(t,x) x < 10;

%[t1,x1] = RK2(f,[0,5],x0,0.01);
%[t1,x1] = RK2(f,{0,C},x0,0.01);

%[t1,x1] = RK4(f,[0,5],x0,0.06);
%[t1,x1] = RK4(f,{0,C},x0,0.01);

%[t1,x1] = ABM8(f,[0,5],x0,0.01);
%[t1,x1] = ABM8(f,{0,C},x0,0.01);

%[t2,x2] = ode45(f,[0,5],x0);

tic
[t1,x1] = ABM8(f,[0,5],x0,0.01);
toc
x1(end,:)

tic
[t1,x1] = RK4(f,[0,5],x0,0.01);
toc
x1(end,:)

hold on;
plot(t1,x1);
%plot(t2,x2);
hold off;


function f = func(t,y)
    for i = 1:200
        y = y*1.01*sin(pi/2+0.00001);
        y = y+tan(y)^2;
        y = y-tan(y)^2;
    end
    f = y*1.01'/10;
end