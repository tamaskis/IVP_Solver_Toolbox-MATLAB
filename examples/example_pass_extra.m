%% example_pass_extra.m
% IVP Solver Toolbox
%
% Example of passing extra parameters to functions.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-06
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% EXAMPLE #1

% function definition
a = 5; b = 5; c = 5; d = 5;
f = @(x,y) -a*(x-b)^2-c*(y-d)^2;

% call function, update "a", then call function again
f(2,2)
a = 20;
f(2,2)

% NOTE: Both evaluations of "f" yield the same result.



%% EXAMPLE #2

% original function definition
a = 5; b = 5; c = 5; d = 5;
f = @(x,y) -a*(x-b)^2-c*(y-d)^2;
f(2,2)

% update values of constants
a = 10; b = 10; c = 10; d = 10;
f = @(x,y) -a*(x-b)^2-c*(y-d)^2;
f(2,2)

% NOTE: The two evaluations of "f" no longer yield the same result.



%% EXAMPLE #3

% define function where constants can vary as well
f_extra = @(x,y,a,b,c,d) -a*(x-b)^2-c*(y-d)^2;

% define original function by assigning function handle to f_extra
f = @(x,y) f_extra(x,y,5,5,5,5);
f(2,2)

% update values of constants
f = @(x,y) f_extra(x,y,10,10,10,10);
f(2,2)

% NOTE: This example is an alternate solution to Example #2.



%% EXAMPLE #4

% sets values of constants
a = 10; b = 10; c = 10; d = 10;
f = @(x,y) f_extra2(x,y,a,b,c,d);
f(2,2)

% MATLAB function must be declared at end
function f = f_extra2(x,y,a,b,c,d)
    f = -a*(x-b)^2-c*(y-d)^2;
end

% NOTE: This example is an alternate solution to the second part of 
%       Example #2.