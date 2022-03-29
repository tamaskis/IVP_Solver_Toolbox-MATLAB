%==========================================================================
%
% RK2_ralston_step  Propagates the state vector forward one time step using
% Ralston's second-order method (Runge-Kutta second-order method).
%
%   y_next = RK2_ralston_step(f,t,y,h)
%
% See also RK1_euler_step, RK2_step, RK2_heun_step, RK3_step, 
% RK3_heun_step, RK3_ralston_step, SSPRK3_step, RK4_step, RK4_ralston_step,
% RK4_38_step.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-03-29
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
% REFERENCES:
%   [1] https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods
%   [2] https://en.wikipedia.org/wiki/Runge-Kutta_methods
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖ → ℝᵖ) defining ODE
%   t       - (1×1 double) current sample time
%   y       - (p×1 double) state (i.e. solution) at current sample time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   y_next  - (p×1 double) state (i.e. solution) at next sample time, t+h
%
% -----
% NOTE:
% -----
%   --> p = dimension of state vector (for the scalar case, p = 1)
%   --> The documentation above is written specifically for the case of
%       vector-valued ODEs, but this function can also be used for matrix-
%       valued ODEs of the form dM/dt = F(t,M), where 
%       F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ.
%
%==========================================================================
function y_next = RK2_ralston_step(f,t,y,h)
    
    % k terms
    k1 = f(t,y);
    k2 = f(t+2*h/3,y+2*h*k1/3);
    
    % state vector propogated to next sample time
    y_next = y+(h/4)*(k1+3*k2);
    
end