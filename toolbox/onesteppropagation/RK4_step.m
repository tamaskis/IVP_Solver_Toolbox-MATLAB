%==========================================================================
%
% RK4_step  Propagates the state vector forward one time step using the
% (classic) Runge-Kutta fourth-order method.
%
%   y_next = RK4_step(f,t,y,h)
%
% See also RK1_euler_step, RK2_step, RK2_heun_step, RK2_ralston_step,
% RK3_step, RK3_heun_step, RK3_ralston_step, SSPRK3_step, RK4_ralston_step,
% RK4_38_step.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-22
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
% REFERENCES:
%   [1] Burden and Faires, "Numerical Analysis", 9th Ed. (p. 288)
%   [2] Griffiths and Higham, "Numerical Methods for Ordinary Differential 
%       Equations: Initial Value Problems" (p. 131)
%   [3] https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f:R×Rp->Rp) defining ODE
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
%       valued ODEs of the form dM/dt = F(t,M), where F:R×R(p×q)-->R(p×q).
%
%==========================================================================
function y_next = RK4_step(f,t,y,h)
    
    % k terms
    k1 = f(t,y);
    k2 = f(t+h/2,y+h*k1/2);
    k3 = f(t+h/2,y+h*k2/2);
    k4 = f(t+h,y+h*k3);
    
    % state vector propogated to next sample time
    y_next = y+(h/6)*(k1+2*k2+2*k3+k4);
    
end