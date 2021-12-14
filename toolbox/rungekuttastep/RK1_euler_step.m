%==========================================================================
%
% RK1_euler_step  Propagates the state vector forward one time step using
% the Euler (1st-order) method.
%
%   y_next = RK1_euler_step(f,t,y,h)
%
% See also RK2_step, RK2_heun_step, RK2_ralston_step, RK3_step,
% RK3_heun_step, RK3_ralston_step, SSPRK3_step, RK4_step, RK4_ralston_step,
% RK4_38_step.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-14
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
% REFERENCES:
%   [1] https://en.wikipedia.org/wiki/Euler_method
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
%
%==========================================================================
function y_next = RK1_euler_step(f,t,y,h)
    y_next = y+h*f(t,y);
end