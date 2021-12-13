%==========================================================================
%
% RK3_ralston_step  Propagates the state vector forward one time step using 
% Ralston's third-order method (Runge-Kutta third-order method).
%
%   y_next = RK3_ralston_step(f,t,y,h)
%
% See also euler_step, RK2_step, RK2_heun_step, RK2_ralston_step, RK3_step,
% RK3_heun_step, SSPRK3_step, RK4_step, RK4_ralston_step, RK4_38_step.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-12
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
% REFERENCES:
%   [1] https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> function defining
%             ODE (f:R×Rp->Rp)
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
function y_next = RK3_ralston_step(f,t,y,h)
    
    % k terms
    k1 = f(t,y);
    k2 = f(t+h/2,y+h*k1/2);
    k3 = f(t+3*h/4,y+3*h*k2/4);
    
    % state vector propogated to next sample time
    y_next = y+(h/9)*(2*k1+3*k2+4*k3);
    
end