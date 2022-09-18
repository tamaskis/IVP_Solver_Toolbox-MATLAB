%==========================================================================
%
% solve_ivp_matrix  Fixed-step IVP solvers for solving matrix-valued
% initial value problems.
%
%   [t,M] = solve_ivp_matrix(F,[t0,tf],M0,h)
%   [t,M] = solve_ivp_matrix(F,{t0,E},M0,h)
%   [t,M] = solve_ivp_matrix(__,method)
%   [t,M] = solve_ivp_matrix(__,method,wb)
%
% See also solve_ivp.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-09-18
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   F       - (1×1 function_handle) multivariate, matrix-valued function
%             defining matrix-valued ODE, dM/dt = F(t,M) 
%             (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ)
%   I       - defines interval over which to solve the IVP, 2 options:
%               --> [t0,tf] - (1×2 double) initial and final times
%               --> {t0,E}  - (1×2 cell) initial time, t₀, and function 
%                             handle for event function, E(t,M) 
%                             (E : ℝ×ℝᵖˣʳ → B)
%   M0      - (p×r double) initial condition, M₀ = M(t₀)
%   h       - (1×1 double) step size
%   method  - (OPTIONAL) (char) integration method --> 'Euler', 'RK2', 
%             'RK2 Heun', 'RK2 Ralston', 'RK3', 'RK3 Heun', 'RK3 Ralston', 
%             'SSPRK3', 'RK4', 'RK4 Ralston', 'RK4 3/8', 'AB2', 'AB3', 
%             'AB4', 'AB5', 'AB6', 'AB7', 'AB8', 'ABM2', 'ABM3', 'ABM4', 
%             'ABM5', 'ABM6', 'ABM7', 'ABM8' (defaults to 'RK4')
%   wb      - (OPTIONAL) (1×1 logical or char) waitbar parameters (defaults
%             to false)
%               --> input as true if you want waitbar with default message 
%                   displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%
% -------
% OUTPUT:
% -------
%   t       - ((N+1)×1 double) time vector
%   M       - (p×r×(N+1) double) solution array
%
% -----
% NOTE:
% -----
%   --> The nth page of "M" stores the state matrix (i.e. the solution)
%       corresponding to the nth time in "t".
%
%==========================================================================
function [t,M] = solve_ivp_matrix(F,I,M0,h,method,wb)
    
    % defaults optional parameters to empty vector if not input
    if (nargin < 5), method = []; end
    if (nargin < 6), wb = []; end
    
    % number of rows of state matrix
    p = size(M0,1);
    
    % converts matrix-valued ODE to vector-valued ODE
    f = mat2vec_ODE(F,p);
    
    % converts initial condition for matrix-valued IVP into initial 
    % condition for the corresponding vector-valued IVP
    y0 = mat2vec_IC(M0);
    
    % converts event function for matrix-valued IVP into event function for
    % the corresponding vector-valued IVP
    if iscell(I), I(2) = mat2vec_E(I(2),p); end
    
    % solves corresponding vector-valued IVP
    [t,y] = solve_ivp(f,I,y0,h,method,wb);
    
    % transforms solution matrix for vector-valued IVP into solution array
    % for matrix-valued IVP
    M = vec2mat_sol(y,p);
    
end