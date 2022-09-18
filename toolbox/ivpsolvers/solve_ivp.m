%==========================================================================
%
% solve_ivp  Solve vector-valued and matrix-valued initial value problems 
% using fixed-step IVP solvers.
%
%   [t,y] = solve_ivp(f,[t0,tf],y0,h)
%   [t,y] = solve_ivp(f,{t0,E},y0,h)
%   [t,y] = solve_ivp(__,method)
%   [t,y] = solve_ivp(__,method,wb)
%
%   [t,M] = solve_ivp(F,[t0,tf],M0,h)
%   [t,M] = solve_ivp(F,{t0,E},M0,h)
%   [t,M] = solve_ivp(__,method)
%   [t,M] = solve_ivp(__,method,wb)
%
% See also solve_ivp_matrix, solve_ivp_vector.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-09-17
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
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% VECTOR-VALUED INITIAL VALUE PROBLEMS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖ → ℝᵖ) defining ODE
%   I       - defines interval over which to solve the IVP, 2 options:
%               --> [t0,tf] - (1×2 double) initial and final times
%               --> {t0,E}  - (1×2 cell) initial time, t₀, and function 
%                             handle for event function, E(t,y) 
%                             (E : ℝ×ℝᵖ → B)
%   y0      - (p×1 double) initial condition, y₀ = y(t₀)
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
%   y       - ((N+1)×p double) solution matrix
%
% -----
% NOTE:
% -----
%   --> The nth row of "y" stores the TRANSPOSE of the state vector (i.e. 
%       the solution) corresponding to the nth time in "t". This convention
%       is chosen to match the convention used by MATLAB's ODE suite.
%
%--------------------------------------------------------------------------
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% MATRIX-VALUED INITIAL VALUE PROBLEMS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ------
% INPUT:
% ------
%   F       - (1×1 function_handle) dM/dt = F(t,M) --> multivariate, 
%             matrix-valued function (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ) defining
%             matrix-valued ODE
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
function [t,y] = solve_ivp(f,I,y0,h,method,wb)
    
    % defaults optional parameters to empty vector if not input
    if (nargin < 5), method = []; end
    if (nargin < 6), wb = []; end

    % number of columns of the state matrix
    r = size(y0,2);
    
    % determines if the IVP is matrix-valued
    matrix_valued = (size(y0,2) > 1);
    
    % solves IVP
    if matrix_valued
        [t,y] = solve_ivp_matrix(f,I,y0,h,method,wb);
    else
        [t,y] = solve_ivp_vector(f,I,y0,h,method,wb);
    end
    
end