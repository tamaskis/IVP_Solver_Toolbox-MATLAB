%==========================================================================
%
% mat2vec_fun  Transforms a matrix-valued ODE into a vector-valued ODE.
%
%   f = mat2vec_fun(F)
%   f = mat2vec_fun(F,p)
%
% See also mat2vec_IC, mat2vec_C, vec2mat_sol.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
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
%   F       - (1×1 function_handle) dM/dt = F(t,M) --> multivariate, 
%             matrix-valued function (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ) defining
%             matrix-valued ODE
%   p       - (OPTIONAL) (1×1 double) number of rows of state matrix
%
% -------
% OUTPUT:
% -------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖʳ → ℝᵖʳ) defining
%             corresponding vector-valued ODE
%
% -----
% NOTE:
% -----
%   --> If "p" is not input, it is assumed that the state matrix (M) is
%       a square matrix.
%
%==========================================================================
function f = mat2vec_fun(F,p)
    
    % defaults "p" to empty vector if not input
    if (nargin < 2), p = []; end

    % function handle for corresponding vector-valued ODE
    f = @(t,y) state_vector_derivative(F,t,y,p);
    
    %----------------------------------------------------------------------
    % state_vector_derivative  Finds the state vector derivative given a 
    % matrix-valued ODE.
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   F       - (1×1 function_handle) dM/dt = F(t,M) --> multivariate, 
    %             matrix-valued function (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ) defining ODE
    %   t       - (1×1 double) current time
    %   y       - (pr×1 double) state vector at current time
    %   p       - (1×1 double) number of rows of state matrix
    %
    % OUTPUT:
    %   dydt    - (pr×1 double) state vector derivative
    %
    %----------------------------------------------------------------------
    function dydt = state_vector_derivative(F,t,y,p)
        
        % state dimension
        pr = length(y);
        
        % determine "p" if not specified (assuming M is square)
        if isempty(p)
            p = sqrt(pr);
        end
        
        % determines r
        r = pr/p;
        
        % reshapes pr×1 state vector into p×r state matrix
        M = reshape(y,[p,r]);
        
        % evaluates matrix-valued ODE
        dMdt = F(t,M);
        
        % reshapes p×r state matrix derivative into pr×1 state vector
        % derivative
        dydt = dMdt(:);
        
    end
    
end