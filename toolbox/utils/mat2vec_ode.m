%==========================================================================
%
% mat2vec_ODE  Transforms a matrix-valued ODE into a vector-valued ODE.
%
%   f = mat2vec_ODE(F,p)
%
% See also mat2vec_E, mat2vec_IC, vec2mat_sol.
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
%             defining the matrix-valued ODE, dM/dt = F(t,M) 
%             (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ)
%   p       - (1×1 double) number of rows of state matrix
%
% -------
% OUTPUT:
% -------
%   f       - (1×1 function_handle) multivariate, vector-valued function 
%             defining corresponding vector-valued ODE, dy/dt = f(t,y) 
%             (f : ℝ×ℝᵖʳ → ℝᵖʳ)
%
%==========================================================================
function f = mat2vec_ODE(F,p)
    
    % function handle for corresponding vector-valued ODE
    f = @(t,y) state_vector_derivative(F,t,y,p);
    
    %----------------------------------------------------------------------
    % state_vector_derivative  Finds the state vector derivative given a 
    % matrix-valued ODE.
    %----------------------------------------------------------------------
    %
    % ------
    % INPUT:
    % ------
    %   F       - (1×1 function_handle) multivariate, matrix-valued 
    %             function defining matrix-valued ODE, dM/dt = F(t,M) 
    %             (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ)
    %   t       - (1×1 double) current time
    %   y       - (pr×1 double) state vector at current time
    %   p       - (1×1 double) number of rows of state matrix
    %
    % -------
    % OUTPUT:
    % -------
    %   dydt    - (pr×1 double) state vector derivative
    %
    %----------------------------------------------------------------------
    function dydt = state_vector_derivative(F,t,y,p)
        
        % state dimension
        pr = length(y);
        
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