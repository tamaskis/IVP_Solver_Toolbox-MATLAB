%==========================================================================
%
% odefun_mat2vec  Transforms a matrix-valued ODE into a vector-valued ODE.
%
%   f = odefun_mat2vec(F)
%   f = odefun_mat2vec(F,p)
%
% See also odeIC_mat2vec, odesol_vec2mat.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-03-29
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   F       - (1×1 function_handle) dM/dt = F(t,M) --> multivariate, 
%             matrix-valued function (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ) defining
%             matrix-valued ODE
%   p       - (1×1 double) (OPTIONAL) number of rows of state matrix
%
% -------
% OUTPUT:
% -------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖʳ → ℝᵖʳ) defining
%             corresponding vector-valued ODE
%
%==========================================================================
function f = odefun_mat2vec(F,p)
    
    % assigns function handle to subfunction below that produces the state
    % vector derivative, given the matrix-valued ODE
    if nargin == 2
        f = @(t,y) state_vector_derivative(F,t,y,p);
    else
        f = @(t,y) state_vector_derivative(F,t,y);
    end
    
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
    %   p       - (1×1 double) (OPTIONAL) number of rows of state matrix
    %
    % OUTPUT:
    %   dydt    - (pr×1 double) state vector derivative
    %
    %----------------------------------------------------------------------
    function dydt = state_vector_derivative(F,t,y,p)

        % state dimension
        pr = size(y,1);
        
        % determine "p" if not input (assuming M is square)
        if nargin < 4
            p = sqrt(length(y));
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