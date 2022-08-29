%==========================================================================
%
% mat2vec_C  Transforms the condition function for a matrix-valued IVP into
% the condition function for its corresponding vector-valued IVP.
%
%   Cv = mat2vec_C(Cm)
%   Cv = mat2vec_C(Cm,p)
%
% See also mat2vec_ode, mat2vec_IC, vec2mat_sol.
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
%   Cm      - (1×1 function_handle) condition function for matrix-valued 
%             IVP, Cₘ(t,M) (Cₘ : ℝ×ℝᵖˣʳ → B)
%   p       - (OPTIONAL) (1×1 double) number of rows of state matrix
%
% -------
% OUTPUT:
% -------
%   Cv      - (1×1 function_handle) condition function for corresponding
%             vector-valued IVP, Cᵥ(t,y) (Cᵥ : ℝ×ℝᵖʳ → B)
%
% -----
% NOTE:
% -----
%   --> If "p" is not input, it is assumed that the state matrix (M) is
%       a square matrix.
%
%==========================================================================
function Cv = mat2vec_C(Cm,p)
    
    % defaults "p" to empty vector if not input
    if (nargin < 2), p = []; end
    
    % function handle for condition function for corresponding 
    % vector-valued IVP
    Cv = @(t,y) vector_condition_function(Cm,t,y,p);
    
    %----------------------------------------------------------------------
    % state_vector_derivative  Evaluates the condition function for a
    % matrix-valued IVP given the current time and corresponding state
    % vector.
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   Cm      - (1×1 function_handle) condition function for matrix-
    %              valued IVP, Cₘ(t,M) (Cₘ : ℝ×ℝᵖˣʳ → B)
    %   t       - (1×1 double) current time
    %   y       - (pr×1 double) state vector at current time
    %   p       - (1×1 double) number of rows of state matrix
    %
    % OUTPUT:
    %   C       - (1×1 logical) evaluation of condition function
    %----------------------------------------------------------------------
    function C = vector_condition_function(Cm,t,y,p)
        
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
        
        % evaluates condition function
        C = Cm(t,M);
        
    end
    
end