%==========================================================================
%
% mat2vec_E  Transforms the event function for a matrix-valued IVP into the
% event function for its corresponding vector-valued IVP.
%
%   Ev = mat2vec_E(Em,p)
%
% See also mat2vec_IC, mat2vec_ODE, vec2mat_sol.
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
%   Em      - (1×1 function_handle) event function for matrix-valued 
%             IVP, Eₘ(t,M) (Eₘ : ℝ×ℝᵖˣʳ → B)
%   p       - (1×1 double) number of rows of state matrix
%
% -------
% OUTPUT:
% -------
%   Ev      - (1×1 function_handle) event function for corresponding
%             vector-valued IVP, Eᵥ(t,y) (Eᵥ : ℝ×ℝᵖʳ → B)
%
%==========================================================================
function Ev = mat2vec_E(Em,p)
    
    % function handle for event function for corresponding 
    % vector-valued IVP
    Ev = @(t,y) vector_event_function(Em,t,y,p);
    
    %----------------------------------------------------------------------
    % state_vector_derivative  Evaluates the event function for a
    % matrix-valued IVP given the current time and corresponding state
    % vector.
    %----------------------------------------------------------------------
    %
    % ------
    % INPUT:
    % ------
    %   Em      - (1×1 function_handle) event function for matrix-
    %              valued IVP, Eₘ(t,M) (Eₘ : ℝ×ℝᵖˣʳ → B)
    %   t       - (1×1 double) current time
    %   y       - (pr×1 double) state vector at current time
    %   p       - (1×1 double) number of rows of state matrix
    %
    % -------
    % OUTPUT:
    % -------
    %   E       - (1×1 logical) evaluation of event function
    %----------------------------------------------------------------------
    function E = vector_event_function(Em,t,y,p)
        
        % state dimension
        pr = length(y);
        
        % determines r
        r = pr/p;
        
        % reshapes pr×1 state vector into p×r state matrix
        M = reshape(y,[p,r]);
        
        % evaluates event function
        E = Em(t,M);
        
    end
    
end