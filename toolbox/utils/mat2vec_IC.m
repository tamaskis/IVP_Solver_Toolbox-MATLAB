%==========================================================================
%
% mat2vec_IC  Transforms the initial condition for a matrix-valued IVP into
% the initial condition for its corresponding vector-valued IVP.
%
%   y0 = mat2vec_IC(M0)
%
% See also mat2vec_E, mat2vec_ODE, vec2mat_sol.
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
% ------
% INPUT:
% ------
%   M0      - (p×r double) initial condition for matrix-valued IVP, 
%             M₀ = M(t₀)
%
% -------
% OUTPUT:
% -------
%   y0      - (pr×1 double) initial condition for corresponding
%             vector-valued IVP, y₀ = y(t₀)
%
%==========================================================================
function y0 = mat2vec_IC(M0)
    y0 = M0(:);
end