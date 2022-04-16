%==========================================================================
%
% odeIC_mat2vec  Transforms the initial condition for a matrix-valued ODE
% into the initial condition for the corresponding vector-valued ODE.
%
%   y0 = odeIC_mat2vec(M0)
%
% See also odefun_mat2vec, odesol_vec2mat.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-04-16
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/ODE_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   M0      - (p×r double) initial condition for matrix-valued ODE
%
% -------
% OUTPUT:
% -------
%   y0      - (pr×1 double) initial condition for corresponding
%             vector-valued ODE
%
%==========================================================================
function y0 = odeIC_mat2vec(M0)
    y0 = M0(:);
end