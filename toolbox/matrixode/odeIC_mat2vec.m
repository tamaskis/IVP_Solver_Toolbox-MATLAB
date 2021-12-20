%==========================================================================
%
% odeIC_mat2vec  Transforms the initial condition for a matrix-valued ODE
% into the initial condition for the corresponding vector-valued ODE.
%
%   y0 = odeIC_mat2vec(M0)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-20
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
%   M0      - (p×q double) initial condition for matrix-valued ODE
%
% -------
% OUTPUT:
% -------
%   y0      - (pq×1 double) initial condition for corresponding
%             vector-valued ODE
%
%==========================================================================
function y0 = odeIC_mat2vec(M0)
    y0 = M0(:);
end