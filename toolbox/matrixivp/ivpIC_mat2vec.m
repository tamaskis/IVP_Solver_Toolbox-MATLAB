%==========================================================================
%
% ivpIC_mat2vec  Transforms the initial condition for a matrix-valued IVP 
% into the initial condition for the corresponding vector-valued IVP.
%
%   y0 = ivpIC_mat2vec(M0)
%
% See also ivpsol_vec2mat, odefun_mat2vec.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-07
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   M0      - (p×r double) initial condition for matrix-valued IVP
%
% -------
% OUTPUT:
% -------
%   y0      - (pr×1 double) initial condition for corresponding
%             vector-valued IVP
%
%==========================================================================
function y0 = ivpIC_mat2vec(M0)
    y0 = M0(:);
end