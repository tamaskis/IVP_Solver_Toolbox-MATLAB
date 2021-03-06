%==========================================================================
%
% ivpsol_vec2mat  Transforms the solution matrix for a vector-valued IVP
% into the solution array for the corresponding matrix-valued IVP.
%
%   M = ivpsol_vec2mat(y)
%   M = ivpsol_vec2mat(y,p)
%
% See also ivpIC_mat2vec, odefun_mat2vec.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-07-06
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
%   y       - ((N+1)×pr double) matrix storing time history of state vector
%   p       - (OPTIONAL) (1×1 double) number of rows of state matrix
%
% -------
% OUTPUT:
% -------
%   M       - (p×r×(N+1) double) 3D array storing time history of state 
%             matrix
%
% -----
% NOTE:
% -----
%   --> [p,r] = size of the state matrix, M
%   --> pr = dimension of the state vector, y
%   --> If M is square, then p = √(x).
%   --> N+1 = length of time vector (i.e. number of points in time that we
%       are given y)
%   --> The ith row of "y" is the TRANSPOSE of the state vector (i.e. the
%       solution corresponding to the ith time in "t"). This convention is
%       chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function M = ivpsol_vec2mat(y,p)
    
    % state dimension (pr) and N (IVP solution given at N+1 times)
    pr = size(y,2);
    N = size(y,1)-1;
    
    % determine state dimension if not input (assuming Y is square)
    if nargin < 2
        p = sqrt(pr);
    end
    
    % determines r
    r = pr/p;
    
    % preallocate solution array to store time history of state matrix
    M = zeros(p,r,N+1);
    
    % populate solution array
    for n = 1:(N+1)
        M(:,:,n) = reshape(y(n,:),[p,r]);
    end
    
end