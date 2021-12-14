%==========================================================================
%
% odesol_vec2mat  Transforms the solution matrix for a vector-valued 
% differential equation into the solution array for matrix-valued
% differential equation.
%
%   Y = odesol_vec2mat(y)
%   Y = odesol_vec2mat(y,p)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-13
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
%   y       - ((N+1)×p double) matrix storing time history of state vector
%   p       - (OPTIONAL) (1×1 double) state dimension
%
% -------
% OUTPUT:
% -------
%   Y       - (p×q×(N+1) double) 3D array storing time history of state 
%             matrix
%
% -----
% NOTE:
% -----
%   --> [p,q] = size of the state matrix, Y
%   --> pq = dimension of the state vector, y
%   --> If Y is square, then p = √(x).
%   --> N+1 = length of time vector (i.e. number of points in time that we
%       are given y)
%   --> The ith row of "y" is the TRANSPOSE of the state vector (i.e. the
%       solution corresponding to the ith time in "t"). This convention is
%       chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function Y = odesol_vec2mat(y,p)
    
    % state vector dimension
    pq = size(y,2);

    % determine state dimension if not input (assuming Y is square)
    if nargin < 2
        p = sqrt(pq);
    end

    % determines q, where Y is a p×q matrix
    q = pq/p;
    
    % determines N, where the ODE solution is given at N+1 points in time
    N = size(y,1)-1;
    
    % preallocate array to store time history of state matrix
    Y = zeros(p,q,N+1);

    % populate state matrix time history
    for i = 1:(N+1)
        Y(:,:,i) = reshape(y(i,:),[p,q]);
    end
    
end