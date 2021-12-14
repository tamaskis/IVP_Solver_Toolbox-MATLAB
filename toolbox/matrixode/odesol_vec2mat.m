%==========================================================================
%
% odesol_vec2mat  Transforms the solution matrix for a vector-valued 
% differential equation into the solution array for the corresponding 
% matrix-valued differential equation.
%
%   Y = odesol_vec2mat(y)
%   Y = odesol_vec2mat(y,p)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-14
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
%   M       - (p×q×(N+1) double) 3D array storing time history of state 
%             matrix
%
% -----
% NOTE:
% -----
%   --> [p,q] = size of the state matrix, M
%   --> pq = dimension of the state vector, y
%   --> If M is square, then p = √(x).
%   --> N+1 = length of time vector (i.e. number of points in time that we
%       are given y)
%   --> The ith row of "y" is the TRANSPOSE of the state vector (i.e. the
%       solution corresponding to the ith time in "t"). This convention is
%       chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function M = odesol_vec2mat(y,p)
    
    % state dimension (pq) and N (ODE solution given at N+1 times)
    pq = size(y,2);
    N = size(y,1)-1;

    % determine state dimension if not input (assuming Y is square)
    if nargin < 2
        p = sqrt(pq);
    end

    % determines q
    q = pq/p;
    
    % preallocate solution array to store time history of state matrix
    M = zeros(p,q,N+1);

    % populate solution array
    for n = 1:(N+1)
        M(:,:,n) = reshape(y(n,:),[p,q]);
    end
    
end