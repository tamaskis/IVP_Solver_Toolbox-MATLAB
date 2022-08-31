%==========================================================================
%
% vec2mat_sol  Transforms the solution matrix for a vector-valued IVP into
% the solution array for the corresponding matrix-valued IVP.
%
%   M = vec2mat_sol(y)
%   M = vec2mat_sol(y,p)
%
% See also mat2vec_ode, mat2vec_IC, mat2vec_C.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-30
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
%   --> If "p" is not input, it is assumed that the state matrix (M) is
%       a square matrix.
%   --> M ∈ ℝᵖˣʳ
%   --> y ∈ ℝᵖʳ
%   --> t ∈ ℝᴺ⁺¹
%   --> The nth page of "M" stores the state matrix (i.e. the solution)
%       corresponding to the nth time in "t".
%
%==========================================================================
function M = vec2mat_sol(y,p)
    
    % state dimension (pr) and N (IVP solution given at N+1 times)
    pr = size(y,2);
    N = size(y,1)-1;
    
    % determine state dimension if not input (assuming Y is square)
    if (nargin < 2) || isempty(p)
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