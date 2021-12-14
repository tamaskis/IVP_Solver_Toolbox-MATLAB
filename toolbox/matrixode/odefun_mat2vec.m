%==========================================================================
%
% ode_mat2vec  Transforms a matrix-valued ODE into a vector-valued ODE.
%
%   ydot = odefun_mat2vec(F,t,y)
%   ydot = odefun_mat2vec(F,t,y,n)
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
%   F       - (1×1 function_handle) dM/dt = F(t,M) --> multivariate, 
%             matrix-valued function (f:R×R(p×q)->R(p×q)) defining ODE
%   t       - (1×1 double) current time
%   y       - (pq×1 double) state vector at current time
%   p       - (OPTIONAL) (1×1 double) state dimension
%
% -------
% OUTPUT:
% -------
%   dydt    - (pq×1 double) state vector derivative
%
%==========================================================================
function dydt = odefun_mat2vec(F,t,y,p)
    
    % state dimension
    pq = size(y,2);

    % determine "p" if not input (assuming M is square)
    if nargin < 4
        p = sqrt(length(y));
    end

    % determines q
    q = pq/p;
    
    % reshapes pq×1 state vector into p×q state matrix
    M = reshape(y,[p,q]);

    % evaluates matrix-valued ODE
    dMdt = F(t,M);
    
    % reshapes p×q state matrix derivative into pq×1 state vector deriv.
    dydt = dMdt(:);
    
end