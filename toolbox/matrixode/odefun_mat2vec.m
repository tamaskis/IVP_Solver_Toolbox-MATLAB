%==========================================================================
%
% ode_mat2vec  Transforms a matrix-valued differential equation into a
% vector-valued differential equation.
%
%   ydot = odefun_mat2vec(F,t,y)
%   ydot = odefun_mat2vec(F,t,y,n)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-13
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Riccati_Differential_Equation.pdf
%
% REFERENCES:
%   [1] https://www.mathworks.com/matlabcentral/answers/94722-how-can-i-solve-the-matrix-riccati-differential-equation-within-matlab
%   [2] https://en.wikipedia.org/wiki/Linear%E2%80%93quadratic_regulator#Finite-horizon,_continuous-time_LQR
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   F       - (1×1 function_handle) dY/dt = F(t,Y) --> multivariate, 
%             matrix-valued function (f:R×R(p×q)->R(p×q)) defining ODE
%   t       - (1×1 double) current time
%   y       - (pq×1 double) state vector at current time
%   p       - (OPTIONAL) (1×1 double) state dimension
%
% -------
% OUTPUT:
% -------
%   ydot    - (pq×1 double) state vector derivative
%
%==========================================================================
function ydot = odefun_mat2vec(F,t,y,p)
    
    % determine state dimension if not input (assuming a Y is square)
    if nargin < 4
        p = sqrt(length(y));
    end

    % determines q, where Y is a p×q matrix
    q = length(y)/p;
    
    % reshapes pq×1 state vector into p×q state matrix
    Y = reshape(y,[p,q]);

    % evaluates matrix ODE
    Ydot = F(t,Y);
    
    % reshapes p×q state matrix derivative into pq×1 state vector deriv.
    ydot = Ydot(:);
    
end