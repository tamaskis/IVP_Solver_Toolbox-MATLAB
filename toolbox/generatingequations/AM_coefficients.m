%==========================================================================
%
% AM_coefficients  Coefficients for the mth-order Adams-Moulton corrector.
%
%   [num,gcd] = AM_coefficients(m,print)
%   AM_coefficients(m,'print')
%
% See also AB_coefficients, AB_predictor, AM_corrector, ABM_equations.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-22
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
% REFERENCES:
%   [1] Seidu, "A Matrix System for Computing the Coefficients of the Adams
%       Bashforth-Moulton Predictor-Corrector formulae"
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   m       - (1×1 double) order of Adams-Moulton corrector
%   print   - (OPTIONAL) (char) prints coefficients if input as 'print'
%
% -------
% OUTPUT:
% -------
%   num     - (m×1 double) vector storing the numerators of 
%             b = (b1,...,bm)^T, where b stores the coefficients of the 
%             mth-order Adams-Moulton predictor
%   gcd     - (1×1 double) greatest common denimonator of b
%
%==========================================================================
function [num,gcd] = AM_coefficients(m,print)
    
    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % preallocates A matrix and c vector
    A = zeros(m);
    c = ones(m,1);
    
    % populates A matrix and c vector
    for i = 1:m
        for j = 1:m
            A(i,j) = (j-2)^(i-1);
        end
        c(i) = (-1)^(i-1)/i;
    end
    
    % coefficients
	b = A\c;
    
    % converts coefficients so they all have same denominator
    [num,gcd] = same_denominator(b);
    
    % prints if requested
    if strcmp(print,'print')
        fprintf("%dth-order Adams-Moulton corrector coefficients:",m);
        for i = 1:length(num)
            fprintf("\nb(%d) = %d/%d",i,num(i),gcd);
        end
        fprintf("\n\n\n");
    end
    
end