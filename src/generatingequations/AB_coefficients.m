%==========================================================================
%
% AB_coefficients  Coefficients for the nth-order Adams-Bashforth 
% predictor.
%
%   [num,gcd] = AB_coefficients(n)
%   AB_coefficients(n,'print')
%
% See also AM_coefficients, AB_predictor, AM_corrector, ABM_equations.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-09-06
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   n       - (1×1 double) order of Adams-Bashforth predictor
%   print   - (OPTIONAL) (char) prints coefficients if input as 'print'
%
% -------
% OUTPUT:
% -------
%   num     - (n×1 double) vector storing the numerators of
%             b = (b1,...,bn)^T, where b stores the coefficients of the 
%             nth-order Adams-Bashforth predictor
%   gcd     - (1×1 double) greatest common denimonator of b
%
%==========================================================================
function [num,gcd] = AB_coefficients(n,print)
    
    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % preallocates A matrix and c vector
    A = zeros(n);
    c = ones(n,1);
    
    % populates A matrix and c vector
    for i = 1:n
        for j = 1:n
            A(i,j) = (j-1)^(i-1);
        end
        c(i,1) = (-1)^(i-1)/i;
    end
    
    % coefficients
	b = A\c;
    
    % converts coefficients so they all have same denominator
    [num,gcd] = same_denominator(b);
    
    % prints if requested
    if strcmp(print,'print')
        fprintf("%dth-order Adams-Bashforth predictor coefficients:",n);
        for i = 1:length(num)
            fprintf("\nb(%d) = %d/%d",i,num(i),gcd);
        end
        fprintf("\n\n\n");
    end
    
end