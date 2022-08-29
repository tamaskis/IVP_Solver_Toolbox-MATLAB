%==========================================================================
%
% AB_coefficients  Coefficients for the mth-order Adams-Bashforth 
% predictor.
%
%   [num,gcd] = AB_coefficients(m)
%   AB_coefficients(m,'print')
%
% See also AM_coefficients, AB_predictor, AM_corrector, ABM_equations.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
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
%   m       - (1×1 double) order of Adams-Bashforth predictor
%   print   - (OPTIONAL) (char) prints coefficients if input as 'print'
%
% -------
% OUTPUT:
% -------
%   num     - (m×1 double) vector storing the numerators of
%             b = (b1,...,bm)^T, where b stores the coefficients of the 
%             mth-order Adams-Bashforth predictor
%   gcd     - (1×1 double) greatest common denimonator of b
%
%==========================================================================
function [num,gcd] = AB_coefficients(m,print)
    
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
        fprintf("%dth-order Adams-Bashforth predictor coefficients:",m);
        for i = 1:length(num)
            fprintf("\nb(%d) = %d/%d",i,num(i),gcd);
        end
        fprintf("\n\n\n");
    end
    
end