%==========================================================================
%
% AM_corrector  mth-order Adams-Moulton corrector.
%
%   eqn = AM_corrector(m)
%   AM_corrector(m,'print')
%
% See also AB_coefficients, AM_coefficients, AB_predictor, ABM_equations.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-03-29
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
%   print   - (char) (OPTIONAL) prints coefficients if input as 'print'
%
% -------
% OUTPUT:
% -------
%   eqn     - (string) mth-order Adams-Moulton corrector
%
%==========================================================================
function eqn = AM_corrector(m,print)
    
    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % obtains coefficients of mth-order Adams-Moulton corrector
    [num,gcd] = AM_coefficients(m);
    
    % assemble string storing equation
    if gcd == 1
        eqn = "y(n+1) = y(n) + h(";
    else
        eqn = "y(n+1) = y(n) + (h/"+gcd+")(";
    end
    for n = 1:length(num)
        if num(n) >= 0
            num_str = ""+num(n);
            if strcmp(num_str,"1"), num_str = ""; end
            if n == 1
                eqn = eqn+num_str+"f(n+1)";
            elseif n == 2
                eqn = eqn+" + "+num_str+"f(n)";
            else
                eqn = eqn+" + "+num_str+"f(n"+(2-n)+")";
            end
        else
            num_str = extractAfter(""+num(n),1);
            if strcmp(num_str,"1"), num_str = ""; end
            if n == 1
                eqn = eqn+" -"+num_str+"f(n+1)";
            elseif n == 2
                eqn = eqn+" - "+num_str+"f(n)";
            else
                eqn = eqn+" - "+num_str+"f(n"+(2-n)+")";
            end
        end
    end
    eqn = eqn+")";
    
    % prints if requested
    if strcmp(print,'print')
        if m == 1
            fprintf("1st-order Adams-Moulton corrector:\n");
        elseif m == 2
            fprintf("2nd-order Adams-Moulton corrector:\n");
        elseif m == 3
            fprintf("3rd-order Adams-Moulton corrector:\n");
        else
            fprintf("%dth-order Adams-Moulton corrector:\n",m);
        end
        fprintf(eqn+"\n\n\n");
    end
    
end