%==========================================================================
%
% AB_predictor  Determines and prints the mth-order Adams-Bashforth
% predictor.
%
%   eqn = AB_predictor(m)
%   AB_predictor(m,'print')
%
% See also AB_coefficients, AM_coefficients, AM_corrector, ABM_equations.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-12
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
%   m       - (1×1 double) order of Adams-Bashforth predictor
%   print   - (OPTIONAL) (char) prints coefficients if input as 'print'
%
% -------
% OUTPUT:
% -------
%   eqn     - (string) mth-order Adams-Bashforth predictor
%
%==========================================================================
function eqn = AB_predictor(m,print)
    
    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % obtains coefficients of mth-order Adams-Bashforth predictor
    [num,gcd] = AB_coefficients(m);
    
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
                eqn = eqn+num_str+"f(n)";
            else
                eqn = eqn+" + "+num_str+"f(n"+(1-n)+")";
            end
        else
            num_str = extractAfter(""+num(n),1);
            if strcmp(num_str,"1"), num_str = ""; end
            if n == 1
                eqn = eqn+" -"+num_str+"f(n)";
            else
                eqn = eqn+" - "+num_str+"f(n"+(1-n)+")";
            end
        end
    end
    eqn = eqn+")";
    
    % prints if requested
    if strcmp(print,'print')
        if m == 1
            fprintf("1st-order Adams-Bashforth predictor:\n");
        elseif m == 2
            fprintf("2nd-order Adams-Bashforth predictor:\n");
        elseif m == 3
            fprintf("3rd-order Adams-Bashforth predictor:\n");
        else
            fprintf("%dth-order Adams-Bashforth predictor:\n",m);
        end
        fprintf(eqn+"\n\n\n");
    end
    
end