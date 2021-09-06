%==========================================================================
%
% AB_predictor  Determines and prints the nth-order Adams-Bashforth
% predictor.
%
%   eqn = AB_predictor(n)
%   AB_predictor(n,'print')
%
% See also AB_coefficients, AM_coefficients, AM_corrector, ABM_equations.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-08-23
% Website: tamaskis.github.io
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
%   eqn     - (string) nth-order Adams-Bashforth predictor
%
%==========================================================================
function eqn = AB_predictor(n,print)
    
    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % obtains coefficients of nth-order Adams-Bashforth predictor
    [num,gcd] = AB_coefficients(n);
    
    % assemble string storing equation
    if gcd == 1
        eqn = "y(i+1) = y(i) + h(";
    else
        eqn = "y(i+1) = y(i) + (h/"+gcd+")(";
    end
    for i = 1:length(num)
        if num(i) >= 0
            num_str = ""+num(i);
            if strcmp(num_str,"1"), num_str = ""; end
            if i == 1
                eqn = eqn+num_str+"f(i)";
            else
                eqn = eqn+" + "+num_str+"f(i"+(1-i)+")";
            end
        else
            num_str = extractAfter(""+num(i),1);
            if strcmp(num_str,"1"), num_str = ""; end
            if i == 1
                eqn = eqn+" -"+num_str+"f(i)";
            else
                eqn = eqn+" - "+num_str+"f(i"+(1-i)+")";
            end
        end
    end
    eqn = eqn+")";
    
    % prints if requested
    if strcmp(print,'print')
        if n == 1
            fprintf("1st-order Adams-Bashforth predictor:\n");
        elseif n == 2
            fprintf("2nd-order Adams-Bashforth predictor:\n");
        elseif n == 3
            fprintf("3rd-order Adams-Bashforth predictor:\n");
        else
            fprintf("%dth-order Adams-Bashforth predictor:\n",n);
        end
        fprintf(eqn+"\n\n\n");
    end
    
end