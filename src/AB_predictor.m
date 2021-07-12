%==========================================================================
%
% AB_predictor  Determines and prints the nth-order Adams-Bashforth
% predictor.
%
%   eqn = AB_predictor(n)
%   AB_predictor(n,'print')
%
% Copyright © 2021 Tamas Kis
% Contact: tamas.a.kis@outlook.com
% Last Update: 2021-07-10
%
%--------------------------------------------------------------------------
%
% MATLAB Central File Exchange:
% GitHub: https://github.com/tamaskis/ode_solver_toolbox-MATLAB
%
% See [ADD LINK] for examples and "DOCUMENTATION.pdf" for additional 
% documentation.
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   n       - (1×1) order of Adams-Bashforth predictor
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
        eqn = "y(i+1) = h(";
    else
        eqn = "y(i+1) = (h/"+gcd+")(";
    end
    for i = 1:length(num)
        if num(i) >= 0
            num_str = num(i);
            if strcmp(num_str,"1"), num_str = ""; end
            if i == 1
                eqn = eqn+num_str+"f(t(i),y(i))";
            else
                eqn = eqn+" + "+num_str+"f(t(i"+(1-i)+"),y(i"+(1-i)+"))";
            end
        else
            num_str = extractAfter(""+num(i),1);
            if strcmp(num_str,"1"), num_str = ""; end
            if i == 1
                eqn = eqn+" -"+num_str+"f(t(i),y(i))";
            else
                eqn = eqn+" - "+num_str+"f(t(i"+(1-i)+"),y(i"+(1-i)+"))";
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