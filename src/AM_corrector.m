%==========================================================================
%
% AM_corrector  nth-order Adams-Moulton corrector.
%
%   eqn = AM_corrector(n)
%   AM_corrector(n,'print')
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
%   n       - (1×1) order of Adams-Moulton corrector
%   print   - (OPTIONAL) (char) prints coefficients if input as 'print'
%
% -------
% OUTPUT:
% -------
%   eqn     - (string) nth-order Adams-Moulton corrector
%
%==========================================================================
function eqn = AM_corrector(n,print)
    
    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % obtains coefficients of nth-order Adams-Moulton corrector
    [num,gcd] = AM_coefficients(n);
    
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
                eqn = eqn+num_str+"f(t(i+1),y(i+1))";
            elseif i == 2
                eqn = eqn+" + "+num_str+"f(t(i),y(i))";
            else
                eqn = eqn+" + "+num_str+"f(t(i"+(2-i)+"),y(i"+(2-i)+"))";
            end
        else
            num_str = extractAfter(""+num(i),1);
            if strcmp(num_str,"1"), num_str = ""; end
            if i == 1
                eqn = eqn+" -"+num_str+"f(t(i+1),y(i+1))";
            elseif i == 2
                eqn = eqn+" - "+num_str+"f(t(i),y(i))";
            else
                eqn = eqn+" - "+num_str+"f(t(i"+(2-i)+"),y(i"+(2-i)+"))";
            end
        end
    end
    eqn = eqn+")";
    
    % prints if requested
    if strcmp(print,'print')
        if n == 1
            fprintf("1st-order Adams-Moulton corrector:\n");
        elseif n == 2
            fprintf("2nd-order Adams-Moulton corrector:\n");
        elseif n == 3
            fprintf("3rd-order Adams-Moulton corrector:\n");
        else
            fprintf("%dth-order Adams-Moulton corrector:\n",n);
        end
        fprintf(eqn+"\n\n\n");
    end
    
end