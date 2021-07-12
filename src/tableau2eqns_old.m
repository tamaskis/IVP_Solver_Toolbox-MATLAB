%==========================================================================
%
% tableau2eqns(T)  Determines the propagation equations for an explicit 
% Runge-Kutta method from its Butcher tableau.
%
%   tableau2eqns(T)
%   tableau2eqns(T,name)
%
% Copyright © 2021 Tamas Kis
% Contact: tamas.a.kis@outlook.com
% Last Update: 2021-07-11
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
%   T       - ((s+1)×(s+1)) Butcher tableau
%   name    - (OPTIONAL) (string) name of the method
%
% -----
% NOTE:
% -----
%   --> s = number of stages of the explicit Runge-Kutta method
%   --> Butcher tableau:
%
%                c1 | a11  a12  ...  a1s
%                c2 | a21  a22  ...  a2s
%           T =   : |  :    :         :
%                c3 | as1  as2  ...  ass
%                -----------------------
%                   | b1   b2   ...  bs
%
%==========================================================================
function tableau2eqns(T,name)
    
    % prints name of method
    if nargin == 2
        fprintf(name+"\n");
        for i = 1:strlength(name)
            fprintf("-");
        end
        fprintf("\n");
    end
    
    % number of rows (m) and columns (n)
    m = size(T,1);
    n = size(T,2);
    
    % extracts A, c, and b
    A = T(1:(m-1),2:n);
    c = T(1:(m-1),1);
    b = T(m,2:n);
    
    % determines s
    s = length(c);
    
    % k terms
    for i = 1:s
        
        % left side of f(__,__)
        [nl,dl] = rat(c(i));
        if nl == 0            
            fprintf("k%.0f = f(ti, y(i)",i);
        elseif (nl == 1) && (dl == 1)
            fprintf("k%.0f = f(ti + h, y(i)",i);
        elseif (nl ~= 1) && (dl == 1)
            fprintf("k%.0f = f(ti + %.0fh, y(i)",i,nl);
        elseif (nl == 1) && (dl ~= 1)
            fprintf("k%.0f = f(ti + h/%.0f, y(i)",i,dl);
        elseif (nl ~= 1) && (dl ~= 1)
            fprintf("k%.0f = f(ti + %.0fh/%.0f, y(i)",i,nl,dl);
        end
        
        % right side of f(__,__)
        for j = 1:s
            [nr,dr] = rat(A(i,j));
            if nr > 0
                fprintf(" + ");
            elseif nr < 0
                fprintf(" - ");
            end
            nr = abs(nr);
            if (nr == 1) && (dr == 1)
                fprintf("h*k%.0f",j);
            elseif (nr ~= 1) && (dr == 1) && (nr ~= 0)
                fprintf("%.0fh*k%.0f",nr,j);
            elseif (nr == 1) && (dr ~= 1)
                fprintf("h*k%.0f/%.0f",j,dr);
            elseif (nr ~= 1) && (dr ~= 1)
                fprintf("%.0fh*k%.0f/%.0f",nr,j,dr);
            end
        end
        fprintf(")\n");
        
    end
    
    % converts elements of b so they all have same denominator
    [num,gcd] = same_denominator(b);
    
    % propagation equation
    if gcd == 1
        fprintf("y(i+1) = y(i) + h(");
    else
        fprintf("y(i+1) = y(i) + (h/%.0f)(",gcd);
    end
    dont_print_sign = false;
    for i = 1:length(b)
        n = num(i);
        if i == 1
            if n == 1
                fprintf("k%.0f",i);
            elseif n ~=0
                fprintf("%.0fk%.0f",n,i);
            else
                dont_print_sign = true;
            end
        else
            if n == 1 
                if dont_print_sign
                    fprintf("k%.0f",i);
                else
                    fprintf(" + k%.0f",i);
                end
            elseif (n > 0)
                if dont_print_sign
                    fprintf("%.0fk%.0f",n,i);
                else
                    fprintf(" + %.0fk%.0f",n,i);
                end
            elseif (n < 0)
                if dont_print_sign
                    fprintf("%.0fk%.0f",n,i);
                else
                    fprintf(" - %.0fk%.0f",abs(n),i);
                end
            end
        end
        if i == length(b)
            fprintf(")\n\n\n");
        end
    end
    
end