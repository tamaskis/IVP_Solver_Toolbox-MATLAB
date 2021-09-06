%==========================================================================
%
% tableau2eqns(T)  Determines the propagation equations for an explicit 
% Runge-Kutta method from its Butcher tableau.
%
%   tableau2eqns(T)
%   tableau2eqns(T,name)
%   tableau2eqns(__,'decimal')
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
%   T       - ((s+1)×(s+1) double) Butcher tableau
%   name    - (OPTIONAL) (string) name of the method
%   type    - (OPTIONAL) (char) print coefficients as 'decimal' or 
%             'fraction' (defaults to 'fraction')
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
function tableau2eqns(T,name,type)
    
    % prints name of method
    if nargin >= 2
        fprintf(name+"\n");
        for i = 1:strlength(name)
            fprintf("-");
        end
        fprintf("\n");
    end
    
    % defaults "type" to 'fraction'
    if (nargin < 3) || isempty(type)
        type = 'fraction';
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
    
    % k terms with coefficients as fractions
    if strcmpi(type,'fraction')
        for i = 1:s

            % left side of f(__,__)
            [nl,dl] = rat(c(i));
            if nl == 0            
                fprintf("k%d = f(t(i), y(i)",i);
            elseif (nl == 1) && (dl == 1)
                fprintf("k%d = f(t(i) + h, y(i)",i);
            elseif (nl ~= 1) && (dl == 1)
                fprintf("k%d = f(t(i) + %dh, y(i)",i,nl);
            elseif (nl == 1) && (dl ~= 1)
                fprintf("k%d = f(t(i) + h/%d, y(i)",i,dl);
            elseif (nl ~= 1) && (dl ~= 1)
                fprintf("k%d = f(t(i) + %dh/%d, y(i)",i,nl,dl);
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
                    fprintf("h*k%d",j);
                elseif (nr ~= 1) && (dr == 1) && (nr ~= 0)
                    fprintf("%dh*k%d",nr,j);
                elseif (nr == 1) && (dr ~= 1)
                    fprintf("h*k%d/%d",j,dr);
                elseif (nr ~= 1) && (dr ~= 1)
                    fprintf("%dh*k%d/%d",nr,j,dr);
                end
            end
            fprintf(")\n");
        end
    
    % k terms with coefficients as decimals
    else
        for i = 1:s

            % left side of f(__,__)
            if c(i) == 0            
                fprintf("k%d = f(t(i), y(i)",i);
            elseif (c(i) == 1)
                fprintf("k%d = f(t(i) + h, y(i)",i);
            else
                fprintf("k%d = f(t(i) + %.8fh, y(i)",i,c(i));
            end

            % right side of f(__,__)
            for j = 1:s
                if A(i,j) > 0
                    fprintf(" + ");
                elseif A(i,j) < 0
                    fprintf(" - ");
                    A(i,j) = abs(A(i,j));
                end
                if A(i,j) == 1
                    fprintf("h*k%d",j);
                elseif A(i,j) ~= 0
                    fprintf("%.8fh*k%d",A(i,j),j);
                end
            end
            fprintf(")\n");
        end
        
    end
    
    % converts elements of b so they all have same denominator
    [num,gcd] = same_denominator(b);
    
    % propagation equation using fractions
    if strcmpi(type,'fraction')
        if gcd == 1
            fprintf("y(i+1) = y(i) + h(");
        else
            fprintf("y(i+1) = y(i) + (h/%d)(",gcd);
        end
        dont_print_sign = false;
        for i = 1:length(b)
            n = num(i);
            if i == 1
                if n == 1
                    fprintf("k%d",i);
                elseif n ~=0
                    fprintf("%dk%d",n,i);
                else
                    dont_print_sign = true;
                end
            else
                if n == 1 
                    if dont_print_sign
                        fprintf("k%d",i);
                    else
                        fprintf(" + k%d",i);
                    end
                elseif (n > 0)
                    if dont_print_sign
                        fprintf("%dk%d",n,i);
                    else
                        fprintf(" + %dk%d",n,i);
                    end
                elseif (n < 0)
                    if dont_print_sign
                        fprintf("%dk%d",n,i);
                    else
                        fprintf(" - %dk%d",abs(n),i);
                    end
                end
            end
            if i == length(b)
                fprintf(")\n\n\n");
            end
        end
        
    % propagation equation using decimals
    else
        fprintf("y(i+1) = y(i) + h(");
        dont_print_sign = false;
        for i = 1:length(b)
            if i == 1
                if b(i) == 1
                    fprintf("k%d",i);
                elseif b(i) ~=0
                    fprintf("%.8fk%d",b(i),i);
                else
                    dont_print_sign = true;
                end
            else
                if b(i) == 1 
                    if dont_print_sign
                        fprintf("k%d",i);
                    else
                        fprintf(" + k%d",i);
                    end
                elseif (b(i) > 0)
                    if dont_print_sign
                        fprintf("%.8fk%d",b(i),i);
                    else
                        fprintf(" + %.8fk%d",b(i),i);
                    end
                elseif (b(i) < 0)
                    if dont_print_sign
                        fprintf("%.8fk%d",b(i),i);
                    else
                        fprintf(" - %.8fk%d",abs(b(i)),i);
                    end
                end
            end
            if i == length(b)
                fprintf(")\n\n\n");
            end
        end
    end
    
end