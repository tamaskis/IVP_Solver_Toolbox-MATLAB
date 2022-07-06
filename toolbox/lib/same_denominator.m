%==========================================================================
%
% same_denominator  Scales a set of fractions so they each have the same
% denominator.
%
%   [num,gcd] = same_denominator(X)
%   [num,gcd] = same_denominator(X,'print')
%   same_denominator(X,'print')
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-07-05
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   X       - (double array of any size) array of fractions
%   print   - (OPTIONAL) (char) prints new fractions if input as 'print'
%
% -------
% OUTPUT:
% -------
%   num     - (double array of same size as X) scaled numerators of X
%   gcd     - (1×1 double) greatest common divisor of the elements of X
%
%==========================================================================
function [num,gcd] = same_denominator(X,print)

    % defaults "print" to 'do not print'
    if (nargin < 2) || isempty(print)
        print = 'do not print';
    end
    
    % stores original size of array for later
    s = size(X);
    
    % reshapes array
    N = numel(X);
    Y = reshape(X,[N,1]);
    
    % preallocates arrays
    num = zeros(size(Y));
    den = zeros(size(Y));
    
    % finds the numerator and denominator of each element of X
    for i = 1:length(Y)
        [num(i),den(i)] = rat(Y(i),1e-10);
        format long;
    end
    
    % saves original numerators and denominators for printing later
    num0 = num;
    den0 = den;
    
    % new array of denominators sorted from smallest to greatest
    den_increasing = sort(den);
    
    % initialize the greatest common divisor (i.e. least common multiple of
    % the denominator) to be the smallest denominator
    gcd = den_increasing(1);
    
    % finds the least common multiple of the set of denominators
    for i = 2:length(Y)
        gcd = lcm(gcd,den_increasing(i));
    end
    
    % scale the numerators so the denominator can be replaced with gcd
    for i = 1:length(Y)
        num(i) = num(i)*gcd/den(i);
    end
    
    % prints if requested
    if strcmpi(print,'print')
        fprintf('Updated fractions:');
        for i = 1:length(num)
            fprintf('\n%d/%d   -->   %d/%d',num0(i),den0(i),num(i),gcd);
        end
        fprintf('\n\n\n');
    end
    
    % reshapes result to original size
    num = reshape(num,s);
    
end