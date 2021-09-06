%==========================================================================
%
% ABM_equations  Determines and prints the nth-order
% Adams-Bashforth-Moulton predictor/corrector equations.
%
%   ABM_equations(n)
%
% See also AB_coefficients, AM_coefficients, AB_predictor, AM_corrector.
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
%   n   - (1×1 double) order of Adams-Bashforth-Moulton method
%
%==========================================================================
function ABM_equations(n)
    
    % nth-order Adams-Bashforth predictor
    eqn_AB = AB_predictor(n);
    eqn_AB = strrep(eqn_AB,"y(i+1)","yp(i+1)");
    
    % nth-order Adams-Moulton corrector
    eqn_AM = AM_corrector(n);
    eqn_AM = strrep(eqn_AM,"f(i+1)","f(t(i+1),yp(i+1))");
    
    % prints nth-order Adams-Bashforth-Moulton method
    if n == 1
        fprintf("1st-order Adams-Bashforth-Moulton method:\n");
    elseif n == 2
        fprintf("2nd-order Adams-Bashforth-Moulton method:\n");
    elseif n == 3
        fprintf("3rd-order Adams-Bashforth-Moulton method:\n");
    else
        fprintf("%dth-order Adams-Bashforth-Moulton method:\n",n);
    end
    fprintf(eqn_AB+"\n");
    fprintf(eqn_AM+"\n\n\n");
    
end