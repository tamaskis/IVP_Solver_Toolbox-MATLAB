%==========================================================================
%
% ABM_equations  Determines and prints the nth-order
% Adams-Bashforth-Moulton predictor/corrector equations.
%
%   ABM_equations(n)
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
%   n   - (1×1) order of Adams-Bashforth-Moulton method
%
%==========================================================================
function ABM_equations(n)
    
    % nth-order Adams-Bashforth predictor
    eqn_AB = AB_predictor(n);
    eqn_AB = strrep(eqn_AB,"y(i+1)","yp(i+1)");
    
    % nth-order Adams-Moulton corrector
    eqn_AM = AM_corrector(n);
    eqn_AM = extractBetween(eqn_AM,7,strlength(eqn_AM));
    eqn_AM = strrep(eqn_AM,"y(i+1)","yp(i+1)");
    eqn_AM = " y(i+1)"+eqn_AM;
    
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