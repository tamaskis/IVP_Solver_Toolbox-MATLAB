%==========================================================================
%
% ABM_equations  Determines and prints the mth-order
% Adams-Bashforth-Moulton predictor-corrector equations.
%
%   ABM_equations(m)
%
% See also AB_coefficients, AM_coefficients, AB_predictor, AM_corrector.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   m   - (1×1 double) order of Adams-Bashforth-Moulton method
%
%==========================================================================
function ABM_equations(m)
    
    % mth-order Adams-Bashforth predictor
    eqn_AB = AB_predictor(m);
    eqn_AB = strrep(eqn_AB,"y(n+1)","yp(n+1)");
    
    % mth-order Adams-Moulton corrector
    eqn_AM = AM_corrector(m);
    eqn_AM = strrep(eqn_AM,"f(n+1)","f(t(n+1),yp(n+1))");
    
    % prints mth-order Adams-Bashforth-Moulton method
    if m == 1
        fprintf("1st-order Adams-Bashforth-Moulton method:\n");
    elseif m == 2
        fprintf("2nd-order Adams-Bashforth-Moulton method:\n");
    elseif m == 3
        fprintf("3rd-order Adams-Bashforth-Moulton method:\n");
    else
        fprintf("%dth-order Adams-Bashforth-Moulton method:\n",m);
    end
    fprintf(eqn_AB+"\n");
    fprintf(eqn_AM+"\n\n\n");
    
end