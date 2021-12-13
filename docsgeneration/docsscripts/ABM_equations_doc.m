%% |ABM_equations|
% mth-order Adams-Bashforth-Moulton equations.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   ABM_equations(m)
%% Description
% |ABM_equations(m)| prints the mth-order Adams-Bashforth-Moulton
% predictor-corrector equations to the Command Window.
%% Input/Output Parameters
% <html>
%   <table border=1>
%       <tr>
%           <td></td>
%           <td style="text-align:center"><b>Variable</b></td>
%           <td style="text-align:center"><b>Symbol</b></td>
%           <td style="text-align:center"><b>Description</b></td>
%           <td style="text-align:center"><b>Format</b></td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><TT>m</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;m" title="m" /></td>
%           <td>order of Adams-Bashforth-Moulton method</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%   </table>
% </html>
%% Example #1: 1st-order Adams-Bashforth-Moulton equations.
ABM_equations(1);
%% Example #2: 2nd-order Adams-Bashforth-Moulton equations.
ABM_equations(2);
%% Example #3: 3rd-order Adams-Bashforth-Moulton equations.
ABM_equations(3);
%% Example #4: 4th-order Adams-Bashforth-Moulton equations.
ABM_equations(4);
%% Example #5: 5th-order Adams-Bashforth-Moulton equations.
ABM_equations(5);
%% Example #6: 6th-order Adams-Bashforth-Moulton equations.
ABM_equations(6);
%% Example #7: 7th-order Adams-Bashforth-Moulton equations.
ABM_equations(7);
%% Example #8: 8th-order Adams-Bashforth-Moulton equations.
ABM_equations(8);
%% See also
% <AB_coefficients_doc.html |AB_coefficients|> | 
% <AM_coefficients_doc.html |AM_coefficients|> | 
% <ABM_equations_doc.html |AB_predictor|> | 
% <AM_corrector_doc.html |AM_corrector|>