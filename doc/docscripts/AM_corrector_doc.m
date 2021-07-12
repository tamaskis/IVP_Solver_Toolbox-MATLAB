%% |AM_predictor|
% nth-order Adams-Moulton corrector.
% 
% <ODE_Solver_Toolbox_Contents.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   eqn = AM_corrector(n)
%   AM_corrector(n,'print')
%% Input/Output Parameters
% <html>
%   <table border=1>
%       <tr>
%           <td></td>
%           <td style="text-align:center"><b>Symbol</b></td>
%           <td style="text-align:center"><b>Variable</b></td>
%           <td style="text-align:center"><b>Description</b></td>
%           <td style="text-align:center"><b>Format</b></td>
%       </tr>
%       <tr>
%           <td rowspan="2" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;n" title="n" /></td>
%           <td><TT>n</TT></td>
%           <td>order of Adams-Moulton corrector</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center">-</td>
%           <td><TT>print</TT></td>
%           <td>(OPTIONAL) specify as <TT>'print'</TT> if you want to print the coefficients to the command window</td>
%           <td style="text-align:center">char</td>
%       </tr>
%       <tr>
%           <td rowspan="1" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center">-</td>
%           <td><TT>eqn</TT></td>
%           <td>nth-order Adams-Moulton corrector</td>
%           <td style="text-align:center">string</td>
%       </tr>
%   </table>
% </html>
%% Example #1
% _Return the string storing the 3rd-order Adams-Moulton corrector._
AB3_eqn = AM_corrector(3);
%% Example #2
% _Print the 3rd-order Adams-Moulton corrector._
AM_corrector(3,'print');
%% See also
% <AB_coefficients_doc.html |AB_coefficients|> | 
% <AM_coefficients_doc.html |AM_coefficients|> | 
% <AB_predictor_doc.html |AB_predictor|> | 
% <ABM_equations_doc.html |ABM_equations|>