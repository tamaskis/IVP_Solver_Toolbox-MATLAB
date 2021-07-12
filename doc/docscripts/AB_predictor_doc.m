%% |AB_predictor|
% nth-order Adams-Bashforth predictor.
% 
% <ODE_Solver_Toolbox_Contents.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   eqn = AB_predictor(n)
%   AB_predictor(n,'print')
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
%           <td>order of Adams-Bashforth predictor</td>
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
%           <td>nth-order Adams-Bashforth predictor</td>
%           <td style="text-align:center">string</td>
%       </tr>
%   </table>
% </html>
%% Example #1
% _Return the string storing the 3rd-order Adams-Bashforth predictor._
AB3_eqn = AB_predictor(3);
%% Example #2
% _Print the 3rd-order Adams-Bashforth predictor._
AB_predictor(3,'print');
%% See also
% <AB_coefficients_doc.html |AB_coefficients|> | 
% <AM_coefficients_doc.html |AM_coefficients|> | 
% <AM_corrector_doc.html |AM_corrector|> | 
% <ABM_equations_doc.html |ABM_equations|>