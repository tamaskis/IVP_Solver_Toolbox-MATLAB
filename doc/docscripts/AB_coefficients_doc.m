%% |AB_coefficients|
% Coefficients for the nth-order Adams-Bashforth predictor.
% 
% <ODE_Solver_Toolbox_Contents.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   [num,gcd] = AB_coefficients(n)
%   AB_coefficients(n,'print')
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
%           <td rowspan="2" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center">-</td>
%           <td><TT>num</TT></td>
%           <td>vector storing the numerators of <img src="https://latex.codecogs.com/svg.latex?\mathbf{b}=(b_{1},...,b_{n})^{T}" title="\mathbf{b}=(b_{1},...,b_{n})^{T}" />, where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{b}" title="\mathbf{b}" /> stores the coefficients of the nth-order Adams-Bashforth predictor</td>
%           <td style="text-align:center">n×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center">-</td>
%           <td><TT>gcd</TT></td>
%           <td>greatest common denimonator of <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{b}" title="\mathbf{b}" /></td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%   </table>
% </html>
%
% *NOTE:* This function assumes the Adams-Bashforth predictor is written as
% 
% $$\mathbf{y}_{i+1}=\mathbf{y}_{i}+h\sum_{j=1}^{n}b_{j}\mathbf{f}(t_{i-j+1},\mathbf{y}_{i-j+1})$$
%% Example #1
% _Return the numerators and denominators of the coefficients of the 3rd-order 
% Adams-Bashforth predictor._
%
% |AB_coefficients| returns two quantities; |num| stores the numerators of
% the coefficients, while |gcd| stores the greatest common divisor of the
% coefficients. The numerators in |num| are already scaled to be used with
% |gcd|.
[num,gcd] = AB_coefficients(3)
%% Example #2
% _Print the coefficients of the 3rd-order Adams-Bashforth predictor._
AB_coefficients(3,'print');
%% See also
% <AM_coefficients_doc.html |AM_coefficients|> | 
% <AB_predictor_doc.html |AB_predictor|> | 
% <AM_corrector_doc.html |AM_corrector|> | 
% <ABM_equations_doc.html |ABM_equations|>