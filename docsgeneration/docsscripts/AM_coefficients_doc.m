%% |AM_coefficients|
% Coefficients for the mth-order Adams-Moulton corrector.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   [num,gcd] = AM_coefficients(m)
%   AM_coefficients(m,'print')
%% Description
% |[num,gcd] = AM_coefficients(m)| returns a vector of numerators, |num|,
% and the greatest common denominator, |gcd|, of the vector
% $\mathbf{b}=(b_{1},...,b_{m})^{T}$ storing the coefficients of the
% mth-order Adams-Moulton corrector.
%%
% |AM_coefficients(m,'print')| prints the coefficients $b_{1},...,b_{m}$ of 
% the mth-order Adams-Moulton corrector to the Command Window.
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
%           <td rowspan="2" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><TT>m</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;m" title="m" /></td>
%           <td>order of Adams-Moulton corrector</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>print</TT></td>
%           <td style="text-align:center">-</td>
%           <td>(OPTIONAL) specify as <TT>'print'</TT> if you want to print the coefficients to the command window</td>
%           <td style="text-align:center">char</td>
%       </tr>
%       <tr>
%           <td rowspan="2" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>num</TT></td>
%           <td style="text-align:center">-</td>
%           <td>vector storing the numerators of <img src="https://latex.codecogs.com/svg.latex?\mathbf{b}=(b_{1},...,b_{m})^{T}" title="\mathbf{b}=(b_{1},...,b_{m})^{T}" />, where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{b}" title="\mathbf{b}" /> stores the coefficients of the mth-order Adams-Moulton corrector</td>
%           <td style="text-align:center">n×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>gcd</TT></td>
%           <td style="text-align:center">-</td>
%           <td>greatest common denimonator of <img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{b}" title="\mathbf{b}" /></td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%   </table>
% </html>
%
% *NOTE:* This function assumes the Adams-Moulton corrector is written as
% 
% $$\mathbf{y}_{n+1}=\mathbf{y}_{n}+h\sum_{i=1}^{m}b_{i}\mathbf{f}(t_{n-i+2},\mathbf{y}_{n-i+2})$$
%% Example #1
% _Return the numerators and denominators of the coefficients of the 3rd-order 
% Adams-Moulton corrector._
%
% |AM_coefficients| returns two quantities; |num| stores the numerators of
% the coefficients, while |gcd| stores the greatest common divisor of the
% coefficients. The numerators in |num| are already scaled to be used with
% |gcd|.
[num,gcd] = AM_coefficients(3)
%% Example #2
% _Print the coefficients of the 3rd-order Adams-Moulton corrector._
AM_coefficients(3,'print');
%% See also
% <AB_coefficients_doc.html |AB_coefficients|> | 
% <AB_predictor_doc.html |AB_predictor|> | 
% <AM_corrector_doc.html |AM_corrector|> | 
% <ABM_equations_doc.html |ABM_equations|>