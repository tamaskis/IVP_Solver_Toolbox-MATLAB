%% |odesol_vec2mat|
% Transforms the solution matrix for a vector-valued ODE into the solution
% array for the corresponding matrix-valued ODE.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   M = odesol_vec2mat(y)
%   M = odesol_vec2mat(y,p)
%% Description
% |M = odesol_vec2mat(y)| transforms the solution matrix, $\mathbf{Y}\in\mathrm{R}^{(N+1)\times p^{2}}$
% (|y| in MATLAB), for the vector-valued ODE into the solution array, $\mathbf{M}_{\mathrm{sol}}\in\mathrm{R}^{p\times p\times (N+1)}$
% (|M| in MATLAB) for the corresponding matrix-valued ODE. It is assumed 
% that $\mathbf{M}$ is a square matrix.
%%
% |M = odesol_vec2mat(y,p)| transforms the solution matrix, $\mathbf{Y}\in\mathrm{R}^{(N+1)\times pq}$
% (|y| in MATLAB), for the vector-valued ODE into the solution array, $\mathbf{M}_{\mathrm{sol}}\in\mathrm{R}^{p\times q\times (N+1)}$
% (|M| in MATLAB) for the corresponding matrix-valued ODE. |p| specifies 
% the number of rows of |M|.
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
%           <td style="text-align:center"><TT>y</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Y}" title="\mathbf{Y}" /></td>
%           <td>solution matrix<BR>- the nth row of <TT>y</TT> stores the <i>transpose</i> of the solution corresponding to the nth time in the time vector, <img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="" /></td>
%           <td style="text-align:center">(N+1)×pq<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>p</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;p" title="" /></td>
%           <td>(OPTIONAL) number of rows of state matrix</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="1" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>M</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{M}_{\mathrm{sol}}" title="" /></td>
%           <td>solution array<BR>- the nth layer of <TT>M</TT> stores the solution corresponding to the nth time in the time vector, <img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="" /></td>
%           <td style="text-align:center">p×q×(N+1)<BR>double</td>
%       </tr>
%   </table>
% </html>
%% Example
% Click <Matrix_ODE_Example_doc.html here> for an example.
%% See also
% <odefun_mat2vec_doc.html |odefun_mat2vec|> | 
% <odeIC_mat2vec_doc.html |odeIC_mat2vec|>