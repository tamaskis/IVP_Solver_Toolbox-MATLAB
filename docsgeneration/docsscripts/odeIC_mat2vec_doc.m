%% |odeIC_mat2vec|
% Transforms the initial condition for a matrix-valued ODE into the initial
% condition for the corresponding vector-valued ODE.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   y0 = odeIC_mat2vec(M0)
%% Description
% |y0 = odeIC_mat2vec(M0)| transforms the initial condition $\mathbf{M}_{0}\in\mathrm{R}^{p\times q}$
% for a matrix-valued ODE into the initial condition $\mathbf{y}_{0}\in\mathrm{R}^{pq}$
% for the corresponding vector-valued ODE.
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
%           <td rowspan="1" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><TT>M0</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{M}_{0}" title="" /></td>
%           <td>initial condition for matrix-valued ODE</td>
%           <td style="text-align:center">p×q<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="1" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>y0</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{y}_{0}" title="" /></td>
%           <td>initial condition for corresponding vector-valued ODE</td>
%           <td style="text-align:center">pq×1<BR>double</td>
%       </tr>
%   </table>
% </html>
%% Example
% Click <Matrix_ODE_Example_doc.html here> for an example.
%% See also
% <odefun_mat2vec_doc.html |odefun_mat2vec|> | 
% <odesol_vec2mat_doc.html |odesol_vec2mat|>