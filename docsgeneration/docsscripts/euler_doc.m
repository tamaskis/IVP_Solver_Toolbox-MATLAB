%% |euler|
% Euler (1st-order) method.
% 
% <ODE_Solver_Toolbox_Contents.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   [t,y] = euler(f,[t0,tf],y0,h)
%   [t,y] = euler(f,{t0,C},y0,h)
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
%           <td rowspan="5" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}(t,\mathbf{y})" title="\mathbf{f}(t,\mathbf{y})" /></td>
%           <td><TT>f</TT></td>
%           <td>ordinary differential equation<BR>  - inputs are <TT>t</TT> (1×1 double) and <TT>y</TT> (n×1 double)<BR>  - output is an n×1 double</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{0}" title="t_{0}" /></td>
%           <td><TT>t0</TT></td>
%           <td>initial time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{f}" title="t_{f}" /></td>
%           <td><TT>tf</TT></td>
%           <td>final time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;C(t,\mathbf{y})" title="C(t,\mathbf{y})" /></td>
%           <td><TT>C</TT></td>
%           <td>condition function<BR>  - inputs are <TT>t</TT> (1×1 double) and <TT>y</TT> (n×1 double)<BR>  - output is a 1×1 logical</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{y}_{0}" title="\mathbf{y}_{0}" /></td>
%           <td><TT>y0</TT></td>
%           <td>initial condition</td>
%           <td style="text-align:center">n×1<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="2" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></td>
%           <td><TT>t</TT></td>
%           <td>time vector</td>
%           <td style="text-align:center">(N+1)×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Y}" title="\mathbf{Y}" /></td>
%           <td><TT>y</TT></td>
%           <td>solution matrix<BR>- the ith row of <TT>y</TT> stores the transpose of the solution corresponding the ith time in <TT>t</TT></td>
%           <td style="text-align:center">(N+1)×n<BR>double</td>
%       </tr>
%   </table>
% </html>
%% Examples
% See the <Examples.html examples> page for a list of examples. The
% examples are performed using the <RK4_doc.html |RK4|> function, but the
% |euler| function could be used as well.
%% See also
% <RK2_doc.html |RK2|> | 
% <RK2_heun_doc.html |RK2_heun|> | 
% <RK2_ralston_doc.html |RK2_ralston|> |
% <RK3_doc.html |RK3|> |
% <RK3_heun_doc.html |RK3_heun|> |
% <RK3_ralston_doc.html |RK3_ralston|> |
% <SSPRK3_doc.html |SSPRK3|> |
% <RK4_doc.html |RK4|> |
% <RK4_ralston_doc.html |RK4_ralston|> |
% <RK4_38_doc.html |RK4_38|>