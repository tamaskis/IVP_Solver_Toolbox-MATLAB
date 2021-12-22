%% |odefun_mat2vec|
% Transforms a matrix-valued ODE into a vector-valued ODE.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   f = odefun_mat2vec(F)
%   f = odefun_mat2vec(F,p)
%% Description
% |f = odefun_mat2vec(F)| transforms the matrix-valued ODE
% $d\mathbf{M}/dt=\mathbf{F}(t,\mathbf{M})$ (where 
% $\mathbf{F}:\mathrm{R}\times\mathrm{R}^{p\times p}\to\mathrm{R}^{p\times p}$
% into the vector-valued ODE $d\mathbf{y}/dt=\mathbf{f}(t,\mathbf{y})$
% (where $\mathbf{f}:\mathrm{R}\times\mathrm{R}^{p^{2}}\to\mathrm{R}^{p^{2}}$).
% It is assumed that $\mathbf{M}$ is a square matrix.
%%
% |f = odefun_mat2vec(F,p)| transforms the matrix-valued ODE 
% $d\mathbf{M}/dt=\mathbf{F}(t,\mathbf{M})$ (where 
% $\mathbf{F}:\mathrm{R}\times\mathrm{R}^{p\times q}\to\mathrm{R}^{p\times q}$)
% into the vector-valued ODE $d\mathbf{y}/dt=\mathbf{f}(t,\mathbf{y})$
% (where $\mathbf{f}:\mathrm{R}\times\mathrm{R}^{pq}\to\mathrm{R}^{pq}$).
% |p| specifies the number of rows of |M|.
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
%           <td style="text-align:center"><TT>F</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{F}(t,\mathbf{M})" title="" /></td>
%           <td>multivariate, matrix-valued function (<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{F}:\mathbb{R}\times\mathbb{R}^{p\times q}\rightarrow\mathbb{R}^{p\times q}" title="" />) defining the ordinary differential equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\frac{d\mathbf{M}}{dt}=\mathbf{f}(t,\mathbf{M})" title="" />
%               <BR> - inputs to <TT>F</TT> are the current time (<TT>t</TT>, 1×1 double) and the current state matrix (<TT>M</TT>, p×q double)
%               <BR> - output of <TT>F</TT> is the state matrix derivative (<TT>dMdt</TT>, p×q double) at the current time/state</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>p</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;p" title="" /></td>
%           <td>(OPTIONAL) number of rows of state matrix</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="1" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>f</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}(t,\mathbf{y})" title="" /></td>
%           <td>multivariate, vector-valued function (<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}:\mathbb{R}\times\mathbb{R}^{pq}\rightarrow\mathbb{R}^{pq}" title="" />) defining the ordinary differential equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\frac{d\mathbf{y}}{dt}=\mathbf{f}(t,\mathbf{y})" title="" />
%               <BR> - inputs to <TT>f</TT> are the current time (<TT>t</TT>, 1×1 double) and the current state vector (<TT>y</TT>, pq×1 double)
%               <BR> - output of <TT>f</TT> is the state vector derivative (<TT>dydt</TT>, pq×1 double) at the current time/state</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%   </table>
% </html>
%% Example
% Click <Matrix_ODE_Example_doc.html here> for an example.
%% See also
% <odeIC_mat2vec_doc.html |odeIC_mat2vec|> | 
% <odesol_vec2mat_doc.html |odesol_vec2mat|>