%% |odefun_mat2vec|
% Transforms a matrix-valued ODE into a vector-valued ODE.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   dydt = odefun_mat2vec(F,t,y)
%   dydt = odefun_mat2vec(F,t,y,p)
%% Description
% |dydt = odefun_mat2vec(F,t,y)| transforms the matrix-valued ODE 
% dM/dt = F(t,M) (where F:R×R(p×p)->R(p×p) into the vector-valued ODE
% dy/dt = f(t,y) (where f:R×R(p^2)->R(p^2)). It is assumed that M is a
% square matrix.
%%
% |dydt = odefun_mat2vec(F,t,y,p)| transforms the matrix-valued ODE 
% dM/dt = F(t,M) (where F:R×R(p×q)->R(p×q) into the vector-valued ODE
% dy/dt = f(t,y) (where f:R×Rpq->Rpq). |p| specifies the number of rows of
% |M|.
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
%           <td rowspan="4" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><TT>F</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{F}(t,\mathbf{M})" title="" /></td>
%           <td>multivariate, vector-valued function (<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{F}:\mathbb{R}\times\mathbb{R}^{p\times q}\rightarrow\mathbb{R}^{p\times q}" title="" />) defining the ordinary differential equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\frac{d\mathbf{M}}{dt}=\mathbf{f}(t,\mathbf{M})" title="" />
%               <BR> - inputs to <TT>F</TT> are the current time (<TT>t</TT>, 1×1 double) and the current state matrix (<TT>M</TT>, p×q double)
%               <BR> - output of <TT>F</TT> is the state matrix derivative (<TT>dMdt</TT>, p×q double) at the current time/state</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>t</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t" title="" /></td>
%           <td>current time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>y</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{y}" title="" /></td>
%           <td>state vector at current time</td>
%           <td style="text-align:center">pq×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>p</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;p" title="" /></td>
%           <td>(OPTIONAL) number of rows of state matrix</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="1" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>dydt</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\frac{d\mathbf{y}}{dt}" title="" /></td>
%           <td>state vector derivative</td>
%           <td style="text-align:center">pq×1<BR>double</td>
%       </tr>
%   </table>
% </html>
%% See also
% <odeIC_mat2vec_doc.html |odeIC_mat2vec|> | 
% <odesol_vec2mat_doc.html |odesol_vec2mat|>