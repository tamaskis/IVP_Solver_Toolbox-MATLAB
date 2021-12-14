%% |RK2_ralston_step|
% Propagates the state vector forward one time step using Ralston's 
% second-order method (Runge-Kutta second-order method).
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   y_next = RK2_ralston_step(f,t,y,h)
%% Description
% |y_next = RK2_ralston_step(f,t,y,h)| returns the state vector at the next 
% sample time, |y_next|, given the current state vector |y| at time |t|,
% the function |f(t,y)| defining the ODE 
% $\dot{\mathbf{y}}=\mathbf{f}(t,\mathbf{y})$, and the step size |h|.
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
%           <td style="text-align:center"><TT>f</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}(t,\mathbf{y})" title="" /></td>
%           <td>multivariate, vector-valued function (<img
%           src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}:\mathbb{R}\times\mathbb{R}^{p}\rightarrow\mathbb{R}^{p}"
%           title="" />) defining the ordinary differential equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\frac{d\mathbf{y}}{dt}=\mathbf{f}(t,\mathbf{y})" title="" />
%               <BR> - inputs to <TT>f</TT> are the current time (<TT>t</TT>, 1×1 double) and the current state vector (<TT>y</TT>, p×1 double)
%               <BR> - output of <TT>f</TT> is the state vector derivative (<TT>ydot</TT>, p×1 double) at the current time/state</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>t</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{n}" title="" /></td>
%           <td>current sample time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>y</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{y}_{n}=\mathbf{y}(t_{n})" title="" /></td>
%           <td>state vector (i.e. solution) at the current sample time</td>
%           <td style="text-align:center">p×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>h</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;h" title="h" /></td>
%           <td>step size</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="1" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>y_next</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{y}_{n+1}=\mathbf{y}(t_{n+1})" title="" /></td>
%           <td>state vector (i.e. solution) at the next sample time, <img src="https://latex.codecogs.com/svg.latex?t_{n+1}=t_{n}+h" title="" /></td>
%           <td style="text-align:center">p×1<BR>double</td>
%       </tr>
%   </table>
% </html>
%% Example
% _Consider the initial value problem_
%
% $$\frac{dy}{dt}=y,\quad y(2)=3$$
%
% _Find the solution $y(t)$ until $t=10$ using |RK2_ralston_step|. Then, 
% compare your result to the solution found by |RK2_ralston|._
%
% First, let's define our ODE ($\frac{dy}{dt}=f(t,y)$) and initial
% condition in MATLAB.
f = @(t,y) y;
y2 = 3;
%%
% Let's define a time vector between $t=2$ and $t=10$ with a spacing of
% $h=0.01$.
h = 0.01;
t = (2:h:10)';
%%
% Solving for $y(t)$ using |RK2_ralston_step| and comparing the result to
% the result obtained using |RK2_ralston|,

% preallocate vector to store solution
y = zeros(size(t));

% store initial condition
y(1) = y2;

% solving using "RK2_ralston_step"
for i = 1:(length(t)-1)
    y(i+1) = RK2_ralston_step(f,t(i),y(i),h);
end

% solving using "RK2_ralston"
[t_RK2_ralston,y_RK2_ralston] = RK2_ralston(f,[2,10],y2,h);

% maximum absolute error between the two results
max(abs(y_RK2_ralston-y))
%%
% As expected, the two methods obtain identical results.
%% See also
% <RK1_euler_step_doc.html |RK1_euler_step|> | 
% <RK2_step_doc.html |RK2_step|> | 
% <RK2_heun_step_doc.html |RK2_heun_step|> | 
% <RK3_step_doc.html |RK3_step|> |
% <RK3_heun_step_doc.html |RK3_heun_step|> |
% <RK3_ralston_step_doc.html |RK3_ralston_step|> |
% <SSPRK3_step_doc.html |SSPRK3_step|> |
% <RK4_step_doc.html |RK4_step|> |
% <RK4_ralston_step_doc.html |RK4_ralston_step|> |
% <RK4_38_step_doc.html |RK4_38_step|>