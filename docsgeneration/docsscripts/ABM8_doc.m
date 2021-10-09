%% |ABM8|
% Adams-Bashforth-Moulton 8th-order method.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   [t,y] = ABM8(f,[t0,tf],y0,h)
%   [t,y] = ABM8(f,{t0,C},y0,h)
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
%           <td rowspan="5" style="text-align:center"><b>Input</b></td>
%           <td><TT>f</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}(t,\mathbf{y})" title="\mathbf{f}(t,\mathbf{y})" /></td>
%           <td>multivariate, vector-valued function (<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}:\mathbb{R}^{n+1}\rightarrow\mathbb{R}^{n}" title="" />) defining the ordinary differential equation <BR>  - inputs to <TT>f</TT> are the current time (1×1 double) and the current state vector (n×1 double)<BR>  - output of <TT>f</TT> is the state vector derivative at the current time/state (n×1 double)</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td><TT>t0</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{0}" title="t_{0}" /></td>
%           <td>initial time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td><TT>tf</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{f}" title="t_{f}" /></td>
%           <td>final time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td><TT>C</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;C(t,\mathbf{y})" title="C(t,\mathbf{y})" /></td>
%           <td>condition function<BR>  - inputs are the current time (1×1 double) and the current state vector (n×1 double)<BR>  - output is a 1×1 logical (<TT>true</TT> if solver should continue running, <TT>false</TT> if solver should terminate)</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td><TT>y0</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{y}_{0}" title="\mathbf{y}_{0}" /></td>
%           <td>initial condition</td>
%           <td style="text-align:center">n×1<BR>double</td>
%       </tr>
%       <tr>
%           <td rowspan="2" style="text-align:center"><b>Output</b></td>
%           <td><TT>t</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></td>
%           <td>time vector</td>
%           <td style="text-align:center">(N+1)×1<BR>double</td>
%       </tr>
%       <tr>
%           <td><TT>y</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Y}" title="\mathbf{Y}" /></td>
%           <td>solution matrix<BR>- the ith row of <TT>y</TT> stores the transpose of the solution corresponding the ith time in <TT>t</TT></td>
%           <td style="text-align:center">(N+1)×n<BR>double</td>
%       </tr>
%   </table>
% </html>
%% Example #1: Time detection.
% _Consider the Lorenz system, with $\rho=28$, $\sigma=10$, and
% $\beta=8/3$:_
%
% $$\dot{x}=\sigma(y-x)$$
%
% $$\dot{y}=x(\rho-z)-y$$
%
% $$\dot{z}=xy-\beta z$$
%
% $$x(0)=0,\quad y(0)=1,\quad z(0)=1.05$$
%
% _Plot the solution of this system for $t$ in the interval $[0,100]$._
%
% First, we can rewrite this system in vector form by letting $x_{1}=x$,
% $x_{2}=y$, $x_{3}=z$, and $\mathbf{x}=(x_{1},x_{2},x_{3})^{T}$.
%
% $$\dot{\mathbf{x}}=\mathbf{f}(t,\mathbf{x})=\pmatrix{\sigma(x_{2}-x_{1}) 
% \cr x_{1}(\rho-x_{3})-x_{2} \cr x_{1}x_{2}-\beta x_{3}}$$
%
% $$\mathbf{x}(0)=\pmatrix{0 \cr 1 \cr 1.05}$$
%
% Defining the system and its initial condition in MATLAB,
%%

% Lorenz parameters
rho = 28;
sigma = 10;
beta = 8/3;

% Lorenz equations in vector form
f = @(t,x) [sigma*(x(2)-x(1));
            x(1)*(rho-x(3))-x(2);
            x(1)*x(2)-beta*x(3)];
        
% initial condition
x0 = [0;
      1;
      1.05];
%%
% Solving the system for $t$ in the interval $[0,100]$ using a time step 
% (i.e. step size) of $\Delta t=0.001$,
[t,x] = ABM8(f,[0,100],x0,0.001);
%%
% Plotting the solution,
figure;
plot3(x(:,1),x(:,2),x(:,3));
view(45,20);
grid on;
xlabel('$x$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
zlabel('$z$','interpreter','latex','fontsize',18);
%% Example #2: Event detection.
% _Consider the initial value problem_
%
% $$\frac{dy}{dx}=y,\quad y(2)=3$$
%
% _Find the solution $y(t)$ until $y=10$._
%
% First, let's define our ODE ($\frac{dy}{dt}=f(t,y)$) and initial
% condition in MATLAB.
f = @(t,y) y;
y2 = 3;
%%
% Next, let's define the condition function, $C(t,y)$. Since we want to
% continue solving until $y=10$, we want the solver to run while $y\leq10$;
% this forms our condition. Therefore,
C = @(t,y) y <= 10;
%%
% Solving for $y$ using a step size of $h=0.01$,
[t,y] = ABM8(f,{2,C},y2,0.01);
%%
% Plotting the solution,
figure;
plot(t,y);
grid on;
xlabel('$t$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
%% Example #3:  Backward integration.
% _Consider the same ODE as in Example #2, but now we know its value at
% $t=20$._
%
% $$\frac{dy}{dx}=y,\quad y(20)=50$$
%
% _Find $y(10)$. Then, confirm your result by solving the same ODE from
% $t=10$ to $t=20$, using $y(10)$ as the initial condition._
%
% In this case, we know $y$ at $t=20$, and want to find $y$ at a _previous_
% time $t=10$. To do this, we need to solve the ODE backwards, _from_
% $t=20$ _to_ $t=10$. First, let's define our ODE ($\frac{dy}{dt}=f(t,y)$) 
% and initial condition in MATLAB.
f = @(t,y) y;
y20 = 50;
%%
% Let's use a step size of magnitude $0.001$ this time. Since we are
% integrating _backwards_, we actually need to use a _negative_ step size
% of $h=-0.001$.
h = -0.001;
%%
% Solving for $y(t)$ from $t=20$ to $t=10$,
[t,y] = ABM8(f,[20,10],y20,h);
%%
% The solution for $y$ corresponding to $t=10$ will be located at the last
% element of the solution matrix, since $t=10$ is stored in the last
% element of the time vector. Therefore, $y(10)$ is
y10 = y(end)
%%
% Confirming our result by solving the same ODE but from $t=10$ to $t=20$
% and using our result for $y(10)$ as the initial condition,
[t,y] = ABM8(f,[10,20],y10,0.001);
y20 = y(end)
%% See also
% <ABM2_doc.html |ABM2|> | 
% <ABM3_doc.html |ABM3|> | 
% <ABM4_doc.html |ABM4|> | 
% <ABM5_doc.html |ABM5|> | 
% <ABM6_doc.html |ABM6|> | 
% <ABM7_doc.html |ABM7|> | 
% <ABM8_doc.html |ABM8|>