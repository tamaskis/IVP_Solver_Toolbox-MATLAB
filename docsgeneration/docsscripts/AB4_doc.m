%% |AB4|
% Adams-Bashforth 4th-order method.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   [t,y] = AB4(f,[t0,tf],y0,h)
%   [t,y] = AB4(f,{t0,C},y0,h)
%   [t,y] = AB4(__,wb)
%% Description
% |[t,y] = AB4(f,[t0,tf],y0,h)| solves the ODE defined by |f(t,y)| from 
% |t0| until |tf| using the Adams-Bashforth 4th-order method with an
% initial condition |y0| and step size |h|.
%%
% |[t,y] = AB4(f,{t0,C},y0,h)| does the same as the syntax above, but
% instead of terminating at a final time |tf|, the solver terminates once
% the condition function |C(t,y)| is no longer satisfied.
%%
% |[t,y] = AB4(...,wb)| can be used with either of the syntaxes above to
% define a waitbar. If |wb| is input as |true|, then a waitbar is displayed
% with the default message 'Solving ODE...'. To specify a custom waitbar
% message, input |wb| as a char array storing the desired message.
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
%           <td rowspan="7" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><TT>f</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}(t,\mathbf{y})" title="" /></td>
%           <td>multivariate, vector-valued function (<img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{f}:\mathbb{R}\times\mathbb{R}^{p}\rightarrow\mathbb{R}^{p}" title="" />) defining the ordinary differential equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\frac{d\mathbf{y}}{dt}=\mathbf{f}(t,\mathbf{y})" title="" />
%               <BR> - inputs to <TT>f</TT> are the current time (<TT>t</TT>, 1×1 double) and the current state vector (<TT>y</TT>, p×1 double)
%               <BR> - output of <TT>f</TT> is the state vector derivative (<TT>dydt</TT>, p×1 double) at the current time/state</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>t0</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{0}" title="t_{0}" /></td>
%           <td>initial time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>tf</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;t_{f}" title="t_{f}" /></td>
%           <td>final time</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>C</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;C(t,\mathbf{y})" title="C(t,\mathbf{y})" /></td>
%           <td>condition function (<img src="https://latex.codecogs.com/svg.latex?\inline&space;C:\mathbb{R}\times\mathbb{R}^{p}\rightarrow\mathbb{B}" title="" />)
%               <BR> - inputs are the current time (<TT>t</TT>, 1×1 double) and the current state vector (<TT>y</TT>, p×1 double)
%               <BR> - output is a 1×1 logical (<TT>true</TT> if solver should continue running, <TT>false</TT> if solver should terminate)</td>
%           <td style="text-align:center">1×1<BR>function_handle</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>y0</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\mathbf{y}_{0}" title="\mathbf{y}_{0}" /></td>
%           <td>initial condition</td>
%           <td style="text-align:center">p×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>h</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\inline&space;h" title="h" /></td>
%           <td>step size</td>
%           <td style="text-align:center">1×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>wb</TT></td>
%           <td style="text-align:center">-</td>
%           <td>waitbar parameters<BR>  - input as "<TT>true</TT>" if you want waitbar with default message displayed<BR>  - input as a char array storing a message if you want a custom message displayed on the waitbar</td>
%           <td style="text-align:center">char array<BR><i>or</i><BR>1×1 logical</td>
%       </tr>
%       <tr>
%           <td rowspan="2" style="text-align:center"><b>Output</b></td>
%           <td style="text-align:center"><TT>t</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></td>
%           <td>time vector</td>
%           <td style="text-align:center">(N+1)×1<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>y</TT></td>
%           <td style="text-align:center"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Y}" title="\mathbf{Y}" /></td>
%           <td>solution matrix<BR>- the ith row of <TT>y</TT> stores the <i>transpose</i> of the solution corresponding to the ith time in <TT>t</TT> (see below)</td>
%           <td style="text-align:center">(N+1)×p<BR>double</td>
%       </tr>
%   </table>
% </html>
%
% *Time vector and solution matrix:*
% 
% $$\mathbf{t}=\pmatrix{t_{0} \cr\vdots\cr t_{N}},\quad\quad\mathbf{Y}=\pmatrix{\mathbf{y}(t_{0})^{T} \cr\vdots\cr\mathbf{y}(t_{N})^{T}}$$
%% Example #1: Time detection.
% _Consider the <https://en.wikipedia.org/wiki/Lorenz_system Lorenz system>, with $\rho=28$, $\sigma=10$, and
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
[t,x] = AB4(f,[0,100],x0,0.001);
%%
% Plotting the solution,
figure;
plot3(x(:,1),x(:,2),x(:,3));
view(45,20);
grid on;
xlabel('$x$','Interpreter','latex','FontSize',18);
ylabel('$y$','Interpreter','latex','FontSize',18);
zlabel('$z$','Interpreter','latex','FontSize',18);
%% Example #2: Event detection.
% _Consider the initial value problem_
%
% $$\frac{dy}{dt}=y,\quad y(2)=3$$
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
[t,y] = AB4(f,{2,C},y2,0.01);
%%
% Plotting the solution,
figure;
plot(t,y,'LineWidth',1.5);
grid on;
xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$y$','Interpreter','latex','FontSize',18);
%% Example #3: Backward integration (time detection case).
% _Consider the same ODE as in Example #2, but we are now given its value
% at $t=20$._
%
% $$\frac{dy}{dt}=y,\quad y(20)=50$$
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
% Solving for $y(t)$ from $t=20$ to $t=10$ using a step size of $h=0.001$,
[t,y] = AB4(f,[20,10],y20,0.001);
%%
% The solution for $y$ corresponding to $t=10$ will be located at the last
% element of the solution matrix, since $t=10$ is stored in the last
% element of the time vector. Therefore, $y(10)$ is
y10 = y(end)
%%
% Confirming our result by solving the same ODE but from $t=10$ to $t=20$
% and using our result for $y(10)$ as the initial condition,
[t,y] = AB4(f,[10,20],y10,0.001);
y20 = y(end)
%% Example #4: Backward integration (event-detection case).
% _Once again, consider_ 
%
% $$\frac{dy}{dt}=y,\quad y(20)=50$$
%
% _Find $y(10)$ using an ODE solver with a condition function (i.e. event-detection)._
%%
% First, let's define our ODE ($\frac{dy}{dt}=f(t,y)$) and initial 
% condition in MATLAB.
f = @(t,y) y;
y20 = 50;
%%
% The event that terminates the solver is when $t=10$. Therefore, we define
% the condition function as
C = @(t,y) t > 10;
%%
% In the time-detection case (Example #3) we input |[t0,tf] = [20,10]|, so 
% the solver knew to integrate backwards in time since |t0 > tf|. 
% Consequently, internally, the solver made the step size negative. 
% However, for the event-detection case, just given |t| and |C(t,y)|, the 
% solver won't know to use a negative step size to integrate backwards in
% time. Therefore, we must manually specify a negative step size. Solving
% for $y(10)$,
[t,y] = AB4(f,{20,C},y20,-0.001);
y10 = y(end)
%% 
% Note that this is the same result we obtained earlier in Example #3.
%% See also
% <AB2_doc.html |AB2|> | 
% <AB3_doc.html |AB3|> | 
% <AB5_doc.html |AB5|> | 
% <AB6_doc.html |AB6|> | 
% <AB7_doc.html |AB7|> | 
% <AB8_doc.html |AB8|>