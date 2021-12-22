%% Matrix ODE Example
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Problem Statement
% Consider the Riccati differential equation for the finite-horizon linear
% quadratic regulator problem:
%
% $$\dot{\mathbf{P}}=-\left[\mathbf{A}^{T}\mathbf{P}+\mathbf{P}\mathbf{A}-(\mathbf{P}\mathbf{B}+\mathbf{N})\mathbf{R}^{-1}(\mathbf{B}^{T}\mathbf{P}+\mathbf{N}^{T})+\mathbf{Q}\right]$$
%
% Find $\mathbf{P}(0)=\mathbf{P}_{0}$ given that
%
% $$\mathbf{P}(T)=\mathbf{P}_{T}=\pmatrix{1&1\cr1&1}$$
%
% where $T=5$. The state ($\mathbf{A}$) and input ($\mathbf{B}$) matrices
% are
%
% $$\mathbf{A}=\pmatrix{1&1\cr2&1},\quad\quad\mathbf{B}=\pmatrix{1\cr1}$$
%
% The cross-coupling ($\mathbf{N}$), state ($\mathbf{Q}$), and input
% ($\mathbf{R}$) weighting matrices are
%
% $$\mathbf{N}=\pmatrix{0&0\cr0&0},\quad\quad\mathbf{Q}=\pmatrix{2&1\cr1&1},\quad\quad\mathbf{R}=\pmatrix{1&0\cr0&1}$$

%% Matrices

% state matrix
A = [1   1;
     2   1];

% input matrix
B = [1;
     1];

% cross-coupling weighting matrix
N = [0  0;
     0  0];

% state weighting matrix
Q = [2   1;
     1   1];

% input weighting matrix
R = [1   0;
     0   1];
%% Final Condition

% final condition
PT = [1   1;
      1   1];

% final time
T = 5;
%% |odefun_mat2vec|
% Defining the Riccati differential equation (a matrix-valued ODE),
F = @(t,P) -(A.'*P+P*A-(P*B+N)/R*(B.'*P+N.')+Q);
%%
% Converting this matrix-valued ODE to a vector-valued ODE,
f = odefun_mat2vec(F);
%%
% Note that since $\mathbf{P}$ is a square matrix, we did not have to
% specify its number of rows for the |odefun_mat2vec| function.
%% |odeIC_mat2vec|
% In this case, we have a final condition at time $t=T=5$, and want to find
% the initial condition at time $t=t_{0}=0$. To do this, we will integrate
% backwards using an ODE solver, so the final condition $\mathbf{P}_{T}$
% will actually form the initial condition. Therefore, finding the final
% condition for the corresponding vector-valued ODE,
yT = odeIC_mat2vec(PT);
%%
% Note that since $\mathbf{P}$ is a square matrix, we did not have to
% specify its number of rows for the |odeIC_mat2vec| function.
%% |odesol_vec2mat|
% First, let's solve the vector-valued ODE for $\mathbf{y}(t)$. Solving the
% ODE from $t=T$ to $t=0$,
[t,y] = ode45(f,[T,0],yT);
%%
% Transforming the solution matrix (|y|) for the vector-valued ODE into the
% solution _array_ (|P|) for the matrix-valued ODE,
P = odesol_vec2mat(y);
%%
% Note that since $\mathbf{P}$ is a square matrix, we did not have to
% specify its number of rows for the |odesol_vec2mat| function.
%% Solution for $\mathbf{P}_{0}$
% Our original goal was to find $\mathbf{P}_{0}$. Extracting this matrix
% from the solution array (noting that the solution corresponding to $t=0$
% will be stored at the end of the solution array since we integrated
% backwards in time), we find
P0 = P(:,:,end)
%% References
% This example is adapted from the following sources:
%
% * https://www.mathworks.com/matlabcentral/answers/94722-how-can-i-solve-the-matrix-riccati-differential-equation-within-matlab
% * https://en.wikipedia.org/wiki/Linear-quadratic_regulator