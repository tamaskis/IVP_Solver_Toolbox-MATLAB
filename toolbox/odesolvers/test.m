clear;clc;close all;


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

% Solving the system for $t$ in the interval $[0,100]$ using the AB8 method 
% with a time step (i.e. step size) of $\Delta t=0.001$,
%tic
[t,x] = odeRK(f,[0,100],x0,0.001,'RK1_euler',true);
%toc

tic
C = @(t,x) t <= 100;
toc

%tic
[t,x] = odeRK(f,{0,C},x0,0.001);
%toc

% Plotting the solution,
figure;
plot3(x(:,1),x(:,2),x(:,3));
view(45,20);
grid on;
xlabel('$x$','Interpreter','latex','FontSize',18);
ylabel('$y$','Interpreter','latex','FontSize',18);
zlabel('$z$','Interpreter','latex','FontSize',18);
