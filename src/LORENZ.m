clear;clc;

% Lorenz parameters
rho = 28;
sigma = 10;
beta = 8/3;

% Lorenz equations in vector form
f = @(t,x) [sigma*(x(2)-x(1));
            x(1)*(rho-x(3))-x(2);
            x(1)*x(2)-beta*x(3)];

% solution
[t,y] = ABM8(f,[0,100],[0;1;1.05],0.001);


figure;
plot3(y(:,1),y(:,2),y(:,3));
view(45,20);
axis off;
