%==========================================================================
%
% RK4_ralston  Ralston's fourth-order method (Runge-Kutta fourth-order 
% method).
%
%   [t,y] = RK4_ralston(f,[t0,tf],y0,h)
%   [t,y] = RK4_ralston(f,{t0,C},y0,h)
%
% See also EDIT.
%
% Copyright © 2021 Tamas Kis
% Contact: tamas.a.kis@outlook.com
% Last Update: 2021-07-11
%
%--------------------------------------------------------------------------
%
% MATLAB Central File Exchange:
% GitHub: https://github.com/tamaskis/ode_solver_toolbox-MATLAB
%
% See [ADD LINK] for examples and "DOCUMENTATION.pdf" for additional 
% documentation.
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f           - (function handle) function defining ODE dy/dt = f(t,y)
%   interval    - defines interval over which to solve the ODE, 2 options:
%                   --> [t0,tf] - (1×2) initial and final times
%                   --> {t0,C}  - (1×2 cell) initial time and function 
%                                 handle for condition function C(t,y)
%   y0          - (n×1) initial condition
% 	h           - (1×1) step size
%
% -------
% OUTPUT:
% -------
%   t       (m×1) time vector
%   y       (m×n) matrix storing time history of state vector
%
% -----
% NOTE:
% -----
%   --> n = dimension of state vector
%   --> m = length of time vector
%   --> The ith row of "y" is the TRANSPOSE of the state vector (i.e. the
%       solution corresponding to the ith time in "t"). This convention is
%       chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function [t,y] = RK4_ralston(f,domain,y0,h)
    
    % determines type of implementation based on passed parameters
    if iscell(domain)
        t0 = domain{1};
        C = domain{2};
        implementation = 'event';
    else
        t0 = domain(1);
        tf = domain(2);
        implementation = 'time';
    end
        
    % time detection implementation (solves until final time)
    if strcmp(implementation,'time')

        % number of subintervals between iterations
        N = ceil((tf-t0)/h);
        
        % last element of the time vector
        tN = t0+N*h;
        
        % defines time vector and preallocates solution matrix
        t = (t0:h:tN)';
        y = zeros(length(y0),length(t));

        % stores initial condition in solution matrix
        y(:,1) = y0;
        
        % propagating state vector
        for i = 1:N

            % solution at current iteration
            ti = t(i);
            yi = y(:,i);

            % k terms
            k1 = f(ti,yi);
            k2 = f(ti+0.4*h,yi+0.4*h*k1);
            k3 = f(ti+0.45573725*h,yi+0.29697761*h*k1+0.15875964*h*k2);
            k4 = f(ti+h,yi+0.21810040*h*k1-3.05096516*h*k2+3.83286476*h*...
                k3);

            % state vector propogated to next iteration
            y(:,i+1) = yi+h*(0.17476028*k1-0.55148066*k2+1.20553560*k3+...
                0.17118478*k4);

        end
        
        % linearly interpolates for solution at tf
        y(:,N+1) = y(:,N)+((y(:,N+1)-y(:,N))/(t(N+1)-t(N)))*(tf-t(N));
        
        % replaces last element of "t" with "tf"
        t(N+1) = tf;
    
    % event detection implementation (solves while condition is satisfied)
    else

        % preallocates time vector and solution matrix
        t = zeros(10000,1);
        y = zeros(length(y0),length(t));
        
        % stores initial condition in solution matrix
        t(1) = t0;
        y(:,1) = y0;

        % state vector propagation while condition is satisfied
        i = 1;
        while C(t(i),y(:,i))

            % expands t and y if needed
            if (i+1) > length(t)
                t = [t;zeros(size(t))];
                y = [y,zeros(size(y))];
            end
            
            % solution at current iteration
            ti = t(i);
            yi = y(:,i);

            % k terms
            k1 = f(ti,yi);
            k2 = f(ti+0.4*h,yi+0.4*h*k1);
            k3 = f(ti+0.45573725*h,yi+0.29697761*h*k1+0.15875964*h*k2);
            k4 = f(ti+h,yi+0.21810040*h*k1-3.05096516*h*k2+3.83286476*h*...
                k3);

            % state vector propogated to next iteration
            y(:,i+1) = yi+h*(0.17476028*k1-0.55148066*k2+1.20553560*k3+...
                0.17118478*k4);
            
            % increments time and loop index
            t(i+1) = t(i)+h;
            i = i+1;

        end

        % trims arrays
        y = y(:,1:(i-1));
        t = t(1:(i-1));
        
    end
    
    % transposes solution array so it is returned in "standard form"
    y = y';
    
end