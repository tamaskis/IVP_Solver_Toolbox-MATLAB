%==========================================================================
%
% RK2_heun  Heun's method (Runge-Kutta second-order method).
%
%   [t,y] = RK2_heun(f,[t0,tf],y0,h)
%   [t,y] = RK2_heun(f,{t0,C},y0,h)
%
% See also euler, RK2, RK2_ralston, RK3, RK3_heun, RK3_ralston, SSPRK3, 
% RK4, RK4_ralston, RK4_38.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-09-06
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) function defining ODE dy/dt = f(t,y)
%   I       - defines interval over which to solve the ODE, 2 options:
%               --> [t0,tf] - (1×2 double) initial and final times
%               --> {t0,C}  - (1×2 cell) initial time and function handle
%                             for condition function C(t,y)
%   y0      - (n×1 double) initial condition
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   t       - (m×1 double) time vector
%   y       - (m×n double) matrix storing time history of state vector
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
function [t,y] = RK2_heun(f,I,y0,h)
    
    % ---------------------------------------
    % Determines which implementation to use.
    % ---------------------------------------
    
    % event detection implementation
    if iscell(I)
        t0 = I{1};
        C = I{2};
        implementation = 'event';
        
    % time detection implementation
    else
        t0 = I(1);
        tf = I(2);
        implementation = 'time';
    end
        
    % --------------------------------------------------------
    % Time detection implementation (solves until final time).
    % --------------------------------------------------------
    
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
            k2 = f(ti+h,yi+h*k1);

            % state vector propogated to next iteration
            y(:,i+1) = yi+(h/2)*(k1+k2);

        end
        
        % linearly interpolates for solution at tf
        y(:,N+1) = y(:,N)+((y(:,N+1)-y(:,N))/(t(N+1)-t(N)))*(tf-t(N));
        
        % replaces last element of "t" with "tf"
        t(N+1) = tf;
    
    % ---------------------------------------------------------------------
    % Event detection implementation (solves while condition is satisfied).
    % ---------------------------------------------------------------------
    
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
            k2 = f(ti+h,yi+h*k1);

            % state vector propogated to next iteration
            y(:,i+1) = yi+(h/2)*(k1+k2);
            
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