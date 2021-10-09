%==========================================================================
%
% ABM8  Adams-Bashforth-Moulton 8th-order method.
%
%   [t,y] = ABM8(f,[t0,tf],y0,h)
%   [t,y] = ABM8(f,{t0,C},y0,h)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-10-09
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) multivariate, vector-valued function 
%             (f:R(n+1)->Rn) defining the ODE dy/dt = f(t,y)
%               --> inputs to f are the current time (1×1 double) and the
%                   current state vector (n×1 double)
%               --> output of f is the state vector derivative at the
%                   current time/state (n×1 double)
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
%   t       - ((N+1)×1 double) time vector
%   y       - ((N+1)×n double) matrix storing time history of state vector
%
% -----
% NOTE:
% -----
%   --> n = dimension of state vector
%   --> N+1 = length of time vector
%   --> The ith row of "y" is the TRANSPOSE of the state vector (i.e. the
%       solution corresponding to the ith time in "t"). This convention is
%       chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function [t,y] = ABM8(f,I,y0,h)
    
    % ---------------------------------------
    % Determines which implementation to use.
    % ---------------------------------------
    
    if iscell(I)
        t0 = I{1};
        C = I{2};
        implementation = 'event';
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
        
        % defines time vector, preallocates solution matrix and array to 
        % store last 8 function evaluations
        t = (t0:h:tN)';
        y = zeros(length(y0),length(t));
        f8 = zeros(length(y0),8);

        % stores initial condition in solution matrix
        y(:,1) = y0;
        
        % propagating state vector using RK4 for first 8 iterations
        for i = 1:8

            % solution at current iteration
            ti = t(i);
            yi = y(:,i);

            % k terms
            k1 = f(ti,yi);
            k2 = f(ti+h/2,yi+h*k1/2);
            k3 = f(ti+h/2,yi+h*k2/2);
            k4 = f(ti+h,yi+h*k3);

            % dependent variable propogated to next time level
            y(:,i+1) = yi+(h/6)*(k1+2*k2+2*k3+k4);

        end
        
        % stores function evaluations for first 8 iterations
        for i = 1:8
            f8(:,i) = f(t(i),y(:,i));
        end
        
        % propagating state vector using ABM8
        i = 8;
        for j = 9:N
            
            % updates stored function evaluations
            f8 = [f8(:,2:end),f(t(j),y(:,j))];
            
            % predictor step
            yp = y(:,j)+(h/120960)*(434241*f8(:,i)-1152169*f8(:,i-1)+...
                2183877*f8(:,i-2)-2664477*f8(:,i-3)+2102243*f8(:,i-4)-...
                1041723*f8(:,i-5)+295767*f8(:,i-6)-36799*f8(:,i-7));
            
            % corrector step (state vector propagated to next iteration)
            y(:,j+1) = y(:,j)+(h/120960)*(36799*f(t(j+1),yp)+139849*...
                f8(:,i)-121797*f8(:,i-1)+123133*f8(:,i-2)-88547*f8(:,i-...
                3)+41499*f8(:,i-4)-11351*f8(:,i-5)+1375*f8(:,i-6));
            
        end
        
        % linearly interpolates for solution at tf
        y(:,N+1) = y(:,N)+((y(:,N+1)-y(:,N))/(t(N+1)-t(N)))*(tf-t(N));
        
        % replaces last element of "t" with tf
        t(N+1) = tf;
    
    % ---------------------------------------------------------------------
    % Event detection implementation (solves while condition is satisfied).
    % ---------------------------------------------------------------------
    
    else

        % preallocates time vector, solution matrix, and array to store 
        % last 8 function evaluations
        t = zeros(10000,1);
        y = zeros(length(y0),length(t));
        f8 = zeros(length(y0),8);
        
        % time vector for first 8 iterations
        t(1:8) = (t0:h:(t0+7*h))';
        
        % stores initial condition in solution matrix
        t(1) = t0;
        y(:,1) = y0;
        
        % propagating state vector using RK4 for first 8 iterations
        for i = 1:8

            % solution at current iteration
            ti = t(i);
            yi = y(:,i);

            % k terms
            k1 = f(ti,yi);
            k2 = f(ti+h/2,yi+k1*h/2);
            k3 = f(ti+h/2,yi+k2*h/2);
            k4 = f(ti+h,yi+k3*h);

            % dependent variable propogated to next time level
            y(:,i+1) = yi+(h/6)*(k1+2*k2+2*k3+k4);

        end

        % state vector propagation while condition is satisfied
        i = 8;
        j = 8;
        while C(t(j),y(:,j))
        
            % expands t and y if needed
            if (j+1) > length(t)
                [t,y] = expand_solution_arrays(t,y);
            end
            
            % updates stored function evaluations
            f8 = [f8(:,2:end),f(t(j),y(:,j))];
            
            % predictor step
            yp = y(:,j)+(h/120960)*(434241*f8(:,i)-1152169*f8(:,i-1)+...
                2183877*f8(:,i-2)-2664477*f8(:,i-3)+2102243*f8(:,i-4)-...
                1041723*f8(:,i-5)+295767*f8(:,i-6)-36799*f8(:,i-7));
            
            % corrector step (state vector propagated to next iteration)
            y(:,j+1) = y(:,j)+(h/120960)*(36799*f(t(j+1),yp)+139849*...
                f8(:,i)-121797*f8(:,i-1)+123133*f8(:,i-2)-88547*f8(:,i-...
                3)+41499*f8(:,i-4)-11351*f8(:,i-5)+1375*f8(:,i-6));
            
            % increments time and loop index
            t(j+1) = t(j)+h;
            j = j+1;

        end

        % trims arrays
        y = y(:,1:(j-1));
        t = t(1:(j-1));
        
    end
    
    % transposes solution array so it is returned in "standard form"
    y = y';
    
    % -------------
    % Subfunctions.
    % -------------
    
    %----------------------------------------------------------------------
    % expand_solution_arrays
    %
    % Expands the arrays storing the ODE solution. This function is used to
    % preallocate more space for the solution arrays (i.e. the time vector
    % and the solution matrix) once they have become full.
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   t       - ((N+1)×1 double) time vector
    %  	y       - (n×(N+1) double) solution matrix
    %
    % OUTPUT:
    %   t_new   - (2(N+1)×1 double) expanded time vector
    %  	y_new   - (n×2(N+1) double) expanded solution matrix
    %
    % NOTE:
    %   --> n = dimension of state vector
    %   --> N+1 = length of time vector
    %
    %----------------------------------------------------------------------
    function [t_new,y_new] = expand_solution_arrays(t,y)
        t_new = [t;zeros(2*length(t),1)];
        y_new = [y,zeros(size(y,1),2*length(t))];
    end

end