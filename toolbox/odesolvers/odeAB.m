%==========================================================================
%
% odeAB  Fixed-step ODE solver using the Adams-Bashforth methods.
%
%   [t,y] = odeAB(f,[t0,tf],y0,h)
%   [t,y] = odeAB(f,{t0,C},y0,h)
%   [t,y] = odeAB(__,method)
%   [t,y] = odeAB(__,method,wb)
%
% See also odeRK, odeABM.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-04
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/ODE_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖ → ℝᵖ) defining ODE
%   I       - defines interval over which to solve the ODE, 2 options:
%               --> [t0,tf] - (1×2 double) initial and final times
%               --> {t0,C}  - (1×2 cell) initial time, t₀, and function 
%                             handle for condition function, C(t,y) 
%                             (C : ℝ×ℝᵖ → 𝔹)
%   y0      - (p×1 double) initial condition, y₀ = y(t₀)
%   h       - (1×1 double) step size
%   method  - (char) (OPTIONAL) Adams-Bashforth method --> 'AB2', 'AB3', 
%             'AB4', 'AB5', 'AB6', 'AB7', 'AB8' (defaults to 'AB8')
%   wb      - (1×1 logical or char) (OPTIONAL) waitbar parameters
%               --> input as "true" if you want waitbar with default 
%                   message displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%
% -------
% OUTPUT:
% -------
%   t       - ((N+1)×1 double) time vector
%   y       - ((N+1)×p double) matrix storing time history of state vector
%
% -----
% NOTE:
% -----
%   --> p = dimension of state vector (for the scalar case, p = 1)
%   --> N+1 = length of time vector
%   --> The ith row of "y" is the TRANSPOSE of the state vector (i.e. the
%       solution corresponding to the ith time in "t"). This convention is
%       chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function [t,y] = odeAB(f,I,y0,h,method,wb)
    
    % -------------------
    % Setting up waitbar.
    % -------------------
    
    % turns waitbar on with custom message
    if (nargin == 6) && ischar(wb)
        display_waitbar = true;
        [wb,prop] = initialize_waitbar(msg);
        
    % turns waitbar on with default message
    elseif (nargin == 6) && wb
        display_waitbar = true;
        [wb,prop] = initialize_waitbar('Solving ODE...');
        
    % turns waitbar off otherwise
    else
        display_waitbar = false;
        
    end
    
    % ---------------------------------------
    % Determines which implementation to use.
    % ---------------------------------------
    
    % event detection implementation
    if iscell(I)
        t0 = I{1};
        C = I{2};
        time_detection_implementation = false;
        
    % time detection implementation
    else
        t0 = I(1);
        tf = I(2);
        time_detection_implementation = true;
        
    end
    
    % --------------------
    % Sets up integration.
    % --------------------
    
    % defaults integration method to AB8
    if (nargin < 5) || isempty(method)
        method = 'AB8';
    end
    
    % makes step size negative if t0 > tf
    if time_detection_implementation && (t0 > tf)
        h = -h;
    end
    
    % sets propagation function and order of the Adams-Bashforth method
    if strcmpi(method,'AB8')
        propagate = @(t,F) AB8(f,t,F,h);
        m = 8;
    elseif strcmpi(method,'AB2')
        propagate = @(t,F) AB2(f,t,F,h);
        m = 2;
    elseif strcmpi(method,'AB3')
        propagate = @(t,F) AB3(f,t,F,h);
        m = 3;
    elseif strcmpi(method,'AB4')
        propagate = @(t,F) AB4(f,t,F,h);
        m = 4;
    elseif strcmpi(method,'AB5')
        propagate = @(t,F) AB5(f,t,F,h);
        m = 5;
    elseif strcmpi(method,'AB6')
        propagate = @(t,F) AB6(f,t,F,h);
        m = 6;
    elseif strcmpi(method,'AB7')
        propagate = @(t,F) AB7(f,t,F,h);
        m = 7;
    end
    
    % --------------------------------------------------------
    % Time detection implementation (solves until final time).
    % --------------------------------------------------------
    
    if time_detection_implementation
        
        % number of subintervals between sample times
        N = ceil((tf-t0)/h);
        
        % last element of the time vector
        tN = t0+N*h;
        
        % defines time vector and preallocates solution matrix
        t = (t0:h:tN)';
        y = zeros(length(y0),length(t));
        
        % stores initial condition in solution matrix
        y(:,1) = y0;
        
        % propagating state vector using RK4 for first "m" sample times
        for n = 1:m
            
            % state vector propagated to next sample time
            y(:,n+1) = RK4(f,t(n),y(:,n),h);
            
            % updates waitbar
            if display_waitbar, prop = update_waitbar(n,N,wb,prop); end
            
        end
        
        % initializes F matrix (stores function evaluations for first m 
        % sample times in first m columns, state vector at (m+1)th sample
        % time in (m+1)th column)
        F = zeros(length(y0),m+1);
        for n = 1:m
            F(:,n) = f(t(n),y(:,n));
        end
        F(:,n+1) = y(:,m+1);
        
        % propagating state vector using an Adams-Bashforth method
        for n = (m+1):N
            
            % updates F matrix
            F = propagate(t(n),F);
            
            % extracts/stores state vector propagated to next sample time
            y(:,n+1) = F(:,m+1);
            
            % updates waitbar
            if display_waitbar, prop = update_waitbar(n,N,wb,prop); end
            
        end
        
        % linearly interpolates for solution at tf
        y(:,N+1) = y(:,N)+((y(:,N+1)-y(:,N))/(t(N+1)-t(N)))*(tf-t(N));
        
        % replaces last element of "t" with tf
        t(N+1) = tf;
        
        % closes waitbar
        if display_waitbar, close(wb); end
        
    % ---------------------------------------------------------------------
    % Event detection implementation (solves while condition is satisfied).
    % ---------------------------------------------------------------------
    
    else
        
        % preallocates time vector and solution matrix
        t = zeros(10000,1);
        y = zeros(length(y0),length(t));
        
        % time vector for first m sample times
        t(1:m) = (t0:h:(t0+(m-1)*h))';
        
        % stores initial condition in solution matrix
        t(1) = t0;
        y(:,1) = y0;
        
        % propagating state vector using RK4 for first "m" sample times
        for n = 1:m
            y(:,n+1) = RK4(f,t(n),y(:,n),h);
        end
        
        % initializes F matrix (stores function evaluations for first m 
        % sample times in first m columns, state vector at (m+1)th sample
        % time in (m+1)th column)
        F = zeros(length(y0),m+1);
        for n = 1:m
            F(:,n) = f(t(n),y(:,n));
        end
        F(:,n+1) = y(:,m+1);
        
        % state vector propagation while condition is satisfied
        n = 1;
        while C(t(n),y(:,n))
            
            % expands t and y if needed
            if (n+1) > length(t)
                [t,y] = expand_solution_arrays(t,y);
            end
            
            % updates F matrix
            F = propagate(t(n),F);
            
            % extracts/stores state vector propagated to next sample time
            y(:,n+1) = F(:,m+1);
            
            % increments time and loop index
            t(n+1) = t(n)+h;
            n = n+1;
            
        end
        
        % trims arrays
        y = y(:,1:(n-1));
        t = t(1:(n-1));
        
    end
    
    % transposes solution matrix so it is returned in "standard form"
    y = y';
    
end