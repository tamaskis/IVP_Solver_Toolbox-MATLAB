%==========================================================================
%
% solve_ivp  Fixed-step IVP solvers for solving vector-valued initial value 
% problems.
%
%   [t,y] = solve_ivp(f,[t0,tf],y0,h)
%   [t,y] = solve_ivp(f,{t0,E},y0,h)
%   [t,y] = solve_ivp(__,method)
%   [t,y] = solve_ivp(__,method,wb)
%
% See also solve_ivp_matrix.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-09-18
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) multivariate, vector-valued function 
%             defining vector-valued ODE, dy/dt = f(t,y) (f : ℝ×ℝᵖ → ℝᵖ)
%   I       - defines interval over which to solve the IVP, 2 options:
%               --> [t0,tf] - (1×2 double) initial and final times
%               --> {t0,E}  - (1×2 cell) initial time, t₀, and function 
%                             handle for event function, E(t,y) 
%                             (E : ℝ×ℝᵖ → B)
%   y0      - (p×1 double) initial condition, y₀ = y(t₀)
%   h       - (1×1 double) step size
%   method  - (OPTIONAL) (char) integration method --> 'Euler', 'RK2', 
%             'RK2 Heun', 'RK2 Ralston', 'RK3', 'RK3 Heun', 'RK3 Ralston', 
%             'SSPRK3', 'RK4', 'RK4 Ralston', 'RK4 3/8', 'AB2', 'AB3', 
%             'AB4', 'AB5', 'AB6', 'AB7', 'AB8', 'ABM2', 'ABM3', 'ABM4', 
%             'ABM5', 'ABM6', 'ABM7', 'ABM8' (defaults to 'RK4')
%   wb      - (OPTIONAL) (1×1 logical or char) waitbar parameters (defaults
%             to false)
%               --> input as true if you want waitbar with default message
%                   displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%
% -------
% OUTPUT:
% -------
%   t       - ((N+1)×1 double) time vector
%   y       - ((N+1)×p double) solution matrix
%
% -----
% NOTE:
% -----
%   --> The nth row of "y" stores the TRANSPOSE of the state vector (i.e. 
%       the solution) corresponding to the nth time in "t". This convention
%       is chosen to match the convention used by MATLAB's ODE suite.
%
%==========================================================================
function [t,y] = solve_ivp(f,I,y0,h,method,wb)
    
    % --------------------
    % Time detection mode.
    % --------------------
    
    if ~iscell(I)
        
        % extracts initial and final times
        t0 = I(1);
        tf = I(2);
        
        % makes step size negative if t0 > tf
        if (t0 > tf)
            h = -h;
        end
        
        % defines event function
        if (tf > t0)
            E = @(t,y) t-h <= tf;
        else
            E = @(t,y) t-h >= tf;
        end
        
        % indicates that final time is known
        final_time_known = true;
        
        % number of subintervals between sample times
        N = ceil((tf-t0)/h);
        
        % turns waitbar on with custom message
        if (nargin == 6) && ~isempty(wb) && ischar(wb)
            display_waitbar = true;
            [wb,prop] = initialize_waitbar(wb);
            
        % turns waitbar on with default message
        elseif (nargin == 6) && ~isempty(wb) && wb
            display_waitbar = true;
            [wb,prop] = initialize_waitbar('Solving IVP...');
            
        % turns waitbar off otherwise
        else
            display_waitbar = false;
            
        end
        
    % ---------------------
    % Event detection mode.
    % ---------------------
    
    else
        
        % extracts initial time and event function
        t0 = I{1};
        E = I{2};
        
        % indicates that final time is unknown
        final_time_known = false;
        
    end
    
    % -------------------
    % Integration method.
    % -------------------
    
    % defaults integration method to RK4
    if (nargin < 5) || isempty(method)
        method = 'RK4';
    end
    
    % sets propagation function for single-step methods
    if strcmpi(method,'Euler')
        g = @(t,y) RK1_euler(f,t,y,h);
    elseif strcmpi(method,'RK2')
        g = @(t,y) RK2(f,t,y,h);
    elseif strcmpi(method,'RK2 Heun')
        g = @(t,y) RK2_heun(f,t,y,h);
    elseif strcmpi(method,'RK2 Ralston')
        g = @(t,y) RK2_ralston(f,t,y,h);
    elseif strcmpi(method,'RK3')
        g = @(t,y) RK3(f,t,y,h);
    elseif strcmpi(method,'RK3 Heun')
        g = @(t,y) RK3_heun(f,t,y,h);
    elseif strcmpi(method,'RK3 Ralston')
        g = @(t,y) RK3_ralston(f,t,y,h);
    elseif strcmpi(method,'SSPRK3')
        g = @(t,y) SSPRK3(f,t,y,h);
    elseif strcmpi(method,'RK4')
        g = @(t,y) RK4(f,t,y,h);
    elseif strcmpi(method,'RK4 3/8')
        g = @(t,y) RK4_38(f,t,y,h);
    elseif strcmpi(method,'RK4 Ralston')
        g = @(t,y) RK4_ralston(f,t,y,h);
    end
    
    % sets propagation function for multi-step predictor methods
    if strcmpi(method,'AB2')
        g = @(t,F) AB2(f,t,F,h);
    elseif strcmpi(method,'AB3')
        g = @(t,F) AB3(f,t,F,h);
    elseif strcmpi(method,'AB4')
        g = @(t,F) AB4(f,t,F,h);
    elseif strcmpi(method,'AB5')
        g = @(t,F) AB5(f,t,F,h);
    elseif strcmpi(method,'AB6')
        g = @(t,F) AB6(f,t,F,h);
    elseif strcmpi(method,'AB7')
        g = @(t,F) AB7(f,t,F,h);
    elseif strcmpi(method,'AB8')
        g = @(t,F) AB8(f,t,F,h);
    end
    
    % sets propagation function for multi-step predictor-corrector methods
    if strcmpi(method,'ABM2')
        g = @(t,F) ABM2(f,t,F,h);
    elseif strcmpi(method,'ABM3')
        g = @(t,F) ABM3(f,t,F,h);
    elseif strcmpi(method,'ABM4')
        g = @(t,F) ABM4(f,t,F,h);
    elseif strcmpi(method,'ABM5')
        g = @(t,F) ABM5(f,t,F,h);
    elseif strcmpi(method,'ABM6')
        g = @(t,F) ABM6(f,t,F,h);
    elseif strcmpi(method,'ABM7')
        g = @(t,F) ABM7(f,t,F,h);
    elseif strcmpi(method,'ABM8')
        g = @(t,F) ABM8(f,t,F,h);
    end
    
    % determines if the method is a single-step method
    if ~strcmpi(method(1:2),'AB')
        single_step = true;
    else
        single_step = false;
    end
    
    % determines order for multistep methods
    if ~single_step
        m = str2double(method(end));
    end
    
    % --------------------------------------------------
    % Preallocates arrays and stores initial conditions.
    % --------------------------------------------------
    
    % state dimension
    p = length(y0);
    
    % preallocates time vector and solution matrix
    t = zeros(10000,1);
    y = zeros(p,10000);
    
    % stores initial conditions
    t(1) = t0;
    y(:,1) = y0;
    
    % ---------------------------------------
    % Solves IVP (using single-step methods).
    % ---------------------------------------
    
    if single_step
        
        % state vector propagation while event function is satisfied
        n = 1;
        while E(t(n),y(:,n))
            
            % expands t and y if needed
            if (n+1) > length(t)
                [t,y] = expand_ivp_arrays(t,y);
            end
            
            % propagates state vector to next sample time
            y(:,n+1) = g(t(n),y(:,n));
            
            % increments time and loop index
            t(n+1) = t(n)+h;
            n = n+1;
            
            % updates waitbar
            if final_time_known && display_waitbar
                prop = update_waitbar(n,N,wb,prop);
            end
            
        end
        
    % --------------------------------------
    % Solves IVP (using multi-step methods).
    % --------------------------------------
    
    else
        
        % time vector for first m+1 sample times
        t(1:(m+1)) = (t0:h:(t0+m*h)).';
        
        % propagating state vector using RK4 for first "m" sample times
        for n = 1:m
            y(:,n+1) = RK4(f,t(n),y(:,n),h);
        end
        
        % initializes F matrix (stores function evaluations for first m 
        % sample times in first m columns, state vector at (m+1)th sample
        % time in (m+1)th column)
        F = zeros(p,m+1);
        for n = 1:m
            F(:,n) = f(t(n),y(:,n));
        end
        F(:,n+1) = y(:,m+1);
        
        % state vector propagation while event function is satisfied
        n = m+1;
        while E(t(n),y(:,n))
            
            % expands t and y if needed
            if (n+1) > length(t)
                [t,y] = expand_ivp_arrays(t,y);
            end
            
            % updates F matrix (propagates to next sample time)
            F = g(t(n),F);
            
            % extracts/stores state vector at next sample time
            y(:,n+1) = F(:,m+1);
            
            % increments time and loop index
            t(n+1) = t(n)+h;
            n = n+1;
            
            % updates waitbar
            if final_time_known && display_waitbar
                prop = update_waitbar(n,N,wb,prop);
            end
            
        end
        
    end
    
    % -----------------
    % Final formatting.
    % -----------------
    
    % trims arrays
    y = y(:,1:(n-1));
    t = t(1:(n-1));
    
    % number of subintervals
    N = length(t)-1;
    
    % linearly interpolates to find solution at desired final time
    if final_time_known
        
        % linearly interpolates for solution at tf
        y(:,N+1) = y(:,N)+((y(:,N+1)-y(:,N))/(t(N+1)-t(N)))*(tf-t(N));
        
        % replaces last element of "t" with tf
        t(N+1) = tf;
        
    end
    
    % deletes penultimate solution if at same time as last solution
    if (abs(t(N+1)-t(N)) < 1e-10)
        t(N) = [];
        y(:,N) = [];
    end
    
    % transposes solution matrix so it is returned in "standard form"
    y = y.';
    
    % closes waitbar
    if final_time_known && display_waitbar
        close(wb);
    end
    
end