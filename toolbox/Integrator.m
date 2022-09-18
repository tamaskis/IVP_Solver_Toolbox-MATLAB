%==========================================================================
%
% Integrator  Class defining an integrator.
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-17
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%==========================================================================

classdef Integrator < handle
    
    % -----------
    % Properties.
    % -----------
    
    properties
        dt          % (1×1 double) macro step size
        dtau        % (1×1 double) micro step size
        macro_step  % (1×1 function_handle) obtains state vector or state matrix at next sample time using a single macro step
        method      % (char) integration method (defaults to 'RK4')
    end
    
    % ---------------
    % Public methods.
    % ---------------
    
    methods (Access = public)
        
        function obj = Integrator(dt,method,dtau)
            % obj = Integrator(dt,method,dtau)
            %
            % Constructor.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   dt      - (1×1 double) major time step, Δt
            %   method  - (OPTIONAL) (char) integration method  --> 
            %             'Euler', 'RK2', 'RK2 Heun', 'RK2 Ralston', 'RK3',
            %             'RK3 Heun', 'RK3 Ralston', 'SSPRK3', 'RK4', 
            %             'RK4 Ralston', 'RK4 3/8', 'AB2', 'AB3', 'AB4', 
            %             'AB5', 'AB6', 'AB7', 'AB8', 'ABM2', 'ABM3', 
            %             'ABM4', 'ABM5', 'ABM6', 'ABM7', 'ABM8' (defaults 
            %             to 'RK4')
            %   dtau    - (OPTIONAL) (1×1 double) minor time step, Δτ
            %             (defaults to Δt/10)
            %   
            %
            % -------
            % OUTPUT:
            % -------
            %   obj     - (1×1 Integrator) integrator object
            %
            %--------------------------------------------------------------
            
            % defaults "method" to 'RK4' if not input
            if (nargin < 2) || isempty(method)
                obj.method = 'RK4';
            else
                obj.method = method;
            end
            
            % sets time step
            obj.dt = dt;
            
            % sets micro time step (defaults to Δt/10)
            if (nargin < 3) || isempty(dtau)
                obj.dtau = dt/10;
            else
                obj.dtau = dtau;
            end
            
            % sets up macro step function (NOTE: if the integration method 
            % is a multistep method, the macro step function will default 
            % to using RK4)
            if strcmpi(obj.method,'Euler')
                obj.macro_step = @(f,t,y) RK1_euler(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK2')
                obj.macro_step = @(f,t,y) RK2(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK2 Heun')
                obj.macro_step = @(f,t,y) RK2_heun(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK2 Ralston')
                obj.macro_step = @(f,t,y) RK2_ralston(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK3')
                obj.macro_step = @(f,t,y) RK3(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK3 Heun')
                obj.macro_step = @(f,t,y) RK3_heun(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK3 Ralston')
                obj.macro_step = @(f,t,y) RK3_ralston(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'SSPRK3')
                obj.macro_step = @(f,t,y) SSPRK3(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK4')
                obj.macro_step = @(f,t,t) RK4(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK4 Ralston')
                obj.macro_step = @(f,t,y) RK4_ralston(f,t,xk,obj.dt);
            elseif strcmpi(obj.method,'RK4 3/8')
                obj.macro_step = @(f,t,y) RK4_38(f,t,xk,obj.dt);
            else
                obj.macro_step = @(f,t,y) RK4(f,t,xk,obj.dt);
            end
            
        end
        
        function y_next = micro_step_vector(obj,f,t,y)
            % y_next = Integrator.micro_step_vector(f,t,y)
            %
            % Obtains state vector at next sample time using multiple micro
            % steps.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   f       - (1×1 function_handle) dy/dt = f(t,y) --> 
            %             multivariate, vector-valued function 
            %             (f : ℝ×ℝᵖ → ℝᵖ) defining ODE
            %   t       - (1×1 double) current time
            %   y       - (p×1 double) state vector at current sample time
            %
            % -------
            % OUTPUT:
            % -------
            %   y_next  - (p×1 double) state vector at next sample time
            %
            %--------------------------------------------------------------
            
            % solves IVP from t to t + Δt with a step size of Δτ
            [~,y_array] = solve_ivp(f,[t,t+obj.dt],y,obj.dtau,obj.method);
            
            % only keeps solution at t + Δt
            y_next = y_array(end,:).';
            
        end
        
        function M_next = micro_step_matrix(obj,F,t,M)
            % M_next = Integrator.micro_step_matrix(f,t,y)
            %
            % Obtains state matrix at next sample time using multiple micro
            % steps.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   F       - (1×1 function_handle) dM/dt = F(t,M) --> 
            %             multivariate, matrix-valued function 
            %             (F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ) defining matrix-valued ODE
            %   t       - (1×1 double) current time
            %   M       - (p×r double) state matrix at current sample time
            %
            % -------
            % OUTPUT:
            % -------
            %   M_next  - (p×r double) state matrix at next sample time
            %
            %--------------------------------------------------------------
            
            % solves IVP from t to t + Δt with a step size of Δτ
            [~,M_array] = solve_ivp_matrix(F,[t,t+obj.dt],M,obj.dtau,...
                obj.method);
            
            % only keeps solution at t + Δt
            M_next = M_array(:,:,end);
            
        end
        
    end
    
end